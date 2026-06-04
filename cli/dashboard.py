"""Static-HTML status / coverage / missingness / provenance dashboard.

`gvf dashboard` reads the consolidated corpus (`corpus/INDEX.csv` + the per-paper
source under `corpus/<GENE>/<PMID>/`) and the scored SQLite DB(s), and writes an
offline static dashboard under `--out` (default `corpus/dashboard/`). Figures and
supplements are referenced relatively into the corpus (no duplication), so zip or
serve the whole `corpus/` to share it:

* `index.html` — per-gene cards: source coverage (usable vs stub), extraction
  funnel (in-corpus -> extracted -> has-variants), variant/unique counts.
* `<GENE>.html` — funnel, coverage-by-extraction-method facet, a
  provenance-completeness audit, the "what's left" (stub / no-variant) list, and
  a sortable/filterable paper table linking to per-paper adjudication pages.
* `paper_<PMID>.html` — the ADJUDICATION view: paper header with one-click links
  to PubMed / DOI / PMC, the extracted records (variant provenance +
  per-patient characteristics) on the left, and the EXACT on-disk full text the
  LLM parsed rendered on the right. Clicking a record jumps to and highlights
  the sentence/section it was extracted from, so you can verify it yourself.

Nothing here mutates the DB or hits the network; it is a read-only view.
"""

from __future__ import annotations

import csv
import html
import json
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Links / small helpers
# ---------------------------------------------------------------------------
def pubmed_url(pmid: str) -> str:
    return f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"


def pmc_url(pmcid: str) -> str:
    pmcid = pmcid.strip()
    if not pmcid.upper().startswith("PMC"):
        pmcid = "PMC" + pmcid
    return f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"


def doi_url(doi: str) -> str:
    return f"https://doi.org/{doi.strip()}"


def esc(x) -> str:
    return html.escape("" if x is None else str(x))


def jq(x) -> str:
    """A JS-string argument safe to embed inside a double-quoted onclick="...".

    json.dumps gives a valid JS string literal; html.escape then turns its
    double quotes into &quot; so they don't terminate the onclick attribute
    (the browser decodes &quot; back to " when it parses the handler).
    """
    return html.escape(json.dumps("" if x is None else str(x)))


def pct(num: int, den: int) -> float:
    return (100.0 * num / den) if den else 0.0


def _to_int(x) -> int:
    """Coerce a possibly-TEXT/None DB count to int (0 on junk)."""
    try:
        return int(float(str(x).strip()))
    except (TypeError, ValueError):
        return 0


_BOILERPLATE = (
    "skip to main content",
    "log in",
    "create a free account",
    "your email address",
    "forgot password",
    "advertisement",
    "cookie",
    "subscribe",
    "sign in",
    "institutional access",
    "get access",
    "add to mendeley",
    "reset your password",
)


def boilerplate_warning(text: str) -> str:
    """Flag sources that look like a captured page shell (login/nav/ads, thin
    body) rather than article text — i.e. acquisition, not extraction, failed."""
    low = text.lower()
    hits = sum(1 for m in _BOILERPLATE if m in low)
    body_chars = sum(len(ln) for ln in text.splitlines() if len(ln.strip()) > 40)
    if hits >= 4 and body_chars < 6000:
        return (
            "<div class='note'>⚠ This source looks boilerplate-heavy "
            f"({hits} login/nav/ad markers, ~{body_chars} body chars) — acquisition may "
            "have captured the page shell, not the article. Records below may be "
            "unverifiable in this text; this is a fetch problem, not an extraction one.</div>"
        )
    return ""


# ---------------------------------------------------------------------------
# Minimal, dependency-free Markdown -> HTML (faithful: every block is taggable
# so the adjudication JS can jump to and highlight the source sentence/row).
# ---------------------------------------------------------------------------
def md_to_html(text: str) -> str:
    lines = text.replace("\r\n", "\n").split("\n")
    out: list[str] = []
    i, n = 0, len(lines)
    para: list[str] = []
    hid = 0

    def flush_para():
        if para:
            joined = " ".join(s.strip() for s in para if s.strip())
            if joined:
                out.append(f'<p class="blk">{esc(joined)}</p>')
            para.clear()

    while i < n:
        line = lines[i]
        stripped = line.strip()
        # fenced code
        if stripped.startswith("```"):
            flush_para()
            i += 1
            buf = []
            while i < n and not lines[i].strip().startswith("```"):
                buf.append(lines[i])
                i += 1
            i += 1
            out.append(f'<pre class="blk">{esc(chr(10).join(buf))}</pre>')
            continue
        # ATX heading
        m = re.match(r"^(#{1,6})\s+(.*)$", stripped)
        if m:
            flush_para()
            level = min(len(m.group(1)) + 1, 6)
            txt = m.group(2)
            out.append(f'<h{level} class="blk hd" id="sec{hid}">{esc(txt)}</h{level}>')
            hid += 1
            i += 1
            continue

        # pipe table — detect a run of >=2 consecutive pipe rows, WITH or WITHOUT
        # a |---|---| separator (many converted sources omit the separator, which
        # is why tables were silently dropped before).
        def _is_sep(s: str) -> bool:
            return bool(re.match(r"^\s*\|?[\s:|-]+\|?\s*$", s)) and "-" in s

        def _cells(s: str) -> list[str]:
            return [c.strip() for c in s.strip().strip("|").split("|")]

        if (
            line.count("|") >= 2
            and i + 1 < n
            and (lines[i + 1].count("|") >= 2 or _is_sep(lines[i + 1]))
        ):
            flush_para()
            header = _cells(line)
            i += 1
            if i < n and _is_sep(lines[i]):
                i += 1  # skip the separator row if present
            rows = []
            while i < n and lines[i].count("|") >= 2 and lines[i].strip():
                if _is_sep(lines[i]):
                    i += 1
                    continue
                rows.append(_cells(lines[i]))
                i += 1
            thead = "".join(f"<th>{esc(c)}</th>" for c in header)
            bodyrows = "".join(
                "<tr class='blk'>" + "".join(f"<td>{esc(c)}</td>" for c in r) + "</tr>"
                for r in rows
            )
            out.append(
                f'<div class="tblwrap"><table class="srctbl"><thead><tr>{thead}</tr>'
                f"</thead><tbody>{bodyrows}</tbody></table></div>"
            )
            continue
        # inline figure marker (`_image_: <url>` or markdown ![alt](url)) — render
        # as a clickable figure chip (the downloaded figure gallery is appended
        # separately; remote URLs open in a new tab).
        mimg = re.match(r"^_image_:\s*(\S+)", stripped) or re.match(
            r"^!\[[^\]]*\]\((\S+?)\)", stripped
        )
        if mimg:
            flush_para()
            url = mimg.group(1)
            name = url.rstrip("/").rsplit("/", 1)[-1][:60]
            out.append(
                f"<div class='blk'><a class='figref' href='{esc(url)}' target='_blank'>"
                f"🖼 figure: {esc(name)}</a></div>"
            )
            i += 1
            continue
        # blank line
        if not stripped:
            flush_para()
            i += 1
            continue
        para.append(line)
        i += 1
    flush_para()
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_corpus_index(corpus_dir: Path) -> dict[str, dict[str, dict]]:
    """gene -> pmid -> {full_text_status, fulltext, n_figures, n_supplement_files, ...}."""
    idx_csv = corpus_dir / "INDEX.csv"
    by_gene: dict[str, dict[str, dict]] = defaultdict(dict)
    if not idx_csv.exists():
        return by_gene
    with idx_csv.open(newline="") as f:
        for row in csv.DictReader(f):
            by_gene[row["gene"]][row["pmid"]] = row
    return by_gene


def _table_cols(con: sqlite3.Connection, table: str) -> set[str]:
    try:
        return {r[1] for r in con.execute(f"PRAGMA table_info({table})")}
    except sqlite3.Error:
        return set()


def load_db(db_path: Path) -> dict:
    """Pull the provenance/coverage surface from a scored DB. Degrades gracefully."""
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    data: dict = {
        "papers": {},
        "by_pmid": defaultdict(lambda: {"variants": [], "patients": [], "tables": []}),
        "source_method": Counter(),
        "prov": Counter(),
        "variants_agg": defaultdict(
            lambda: {
                "sig": "",
                "pmids": set(),
                "affected": 0,
                "unaffected": 0,
                "carriers": 0,
            }
        ),
        # cohort counts from penetrance_data (the real patient/carrier numbers;
        # individual_records is per-person and far sparser)
        "pen_by_pmid": defaultdict(
            lambda: {"carriers": 0, "affected": 0, "unaffected": 0}
        ),
    }
    try:
        pcols = _table_cols(con, "papers")
        bib = [
            c
            for c in ("first_author", "journal", "publication_date", "doi", "pmc_id")
            if c in pcols
        ]
        for r in con.execute(
            "SELECT pmid, title" + ("," + ",".join(bib) if bib else "") + " FROM papers"
        ):
            rd = dict(r)
            data["papers"][str(rd["pmid"])] = {
                "title": rd.get("title"),
                "meta": {},
                "bib": {
                    k: rd.get(k)
                    for k in (
                        "first_author",
                        "journal",
                        "publication_date",
                        "doi",
                        "pmc_id",
                    )
                },
            }
        if _table_cols(con, "extraction_metadata"):
            for r in con.execute(
                "SELECT pmid, model_used, source_file, source_type, extraction_confidence, "
                "total_variants_found, notes, abstract_only FROM extraction_metadata"
            ):
                data["papers"].setdefault(
                    str(r["pmid"]), {"title": None, "meta": {}, "bib": {}}
                )
                data["papers"][str(r["pmid"])]["meta"] = dict(r)
        # variant provenance (count_provenance column added 2026-06-04; guard for old DBs)
        has_cp = "count_provenance" in _table_cols(con, "variant_papers")
        q = (
            "SELECT vp.pmid, v.protein_notation, v.cdna_notation, v.genomic_position, "
            "v.clinical_significance, vp.source_location, vp.additional_notes, vp.key_quotes"
            + (", vp.count_provenance" if has_cp else "")
            + " FROM variant_papers vp JOIN variants v ON v.variant_id = vp.variant_id"
        )
        for r in con.execute(q):
            pmid = str(r["pmid"])
            loc = r["source_location"] or ""
            vname = (
                r["protein_notation"]
                or r["cdna_notation"]
                or r["genomic_position"]
                or "?"
            )
            data["source_method"][_method_bucket(loc)] += 1
            data["prov"]["total"] += 1
            if loc.strip():
                data["prov"]["with_location"] += 1
            if (r["additional_notes"] or "").strip():
                data["prov"]["with_notes"] += 1
            quotes = _parse_quotes(r["key_quotes"])
            if quotes:
                data["prov"]["with_quotes"] += 1
            cprov = _parse_count_provenance(r["count_provenance"]) if has_cp else ""
            if cprov:
                data["prov"]["with_count_prov"] += 1
            agg = data["variants_agg"][vname]
            agg["sig"] = agg["sig"] or (r["clinical_significance"] or "")
            agg["pmids"].add(pmid)
            data["by_pmid"][pmid]["variants"].append(
                {
                    "variant": vname,
                    "cdna": r["cdna_notation"],
                    "sig": r["clinical_significance"],
                    "location": loc,
                    "notes": r["additional_notes"] or "",
                    "quotes": quotes,
                    "count_prov": cprov,
                }
            )
        # patient characteristics (ethnicity/geographic_origin added 2026-06-04)
        ir_cols = _table_cols(con, "individual_records")
        if ir_cols:
            extra = [c for c in ("ethnicity", "geographic_origin") if c in ir_cols]
            sel = (
                "ir.pmid, v.protein_notation, ir.individual_id, ir.age_at_evaluation, "
                "ir.age_at_onset, ir.sex, ir.affected_status, ir.phenotype_details, ir.evidence_sentence"
                + ("," + ",".join("ir." + c for c in extra) if extra else "")
            )
            for r in con.execute(
                f"SELECT {sel} FROM individual_records ir "
                "LEFT JOIN variants v ON v.variant_id = ir.variant_id"
            ):
                rd = dict(r)
                vname = rd.get("protein_notation") or ""
                data["by_pmid"][str(rd["pmid"])]["patients"].append(
                    {
                        "variant": vname,
                        "id": rd.get("individual_id") or "",
                        "age": rd.get("age_at_evaluation")
                        if rd.get("age_at_evaluation") is not None
                        else rd.get("age_at_onset"),
                        "sex": rd.get("sex") or "",
                        "status": rd.get("affected_status") or "",
                        "phenotype": rd.get("phenotype_details") or "",
                        "ethnicity": rd.get("ethnicity") or "",
                        "origin": rd.get("geographic_origin") or "",
                        "evidence": rd.get("evidence_sentence") or "",
                    }
                )
        # cohort counts (penetrance_data) — the substantive carrier/affected/
        # unaffected numbers per paper and per variant.
        if _table_cols(con, "penetrance_data"):
            for r in con.execute(
                "SELECT pd.pmid, v.protein_notation, "
                "COALESCE(pd.total_carriers_observed,0) c, COALESCE(pd.affected_count,0) a, "
                "COALESCE(pd.unaffected_count,0) u "
                "FROM penetrance_data pd LEFT JOIN variants v ON v.variant_id = pd.variant_id"
            ):
                pmid = str(r["pmid"])
                c, a, u = _to_int(r["c"]), _to_int(r["a"]), _to_int(r["u"])
                data["pen_by_pmid"][pmid]["carriers"] += c
                data["pen_by_pmid"][pmid]["affected"] += a
                data["pen_by_pmid"][pmid]["unaffected"] += u
                vname = r["protein_notation"]
                if vname:
                    pv = data.setdefault("pen_pv", {}).setdefault(
                        (pmid, vname), {"carriers": 0, "affected": 0, "unaffected": 0}
                    )
                    pv["carriers"] += c
                    pv["affected"] += a
                    pv["unaffected"] += u
                    if vname in data["variants_agg"]:
                        agg = data["variants_agg"][vname]
                        agg["carriers"] += c
                        agg["affected"] += a
                        agg["unaffected"] += u
        if _table_cols(con, "tables_processed"):
            for r in con.execute(
                "SELECT pmid, table_name, table_caption, variants_extracted FROM tables_processed"
            ):
                data["by_pmid"][str(r["pmid"])]["tables"].append(dict(r))
    finally:
        con.close()
    return data


def _method_bucket(loc: str) -> str:
    low = loc.lower()
    if "clinvar" in low:
        return "ClinVar"
    if "pubtator" in low:
        return "PubTator"
    if "figure" in low:
        return "Figure"
    if "supplement" in low or "suppl" in low:
        return "Supplement"
    if "abstract" in low:
        return "Abstract"
    if "table" in low:
        return "Table"
    if "text scan" in low or "text" in low:
        return "Text"
    return "Other"


def _parse_quotes(raw) -> list[str]:
    if not raw:
        return []
    try:
        v = json.loads(raw)
        if isinstance(v, list):
            return [str(x) for x in v if str(x).strip()]
        if isinstance(v, str) and v.strip():
            return [v]
    except (json.JSONDecodeError, TypeError):
        if str(raw).strip():
            return [str(raw)]
    return []


def _parse_count_provenance(raw) -> str:
    """Render count_provenance JSON as a short 'why this count' string."""
    if not raw:
        return ""
    try:
        d = json.loads(raw) if isinstance(raw, str) else raw
    except (json.JSONDecodeError, TypeError):
        return str(raw)
    if not isinstance(d, dict):
        return str(raw)
    bits = []
    for kind in ("carriers", "affected", "unaffected"):
        label = d.get(f"{kind}_column_label")
        ctype = d.get(f"{kind}_count_type")
        if label or (ctype and ctype not in ("unknown", None)):
            piece = kind
            if label:
                piece += f" from “{label}”"
            if ctype and ctype != "unknown":
                piece += f" ({ctype})"
            bits.append(piece)
    return "; ".join(bits)


def _artifacts_links(corpus_dir: Path, gene: str, pmid: str) -> dict:
    art = corpus_dir / gene / pmid / f"{pmid}_artifacts.json"
    out = {"doi": None, "pmcid": None}
    if art.exists():
        try:
            d = json.loads(art.read_text())
            out["doi"] = d.get("doi")
            out["pmcid"] = d.get("pmcid")
        except (json.JSONDecodeError, OSError):
            pass
    return out


# ---------------------------------------------------------------------------
# HTML scaffolding (self-contained: inline CSS + JS, no network assets)
# ---------------------------------------------------------------------------
CSS = """
:root{--bg:#eef1f6;--card:#ffffff;--ink:#1b2430;--mut:#5b6776;--ok:#15a34a;--warn:#b45309;--bad:#dc2626;--acc:#2563eb;--acc2:#0e7490;--line:#dde3ec;--hl:#fff3b0}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.55 -apple-system,Segoe UI,Roboto,Helvetica,sans-serif}
a{color:var(--acc);text-decoration:none}a:hover{text-decoration:underline}
header{padding:18px 24px;background:linear-gradient(180deg,#fff,#f3f6fb);border-bottom:1px solid var(--line)}
h1{margin:0;font-size:21px;color:#0f172a}.sub{color:var(--mut);font-size:13px}
.wrap{padding:20px 24px;max-width:1500px;margin:0 auto}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(300px,1fr));gap:14px}
.card{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:16px;box-shadow:0 1px 3px rgba(15,23,42,.06)}
.card h2{margin:0 0 8px;font-size:16px;color:#0f172a}.kv{display:flex;justify-content:space-between;color:var(--mut);font-size:13px;padding:3px 0}.kv b{color:var(--ink)}
.bar{height:9px;background:#e9edf3;border-radius:6px;overflow:hidden;margin:7px 0}.bar>span{display:block;height:100%}
.bar .ok{background:linear-gradient(90deg,#22c55e,#15a34a)}.bar .warn{background:var(--warn)}
.tag{display:inline-block;font-size:11px;padding:2px 8px;border-radius:20px;border:1px solid var(--line);color:var(--mut);margin:2px 3px 0 0;background:#f8fafc}
.tag.ok{color:#fff;background:var(--ok);border-color:var(--ok)}.tag.bad{color:#fff;background:var(--bad);border-color:var(--bad)}.tag.warn{color:#fff;background:var(--warn);border-color:var(--warn)}
table{width:100%;border-collapse:collapse;font-size:13px;background:#fff}
th,td{text-align:left;padding:7px 9px;border-bottom:1px solid var(--line);vertical-align:top}
th{position:sticky;top:0;background:#eef2f8;color:#0f172a;cursor:pointer;user-select:none;font-weight:600}
tbody tr:nth-child(even) td{background:#f8fafc}tr:hover td{background:#eef5ff}
.search{width:100%;padding:9px 11px;margin:8px 0;background:#fff;border:1px solid var(--line);border-radius:9px;color:var(--ink)}
.split{display:grid;grid-template-columns:minmax(430px,1fr) 1.25fr;gap:16px;align-items:start}
.pane{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:14px;max-height:84vh;overflow:auto;box-shadow:0 1px 3px rgba(15,23,42,.05)}
.rec{border:1px solid var(--line);border-radius:9px;padding:8px 10px;margin:8px 0;background:#fbfcfe}.rec .vn{font-weight:700;color:var(--acc)}
.btn{cursor:pointer;font-size:12px;background:#eff4ff;border:1px solid #c7d7f5;color:var(--acc);border-radius:7px;padding:3px 9px;margin:2px 4px 0 0;display:inline-block}
.btn:hover{background:var(--acc);color:#fff}.mut{color:var(--mut)}.q{color:#334155;font-style:italic}
#src{font-size:13.5px;color:#243040}#src .blk{padding:3px 6px;border-radius:5px}
#src .hd{color:var(--acc2);font-weight:700;margin:14px 0 4px;border-top:1px solid var(--line);padding-top:8px;scroll-margin-top:60px}
#src .hl{background:var(--hl);outline:2px solid #f59e0b;outline-offset:1px}
.toc{position:sticky;top:0;z-index:5;background:#fff;border:1px solid var(--line);border-radius:9px;padding:8px 10px;margin-bottom:10px;max-height:150px;overflow:auto}
.toc a{display:inline-block;font-size:12px;padding:2px 9px;margin:2px;border-radius:20px;background:#eef2f8;color:#334155;cursor:pointer}
.toc a:hover{background:var(--acc);color:#fff;text-decoration:none}.toc .h3{font-weight:600}
.tblwrap{overflow-x:auto;margin:8px 0;border:1px solid var(--line);border-radius:8px}
.srctbl{margin:0;font-size:12px;min-width:100%}.srctbl th{position:static;background:#eef2f8}.srctbl th,.srctbl td{border:1px solid var(--line);padding:4px 7px;white-space:nowrap}
.links a{margin-right:14px;font-weight:600}.note{background:#fff7ed;border-left:3px solid var(--warn);padding:8px 12px;border-radius:6px;margin:10px 0;color:#7c4a12}
.big{font-size:28px;font-weight:800;color:#0f172a}.flex{display:flex;gap:24px;flex-wrap:wrap}
.gold{background:linear-gradient(180deg,#fffbeb,#fff);border-color:#fde68a}
.figref{display:inline-block;font-size:12px;padding:2px 9px;border-radius:7px;background:#eef2f8;color:var(--acc2)}
.flag{color:var(--warn);cursor:help;font-weight:700}
"""

SORT_JS = """
function filt(id){var q=document.getElementById(id+'-s').value.toLowerCase();
 document.querySelectorAll('#'+id+' tbody tr').forEach(function(r){r.style.display=r.textContent.toLowerCase().includes(q)?'':'none'});}
function srt(id,n){var t=document.getElementById(id),tb=t.tBodies[0],rs=[].slice.call(tb.rows);
 var asc=t.getAttribute('data-sc')!=n+'a';t.setAttribute('data-sc',asc?n+'a':n+'d');
 rs.sort(function(a,b){var x=a.cells[n].getAttribute('data-v')||a.cells[n].textContent,y=b.cells[n].getAttribute('data-v')||b.cells[n].textContent;
  var nx=parseFloat(x),ny=parseFloat(y);if(!isNaN(nx)&&!isNaN(ny)){x=nx;y=ny}return (x>y?1:x<y?-1:0)*(asc?1:-1)});
 rs.forEach(function(r){tb.appendChild(r)})}
"""

JUMP_JS = r"""
function nrm(s){return (s||'').toLowerCase().replace(/\s+/g,' ').trim();}
function _blocks(){return [].slice.call(document.querySelectorAll('#src .blk'));}
// Find the source block that best matches a phrase: exact substring, then a
// sliding word-window (8..3 words) so partial / lightly-reworded text still hits.
function _find(s){s=nrm(s);if(s.length<3)return null;var blks=_blocks();
 var hit=blks.find(function(b){return nrm(b.textContent).indexOf(s)>=0});if(hit)return hit;
 var w=s.split(' ');for(var len=Math.min(8,w.length);len>=3;len--){
  for(var i=0;i+len<=w.length;i++){var win=w.slice(i,i+len).join(' ');
   var h=blks.find(function(b){return nrm(b.textContent).indexOf(win)>=0});if(h)return h;}}
 return null;}
var AA3={ala:'A',arg:'R',asn:'N',asp:'D',cys:'C',gln:'Q',glu:'E',gly:'G',his:'H',ile:'I',leu:'L',lys:'K',met:'M',phe:'F',pro:'P',ser:'S',thr:'T',trp:'W',tyr:'Y',val:'V',ter:'*'};
// p.Arg100Trp -> R100W (so a DB 3-letter notation also matches source 1-letter)
function shortVar(s){return s.replace(/p\.?\s*([A-Za-z]{3})(\d{1,4})([A-Za-z]{3}|Ter|\*)?/g,function(m,a,num,b){
 var x=AA3[a.toLowerCase()];if(!x)return m;var y=b?(AA3[b.toLowerCase()]||(b==='*'?'*':'')):'';return x+num+y;});}
// source_location is often a structured pointer ("Results - Family 2; Table 1;
// Discussion - Variant p.V822L section"), not verbatim text. Try the VARIANT
// notation first (most specific -> lands on the exact row/sentence), in both
// 3- and 1-letter forms, then quotes, then Table/Figure refs, then phrases.
function jump(q){q=q||'';
 document.querySelectorAll('#src .hl').forEach(function(e){e.classList.remove('hl')});
 var vars=(q.match(/p\.[A-Za-z]{3}\d{1,4}(?:[A-Za-z]{3}|Ter|\*)?|[pc]\.[A-Za-z0-9_>*()+-]+|\b[A-Z]\d{1,4}[A-Z*]\b/g)||[]);
 var shorts=[];vars.forEach(function(v){var s=shortVar(v);if(s!==v&&s.length>2)shorts.push(s)});
 var refs=(q.match(/Table\s*\d+|Figure\s*\d+|Fig\.?\s*\d+|Family\s+[A-Za-z0-9.\s]{1,18}/g)||[]);
 var phrases=q.split(/[;|·]| - |, /).map(nrm).filter(function(x){return x.length>3});
 phrases.sort(function(a,b){return b.length-a.length});
 var order=vars.concat(shorts).concat([q]).concat(refs).concat(phrases);
 var hit=null;for(var i=0;i<order.length;i++){hit=_find(order[i]);if(hit)break;}
 if(hit){hit.classList.add('hl');hit.scrollIntoView({behavior:'smooth',block:'center'})}
 else{alert('Could not locate this in the rendered full text — it may be in a supplement/figure, or reworded by the model.');}}
function srcfind(){var q=nrm(document.getElementById('src-s').value);if(q)jump(q);}
function goto(id){var e=document.getElementById(id);if(!e)return;
 document.querySelectorAll('#src .hl').forEach(function(x){x.classList.remove('hl')});
 e.classList.add('hl');e.scrollIntoView({behavior:'smooth',block:'start'});}
"""


def _page(title: str, body: str, extra_js: str = "") -> str:
    return (
        f"<!doctype html><html><head><meta charset='utf-8'><title>{esc(title)}</title>"
        f"<style>{CSS}</style></head><body>{body}<script>{SORT_JS}{extra_js}</script></body></html>"
    )


def _barbar(ok: int, total: int) -> str:
    p = pct(ok, total)
    return f"<div class='bar'><span class='ok' style='width:{p:.0f}%'></span></div>"


import os  # noqa: E402

REPO = Path(__file__).resolve().parents[1]


def _rel(target: Path, start: Path) -> str:
    return os.path.relpath(target, start)


def find_latest_db(gene: str) -> Optional[Path]:
    cands = []
    for root in ("results", "validation_runs"):
        base = REPO / root
        if base.exists():
            cands += list(base.rglob(f"{gene}.db"))
    return max(cands, key=lambda p: p.stat().st_mtime) if cands else None


def score_genes(genes: list[str], db_map: dict[str, Path]) -> dict[str, dict]:
    """Run the recall suite for genes that have a gold standard + a DB, and
    return {gene: {pmids|variant_rows|unique_variants: (matched, gold, recall)}}.
    Reuses scripts/run_recall_suite.py so the numbers match the official scorer."""
    gold_dir = REPO / "gene_variant_fetcher_gold_standard"
    scoreable: dict[str, Path] = {}
    for g in genes:
        gold = gold_dir / "normalized" / f"{g}_recall_input.csv"
        db = db_map.get(g) or find_latest_db(g)
        if gold.exists() and db and db.exists():
            scoreable[g] = db
    if not scoreable:
        return {}
    out = Path(tempfile.mkdtemp(prefix="gvf_dash_score_"))
    cmd = [
        sys.executable,
        str(REPO / "scripts" / "run_recall_suite.py"),
        "--score",
        "--genes",
        ",".join(scoreable.keys()),
        "--gold-dir",
        str(gold_dir),
        "--outdir",
        str(out),
    ]
    for g, db in scoreable.items():
        cmd += ["--db", f"{g}={db}"]
    res: dict[str, dict] = {}
    try:
        subprocess.run(cmd, capture_output=True, text=True, timeout=2400)
        data = json.loads((out / "summary.json").read_text())
        for item in data.get("gene_results", []):
            if item.get("status") != "scored":
                continue
            g = item["gene"]
            summ = item.get("summary", {})
            rc = summ.get("recall", {})
            agg = {
                k: (
                    rc.get(k, {}).get("matched"),
                    rc.get(k, {}).get("gold"),
                    rc.get(k, {}).get("recall"),
                )
                for k in (
                    "pmids",
                    "variant_rows",
                    "unique_variants",
                    "patients",
                    "affected",
                    "unaffected",
                )
            }
            mae = {
                k: summ.get("mae", {}).get(k, {}).get("mae")
                for k in ("carriers", "affected", "unaffected")
            }
            # Per-paper: gold variant count (from the gold CSV) + missing count
            # (gold variants not matched in the DB, from missing_in_sqlite.csv).
            gold_n: Counter = Counter()
            gold_csv = gold_dir / "normalized" / f"{g}_recall_input.csv"
            seen: set = set()
            if gold_csv.exists():
                with gold_csv.open(newline="") as f:
                    for row in csv.DictReader(f):
                        key = (row.get("pmid", ""), row.get("variant", ""))
                        if row.get("pmid") and key not in seen:
                            seen.add(key)
                            gold_n[row["pmid"]] += 1
            miss_n: Counter = Counter()
            miss_csv = out / g / "missing_in_sqlite.csv"
            if miss_csv.exists():
                with miss_csv.open(newline="") as f:
                    for row in csv.DictReader(f):
                        if row.get("pmid"):
                            miss_n[row["pmid"]] += 1
            per_pmid = {
                p: {"gold": gold_n[p], "missing": int(miss_n.get(p, 0))} for p in gold_n
            }
            res[g] = {"agg": agg, "mae": mae, "per_pmid": per_pmid}
    except Exception:  # noqa: BLE001 - scoring is best-effort; dashboard still renders
        return res
    finally:
        shutil.rmtree(out, ignore_errors=True)
    return res


# ---------------------------------------------------------------------------
# Per-paper ADJUDICATION page
# ---------------------------------------------------------------------------
def render_paper_page(
    gene: str, pmid: str, db: dict, corpus_dir: Path, out_dir: Path
) -> str:
    paper = db["papers"].get(pmid, {"title": None, "meta": {}, "bib": {}})
    rec = db["by_pmid"].get(pmid, {"variants": [], "patients": [], "tables": []})
    links = _artifacts_links(corpus_dir, gene, pmid)
    bib = paper.get("bib") or {}
    doi = bib.get("doi") or links["doi"]
    pmcid = bib.get("pmc_id") or links["pmcid"]
    paper_root = corpus_dir / gene / pmid
    ft = paper_root / f"{pmid}_FULL_CONTEXT.md"

    link_html = [f"<a href='{pubmed_url(pmid)}' target='_blank'>PubMed ↗</a>"]
    if doi:
        link_html.append(f"<a href='{doi_url(doi)}' target='_blank'>DOI ↗</a>")
    if pmcid:
        link_html.append(f"<a href='{pmc_url(pmcid)}' target='_blank'>PMC ↗</a>")
    link_html.append("<a href='../variants.html'>variants</a>")
    link_html.append(f"<a href='../index.html'>← {esc(gene)}</a>")
    link_html.append("<a href='../../index.html'>overview</a>")

    cite = " · ".join(
        x
        for x in [
            esc(bib.get("first_author")) + " et al." if bib.get("first_author") else "",
            f"<i>{esc(bib.get('journal'))}</i>" if bib.get("journal") else "",
            esc(bib.get("publication_date")) if bib.get("publication_date") else "",
        ]
        if x
    )

    meta = paper["meta"]
    meta_bits = []
    for k in (
        "model_used",
        "source_type",
        "source_file",
        "extraction_confidence",
        "total_variants_found",
    ):
        if meta.get(k) not in (None, ""):
            meta_bits.append(f"<span class='tag'>{esc(k)}: {esc(meta[k])}</span>")
    pen_tot = db.get("pen_by_pmid", {}).get(pmid)
    if pen_tot and (
        pen_tot["carriers"] or pen_tot["affected"] or pen_tot["unaffected"]
    ):
        meta_bits.append(
            f"<span class='tag ok'>cohort: {pen_tot['carriers']} carriers · "
            f"{pen_tot['affected']} aff · {pen_tot['unaffected']} unaff</span>"
        )

    # LEFT pane: extracted records with jump-to-source actions
    left = ["<h3>Extracted variants ({}):</h3>".format(len(rec["variants"]))]
    for v in rec["variants"]:
        # Primary, most-precise jump: the variant notation itself (jump() also
        # tries its 1-letter form), so it lands on the exact row/sentence.
        actions = [
            f"<span class='btn' onclick=\"jump({jq(v['variant'])})\">🎯 find “{esc(v['variant'])}”</span>"
        ]
        if v["location"]:
            actions.append(
                f"<span class='btn' onclick=\"jump({jq(v['location'])})\">📍 {esc(v['location'])}</span>"
            )
        for q in v["quotes"][:2]:
            actions.append(
                f"<span class='btn' onclick=\"jump({jq(q)})\">❝ quote</span>"
            )
        if not v["quotes"] and v["notes"]:
            actions.append(
                f"<span class='btn' onclick=\"jump({jq(v['notes'])})\">🔎 notes</span>"
            )
        notes = f"<div class='mut'>{esc(v['notes'])}</div>" if v["notes"] else ""
        cprov = (
            f"<div class='mut'>⚖️ count basis: {esc(v['count_prov'])}</div>"
            if v.get("count_prov")
            else ""
        )
        quotes = "".join(f"<div class='q'>❝ {esc(q)}</div>" for q in v["quotes"][:2])
        pv = db.get("pen_pv", {}).get((pmid, v["variant"]))
        cohort = (
            f"<div class='mut'>👥 cohort: {pv['carriers']} carriers · "
            f"{pv['affected']} affected · {pv['unaffected']} unaffected</div>"
            if pv and (pv["carriers"] or pv["affected"] or pv["unaffected"])
            else ""
        )
        left.append(
            f"<div class='rec'><span class='vn'>{esc(v['variant'])}</span>"
            f" <span class='mut'>{esc(v['sig'] or '')}</span><div>{''.join(actions)}</div>"
            f"{cohort}{quotes}{notes}{cprov}</div>"
        )
    if rec["patients"]:
        left.append(f"<h3>Patient records ({len(rec['patients'])}):</h3>")
        for p in rec["patients"]:
            st = p["status"]
            cls = "ok" if st == "affected" else "warn" if st == "unaffected" else ""
            chars = " · ".join(
                x
                for x in [
                    f"{esc(p['variant'])}" if p["variant"] else "",
                    f"age {esc(p['age'])}" if p["age"] not in (None, "") else "",
                    f"sex {esc(p['sex'])}" if p["sex"] else "",
                    f"🌍 {esc(p['ethnicity'])}" if p.get("ethnicity") else "",
                    f"📍 {esc(p['origin'])}" if p.get("origin") else "",
                    esc(p["phenotype"][:80]) if p["phenotype"] else "",
                ]
                if x
            )
            ev = (
                (
                    f"<span class='btn' onclick=\"jump({jq(p['evidence'])})\">🔎 evidence sentence</span>"
                    f"<div class='q'>❝ {esc(p['evidence'])}</div>"
                )
                if p["evidence"]
                else ""
            )
            left.append(
                f"<div class='rec'><span class='tag {cls}'>{esc(st or '?')}</span> {chars}{ev}</div>"
            )
    if not rec["variants"] and not rec["patients"]:
        left.append(
            "<p class='mut'>No extracted variant/patient records in the DB for this PMID "
            "(it may be abstract-only, a stub, or filtered out).</p>"
        )

    # RIGHT pane: the exact on-disk source + a section table-of-contents
    toc = ""
    bp_warn = ""
    if ft.exists():
        raw = ft.read_text(encoding="utf-8", errors="replace")
        src_html = md_to_html(raw)
        bp_warn = boilerplate_warning(raw)
        src_note = f"<div class='mut'>Exact on-disk source the LLM parsed: <code>{esc(ft.name)}</code></div>"
        heads = re.findall(r'<h(\d) class="blk hd" id="(sec\d+)">(.*?)</h\1>', src_html)
        if heads:
            links = "".join(
                f"<a class='h{lvl}' onclick=\"goto('{hid}')\">{txt[:48]}</a>"
                for lvl, hid, txt in heads
            )
            toc = f"<div class='toc'><b class='mut'>Jump to section:</b> {links}</div>"
    else:
        src_html = "<p class='note'>No full-text source on disk for this PMID (abstract-only or not fetched).</p>"
        src_note = ""

    # Reference figures/supplements RELATIVELY into the corpus (no duplication —
    # copying for thousands of papers bloated the bundle to GBs). These resolve
    # locally and when corpus/ is served as the web root (the dashboard lives at
    # corpus/dashboard/, so ../<gene>/<pmid>/ -> corpus/<gene>/<pmid>/). Zip/serve
    # the whole corpus/ to share it.
    figs = ""
    figdir = paper_root / f"{pmid}_figures"
    if figdir.is_dir():
        imgs = [
            f
            for f in sorted(figdir.iterdir())
            if f.suffix.lower() in {".png", ".jpg", ".jpeg", ".gif"}
        ][:30]
        if imgs:
            items = "".join(
                f"<a href='{esc(_rel(im, out_dir))}' target='_blank'>"
                f"<img src='{esc(_rel(im, out_dir))}' loading='lazy' "
                f"style='max-width:200px;margin:5px;border:1px solid var(--line);border-radius:6px'></a>"
                for im in imgs
            )
            figs = f"<h3>Figures ({len(imgs)})</h3><div>{items}</div>"
    sups = ""
    supdir = paper_root / f"{pmid}_supplements"
    if supdir.is_dir():
        files = [f for f in sorted(supdir.rglob("*")) if f.is_file()][:40]
        if files:
            items = "".join(
                f"<div><a href='{esc(_rel(f, out_dir))}' target='_blank'>📎 {esc(f.name)}</a></div>"
                for f in files
            )
            sups = f"<h3>Supplements ({len(files)})</h3>{items}"

    title = paper["title"] or f"PMID {pmid}"
    cite_line = f"<div class='sub'>{cite}</div>" if cite else ""
    body = (
        f"<header><h1>{esc(title)}</h1>{cite_line}"
        f"<div class='links sub'>PMID {esc(pmid)} · {' '.join(link_html)}</div>"
        f"<div>{''.join(meta_bits)}</div></header>"
        f"<div class='wrap'><div class='split'>"
        f"<div class='pane'>{''.join(left)}</div>"
        f"<div class='pane'><input class='search' id='src-s' placeholder='search the full text…' "
        f"onkeydown=\"if(event.key==='Enter')srcfind()\">{src_note}{bp_warn}{toc}"
        f"<div id='src'>{src_html}</div>{figs}{sups}</div>"
        f"</div></div>"
    )
    return _page(f"{gene} · PMID {pmid}", body, JUMP_JS)


# ---------------------------------------------------------------------------
# Per-gene page + index
# ---------------------------------------------------------------------------
def gene_stats(gene: str, corpus_rows: dict, db: Optional[dict]) -> dict:
    usable = sum(1 for r in corpus_rows.values() if r.get("full_text_status") == "ok")
    with_figs = sum(1 for r in corpus_rows.values() if int(r.get("n_figures") or 0) > 0)
    with_sup = sum(
        1 for r in corpus_rows.values() if int(r.get("n_supplement_files") or 0) > 0
    )
    s = {
        "in_corpus": len(corpus_rows),
        "usable": usable,
        "stub": len(corpus_rows) - usable,
        "with_figs": with_figs,
        "with_sup": with_sup,
        "extracted": 0,
        "with_variants": 0,
        "variant_rows": 0,
        "unique": 0,
    }
    if db:
        s["extracted"] = len(db["papers"])
        s["with_variants"] = sum(1 for v in db["by_pmid"].values() if v["variants"])
        s["variant_rows"] = db["prov"]["total"]
        s["unique"] = len(
            {x["variant"] for p in db["by_pmid"].values() for x in p["variants"]}
        )
    return s


def _gold_card(recall: Optional[dict], mae: Optional[dict] = None) -> str:
    """Gold-standard recall card: all 6 recall metrics + carrier/affected/
    unaffected MAE (only genes with a gold standard get one)."""
    if not recall:
        return ""
    import math

    def line(label: str, key: str) -> str:
        m, g, r = recall.get(key, (None, None, None))
        if g is None:
            return ""
        cls = "ok" if (r or 0) >= 0.9 else ""
        return (
            f"<div class='kv'><span>{label}</span>"
            f"<b class='{cls}'>{m}/{g} ({(r or 0) * 100:.1f}%)</b></div>"
        )

    m, g, r = recall.get("unique_variants", (None, None, None))
    bar = (
        f"<div class='bar'><span class='ok' style='width:{(r or 0) * 100:.0f}%'></span></div>"
        if g
        else ""
    )
    gap = max(0, math.ceil((g or 0) * 0.9) - (m or 0)) if g else 0
    mae = mae or {}
    mae_html = ""
    if any(v is not None for v in mae.values()):
        mae_html = (
            "<div class='kv' style='border-top:1px solid var(--line);margin-top:6px;padding-top:6px'>"
            "<span>count MAE (carriers/aff/unaff)</span><b>"
            + " / ".join(
                f"{mae[k]:.2f}" if mae.get(k) is not None else "—"
                for k in ("carriers", "affected", "unaffected")
            )
            + "</b></div>"
        )
    return (
        "<div class='card gold'><h2>⭐ Gold-standard recall</h2>"
        f"{line('PMIDs', 'pmids')}{line('variant rows', 'variant_rows')}"
        f"{line('unique variants', 'unique_variants')}{bar}"
        f"{line('patients/carriers', 'patients')}{line('affected', 'affected')}"
        f"{line('unaffected', 'unaffected')}"
        f"<div class='kv'><span>unique-variant gap to 90%</span><b>+{gap}</b></div>"
        f"{mae_html}"
        "<div class='mut' style='font-size:11px;margin-top:4px'>recall &amp; MAE via run_recall_suite vs the curated gold set</div></div>"
    )


def gold_coverage_stats(gold_pmids: set, corpus_rows: dict) -> dict:
    """Of the gold-standard PMIDs, how many have usable full text / figures /
    supplements / all three on disk (per corpus INDEX)."""
    n = len(gold_pmids)
    ft = figs = sup = allthree = 0
    for p in gold_pmids:
        r = corpus_rows.get(p)
        has_ft = bool(r) and r.get("full_text_status") == "ok"
        has_fig = bool(r) and int(r.get("n_figures") or 0) > 0
        has_sup = bool(r) and int(r.get("n_supplement_files") or 0) > 0
        ft += has_ft
        figs += has_fig
        sup += has_sup
        allthree += has_ft and has_fig and has_sup
    return {"n": n, "ft": ft, "figs": figs, "sup": sup, "all3": allthree}


def _gold_source_card(gold_pmids: set, corpus_rows: dict) -> str:
    if not gold_pmids:
        return ""
    c = gold_coverage_stats(gold_pmids, corpus_rows)
    n = c["n"]
    ft, figs, sup, allthree = c["ft"], c["figs"], c["sup"], c["all3"]

    def row(label: str, k: int) -> str:
        return f"<div class='kv'><span>{label}</span><b>{k}/{n} ({pct(k, n):.0f}%)</b></div>"

    return (
        "<div class='card'><h2>Gold-PMID source coverage</h2>"
        f"<div class='kv'><span>gold-standard PMIDs</span><b>{n}</b></div>"
        f"{_barbar(ft, n)}"
        f"{row('usable full text', ft)}{row('figures on disk', figs)}"
        f"{row('supplements on disk', sup)}{row('full text + figures + supplements', allthree)}"
        "<div class='mut' style='font-size:11px;margin-top:4px'>over the curated gold-standard PMID set</div></div>"
    )


def render_gene_page(
    gene: str,
    corpus_rows: dict,
    db: Optional[dict],
    out_dir: Path,
    corpus_dir: Path,
    max_papers: int,
    recall: Optional[dict] = None,
    per_pmid: Optional[dict] = None,
    mae: Optional[dict] = None,
) -> tuple[str, dict, list[str]]:
    s = gene_stats(gene, corpus_rows, db)
    per_pmid = per_pmid or {}
    pen_by_pmid = (db or {}).get("pen_by_pmid", {})
    prov = (db or {}).get("prov", Counter())
    method = (db or {}).get("source_method", Counter())

    # provenance-completeness audit (the "test experiment" surface)
    audit = []
    if prov.get("total"):
        audit = [
            ("source_location populated", prov.get("with_location", 0), prov["total"]),
            ("additional_notes populated", prov.get("with_notes", 0), prov["total"]),
            ("key_quotes populated", prov.get("with_quotes", 0), prov["total"]),
        ]

    # which papers still need work (no usable full text, or no variants)
    stub_pmids = [
        p for p, r in corpus_rows.items() if r.get("full_text_status") != "ok"
    ]

    # paper table
    rows = []
    paper_pmids = sorted(
        set(corpus_rows) | set((db or {}).get("papers", {})),
        key=lambda x: int(x) if x.isdigit() else 0,
    )
    written: list[str] = []

    def _rec_count(pm: str) -> int:
        e = (db or {}).get("by_pmid", {}).get(pm)
        return (len(e["variants"]) + len(e["patients"])) if e else 0

    # Papers we can render a page for (have DB rows or on-disk source). EVERY
    # paper with extracted records gets a page (uncapped) so a paper with any
    # variants/patients is never a mere PubMed link-out; the cap applies only to
    # source-only papers (nothing extracted to adjudicate).
    pageable = [
        p
        for p in paper_pmids
        if db is not None
        and (p in (db or {}).get("papers", {}) or (corpus_dir / gene / p).is_dir())
    ]
    with_rec = [p for p in pageable if _rec_count(p) > 0]
    source_only = [p for p in pageable if _rec_count(p) == 0]
    source_only.sort(
        key=lambda p: int(corpus_rows.get(p, {}).get("fulltext_bytes") or 0),
        reverse=True,
    )
    cap_set = set(with_rec) | set(
        source_only if max_papers <= 0 else source_only[:max_papers]
    )

    has_gold = bool(per_pmid)
    for pmid in paper_pmids:
        cr = corpus_rows.get(pmid, {})
        status = cr.get("full_text_status", "—")
        by = (db or {}).get("by_pmid", {}).get(pmid, {})
        nvar = len(by.get("variants", [])) if db else 0
        npat = len(by.get("patients", [])) if db else 0
        title = (db or {}).get("papers", {}).get(pmid, {}).get("title") or ""
        badge = "ok" if status == "ok" else "bad" if status == "stub" else "warn"
        has_page = pmid in cap_set
        pmid_cell = (
            f"<a href='papers/{esc(pmid)}.html'>{esc(pmid)}</a>"
            if has_page
            else f"<a href='{pubmed_url(pmid)}' target='_blank'>{esc(pmid)}↗</a>"
        )
        if has_page:
            written.append(pmid)
        gold_cells = ""
        if has_gold:
            pp = per_pmid.get(pmid)
            if pp and pp["gold"]:
                g, miss = pp["gold"], pp["missing"]
                matched = max(0, g - miss)
                mcls = "ok" if miss == 0 else "bad"
                gold_cells = (
                    f"<td data-v='{g}'>{g}</td>"
                    f"<td data-v='{matched}'>{matched}</td>"
                    f"<td data-v='{miss}'><span class='tag {mcls}'>{'✓' if miss == 0 else '−' + str(miss)}</span></td>"
                )
            else:
                gold_cells = (
                    "<td class='mut'>—</td><td class='mut'>—</td><td class='mut'>—</td>"
                )
        pen = pen_by_pmid.get(pmid, {})
        carriers = pen.get("carriers", 0)
        aff, unaff = pen.get("affected", 0), pen.get("unaffected", 0)
        rows.append(
            f"<tr><td data-v='{esc(pmid)}'>{pmid_cell}</td><td>{esc(title[:80])}</td>"
            f"<td><span class='tag {badge}'>{esc(status)}</span>"
            + (
                " <span class='flag' title='Usable full text but 0 variants extracted — possible boilerplate scrape or extraction miss; check the source'>⚠</span>"
                if status == "ok" and nvar == 0 and db is not None
                else ""
            )
            + "</td>"
            f"<td data-v='{nvar}'>{nvar}</td>{gold_cells}"
            f"<td data-v='{carriers}'>{carriers}</td>"
            f"<td data-v='{aff}'>{aff}</td><td data-v='{unaff}'>{unaff}</td>"
            f"<td data-v='{npat}' class='mut'>{npat}</td></tr>"
        )

    # Build the header (gold genes get gold / matched / missing-vs-gold columns).
    # Carriers/Aff/Unaff are cohort counts from penetrance_data (the real patient
    # numbers); "records" is the sparser per-person individual_records count.
    base_cols = ["PMID", "Title", "Source", "Variants"]
    gold_cols = ["Gold", "Matched", "Δ vs gold"] if has_gold else []
    all_cols = base_cols + gold_cols + ["Carriers", "Aff", "Unaff", "records"]
    thead_html = "".join(
        f"<th onclick=\"srt('pt',{i})\">{esc(c)}</th>" for i, c in enumerate(all_cols)
    )

    audit_html = "".join(
        f"<div class='kv'><span>{esc(name)}</span><b>{n}/{tot} ({pct(n, tot):.0f}%)</b></div>"
        for name, n, tot in audit
    )
    method_html = "".join(
        f"<span class='tag'>{esc(k)}: {v}</span>" for k, v in method.most_common()
    )
    gap_note = (
        ""
        if db
        else "<div class='note'>No scored DB found for this gene — showing source coverage only. "
        "Pass <code>--db {}=path/to/{}.db</code> for variant provenance.</div>".format(
            gene, gene
        )
    )

    body = (
        f"<header><h1>{esc(gene)}</h1><div class='sub'>"
        f"<a href='../index.html'>← overview</a> · <a href='variants.html'>variants ({s['unique']})</a>"
        f"</div></header>"
        f"<div class='wrap'>{gap_note}"
        f"<div class='flex'>"
        f"<div class='card'><h2>Source coverage</h2>"
        f"<div class='kv'><span>papers in corpus</span><b>{s['in_corpus']}</b></div>"
        f"<div class='kv'><span>usable full text</span><b>{s['usable']} ({pct(s['usable'], s['in_corpus']):.0f}%)</b></div>"
        f"{_barbar(s['usable'], s['in_corpus'])}"
        f"<div class='kv'><span>stub / paywalled (to fetch)</span><b>{s['stub']}</b></div>"
        f"<div class='kv'><span>with figures</span><b>{s['with_figs']}</b></div>"
        f"<div class='kv'><span>with supplements</span><b>{s['with_sup']}</b></div></div>"
        f"<div class='card'><h2>Extraction funnel</h2>"
        f"<div class='kv'><span>extracted (in DB)</span><b>{s['extracted']}</b></div>"
        f"<div class='kv'><span>papers with variants</span><b>{s['with_variants']}</b></div>"
        f"<div class='kv'><span>variant rows</span><b>{s['variant_rows']}</b></div>"
        f"<div class='kv'><span>unique variants</span><b>{s['unique']}</b></div></div>"
        f"{_gold_card(recall, mae)}"
        f"{_gold_source_card(set(per_pmid), corpus_rows)}"
        f"<div class='card'><h2>Provenance completeness</h2>{audit_html or '<span class=mut>no DB</span>'}"
        f"<div style='margin-top:8px'>{method_html}</div></div>"
        f"<div class='card'><h2>What's left</h2>"
        f"<div class='kv'><span>papers without usable full text</span><b>{len(stub_pmids)}</b></div>"
        f"<div class='mut' style='font-size:12px;margin-top:6px'>re-run <code>gvf gvf-run {esc(gene)} …</code> "
        f"(corpus cache skips what's done; a new publisher key re-fetches stubs)</div></div>"
        f"</div>"
        f"<h2 style='margin-top:20px'>Papers</h2>"
        f"<div class='mut' style='font-size:12px;margin-bottom:6px'>"
        f"Variants = extracted; Gold/Matched/Δ vs the gold standard (where available); "
        f"Carriers/Aff/Unaff = cohort counts from penetrance_data; records = per-person rows. "
        f"<span class='flag'>⚠</span> = usable full text on disk but 0 variants extracted "
        f"(possible boilerplate/page-shell scrape or an extraction miss — open it to check).</div>"
        f"<input class='search' id='pt-s' placeholder='filter papers…' onkeyup=\"filt('pt')\">"
        f"<table id='pt' data-sc=''><thead><tr>{thead_html}</tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table>"
        f"</div>"
    )
    return _page(f"{gene} — GVF dashboard", body), s, written


def render_variants_page(gene: str, db: dict, page_pmids: set) -> str:
    """Variant-centric view (website-style): one row per unique variant,
    aggregating papers + affected/unaffected patient counts, with clickthrough
    to each paper's adjudication page (or PubMed)."""
    agg = db.get("variants_agg", {})
    rows = []
    for vname, a in sorted(
        agg.items(), key=lambda kv: len(kv[1]["pmids"]), reverse=True
    ):
        pmids = sorted(a["pmids"], key=lambda x: int(x) if x.isdigit() else 0)
        links = [
            (
                f"<a href='papers/{esc(pm)}.html'>{esc(pm)}</a>"
                if pm in page_pmids
                else f"<a href='{pubmed_url(pm)}' target='_blank'>{esc(pm)}↗</a>"
            )
            for pm in pmids[:12]
        ]
        more = (
            f" <span class='mut'>+{len(pmids) - 12}</span>" if len(pmids) > 12 else ""
        )
        rows.append(
            f"<tr><td>{esc(vname)}</td><td>{esc(a['sig'])}</td>"
            f"<td data-v='{len(pmids)}'>{len(pmids)}</td>"
            f"<td data-v='{a.get('carriers', 0)}'>{a.get('carriers', 0)}</td>"
            f"<td data-v='{a['affected']}'>{a['affected']}</td>"
            f"<td data-v='{a['unaffected']}'>{a['unaffected']}</td>"
            f"<td>{' '.join(links)}{more}</td></tr>"
        )
    body = (
        f"<header><h1>{esc(gene)} — variants</h1>"
        f"<div class='sub'><a href='index.html'>← {esc(gene)}</a> · "
        f"<a href='../index.html'>overview</a> · {len(agg)} unique variants · "
        f"carriers/affected/unaffected are cohort counts (penetrance_data)</div></header>"
        f"<div class='wrap'>"
        f"<input class='search' id='vt-s' placeholder='filter variants…' onkeyup=\"filt('vt')\">"
        f"<table id='vt' data-sc=''><thead><tr>"
        f"<th onclick=\"srt('vt',0)\">Variant</th><th onclick=\"srt('vt',1)\">Clinical sig</th>"
        f"<th onclick=\"srt('vt',2)\">Papers</th><th onclick=\"srt('vt',3)\">Carriers</th>"
        f"<th onclick=\"srt('vt',4)\">Affected</th>"
        f"<th onclick=\"srt('vt',5)\">Unaffected</th><th>Papers (→ adjudicate / PubMed)</th>"
        f"</tr></thead><tbody>{''.join(rows)}</tbody></table></div>"
    )
    return _page(f"{gene} variants — GVF", body)


def render_index(summaries: dict[str, dict], generated: str) -> str:
    cards = []
    for gene, s in sorted(summaries.items()):
        rc = s.get("recall") or {}
        gold = ""
        if rc:
            uv = rc.get("unique_variants", (None, None, None))
            pm = rc.get("pmids", (None, None, None))
            gold = (
                "<div class='kv' style='border-top:1px solid var(--line);margin-top:6px;padding-top:6px'>"
                f"<span>⭐ gold recall (uniqV / PMID)</span>"
                f"<b>{(uv[2] or 0) * 100:.0f}% / {(pm[2] or 0) * 100:.0f}%</b></div>"
            )
        gc = s.get("gold_cov")
        if gc and gc["n"]:
            gold += (
                "<div class='kv'><span>gold source (txt/fig/supp/all)</span><b>"
                f"{pct(gc['ft'], gc['n']):.0f}% / {pct(gc['figs'], gc['n']):.0f}% / "
                f"{pct(gc['sup'], gc['n']):.0f}% / {pct(gc['all3'], gc['n']):.0f}%</b></div>"
            )
        cards.append(
            f"<div class='card'><h2><a href='{esc(gene)}/index.html'>{esc(gene)}</a></h2>"
            f"<div class='kv'><span>papers in corpus</span><b>{s['in_corpus']}</b></div>"
            f"{_barbar(s['usable'], s['in_corpus'])}"
            f"<div class='kv'><span>usable full text</span><b>{s['usable']} ({pct(s['usable'], s['in_corpus']):.0f}%)</b></div>"
            f"<div class='kv'><span>extracted</span><b>{s['extracted']}</b></div>"
            f"<div class='kv'><span>papers w/ variants</span><b>{s['with_variants']}</b></div>"
            f"<div class='kv'><span>unique variants</span><b>{s['unique']}</b></div>{gold}"
            f"<div style='margin-top:8px'><a href='{esc(gene)}/index.html'>status</a> · "
            f"<a href='{esc(gene)}/variants.html'>variants</a></div></div>"
        )
    tot = {
        k: sum(s[k] for s in summaries.values())
        for k in ("in_corpus", "usable", "extracted", "variant_rows")
    }
    body = (
        f"<header><h1>GeneVariantFetcher — status & provenance dashboard</h1>"
        f"<div class='sub'>generated {esc(generated)} · single source of truth: docs/RECALL_STATUS.md · "
        f"source corpus: corpus/INDEX.csv</div></header>"
        f"<div class='wrap'><div class='flex' style='margin-bottom:16px'>"
        f"<div><div class='big'>{tot['in_corpus']}</div><div class='mut'>papers in corpus</div></div>"
        f"<div><div class='big'>{tot['usable']}</div><div class='mut'>usable full text</div></div>"
        f"<div><div class='big'>{tot['extracted']}</div><div class='mut'>extracted</div></div>"
        f"<div><div class='big'>{tot['variant_rows']}</div><div class='mut'>variant rows</div></div></div>"
        f"<div class='grid'>{''.join(cards)}</div></div>"
    )
    return _page("GVF dashboard", body)


def generate_dashboard(
    out_dir: Path,
    corpus_dir: Path,
    db_map: dict[str, Path],
    genes: Optional[list[str]],
    max_papers: int,
    generated: str,
    score: bool = True,
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    corpus_idx = load_corpus_index(corpus_dir)
    target_genes = genes or sorted(corpus_idx.keys())
    score_map = score_genes(target_genes, db_map) if score else {}
    summaries: dict[str, dict] = {}
    stats = {"genes": 0, "paper_pages": 0, "scored": len(score_map)}
    for gene in target_genes:
        corpus_rows = corpus_idx.get(gene, {})
        if not corpus_rows and gene not in db_map:
            continue
        db_path = db_map.get(gene) or find_latest_db(gene)
        db = load_db(db_path) if db_path and db_path.exists() else None
        sm = score_map.get(gene) or {}
        recall, per_pmid = sm.get("agg"), sm.get("per_pmid")
        # Per-gene subdir: <out>/<GENE>/{index.html, variants.html, papers/<PMID>.html}
        gene_dir = out_dir / gene
        gene_dir.mkdir(parents=True, exist_ok=True)
        page, s, paper_pmids = render_gene_page(
            gene,
            corpus_rows,
            db,
            gene_dir,
            corpus_dir,
            max_papers,
            recall=recall,
            per_pmid=per_pmid,
            mae=sm.get("mae"),
        )
        (gene_dir / "index.html").write_text(page, encoding="utf-8")
        s["recall"] = recall
        s["gold_cov"] = (
            gold_coverage_stats(set(per_pmid), corpus_rows) if per_pmid else None
        )
        summaries[gene] = s
        stats["genes"] += 1
        if db:
            papers_dir = gene_dir / "papers"
            papers_dir.mkdir(parents=True, exist_ok=True)
            for pmid in paper_pmids:
                (papers_dir / f"{pmid}.html").write_text(
                    render_paper_page(gene, pmid, db, corpus_dir, papers_dir),
                    encoding="utf-8",
                )
                stats["paper_pages"] += 1
            (gene_dir / "variants.html").write_text(
                render_variants_page(gene, db, set(paper_pmids)), encoding="utf-8"
            )
    (out_dir / "index.html").write_text(
        render_index(summaries, generated), encoding="utf-8"
    )
    stats["summaries"] = summaries
    return stats
