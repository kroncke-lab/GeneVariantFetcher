"""Static-HTML status / coverage / missingness / provenance dashboard.

`gvf dashboard` reads the consolidated corpus (`corpus/INDEX.csv` + the per-paper
source under `corpus/<GENE>/<PMID>/`) and the scored SQLite DB(s), and writes a
self-contained, offline, zip-transferable dashboard under `--out`:

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
import sqlite3
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


def pct(num: int, den: int) -> float:
    return (100.0 * num / den) if den else 0.0


# ---------------------------------------------------------------------------
# Minimal, dependency-free Markdown -> HTML (faithful: every block is taggable
# so the adjudication JS can jump to and highlight the source sentence/row).
# ---------------------------------------------------------------------------
def md_to_html(text: str) -> str:
    lines = text.replace("\r\n", "\n").split("\n")
    out: list[str] = []
    i, n = 0, len(lines)
    para: list[str] = []

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
            out.append(f'<h{level} class="blk hd">{esc(txt)}</h{level}>')
            i += 1
            continue
        # pipe table: header row + separator row of ---|---
        if (
            "|" in line
            and i + 1 < n
            and re.match(r"^\s*\|?[\s:|-]+\|?\s*$", lines[i + 1])
            and "-" in lines[i + 1]
        ):
            flush_para()
            header = [c.strip() for c in stripped.strip("|").split("|")]
            i += 2
            rows = []
            while i < n and "|" in lines[i] and lines[i].strip():
                rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
                i += 1
            thead = "".join(f"<th>{esc(c)}</th>" for c in header)
            body = "".join(
                "<tr class='blk'>" + "".join(f"<td>{esc(c)}</td>" for c in r) + "</tr>"
                for r in rows
            )
            out.append(
                f'<table class="srctbl"><thead><tr>{thead}</tr></thead><tbody>{body}</tbody></table>'
            )
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
            lambda: {"sig": "", "pmids": set(), "affected": 0, "unaffected": 0}
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
                status = rd.get("affected_status") or ""
                if vname and vname in data["variants_agg"]:
                    if status == "affected":
                        data["variants_agg"][vname]["affected"] += 1
                    elif status == "unaffected":
                        data["variants_agg"][vname]["unaffected"] += 1
                data["by_pmid"][str(rd["pmid"])]["patients"].append(
                    {
                        "variant": vname,
                        "id": rd.get("individual_id") or "",
                        "age": rd.get("age_at_evaluation")
                        if rd.get("age_at_evaluation") is not None
                        else rd.get("age_at_onset"),
                        "sex": rd.get("sex") or "",
                        "status": status,
                        "phenotype": rd.get("phenotype_details") or "",
                        "ethnicity": rd.get("ethnicity") or "",
                        "origin": rd.get("geographic_origin") or "",
                        "evidence": rd.get("evidence_sentence") or "",
                    }
                )
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
:root{--bg:#0f1419;--card:#1a2230;--ink:#e6edf3;--mut:#9aa7b4;--ok:#3fb950;--warn:#d29922;--bad:#f85149;--acc:#58a6ff;--line:#2d3748}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.5 -apple-system,Segoe UI,Roboto,sans-serif}
a{color:var(--acc);text-decoration:none}a:hover{text-decoration:underline}
header{padding:18px 24px;border-bottom:1px solid var(--line)}h1{margin:0;font-size:20px}.sub{color:var(--mut);font-size:13px}
.wrap{padding:20px 24px;max-width:1400px;margin:0 auto}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(280px,1fr));gap:14px}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:16px}
.card h2{margin:0 0 8px;font-size:17px}.kv{display:flex;justify-content:space-between;color:var(--mut);font-size:13px;padding:2px 0}.kv b{color:var(--ink)}
.bar{height:8px;background:#222c3a;border-radius:6px;overflow:hidden;margin:6px 0}.bar>span{display:block;height:100%}
.bar .ok{background:var(--ok)}.bar .warn{background:var(--warn)}
.tag{display:inline-block;font-size:11px;padding:2px 7px;border-radius:10px;border:1px solid var(--line);color:var(--mut);margin:2px 3px 0 0}
.tag.ok{color:var(--ok);border-color:var(--ok)}.tag.bad{color:var(--bad);border-color:var(--bad)}.tag.warn{color:var(--warn);border-color:var(--warn)}
table{width:100%;border-collapse:collapse;font-size:13px}th,td{text-align:left;padding:6px 8px;border-bottom:1px solid var(--line);vertical-align:top}
th{position:sticky;top:0;background:var(--card);cursor:pointer;user-select:none}tr:hover td{background:#161d29}
.search{width:100%;padding:8px 10px;margin:8px 0;background:#0b0f14;border:1px solid var(--line);border-radius:8px;color:var(--ink)}
.split{display:grid;grid-template-columns:minmax(420px,1fr) 1.2fr;gap:16px;align-items:start}
.pane{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px;max-height:82vh;overflow:auto}
.rec{border:1px solid var(--line);border-radius:8px;padding:8px 10px;margin:8px 0}.rec .vn{font-weight:600;color:var(--acc)}
.btn{cursor:pointer;font-size:12px;background:#222c3a;border:1px solid var(--line);color:var(--ink);border-radius:6px;padding:2px 8px;margin:2px 4px 0 0}
.btn:hover{border-color:var(--acc)}.mut{color:var(--mut)}.q{color:var(--ink);font-style:italic}
#src .blk{padding:2px 4px;border-radius:4px}#src .hd{color:var(--acc);margin:10px 0 4px;border-top:1px solid var(--line);padding-top:8px}
#src .hl{background:rgba(210,153,34,.35);outline:2px solid var(--warn)}
.srctbl{margin:8px 0;font-size:12px}.srctbl th,.srctbl td{border:1px solid var(--line);padding:3px 6px}
.links a{margin-right:12px}.note{background:#1d2533;border-left:3px solid var(--warn);padding:8px 12px;border-radius:6px;margin:10px 0;color:var(--mut)}
.big{font-size:26px;font-weight:700}.flex{display:flex;gap:24px;flex-wrap:wrap}
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

JUMP_JS = """
function nrm(s){return (s||'').toLowerCase().replace(/\\s+/g,' ').trim();}
function jump(q){q=nrm(q);if(!q)return;
 document.querySelectorAll('#src .hl').forEach(function(e){e.classList.remove('hl')});
 var blks=[].slice.call(document.querySelectorAll('#src .blk'));
 var hit=blks.find(function(b){var t=nrm(b.textContent);return t.includes(q)||(q.length>40&&t.length>20&&q.includes(t))});
 if(!hit){var key=q.split(' ').slice(0,6).join(' ');hit=blks.find(function(b){return nrm(b.textContent).includes(key)})}
 if(hit){hit.classList.add('hl');hit.scrollIntoView({behavior:'smooth',block:'center'})}
 else{alert('Could not locate this text in the rendered full text. It may live in a supplement/figure, or the model paraphrased it.');}}
function srcfind(){var q=nrm(document.getElementById('src-s').value);if(q)jump(q);}
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
    link_html.append(f"<a href='{esc(gene)}_variants.html'>variants</a>")
    link_html.append(f"<a href='{esc(gene)}.html'>← {esc(gene)}</a>")

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

    # LEFT pane: extracted records with jump-to-source actions
    left = ["<h3>Extracted variants ({}):</h3>".format(len(rec["variants"]))]
    for v in rec["variants"]:
        actions = []
        if v["location"]:
            actions.append(
                f"<span class='btn' onclick=\"jump({json.dumps(v['location'])})\">📍 {esc(v['location'])}</span>"
            )
        for q in v["quotes"][:2]:
            actions.append(
                f"<span class='btn' onclick=\"jump({json.dumps(q)})\">❝ quote</span>"
            )
        if not v["quotes"] and v["notes"]:
            actions.append(
                f"<span class='btn' onclick=\"jump({json.dumps(v['notes'])})\">🔎 notes</span>"
            )
        notes = f"<div class='mut'>{esc(v['notes'])}</div>" if v["notes"] else ""
        cprov = (
            f"<div class='mut'>⚖️ count basis: {esc(v['count_prov'])}</div>"
            if v.get("count_prov")
            else ""
        )
        quotes = "".join(f"<div class='q'>❝ {esc(q)}</div>" for q in v["quotes"][:2])
        left.append(
            f"<div class='rec'><span class='vn'>{esc(v['variant'])}</span>"
            f" <span class='mut'>{esc(v['sig'] or '')}</span><div>{''.join(actions)}</div>"
            f"{quotes}{notes}{cprov}</div>"
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
                    f"<span class='btn' onclick=\"jump({json.dumps(p['evidence'])})\">🔎 evidence sentence</span>"
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

    # RIGHT pane: the exact on-disk source + figures + supplements
    if ft.exists():
        src_html = md_to_html(ft.read_text(encoding="utf-8", errors="replace"))
        src_note = f"<div class='mut'>Rendered from on-disk source the LLM parsed: <code>{esc(_rel(ft, out_dir))}</code></div>"
    else:
        src_html = "<p class='note'>No full-text source on disk for this PMID (abstract-only or not fetched).</p>"
        src_note = ""
    figdir = paper_root / f"{pmid}_figures"
    figs = ""
    if figdir.is_dir():
        imgs = [
            f
            for f in sorted(figdir.iterdir())
            if f.suffix.lower() in {".png", ".jpg", ".jpeg", ".gif"}
        ]
        if imgs:
            figs = "<h3>Figures on disk:</h3>" + "".join(
                f"<a href='{esc(_rel(im, out_dir))}' target='_blank'><img src='{esc(_rel(im, out_dir))}' "
                f"style='max-width:160px;margin:4px;border:1px solid var(--line);border-radius:4px'></a>"
                for im in imgs[:30]
            )
    supdir = paper_root / f"{pmid}_supplements"
    sups = ""
    if supdir.is_dir():
        files = [f for f in sorted(supdir.rglob("*")) if f.is_file()]
        if files:
            sups = "<h3>Supplements on disk:</h3>" + "".join(
                f"<div><a href='{esc(_rel(f, out_dir))}' target='_blank'>📎 {esc(f.name)}</a></div>"
                for f in files[:40]
            )

    title = paper["title"] or f"PMID {pmid}"
    cite_line = f"<div class='sub'>{cite}</div>" if cite else ""
    body = (
        f"<header><h1>{esc(title)}</h1>{cite_line}"
        f"<div class='links sub'>PMID {esc(pmid)} · {' '.join(link_html)}</div>"
        f"<div>{''.join(meta_bits)}</div></header>"
        f"<div class='wrap'><div class='split'>"
        f"<div class='pane'>{''.join(left)}</div>"
        f"<div class='pane'><input class='search' id='src-s' placeholder='search the full text…' "
        f"onkeydown=\"if(event.key==='Enter')srcfind()\">{src_note}"
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


def render_gene_page(
    gene: str,
    corpus_rows: dict,
    db: Optional[dict],
    out_dir: Path,
    corpus_dir: Path,
    max_papers: int,
) -> tuple[str, dict, list[str]]:
    s = gene_stats(gene, corpus_rows, db)
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

    # Only papers we can actually render a page for (have DB rows or on-disk source),
    # ranked by how much there is to adjudicate so the cap keeps the richest papers.
    pageable = [
        p
        for p in paper_pmids
        if db is not None
        and (p in (db or {}).get("papers", {}) or (corpus_dir / gene / p).is_dir())
    ]
    pageable.sort(key=_rec_count, reverse=True)
    cap_set = set(pageable if max_papers <= 0 else pageable[:max_papers])
    for pmid in paper_pmids:
        cr = corpus_rows.get(pmid, {})
        status = cr.get("full_text_status", "—")
        nvar = (
            len((db or {}).get("by_pmid", {}).get(pmid, {}).get("variants", []))
            if db
            else 0
        )
        npat = (
            len((db or {}).get("by_pmid", {}).get(pmid, {}).get("patients", []))
            if db
            else 0
        )
        title = (db or {}).get("papers", {}).get(pmid, {}).get("title") or ""
        badge = "ok" if status == "ok" else "bad" if status == "stub" else "warn"
        has_page = (
            db is not None
            and pmid in cap_set
            and (
                pmid in (db or {}).get("papers", {})
                or (corpus_dir / gene / pmid).is_dir()
            )
        )
        pmid_cell = (
            f"<a href='paper_{esc(gene)}_{esc(pmid)}.html'>{esc(pmid)}</a>"
            if has_page
            else f"<a href='{pubmed_url(pmid)}' target='_blank'>{esc(pmid)}↗</a>"
        )
        if has_page:
            written.append(pmid)
        rows.append(
            f"<tr><td data-v='{esc(pmid)}'>{pmid_cell}</td><td>{esc(title[:90])}</td>"
            f"<td><span class='tag {badge}'>{esc(status)}</span></td>"
            f"<td data-v='{nvar}'>{nvar}</td><td data-v='{npat}'>{npat}</td></tr>"
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
        f"<header><h1>{esc(gene)}</h1><div class='sub'><a href='index.html'>← overview</a></div></header>"
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
        f"<div class='card'><h2>Provenance completeness</h2>{audit_html or '<span class=mut>no DB</span>'}"
        f"<div style='margin-top:8px'>{method_html}</div></div>"
        f"<div class='card'><h2>What's left</h2>"
        f"<div class='kv'><span>papers without usable full text</span><b>{len(stub_pmids)}</b></div>"
        f"<div class='mut' style='font-size:12px;margin-top:6px'>re-run <code>gvf gvf-run {esc(gene)} …</code> "
        f"(corpus cache skips what's done; a new publisher key re-fetches stubs)</div></div>"
        f"</div>"
        f"<h2 style='margin-top:20px'>Papers</h2>"
        f"<input class='search' id='pt-s' placeholder='filter papers…' onkeyup=\"filt('pt')\">"
        f"<table id='pt' data-sc=''><thead><tr>"
        f"<th onclick=\"srt('pt',0)\">PMID</th><th onclick=\"srt('pt',1)\">Title</th>"
        f"<th onclick=\"srt('pt',2)\">Source</th><th onclick=\"srt('pt',3)\">Variants</th>"
        f"<th onclick=\"srt('pt',4)\">Patients</th></tr></thead><tbody>{''.join(rows)}</tbody></table>"
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
                f"<a href='paper_{esc(gene)}_{esc(pm)}.html'>{esc(pm)}</a>"
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
            f"<td data-v='{a['affected']}'>{a['affected']}</td>"
            f"<td data-v='{a['unaffected']}'>{a['unaffected']}</td>"
            f"<td>{' '.join(links)}{more}</td></tr>"
        )
    body = (
        f"<header><h1>{esc(gene)} — variants</h1>"
        f"<div class='sub'><a href='{esc(gene)}.html'>← {esc(gene)}</a> · "
        f"<a href='index.html'>overview</a> · {len(agg)} unique variants</div></header>"
        f"<div class='wrap'>"
        f"<input class='search' id='vt-s' placeholder='filter variants…' onkeyup=\"filt('vt')\">"
        f"<table id='vt' data-sc=''><thead><tr>"
        f"<th onclick=\"srt('vt',0)\">Variant</th><th onclick=\"srt('vt',1)\">Clinical sig</th>"
        f"<th onclick=\"srt('vt',2)\">Papers</th><th onclick=\"srt('vt',3)\">Affected</th>"
        f"<th onclick=\"srt('vt',4)\">Unaffected</th><th>Papers (→ adjudicate / PubMed)</th>"
        f"</tr></thead><tbody>{''.join(rows)}</tbody></table></div>"
    )
    return _page(f"{gene} variants — GVF", body)


def render_index(summaries: dict[str, dict], generated: str) -> str:
    cards = []
    for gene, s in sorted(summaries.items()):
        cards.append(
            f"<div class='card'><h2><a href='{esc(gene)}.html'>{esc(gene)}</a></h2>"
            f"<div class='kv'><span>papers in corpus</span><b>{s['in_corpus']}</b></div>"
            f"{_barbar(s['usable'], s['in_corpus'])}"
            f"<div class='kv'><span>usable full text</span><b>{s['usable']} ({pct(s['usable'], s['in_corpus']):.0f}%)</b></div>"
            f"<div class='kv'><span>extracted</span><b>{s['extracted']}</b></div>"
            f"<div class='kv'><span>papers w/ variants</span><b>{s['with_variants']}</b></div>"
            f"<div class='kv'><span>unique variants</span><b>{s['unique']}</b></div>"
            f"<div style='margin-top:8px'><a href='{esc(gene)}.html'>status</a> · "
            f"<a href='{esc(gene)}_variants.html'>variants</a></div></div>"
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
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    corpus_idx = load_corpus_index(corpus_dir)
    target_genes = genes or sorted(corpus_idx.keys())
    summaries: dict[str, dict] = {}
    stats = {"genes": 0, "paper_pages": 0}
    for gene in target_genes:
        corpus_rows = corpus_idx.get(gene, {})
        if not corpus_rows and gene not in db_map:
            continue
        db_path = db_map.get(gene) or find_latest_db(gene)
        db = load_db(db_path) if db_path and db_path.exists() else None
        page, s, paper_pmids = render_gene_page(
            gene, corpus_rows, db, out_dir, corpus_dir, max_papers
        )
        (out_dir / f"{gene}.html").write_text(page, encoding="utf-8")
        summaries[gene] = s
        stats["genes"] += 1
        if db:
            for pmid in paper_pmids:
                (out_dir / f"paper_{gene}_{pmid}.html").write_text(
                    render_paper_page(gene, pmid, db, corpus_dir, out_dir),
                    encoding="utf-8",
                )
                stats["paper_pages"] += 1
            (out_dir / f"{gene}_variants.html").write_text(
                render_variants_page(gene, db, set(paper_pmids)), encoding="utf-8"
            )
    (out_dir / "index.html").write_text(
        render_index(summaries, generated), encoding="utf-8"
    )
    stats["summaries"] = summaries
    return stats
