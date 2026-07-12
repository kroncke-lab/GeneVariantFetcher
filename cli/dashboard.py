"""Static-HTML status / coverage / missingness / provenance dashboard.

`gvf dashboard` reads the consolidated corpus (`corpus/INDEX.csv` + the per-paper
source under `corpus/<GENE>/<PMID>/`) and the scored SQLite DB(s), and writes an
offline static dashboard under `--out` (default `corpus/dashboard/`). Figures and
supplements are referenced relatively into the corpus (no duplication), so zip or
serve the whole `corpus/` to share it:

* `index.html` — per-gene cards: source coverage (usable vs stub), extraction
  funnel (in-corpus -> extracted -> has-variants), variant/unique counts, plus
  trust / final-check / since-last-run badges (the gold-free progress read).
* `<GENE>/index.html` — a "since last run" delta line, a **Run health** card (gold-free
  confidence: funnel + trust-tier quarantine + per-paper final-check), a **What to
  change next** ranked worklist (acquisition / extraction / trust / review / gap
  levers, each linking to the affected papers), the coverage-by-method facet, a
  provenance-completeness audit, the "what's left" list, and a
  sortable/filterable paper table (with a final-check Check column) linking to
  per-paper adjudication pages.
* `<GENE>/papers/<PMID>.html` — the ADJUDICATION view: paper header with one-click links
  to PubMed / DOI / PMC, the extracted records (variant provenance +
  per-patient characteristics) on the left, and the DB-recorded source file
  (`extraction_metadata.source_file`, falling back to the corpus full text)
  rendered on the right. Clicking a record jumps to and highlights the
  sentence/section it was extracted from, so you can verify it yourself.
* `<GENE>/papers/<PMID>_process.html` — the searchable paper-level process trail:
  discovery/filter decisions, acquisition/source files, QC and refresh evidence,
  extraction versions, the complete selected-DB snapshot, trust/final-check state,
  and scoring/disagreement records with links to the full local artifacts.

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
from datetime import datetime
from pathlib import Path
from typing import Optional

# Reuse the trust-gate reader rather than re-deriving the two-tier split here
# (scripts is a package cli/ already imports from). Guarded so the dashboard
# still renders if the module ever moves.
try:
    from scripts.trust_report import list_quarantined, summarize_trust
except Exception:  # pragma: no cover - dashboard degrades without trust readers
    summarize_trust = None  # type: ignore[assignment]
    list_quarantined = None  # type: ignore[assignment]


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
    except (TypeError, ValueError, OverflowError):
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


def _load_paper_final_check(con: sqlite3.Connection) -> dict:
    """Read the soft per-paper final-check verdicts (pipeline/paper_final_check.py).

    Returns ``{"by_pmid": {pmid: {...}}, "counts": Counter}``; empty for DBs that
    predate the sniff test (no ``paper_final_check`` table). Never mutates.
    """
    out: dict = {"by_pmid": {}, "counts": Counter()}
    cols = _table_cols(con, "paper_final_check")
    if not cols:
        return out
    want = [
        c
        for c in (
            "pmid",
            "verdict",
            "confidence",
            "n_flagged",
            "n_missing",
            "completeness_status",
            "summary",
            "source_grounded",
        )
        if c in cols
    ]
    if not want:  # table exists but has none of the expected columns (schema drift)
        return out
    for r in con.execute("SELECT " + ",".join(want) + " FROM paper_final_check"):
        rd = dict(r)
        pmid = str(rd.get("pmid"))
        # SQLite is dynamically typed — coerce text-ish columns to str so later
        # slicing / HTML-escaping can't TypeError on an int/float cell.
        verdict = str(rd.get("verdict") or "").strip().lower()
        n_flagged = _to_int(rd.get("n_flagged"))
        n_missing = _to_int(rd.get("n_missing"))
        out["by_pmid"][pmid] = {
            "verdict": verdict,
            "confidence": rd.get("confidence"),
            "n_flagged": n_flagged,
            "n_missing": n_missing,
            "completeness_status": str(rd.get("completeness_status") or ""),
            "summary": str(rd.get("summary") or ""),
            "source_grounded": bool(rd.get("source_grounded")),
        }
        out["counts"]["total"] += 1
        out["counts"][verdict or "unknown"] += 1
        out["counts"]["flagged_facts"] += n_flagged
        out["counts"]["missing_carriers"] += n_missing
    return out


def load_db(db_path: Path) -> dict:
    """Pull the provenance/coverage surface from a scored DB. Degrades gracefully."""
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    data: dict = {
        "db_path": db_path,
        "papers": {},
        "by_pmid": defaultdict(lambda: {"variants": [], "patients": [], "tables": []}),
        "facts_by_pmid_variant": defaultdict(list),
        "facts_by_pmid_individual": defaultdict(list),
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
        "pen_rows_by_pmid": defaultdict(list),
        # soft signals (populated below; kept here so renderers can rely on them)
        "final_check": {"by_pmid": {}, "counts": Counter()},
        "trust": {"tiered": False},
        "quarantined": [],
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
                "SELECT pd.*, "
                "COALESCE(v.protein_notation, v.cdna_notation, v.genomic_position, '?') "
                "AS variant_label "
                "FROM penetrance_data pd LEFT JOIN variants v ON v.variant_id = pd.variant_id"
            ):
                rd = dict(r)
                pmid = str(rd["pmid"])
                c = _to_int(rd.get("total_carriers_observed"))
                a = _to_int(rd.get("affected_count"))
                u = _to_int(rd.get("unaffected_count"))
                data["pen_by_pmid"][pmid]["carriers"] += c
                data["pen_by_pmid"][pmid]["affected"] += a
                data["pen_by_pmid"][pmid]["unaffected"] += u
                data["pen_rows_by_pmid"][pmid].append(rd)
                vname = rd.get("variant_label")
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
        if _table_cols(con, "fact_provenance"):
            for r in con.execute(
                """
                SELECT fp.pmid,
                       COALESCE(v.protein_notation, v.cdna_notation, v.genomic_position, '?') AS variant,
                       fp.fact_type, fp.fact_value, fp.individual_id,
                       fp.source_location, fp.source_section, fp.source_paragraph,
                       fp.source_table, fp.source_row, fp.source_column,
                       fp.evidence_quote, fp.count_type, fp.provenance_kind
                FROM fact_provenance fp
                JOIN variants v ON v.variant_id = fp.variant_id
                ORDER BY fp.provenance_id
                """
            ):
                rd = dict(r)
                pmid = str(rd["pmid"])
                variant = rd.get("variant") or "?"
                fact = {
                    "type": rd.get("fact_type") or "",
                    "value": rd.get("fact_value") or "",
                    "individual_id": rd.get("individual_id") or "",
                    "location": rd.get("source_location") or "",
                    "section": rd.get("source_section") or "",
                    "paragraph": rd.get("source_paragraph") or "",
                    "table": rd.get("source_table") or "",
                    "row": rd.get("source_row") or "",
                    "column": rd.get("source_column") or "",
                    "quote": rd.get("evidence_quote") or "",
                    "count_type": rd.get("count_type") or "",
                    "kind": rd.get("provenance_kind") or "",
                }
                data["prov"]["with_fact_provenance"] += 1
                if fact["table"] or fact["row"] or fact["section"] or fact["quote"]:
                    data["prov"]["with_exact_fact_pointer"] += 1
                data["facts_by_pmid_variant"][(pmid, variant)].append(fact)
                if fact["individual_id"]:
                    data["facts_by_pmid_individual"][
                        (pmid, variant, fact["individual_id"])
                    ].append(fact)
        if _table_cols(con, "tables_processed"):
            for r in con.execute(
                "SELECT pmid, table_name, table_caption, variants_extracted FROM tables_processed"
            ):
                data["by_pmid"][str(r["pmid"])]["tables"].append(dict(r))
        # Soft per-paper final-check verdicts (never mutates counts; a signal).
        data["final_check"] = _load_paper_final_check(con)
    finally:
        con.close()

    # Trust-tier split (reuse scripts/trust_report). Its own read-only connection;
    # degrades to untiered on pre-gate DBs or ones without penetrance_data (the
    # summarizer raises sqlite3.Error there, which we swallow).
    data["trust"] = {"tiered": False}
    data["quarantined"] = []
    if summarize_trust is not None:
        try:
            data["trust"] = summarize_trust(db_path)
            if data["trust"].get("tiered") and list_quarantined is not None:
                data["quarantined"] = list_quarantined(db_path, 60)
        # ValueError covers json.JSONDecodeError from list_quarantined on a
        # malformed trust_reasons cell (summarize_trust guards its own json.loads,
        # so it can return tiered=True and still hand off a corrupt row).
        except (sqlite3.Error, OSError, ValueError):  # pragma: no cover - best-effort
            data["trust"] = {"tiered": False}
            data["quarantined"] = []
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


def _fact_pointer(fact: dict) -> str:
    parts = [
        fact.get("table"),
        fact.get("row"),
        fact.get("column"),
        fact.get("section"),
        fact.get("paragraph"),
        fact.get("location"),
    ]
    return " · ".join(str(p) for p in parts if p)


def _render_fact_provenance(facts: list[dict], limit: int = 6) -> str:
    if not facts:
        return ""
    bits = []
    seen = set()
    for fact in facts:
        key = (
            fact.get("type"),
            fact.get("value"),
            _fact_pointer(fact),
            fact.get("quote"),
        )
        if key in seen:
            continue
        seen.add(key)
        pointer = _fact_pointer(fact)
        label = f"{fact.get('type') or 'fact'}={fact.get('value') or '?'}"
        detail = " · ".join(
            x
            for x in [
                pointer,
                fact.get("count_type"),
                fact.get("kind"),
            ]
            if x
        )
        jump_arg = fact.get("quote") or pointer
        action = (
            f" <span class='btn' onclick=\"jump({jq(jump_arg)})\">evidence</span>"
            if jump_arg
            else ""
        )
        bits.append(
            f"<div class='mut'>↳ {esc(label)}"
            + (f" <span>{esc(detail)}</span>" if detail else "")
            + action
            + "</div>"
        )
        if len(bits) >= limit:
            break
    more = len(facts) - len(bits)
    if more > 0:
        bits.append(f"<div class='mut'>↳ +{more} more evidence facts</div>")
    return "".join(bits)


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


def _resolve_extraction_source(
    meta: dict, corpus_dir: Path, gene: str, pmid: str
) -> tuple[Path, str]:
    """Prefer the exact source file recorded in extraction_metadata.

    Older DBs or externally generated dashboards may not have a resolvable
    source_file, so the corpus FULL_CONTEXT copy remains the fallback.
    """
    fallback = corpus_dir / gene / pmid / f"{pmid}_FULL_CONTEXT.md"
    raw = (meta or {}).get("source_file")
    if raw:
        p = Path(str(raw)).expanduser()
        candidates = (
            [p] if p.is_absolute() else [Path.cwd() / p, REPO / p, corpus_dir / p]
        )
        for cand in candidates:
            if cand.is_file():
                return cand, "extraction_metadata.source_file"
    return fallback, "corpus FULL_CONTEXT fallback"


def _display_source_path(path: Path) -> str:
    resolved = path.expanduser().resolve()
    for base in (Path.cwd(), REPO):
        try:
            return resolved.relative_to(base.resolve()).as_posix()
        except ValueError:
            continue
    return resolved.name


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
#art{table-layout:fixed}#art th:nth-child(1){width:12%}#art th:nth-child(2){width:28%}
#art th:nth-child(3){width:8%}#art th:nth-child(4){width:12%}#art th:nth-child(5){width:40%}
#art td,#art code{overflow-wrap:anywhere;word-break:break-word}#art code{white-space:normal}
#art details{max-width:100%}#art pre{white-space:pre-wrap;overflow-wrap:anywhere;max-height:28rem;overflow:auto}
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
.stagecard{background:#fff;border:1px solid var(--line);border-radius:12px;margin:14px 0;padding:16px;box-shadow:0 1px 3px rgba(15,23,42,.05)}
.stagehead{display:flex;align-items:flex-start;gap:12px;border-bottom:1px solid var(--line);padding-bottom:10px;margin-bottom:10px}
.stagehead h2{margin:0;font-size:17px}.stagehead p{margin:2px 0 0;color:var(--mut);font-size:12px}.stagehead>.tag{margin-left:auto}
.stage-num{display:grid;place-items:center;min-width:34px;height:34px;border-radius:50%;background:var(--acc);color:#fff;font-weight:800;font-size:16px}
.artifact-card{border:1px solid var(--line);border-radius:9px;padding:10px 12px;margin:8px 0;background:#fbfcfe;overflow:hidden}
.artifact-card code{overflow-wrap:anywhere;word-break:break-word;white-space:normal}.dump{white-space:pre-wrap;overflow-wrap:anywhere;max-height:48rem;overflow:auto;background:#0f172a;color:#dbeafe;padding:12px;border-radius:8px;font-size:11px}
.muted-empty{color:var(--mut);font-style:italic;padding:8px 2px}
.stage-nav{display:flex;gap:6px;flex-wrap:wrap;position:sticky;top:0;z-index:6;background:rgba(238,241,246,.96);padding:8px 0;margin:4px 0 10px}
.stage-nav a{font-size:11px;padding:4px 9px;border:1px solid var(--line);border-radius:20px;background:#fff}.stagecard{scroll-margin-top:54px}
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

TRAIL_JS = """
function trailFilter(){var q=(document.getElementById('trail-s').value||'').toLowerCase();
 document.querySelectorAll('.trailitem').forEach(function(x){x.style.display=!q||x.textContent.toLowerCase().includes(q)?'':'none'});
 document.querySelectorAll('.stagecard').forEach(function(s){var any=[].slice.call(s.querySelectorAll('.trailitem')).some(function(x){return x.style.display!='none'});s.style.display=any?'':'none'});}
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

_LOCAL_PATH_RE = re.compile(
    r"((?:/[Uu]sers/|[A-Za-z]:[/\\][Uu]sers[/\\])[^\s<>'\"`\]\)]+"
    r"|/private/var/[^\s<>'\"`\]\)]+|/tmp/[^\s<>'\"`\]\)]+)"
)


def _sanitize_local_paths(text: str) -> str:
    # Split on both separators so Windows (backslash) paths are reduced to their
    # basename too — Path().name only splits on "/" on POSIX and would leak them.
    def _basename(match: str) -> str:
        return re.split(r"[/\\]", match.rstrip("/\\"))[-1] or "redacted"

    return _LOCAL_PATH_RE.sub(lambda m: f"[local path: {_basename(m.group(0))}]", text)


def _rel(target: Path, start: Path) -> str:
    return os.path.relpath(target, start)


def _display_db_source(db_path: Path) -> str:
    """A publish-safe DB identifier: repo-relative when possible, else basename."""
    try:
        return str(db_path.resolve().relative_to(REPO.resolve()))
    except (OSError, ValueError):
        return db_path.name


PROCESS_STAGES = {
    0: "Run setup & status",
    1: "Discovery & filtering",
    2: "Source acquisition",
    3: "Source QC & refresh",
    4: "Extraction",
    5: "Normalization, DB & trust",
    6: "Scoring & review",
    7: "Cross-stage / other",
}
_PROCESS_SUFFIXES = {".json", ".jsonl", ".csv", ".tsv", ".md", ".txt", ".log", ".db"}
_PROCESS_FULL_NAMES = {
    "run_status.json",
    "run_report.md",
    "source_completeness.json",
    "source_acquisition_audit.json",
    "refresh_summary.json",
}
_PROCESS_PER_STAGE_LIMIT = 1500


def _artifact_stage(path: Path) -> int:
    """Classify a run artifact into the stage where an operator acts on it."""
    low = path.as_posix().lower()
    name = path.name.lower()
    if name in {"run_status.json", "run_report.md"} or "doctor" in low:
        return 0
    if path.suffix.lower() == ".db":
        return 5
    if any(x in low for x in ("recall", "precision", "adjudicat", "review", "score")):
        return 6
    if any(
        x in low for x in ("trust", "clinvar", "pubtator", "migration", "penetrance")
    ):
        return 5
    if any(x in low for x in ("extraction", "table", "scout", "retry")):
        return 4
    if any(
        x in low for x in ("discover", "abstract", "filter", "pmid_status", "_pmids")
    ):
        return 1
    if any(
        x in low
        for x in (
            "source_acquisition",
            "source_completeness",
            "acquisition_outcome",
            "refresh_",
            "worklist",
            "supplement_fold",
            "source_qc",
        )
    ):
        return 3
    if any(
        x in low
        for x in (
            "pmc_fulltext",
            "full_context",
            "cleaned",
            "data_zones",
            "paywall",
            "supplement",
            "figure",
            "fetch",
        )
    ):
        return 2
    if "recovery" in low:
        return 5
    if any(
        x in low
        for x in (
            "recall",
            "precision",
            "adjudicat",
            "review",
            "score",
            "report",
            "summary",
        )
    ):
        return 6
    if path.suffix.lower() == ".log":
        return 0
    return 7


def _read_artifact_sample(path: Path, size: int) -> tuple[str, str]:
    """Return (text, mode): full important logs/reports, sampled large dumps."""
    important = (
        path.name.lower() in _PROCESS_FULL_NAMES or path.suffix.lower() == ".log"
    )
    if important and size <= 128_000:
        try:
            return _sanitize_local_paths(
                path.read_text(encoding="utf-8", errors="replace")
            ), "full"
        except OSError:
            return "", "unreadable"
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            head = fh.read(700)
        if size <= 700:
            return _sanitize_local_paths(head), "full"
        return _sanitize_local_paths(head), "sample"
    except (OSError, UnicodeError):
        return "", "binary" if path.suffix.lower() == ".db" else "unreadable"


def discover_process_artifacts(
    gene: str, db_path: Optional[Path], corpus_dir: Path
) -> list[dict]:
    """Inventory searchable stage artifacts without copying the underlying dumps."""
    roots: list[Path] = [corpus_dir / gene]
    explicit: list[Path] = []
    if db_path:
        roots.insert(0, db_path.parent)
        explicit.append(db_path)
    metrics_root = REPO / "recall_metrics"
    if metrics_root.is_dir():
        summaries = list(metrics_root.glob("*/summary.json"))
        if summaries:
            latest = max(summaries, key=lambda p: p.stat().st_mtime).parent
            roots.append(latest / gene)
            explicit.extend(
                p
                for p in (
                    latest / "summary.json",
                    latest / "paper_disagreement_report.csv",
                    latest / "precision_summary.json",
                )
                if p.is_file()
            )
    roots.extend([REPO / "docs" / "dashboard"])

    seen: set[Path] = set()
    candidates: list[dict] = []
    paths: list[Path] = list(explicit)
    for root in roots:
        if not root.is_dir():
            continue
        paths.extend(root.rglob("*"))
    for path in paths:
        if not path.is_file() or path.suffix.lower() not in _PROCESS_SUFFIXES:
            continue
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        # Only the selected live DB is useful here; backups and staging DBs
        # are large, binary, and make the process trail look duplicated.
        if path.suffix.lower() == ".db" and (
            db_path is None or resolved != db_path.resolve()
        ):
            continue
        try:
            st = path.stat()
        except OSError:
            continue
        try:
            rel = resolved.relative_to(REPO.resolve()).as_posix()
        except ValueError:
            rel = path.name
        candidates.append(
            {
                "path": path,
                "relative": rel,
                "stage": _artifact_stage(path),
                "size": st.st_size,
                "mtime": st.st_mtime,
                "modified": datetime.fromtimestamp(st.st_mtime).strftime(
                    "%Y-%m-%d %H:%M"
                ),
            }
        )

    # Historical refreshes often contain byte-for-byte-ish copies under new
    # directories. Collapse copies by stage/name/size, keeping the newest path.
    dedup: dict[tuple, dict] = {}
    for artifact in candidates:
        key = (artifact["stage"], artifact["path"].name.lower(), artifact["size"])
        current = dedup.get(key)
        if current is None or artifact["mtime"] > current["mtime"]:
            dedup[key] = artifact

    selected: list[dict] = []
    by_stage: dict[int, list[dict]] = defaultdict(list)
    for artifact in dedup.values():
        by_stage[artifact["stage"]].append(artifact)
    for stage, items in by_stage.items():
        items.sort(
            key=lambda a: (
                a["path"].name.lower() in _PROCESS_FULL_NAMES
                or a["path"].suffix.lower() == ".log",
                a["mtime"],
                a["size"],
            ),
            reverse=True,
        )
        selected.extend(items[:_PROCESS_PER_STAGE_LIMIT])

    for artifact in selected:
        artifact["text"], artifact["mode"] = _read_artifact_sample(
            artifact["path"], artifact["size"]
        )
    return sorted(selected, key=lambda a: (a["stage"], a["relative"]))


def _human_bytes(n: int) -> str:
    value = float(n)
    for unit in ("B", "KB", "MB", "GB"):
        if value < 1024 or unit == "GB":
            return f"{value:.0f} {unit}" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024
    return f"{n} B"


def render_process_page(
    gene: str, artifacts: list[dict], gene_dir: Path, db_path: Optional[Path]
) -> str:
    """Searchable stage explorer with full small logs and sampled large dumps."""
    counts = Counter(a["stage"] for a in artifacts)
    stage_cards = "".join(
        f"<div class='card'><div class='big'>{counts.get(stage, 0)}</div>"
        f"<div class='mut'>Stage {stage}: {esc(label)}</div></div>"
        for stage, label in PROCESS_STAGES.items()
        if counts.get(stage)
    )
    rows = []
    repo_root = REPO.resolve()
    for artifact in artifacts:
        path = artifact["path"].resolve()
        try:
            path.relative_to(repo_root)
            raw_link = f"<a href='{esc(_rel(path, gene_dir))}' target='_blank'>open full local file</a>"
        except ValueError:
            raw_link = "<span class='mut'>outside repo; preview only</span>"
        mode = artifact["mode"]
        label = {
            "full": "full inline",
            "sample": "sample from large dump",
            "binary": "binary metadata only",
            "unreadable": "unreadable",
        }.get(mode, mode)
        preview = (
            f"<details><summary>{esc(label)}</summary><pre>{esc(artifact['text'])}</pre></details>"
            if artifact["text"]
            else f"<span class='mut'>{esc(label)}</span>"
        )
        stage = artifact["stage"]
        rows.append(
            f"<tr><td data-v='{stage}'><span class='tag'>Stage {stage}</span><br>"
            f"<span class='mut'>{esc(PROCESS_STAGES[stage])}</span></td>"
            f"<td><code>{esc(artifact['relative'])}</code><div>{raw_link}</div></td>"
            f"<td data-v='{artifact['size']}'>{esc(_human_bytes(artifact['size']))}</td>"
            f"<td>{esc(artifact['modified'])}</td><td>{preview}</td></tr>"
        )
    db_note = _display_db_source(db_path) if db_path else "none selected"
    body = (
        f"<header><h1>{esc(gene)} — Process Explorer</h1>"
        f"<div class='sub'><a href='index.html'>← gene status</a> · database: "
        f"<code>{esc(db_note)}</code></div></header><div class='wrap'>"
        "<div class='note'>Search paths and indexed content across every pipeline stage. "
        "Small logs/status reports are shown in full; large JSON/CSV/text dumps show a sample. "
        "Repeat copies from historical refreshes are collapsed and the most recent 1,500 artifacts "
        "per stage are indexed. The full-file links are local-only because run artifacts are intentionally not published.</div>"
        f"<div class='grid' style='margin-bottom:16px'>{stage_cards}</div>"
        f"<h2>Artifacts ({len(artifacts)})</h2>"
        "<input class='search' id='art-s' placeholder='search stage, PMID, filename, error, model, status, variant…' onkeyup=\"filt('art')\">"
        "<table id='art' data-sc=''><thead><tr>"
        "<th onclick=\"srt('art',0)\">Stage</th><th onclick=\"srt('art',1)\">Artifact</th>"
        "<th onclick=\"srt('art',2)\">Size</th><th onclick=\"srt('art',3)\">Modified</th>"
        f"<th>Data / log preview</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div>"
    )
    return _page(f"{gene} — process explorer", body)


PAPER_STAGE_DESCRIPTIONS = {
    0: "Run configuration, orchestration status, and workflow-log evidence for this PMID.",
    1: "How the paper was discovered, represented by its abstract, and passed or failed filtering.",
    2: "What source text, figures, supplements, and publisher-acquisition artifacts landed on disk.",
    3: "Source-quality audits, recovery worklists, supplement folding, and targeted refresh decisions.",
    4: "Every captured extraction input/output version, including staged re-extractions and parser products.",
    5: "The selected SQLite result: paper metadata, extracted facts, counts, provenance, and trust decisions.",
    6: "Gold comparison, disagreement rows, final-check verdicts, and review/adjudication evidence.",
    7: "Additional paper-specific artifacts that do not map cleanly to one pipeline stage.",
}
_PAPER_TEXT_SUFFIXES = {
    ".json",
    ".jsonl",
    ".csv",
    ".tsv",
    ".md",
    ".txt",
    ".log",
    ".xml",
    ".html",
}
_PAPER_SHARED_MARKERS = (
    "audit",
    "candidate",
    "discrep",
    "failure",
    "filter",
    "manifest",
    "missing",
    "outcome",
    "pmids",
    "progress",
    "report",
    "refresh",
    "status",
    "summary",
    "worklist",
)
_PMID_TOKEN_RE = re.compile(r"(?<!\d)(\d{6,9})(?!\d)")


def _latest_metrics_dir() -> Optional[Path]:
    root = REPO / "recall_metrics"
    if not root.is_dir():
        return None
    summaries = list(root.glob("*/summary.json"))
    return max(summaries, key=lambda p: p.stat().st_mtime).parent if summaries else None


def _paper_scan_roots(
    gene: str, db_path: Optional[Path], corpus_dir: Path
) -> list[Path]:
    roots = [corpus_dir / gene]
    if db_path:
        roots.insert(0, db_path.parent)
    latest = _latest_metrics_dir()
    if latest:
        roots.extend([latest / gene, latest])
    return [r for r in roots if r.is_dir()]


def _paper_path_record(
    path: Path, kind: str = "direct", lines: Optional[list] = None
) -> dict:
    st = path.stat()
    try:
        relative = path.resolve().relative_to(REPO.resolve()).as_posix()
    except ValueError:
        relative = path.name
    return {
        "path": path,
        "relative": relative,
        "stage": _artifact_stage(path),
        "size": st.st_size,
        "mtime": st.st_mtime,
        "modified": datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M"),
        "kind": kind,
        "lines": lines or [],
    }


def _line_excerpt(line: str, pmid: str, radius: int = 450) -> str:
    low = line.find(pmid)
    if low < 0:
        return line[: radius * 2]
    start = max(0, low - radius)
    end = min(len(line), low + len(pmid) + radius)
    prefix = "…" if start else ""
    suffix = "…" if end < len(line) else ""
    return prefix + line[start:end].strip() + suffix


def build_paper_process_index(
    gene: str,
    db_path: Optional[Path],
    corpus_dir: Path,
    pmids: set[str],
    selected_artifacts: list[dict],
) -> dict[str, list[dict]]:
    """Build one reusable PMID→artifact map for all paper process pages.

    Direct paper files are indexed by PMID in their path. Shared operational
    logs/reports are scanned once and contribute only the lines around each PMID.
    """
    out: dict[str, list[dict]] = defaultdict(list)
    seen_direct: set[Path] = set()
    for root in _paper_scan_roots(gene, db_path, corpus_dir):
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            resolved = path.resolve()
            if resolved in seen_direct:
                continue
            path_tokens = set(path.parts) & pmids
            name_tokens = set(_PMID_TOKEN_RE.findall(path.name)) & pmids
            matched = path_tokens | name_tokens
            if not matched:
                continue
            seen_direct.add(resolved)
            try:
                record = _paper_path_record(path)
            except OSError:
                continue
            for pmid in matched:
                out[pmid].append(record)

    # Collapse repeated shared reports before scanning their contents. The
    # 24-MB workflow log is intentionally retained: its matching lines explain
    # orchestration behavior that per-paper JSON alone cannot show.
    shared: dict[tuple, dict] = {}
    for artifact in selected_artifacts:
        path = artifact["path"]
        name = path.name.lower()
        if path.suffix.lower() not in _PAPER_TEXT_SUFFIXES:
            continue
        if path.suffix.lower() != ".log" and not any(
            x in name for x in _PAPER_SHARED_MARKERS
        ):
            continue
        if artifact["size"] > 40_000_000:
            continue
        key = (artifact["stage"], name, artifact["size"])
        current = shared.get(key)
        if current is None or artifact["mtime"] > current["mtime"]:
            shared[key] = artifact

    for artifact in sorted(shared.values(), key=lambda a: a["mtime"], reverse=True)[
        :600
    ]:
        path = artifact["path"]
        matches: dict[str, list[str]] = defaultdict(list)
        try:
            with path.open("r", encoding="utf-8", errors="replace") as fh:
                for lineno, line in enumerate(fh, 1):
                    hits = set(_PMID_TOKEN_RE.findall(line)) & pmids
                    for pmid in hits:
                        if len(matches[pmid]) < 24:
                            excerpt = _sanitize_local_paths(_line_excerpt(line, pmid))
                            matches[pmid].append(f"line {lineno}: {excerpt}")
        except (OSError, UnicodeError):
            continue
        for pmid, lines in matches.items():
            try:
                out[pmid].append(_paper_path_record(path, kind="shared", lines=lines))
            except OSError:
                continue

    for pmid, records in out.items():
        # Preserve historical extraction versions but bound pathological trees.
        by_stage: dict[int, list[dict]] = defaultdict(list)
        for record in records:
            by_stage[record["stage"]].append(record)
        kept: list[dict] = []
        for stage, items in by_stage.items():
            items.sort(key=lambda a: (a["mtime"], a["size"]), reverse=True)
            kept.extend(items[:36])
        out[pmid] = sorted(kept, key=lambda a: (a["stage"], -a["mtime"]))
    return out


def _read_paper_artifact(path: Path, size: int, stage: int) -> tuple[str, str]:
    """Inline full paper JSON/abstracts when practical; sample large sources/logs."""
    if path.suffix.lower() not in _PAPER_TEXT_SUFFIXES:
        return "", "binary / external file"
    full_limit = 320_000 if stage in {1, 4, 5, 6} else 120_000
    try:
        if size <= full_limit:
            raw = path.read_bytes().decode("utf-8", errors="replace")
            return _sanitize_local_paths(raw), "full file"
        with path.open("rb") as fh:
            head = fh.read(5000)
            fh.seek(max(0, size - 5000))
            tail = fh.read(5000)
        raw = head.decode("utf-8", errors="replace")
        raw += "\n\n… [middle omitted; open full local file] …\n\n"
        raw += tail.decode("utf-8", errors="replace")
        return _sanitize_local_paths(raw), "head + tail sample"
    except OSError:
        return "", "unreadable"


def _paper_db_snapshot(db: dict, pmid: str) -> dict:
    rec = db.get("by_pmid", {}).get(
        pmid, {"variants": [], "patients": [], "tables": []}
    )
    facts = {
        variant: entries
        for (fact_pmid, variant), entries in db.get("facts_by_pmid_variant", {}).items()
        if fact_pmid == pmid
    }
    return {
        "paper": db.get("papers", {}).get(pmid),
        "extracted_records": rec,
        "cohort_totals": db.get("pen_by_pmid", {}).get(pmid),
        "penetrance_rows": db.get("pen_rows_by_pmid", {}).get(pmid, []),
        "fact_provenance": facts,
        "paper_final_check": db.get("final_check", {}).get("by_pmid", {}).get(pmid),
    }


def _render_trail_artifact(record: dict, papers_dir: Path) -> str:
    path = record["path"].resolve()
    try:
        repo_relative = path.relative_to(REPO.resolve())
        papers_dir.resolve().relative_to(REPO.resolve())
        href = _rel(path, papers_dir)
    except ValueError:
        try:
            repo_relative = path.relative_to(REPO.resolve())
            href = "/" + repo_relative.as_posix()
        except ValueError:
            repo_relative = None
            href = ""
    raw_link = (
        f"<a href='{esc(href)}' target='_blank'>open full local file</a>"
        if repo_relative is not None
        else "<span class='mut'>outside repository; indexed excerpt only</span>"
    )
    if record.get("kind") == "shared":
        text = "\n".join(record.get("lines") or [])
        mode = f"{len(record.get('lines') or [])} PMID-matching line(s)"
    else:
        text, mode = _read_paper_artifact(path, record["size"], record["stage"])
    detail = (
        f"<details{' open' if mode == 'full file' and record['stage'] in {1, 4} else ''}>"
        f"<summary>{esc(mode)}</summary><pre class='dump'>{esc(text)}</pre></details>"
        if text
        else f"<div class='mut'>{esc(mode)}</div>"
    )
    return (
        "<div class='trailitem artifact-card'>"
        f"<div><code>{esc(record['relative'])}</code></div>"
        f"<div class='mut'>{esc(_human_bytes(record['size']))} · {esc(record['modified'])} · {raw_link}</div>"
        f"{detail}</div>"
    )


def render_paper_process_page(
    gene: str,
    pmid: str,
    db: dict,
    artifacts: list[dict],
    papers_dir: Path,
    score: Optional[dict] = None,
) -> str:
    """Expansive, searchable stage-by-stage explanation for one paper."""
    paper = db.get("papers", {}).get(pmid, {})
    title = paper.get("title") or f"PMID {pmid}"
    by_stage: dict[int, list[dict]] = defaultdict(list)
    for record in artifacts:
        by_stage[record["stage"]].append(record)

    db_json = _sanitize_local_paths(
        json.dumps(_paper_db_snapshot(db, pmid), indent=2, default=str)
    )
    by_stage[5].append(
        {
            "synthetic": True,
            "html": "<div class='trailitem artifact-card'><b>Selected database snapshot</b>"
            f"<div class='mut'><code>{esc(_display_db_source(db['db_path']))}</code></div>"
            f"<details open><summary>full paper rows and provenance</summary><pre class='dump'>{esc(db_json)}</pre></details></div>",
        }
    )
    final_check = db.get("final_check", {}).get("by_pmid", {}).get(pmid)
    if score or final_check:
        score_json = _sanitize_local_paths(
            json.dumps(
                {"gold_comparison": score, "paper_final_check": final_check},
                indent=2,
                default=str,
            )
        )
        by_stage[6].append(
            {
                "synthetic": True,
                "html": "<div class='trailitem artifact-card'><b>Current scoring and review state</b>"
                f"<details open><summary>paper result</summary><pre class='dump'>{esc(score_json)}</pre></details></div>",
            }
        )

    sections = []
    for stage, label in PROCESS_STAGES.items():
        items = by_stage.get(stage, [])
        rendered = "".join(
            item["html"]
            if item.get("synthetic")
            else _render_trail_artifact(item, papers_dir)
            for item in items
        )
        if not rendered:
            rendered = "<div class='trailitem muted-empty'>No paper-specific artifact captured for this stage.</div>"
        sections.append(
            f"<section class='stagecard' id='stage-{stage}'><div class='stagehead'>"
            f"<span class='stage-num'>{stage}</span><div><h2>{esc(label)}</h2>"
            f"<p>{esc(PAPER_STAGE_DESCRIPTIONS[stage])}</p></div>"
            f"<span class='tag'>{len(items)} item(s)</span></div>{rendered}</section>"
        )

    body = (
        f"<header><h1>{esc(title)}</h1><div class='sub'>PMID {esc(pmid)} · "
        f"<a href='{esc(pmid)}.html'>evidence &amp; full source</a> · "
        f"<a href='../index.html'>← {esc(gene)}</a> · <a href='../process.html'>gene process explorer</a>"
        "</div></header><div class='wrap'>"
        "<div class='note'>This page reconstructs what the pipeline did to this paper. "
        "Search covers artifact paths, full small files, targeted shared-log excerpts, "
        "and the complete selected-DB snapshot.</div>"
        "<input class='search' id='trail-s' placeholder='search this paper process: error, model, variant, count, trust, source…' "
        "oninput='trailFilter()'>"
        "<div class='stage-nav'>"
        + "".join(
            f"<a href='#stage-{stage}'>Stage {stage}: {esc(label)}</a>"
            for stage, label in PROCESS_STAGES.items()
        )
        + "</div>"
        f"{''.join(sections)}</div>"
    )
    return _page(f"{gene} · PMID {pmid} · process", body, TRAIL_JS)


def find_latest_db(gene: str) -> Optional[Path]:
    """Return the newest pipeline DB, excluding review-publication copies.

    ``review_staging_test`` databases are intentionally refreshed when a review
    packet is staged, so their mtimes are often newer than the underlying run.
    Treating those copies as the latest run makes the dashboard silently score
    and summarize the wrong artifact.
    """
    cands = []
    for root in ("results", "validation_runs"):
        base = REPO / root
        if base.exists():
            cands += list(base.rglob(f"{gene}.db"))
    cands = [p for p in cands if "review_staging_test" not in p.parts]
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
    meta = paper["meta"]
    paper_root = corpus_dir / gene / pmid
    ft, source_basis = _resolve_extraction_source(meta, corpus_dir, gene, pmid)

    link_html = [f"<a href='{pubmed_url(pmid)}' target='_blank'>PubMed ↗</a>"]
    if doi:
        link_html.append(f"<a href='{doi_url(doi)}' target='_blank'>DOI ↗</a>")
    if pmcid:
        link_html.append(f"<a href='{pmc_url(pmcid)}' target='_blank'>PMC ↗</a>")
    link_html.append(f"<a href='{esc(pmid)}_process.html'>full process trail</a>")
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

    meta_bits = []
    for k in (
        "model_used",
        "source_type",
        "source_file",
        "extraction_confidence",
        "total_variants_found",
    ):
        if meta.get(k) not in (None, ""):
            value = _display_source_path(ft) if k == "source_file" else meta[k]
            meta_bits.append(f"<span class='tag'>{esc(k)}: {esc(value)}</span>")
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
        facts_html = _render_fact_provenance(
            db.get("facts_by_pmid_variant", {}).get((pmid, v["variant"]), [])
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
            f"{cohort}{quotes}{notes}{cprov}{facts_html}</div>"
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
            facts_html = _render_fact_provenance(
                db.get("facts_by_pmid_individual", {}).get(
                    (pmid, p["variant"], p["id"]), []
                ),
                limit=3,
            )
            left.append(
                f"<div class='rec'><span class='tag {cls}'>{esc(st or '?')}</span> {chars}{ev}{facts_html}</div>"
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
        raw = _sanitize_local_paths(raw)
        src_html = md_to_html(raw)
        bp_warn = boilerplate_warning(raw)
        src_note = (
            "<div class='mut'>Source shown: "
            f"<code>{esc(_display_source_path(ft))}</code> "
            f"<span class='tag'>{esc(source_basis)}</span></div>"
        )
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


# ---------------------------------------------------------------------------
# Trust / final-check / worklist / delta surfaces (the "progress & what to
# change" views). All read-only, all degrade to nothing on a bare DB, and all
# work WITHOUT a gold standard — the real (no-gold) workload's only signal.
# ---------------------------------------------------------------------------
def _paper_links(pmids: list, page_pmids: set, limit: int = 8) -> str:
    shown = pmids[:limit]
    bits = []
    for pm in shown:
        if pm in page_pmids:
            bits.append(f"<a href='papers/{esc(pm)}.html'>{esc(pm)}</a>")
        else:
            # esc the URL too: worklist PMIDs come from the corpus/trust/final-check
            # tables unfiltered, so a stray non-numeric id can't break the attribute.
            bits.append(
                f"<a href='{esc(pubmed_url(pm))}' target='_blank'>{esc(pm)}↗</a>"
            )
    more = len(pmids) - len(shown)
    tail = f" <span class='mut'>+{more}</span>" if more > 0 else ""
    return " ".join(bits) + tail


def _health_signals(
    s: dict,
    db: Optional[dict],
    zero_var_extracted: int,
    usable_unextracted: int,
    usable_extracted: int,
) -> dict:
    """The component health rates (each None when there's no basis to judge).

    ``zero_rate`` is the genuine extraction-miss rate — extracted-but-empty over
    the papers that had USABLE full text and were extracted (not over all DB
    papers, which would fold in abstract-only rows that can never be in the
    numerator, nor over the not-yet-extracted backlog, which is a throughput
    gap surfaced separately so a partial run doesn't read as broken quality)."""
    trust = (db or {}).get("trust", {}) or {}
    fcc = ((db or {}).get("final_check", {}) or {}).get("counts", {}) or {}
    # trust is only "gated" if the gate actually tiered facts into trusted /
    # quarantine. A present trust_tier COLUMN with only untiered (NULL) rows means
    # the gate never ran (or failed) — that must not read as a healthy 0%.
    trust_gated = (
        bool(trust.get("tiered"))
        and ((trust.get("trusted", 0) or 0) + (trust.get("quarantine", 0) or 0)) > 0
    )
    # final-check flag rate is over REVIEWABLE papers (ok + flag), not the raw
    # row count — "skipped"/"error" rows would otherwise dilute the rate and make
    # the gene read healthier than it is.
    fc_reviewable = fcc.get("ok", 0) + fcc.get("flag", 0)
    return {
        "zero_var": zero_var_extracted,
        "zero_rate": (zero_var_extracted / usable_extracted)
        if usable_extracted
        else None,
        "unextracted": usable_unextracted,
        "tiered": bool(trust.get("tiered")),
        "trust_gated": trust_gated,
        "trusted": trust.get("trusted", 0),
        "quarantine": trust.get("quarantine", 0),
        "quar_rate": trust.get("quarantine_rate") if trust_gated else None,
        "quar_reasons": list((trust.get("by_reason") or {}).items())[:3],
        "fc_total": fcc.get("total", 0),
        "fc_reviewable": fc_reviewable,
        "fc_ok": fcc.get("ok", 0),
        "fc_flag": fcc.get("flag", 0),
        "fc_flag_rate": (fcc.get("flag", 0) / fc_reviewable) if fc_reviewable else None,
        "fc_flagged_facts": fcc.get("flagged_facts", 0),
        "fc_missing_carriers": fcc.get("missing_carriers", 0),
    }


def _health_band(sig: dict) -> tuple[str, str]:
    """Coarse gold-free confidence read from the worst component rate.
    Returns (css_class, label). 'unknown' when nothing has been gated yet."""
    rates = [
        r
        for r in (sig["zero_rate"], sig["quar_rate"], sig["fc_flag_rate"])
        if r is not None
    ]
    if not rates:
        return "", "not yet gated"
    worst = max(rates)
    if worst >= 0.35:
        return "bad", "needs attention"
    if worst >= 0.15:
        return "warn", "watch"
    return "ok", "healthy"


def _health_card(
    s: dict,
    db: Optional[dict],
    zero_var_extracted: int,
    usable_unextracted: int,
    usable_extracted: int,
) -> str:
    """Gold-free 'is this gene going well?' card: extraction-funnel health +
    trust-tier split + final-check verdicts, with an overall confidence chip."""
    if not db:
        return ""
    sig = _health_signals(
        s, db, zero_var_extracted, usable_unextracted, usable_extracted
    )
    band_cls, band_label = _health_band(sig)

    def rate(v: Optional[float]) -> str:
        return f"{v * 100:.0f}%" if v is not None else "—"

    # trust line
    if sig["trust_gated"]:
        reasons = (
            " · ".join(f"{esc(k)} ({v})" for k, v in sig["quar_reasons"])
            if sig["quar_reasons"]
            else "—"
        )
        trust_html = (
            f"<div class='kv'><span>trusted / quarantined counts</span>"
            f"<b>{sig['trusted']} / {sig['quarantine']} ({rate(sig['quar_rate'])})</b></div>"
            f"<div class='mut' style='font-size:11px'>top quarantine reasons: {reasons}</div>"
        )
    elif sig["tiered"]:
        trust_html = (
            "<div class='kv'><span>trust gate</span>"
            "<b class='mut'>ran, 0 facts gated</b></div>"
        )
    else:
        trust_html = (
            "<div class='kv'><span>trust gate</span>"
            "<b class='mut'>not applied</b></div>"
        )
    # final-check line — flag rate is over reviewable (ok+flag) papers; any
    # skipped/error rows are surfaced separately, not folded into the rate.
    non_reviewable = sig["fc_total"] - sig["fc_reviewable"]
    if sig["fc_reviewable"]:
        skipped_note = f" · {non_reviewable} skipped/errored" if non_reviewable else ""
        fc_html = (
            f"<div class='kv'><span>final-check papers ok / flagged</span>"
            f"<b>{sig['fc_ok']} / {sig['fc_flag']} ({rate(sig['fc_flag_rate'])})</b></div>"
            f"<div class='mut' style='font-size:11px'>{sig['fc_flagged_facts']} flagged facts · "
            f"{sig['fc_missing_carriers']} carriers the checker says we missed{skipped_note}</div>"
        )
    elif sig["fc_total"]:
        fc_html = (
            "<div class='kv'><span>per-paper final check</span>"
            f"<b class='mut'>{sig['fc_total']} rows, none reviewable</b></div>"
        )
    else:
        fc_html = (
            "<div class='kv'><span>per-paper final check</span>"
            "<b class='mut'>not run</b></div>"
        )
    zero_flag = (
        " <span class='flag'>⚠</span>" if (sig["zero_rate"] or 0) >= 0.15 else ""
    )
    unextracted_html = (
        f"<div class='kv'><span>usable, not yet extracted</span>"
        f"<b class='mut'>{sig['unextracted']}</b></div>"
        if sig["unextracted"]
        else ""
    )
    return (
        "<div class='card'><h2>Run health "
        f"<span class='tag {band_cls}'>{esc(band_label)}</span></h2>"
        f"<div class='kv'><span>usable full text → extracted → with variants</span>"
        f"<b>{s.get('usable', 0)} → {s.get('extracted', 0)} → {s.get('with_variants', 0)}</b></div>"
        f"{unextracted_html}"
        f"<div class='kv'><span>extracted but 0 variants{zero_flag}</span>"
        f"<b>{sig['zero_var']} ({rate(sig['zero_rate'])})</b></div>"
        f"{trust_html}{fc_html}"
        "<div class='mut' style='font-size:11px;margin-top:4px'>gold-free confidence read: "
        "extraction-miss rate + trust-gate quarantine + per-paper sniff test. "
        "Green/amber/red = worst of the three rates (the not-yet-extracted backlog "
        "is a throughput lever, not a quality signal).</div></div>"
    )


def _worklist_card(
    gene: str,
    s: dict,
    db: Optional[dict],
    stub_pmids: list,
    zero_var_extracted: list,
    usable_unextracted: list,
    page_pmids: set,
    recall: Optional[dict],
    per_pmid: Optional[dict],
) -> str:
    """Ranked 'what to change next' list, built only from data already loaded
    (self-contained, no network, works with no gold). Each lever links to the
    affected papers and names the command to run. The richer gold-based
    disagreement worklist lives in scripts/recall_audit/build_acquisition_worklist.py."""
    import math

    trust = (db or {}).get("trust", {}) or {}
    quarantined = (db or {}).get("quarantined", []) or []
    # true quarantine total (the quarantined LIST is capped at 60 for display)
    quar_total = (
        trust.get("quarantine", 0)
        if (trust.get("tiered") and trust.get("total"))
        else 0
    )
    fc_by = ((db or {}).get("final_check", {}) or {}).get("by_pmid", {}) or {}
    flagged_pmids = [p for p, d in fc_by.items() if d.get("verdict") == "flag"]

    levers: list[dict] = []
    # gap-to-90 (gold genes only) — the headline recall lever
    if recall:
        m, g, _ = recall.get("unique_variants", (None, None, None))
        if g:
            gap = max(0, math.ceil(g * 0.9) - (m or 0))
            if gap > 0:
                worst = sorted(
                    ((p, v) for p, v in (per_pmid or {}).items() if v.get("missing")),
                    key=lambda kv: kv[1]["missing"],
                    reverse=True,
                )
                levers.append(
                    {
                        "n": gap,
                        "label": f"Close the gap to 90% unique-variant recall (+{gap})",
                        "why": "gold variants still missing from the DB — the papers "
                        "with the most missing rows are the highest-yield fixes",
                        "papers": [p for p, _ in worst],
                        "cmd": "scripts/refresh_recall.py / scripts/recall_recovery/run_all_layers.py",
                    }
                )
    # acquisition — stubs (no usable full text)
    if stub_pmids:
        levers.append(
            {
                "n": len(stub_pmids),
                "label": "Acquire missing full text (stub / paywalled)",
                "why": "no usable full text on disk — a publisher key or proxy "
                "fetch is the lever, not prompt tuning",
                "papers": stub_pmids,
                "cmd": f"gvf gvf-run {gene} …  (or scripts/fetch_paywalled.py)",
            }
        )
    # throughput — usable text on disk but never extracted into this DB
    if usable_unextracted:
        levers.append(
            {
                "n": len(usable_unextracted),
                "label": "Extract usable papers not yet in the DB",
                "why": "full text is on disk but these papers were never extracted "
                "(a partial / in-progress run) — the lever is running extraction, "
                "not fetching or prompt tuning",
                "papers": usable_unextracted,
                "cmd": f"gvf gvf-run {gene} …  (corpus cache skips what's done)",
            }
        )
    # extraction quality — extracted but produced 0 variants
    if zero_var_extracted:
        levers.append(
            {
                "n": len(zero_var_extracted),
                "label": "Investigate extraction misses (extracted, 0 variants)",
                "why": "extraction ran on usable text but produced nothing — a "
                "page-shell scrape or a parser/notation miss; open the source to see",
                "papers": zero_var_extracted,
                "cmd": "open the adjudication page → check the rendered source",
            }
        )
    # trust — quarantined counts
    if quar_total or quarantined:
        # dedup so a variant-heavy paper doesn't consume every visible link
        quar_pmids = list(
            dict.fromkeys(str(q.get("pmid")) for q in quarantined if q.get("pmid"))
        )
        levers.append(
            {
                "n": quar_total or len(quar_pmids),
                "label": "Review quarantined counts (trust gate)",
                "why": "counts the trust gate soft-quarantined as likely wrong "
                "(arithmetic / population / outlier) — informational, flagged for "
                "review, not yet auto-excluded from scoring",
                "papers": quar_pmids,
                "cmd": "python scripts/trust_report.py --db <db> --list 20",
            }
        )
    # review — final-check flagged papers
    if flagged_pmids:
        levers.append(
            {
                "n": len(flagged_pmids),
                "label": "Adjudicate final-check flags",
                "why": "papers the per-paper sniff test flagged as likely wrong "
                "(count vs provenance mismatch, missing quote)",
                "papers": flagged_pmids,
                "cmd": "open the adjudication page → verify against source",
            }
        )
    if not levers:
        return (
            "<div class='card'><h2>What to change next</h2>"
            "<div class='mut'>No open levers detected — full text acquired, variants "
            "extracted, no quarantine/flags. (Trust gate + final check may not have "
            "run; those add signal.)</div></div>"
        )
    levers.sort(key=lambda x: x["n"], reverse=True)
    rows = []
    for lev in levers:
        rows.append(
            f"<div class='rec'><b>{esc(lev['label'])}</b> "
            f"<span class='tag warn'>{lev['n']}</span>"
            f"<div class='mut' style='font-size:12px'>{esc(lev['why'])}</div>"
            f"<div style='margin-top:4px'>{_paper_links(lev['papers'], page_pmids)}</div>"
            f"<div class='mut' style='font-size:11px;margin-top:3px'>▶ <code>{esc(lev['cmd'])}</code></div>"
            "</div>"
        )
    return (
        "<div class='card' style='grid-column:1/-1'><h2>What to change next</h2>"
        "<div class='mut' style='font-size:12px;margin-bottom:6px'>ranked by size; "
        "each lever links to the affected papers and names the command.</div>"
        f"{''.join(rows)}</div>"
    )


# --- per-run deltas ("see progress" over time) -----------------------------
def _gene_snapshot(s: dict, db: Optional[dict], recall: Optional[dict]) -> dict:
    """The per-run numbers we diff to show 'what changed since last run'."""
    trust = (db or {}).get("trust", {}) or {}
    fcc = ((db or {}).get("final_check", {}) or {}).get("counts", {}) or {}
    uv = (recall or {}).get("unique_variants", (None, None, None))
    # None means "not measured this run" (so a delta is skipped); use explicit 0
    # defaults for the "measured, happens to be zero" case (Counter.get returns
    # None for an absent key even when total>0 — a subtle all-ok trap).
    # "gated"/"run" require actual signal: facts tiered into trusted/quarantine,
    # and final-check papers that were reviewable (ok/flag) — not merely a present
    # column or a table of all-skipped/errored rows.
    trust_gated = (
        bool(trust.get("tiered"))
        and ((trust.get("trusted", 0) or 0) + (trust.get("quarantine", 0) or 0)) > 0
    )
    fc_run = (fcc.get("ok", 0) + fcc.get("flag", 0)) > 0
    return {
        "usable": s.get("usable"),
        "stub": s.get("stub"),
        "extracted": s.get("extracted"),
        "with_variants": s.get("with_variants"),
        "unique": s.get("unique"),
        "variant_rows": s.get("variant_rows"),
        "quarantine": trust.get("quarantine", 0) if trust_gated else None,
        "flagged_papers": fcc.get("flag", 0) if fc_run else None,
        "missing_carriers": fcc.get("missing_carriers", 0) if fc_run else None,
        "recall_uv_pct": round((uv[2] or 0) * 100, 1) if uv[2] is not None else None,
    }


_SNAPSHOT_FIELDS = (
    ("unique", "unique variants", False),
    ("variant_rows", "variant rows", False),
    ("with_variants", "papers w/ variants", False),
    ("extracted", "extracted papers", False),
    ("usable", "usable full text", False),
    ("stub", "stubs to fetch", True),
    ("quarantine", "quarantined counts", True),
    ("flagged_papers", "flagged papers", True),
    ("missing_carriers", "checker-missed carriers", True),
    ("recall_uv_pct", "uniqV recall %", False),
)


def _history_path(hist_dir: Path, gene: str) -> Path:
    return hist_dir / f"{gene}.jsonl"


def _load_prev_snapshot(hist_dir: Path, gene: str) -> Optional[dict]:
    """The most recent well-formed snapshot dict. Scans backward so a single
    corrupt/partial last line (or a valid-but-non-object line like ``[]``) does
    not discard the whole history or crash the delta consumers."""
    p = _history_path(hist_dir, gene)
    if not p.exists():
        return None
    try:
        # match _append_snapshot's utf-8 write; the platform default (cp1252 on
        # Windows) would raise UnicodeDecodeError (a ValueError, not OSError).
        lines = [ln for ln in p.read_text(encoding="utf-8").splitlines() if ln.strip()]
    except OSError:
        return None
    for ln in reversed(lines):
        try:
            obj = json.loads(ln)
        except json.JSONDecodeError:
            continue
        if isinstance(obj, dict):
            return obj
    return None


def _snapshot_changed(prev: Optional[dict], cur: dict) -> bool:
    if prev is None:
        return True
    return any(prev.get(k) != cur.get(k) for k, _, _ in _SNAPSHOT_FIELDS)


def _append_snapshot(hist_dir: Path, gene: str, snap: dict, generated: str) -> None:
    try:
        hist_dir.mkdir(parents=True, exist_ok=True)
        rec = dict(snap)
        rec["generated"] = generated
        with _history_path(hist_dir, gene).open("a", encoding="utf-8") as f:
            f.write(json.dumps(rec) + "\n")
    except OSError:  # pragma: no cover - history is best-effort
        pass


def _delta_html(prev: Optional[dict], cur: dict) -> str:
    """Compact 'since last run' line. Empty string if there's nothing to say."""
    if prev is None:
        return (
            "<div class='note'>First tracked run — deltas will appear on the next "
            "dashboard build after the numbers change.</div>"
        )
    bits = []
    for key, label, lower_better in _SNAPSHOT_FIELDS:
        a, b = prev.get(key), cur.get(key)
        # only diff like-for-like measured numbers (skip a side that's None =
        # not-measured-this-run, or a corrupt/non-numeric history value)
        if not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
            continue
        if a == b:
            continue
        d = b - a
        good = (d < 0) if lower_better else (d > 0)
        cls = "ok" if good else "bad"
        d_str = f"{d:+.0f}" if float(d).is_integer() else f"{d:+.1f}"
        bits.append(f"<span class='tag {cls}'>{esc(label)} {d_str}</span>")
    prev_when = esc(prev.get("generated") or "last run")
    if not bits:
        return (
            "<div class='mut' style='font-size:12px'>No tracked metric changed "
            f"since {prev_when}.</div>"
        )
    return f"<div class='note'>Since {prev_when}: {' '.join(bits)}</div>"


def _delta_compact(prev: Optional[dict], cur: dict) -> str:
    """One chip for the overview card: unique-variant change since last run."""
    if not prev:
        return ""
    a, b = prev.get("unique"), cur.get("unique")
    # only diff like-for-like measured numbers (a corrupt history value could be
    # a str/None; subtracting would raise or fabricate a delta)
    if not isinstance(a, (int, float)) or not isinstance(b, (int, float)) or a == b:
        return ""
    d = b - a
    cls = "ok" if d > 0 else "bad"
    return f"<span class='tag {cls}'>uniqV {d:+.0f} since last run</span>"


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
    delta_html: str = "",
) -> tuple[str, dict, list[str]]:
    s = gene_stats(gene, corpus_rows, db)
    per_pmid = per_pmid or {}
    pen_by_pmid = (db or {}).get("pen_by_pmid", {})
    prov = (db or {}).get("prov", Counter())
    method = (db or {}).get("source_method", Counter())
    fc_by = (db or {}).get("final_check", {}).get("by_pmid", {}) if db else {}
    has_fc = bool((db or {}).get("final_check", {}).get("counts", {}).get("total"))

    # provenance-completeness audit (the "test experiment" surface)
    audit = []
    if prov.get("total"):
        audit = [
            ("source_location populated", prov.get("with_location", 0), prov["total"]),
            ("additional_notes populated", prov.get("with_notes", 0), prov["total"]),
            ("key_quotes populated", prov.get("with_quotes", 0), prov["total"]),
        ]
        if prov.get("with_fact_provenance"):
            audit.append(
                (
                    "fact provenance exact pointers",
                    prov.get("with_exact_fact_pointer", 0),
                    prov.get("with_fact_provenance", 0),
                )
            )

    # which papers still need work (no usable full text, or no variants)
    stub_pmids = [
        p for p, r in corpus_rows.items() if r.get("full_text_status") != "ok"
    ]

    # paper table. PubMed IDs are numeric; drop any malformed token (e.g. a
    # stray "N/A" extraction row) so it neither links nor tries to write a
    # papers/<pmid>.html path — a "/" in the id would otherwise abort the build.
    rows = []
    paper_pmids = sorted(
        {
            p
            for p in (set(corpus_rows) | set((db or {}).get("papers", {})))
            if p.isdigit()
        },
        key=int,
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
    db_papers = set((db or {}).get("papers", {}))
    # Split usable-but-empty into two very different problems: papers extraction
    # never touched (throughput) vs. papers it processed and got nothing from
    # (a genuine miss). Conflating them makes a partial run look broken.
    zero_var_extracted: list[str] = []  # in DB, extracted, but 0 variants
    usable_unextracted: list[str] = []  # usable text on disk, never extracted
    n_usable_extracted = 0  # usable full text AND in the DB (zero_rate denominator)
    for pmid in paper_pmids:
        cr = corpus_rows.get(pmid, {})
        status = cr.get("full_text_status", "—")
        by = (db or {}).get("by_pmid", {}).get(pmid, {})
        nvar = len(by.get("variants", [])) if db else 0
        npat = len(by.get("patients", [])) if db else 0
        in_db = db is not None and pmid in db_papers
        if status == "ok":
            if not in_db:
                # includes the no-DB case (db is None): usable text on disk but
                # nothing extracted, so surface an "extract these" lever rather
                # than a misleading "no open levers detected".
                usable_unextracted.append(pmid)
            else:
                n_usable_extracted += 1
                if nvar == 0:
                    zero_var_extracted.append(pmid)
        title = (db or {}).get("papers", {}).get(pmid, {}).get("title") or ""
        badge = "ok" if status == "ok" else "bad" if status == "stub" else "warn"
        has_page = pmid in cap_set
        pmid_cell = (
            f"<a href='papers/{esc(pmid)}.html'>{esc(pmid)}</a>"
            f"<div><a class='mut' href='papers/{esc(pmid)}_process.html'>process trail</a></div>"
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
        fc_cell = ""
        if has_fc:
            fc = fc_by.get(pmid)
            if not fc:
                fc_cell = "<td class='mut' data-v='2'>—</td>"
            elif fc.get("verdict") == "flag":
                miss = fc.get("n_missing") or 0
                extra = f" −{miss}" if miss else ""
                title_txt = (fc.get("summary") or "flagged by final check")[:180]
                fc_cell = (
                    f"<td data-v='0'><span class='tag bad' title=\"{esc(title_txt)}\">"
                    f"flag{esc(extra)}</span></td>"
                )
            elif fc.get("verdict") == "ok":
                fc_cell = "<td data-v='1'><span class='tag ok'>ok</span></td>"
            else:
                # skipped / error / unknown — NOT a green pass; show it muted
                v = fc.get("verdict") or "?"
                fc_cell = f"<td class='mut' data-v='2'>{esc(v)}</td>"
        rows.append(
            f"<tr><td data-v='{esc(pmid)}'>{pmid_cell}</td><td>{esc(title[:80])}</td>"
            f"<td><span class='tag {badge}'>{esc(status)}</span>"
            + (
                " <span class='flag' title='Extracted but 0 variants — possible boilerplate scrape or extraction miss; check the source'>⚠</span>"
                if status == "ok" and nvar == 0 and in_db
                else " <span class='mut' title='Usable full text on disk but not yet extracted into this DB'>·</span>"
                if status == "ok" and not in_db and db is not None
                else ""
            )
            + "</td>"
            f"<td data-v='{nvar}'>{nvar}</td>{gold_cells}"
            f"<td data-v='{carriers}'>{carriers}</td>"
            f"<td data-v='{aff}'>{aff}</td><td data-v='{unaff}'>{unaff}</td>"
            f"<td data-v='{npat}' class='mut'>{npat}</td>{fc_cell}</tr>"
        )

    # Build the header (gold genes get gold / matched / missing-vs-gold columns).
    # Carriers/Aff/Unaff are cohort counts from penetrance_data (the real patient
    # numbers); "records" is the sparser per-person individual_records count.
    base_cols = ["PMID", "Title", "Source", "Variants"]
    gold_cols = ["Gold", "Matched", "Δ vs gold"] if has_gold else []
    tail_cols = ["Carriers", "Aff", "Unaff", "records"] + (["Check"] if has_fc else [])
    all_cols = base_cols + gold_cols + tail_cols
    thead_html = "".join(
        f"<th onclick=\"srt('pt',{i})\">{esc(c)}</th>" for i, c in enumerate(all_cols)
    )
    db_source = _display_db_source(db["db_path"]) if db and db.get("db_path") else ""
    db_source_html = (
        f"<div class='sub'>database: <code>{esc(db_source)}</code></div>"
        if db_source
        else ""
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
        f"<a href='../index.html'>← overview</a> · <a href='variants.html'>variants ({s['unique']})</a> · "
        "<a href='process.html'>process explorer</a>"
        f"</div>{db_source_html}</header>"
        f"<div class='wrap'>{gap_note}{delta_html}"
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
        f"{_health_card(s, db, len(zero_var_extracted), len(usable_unextracted), n_usable_extracted)}"
        f"{_gold_card(recall, mae)}"
        f"{_gold_source_card(set(per_pmid), corpus_rows)}"
        f"<div class='card'><h2>Provenance completeness</h2>{audit_html or '<span class=mut>no DB</span>'}"
        f"<div style='margin-top:8px'>{method_html}</div></div>"
        f"<div class='card'><h2>What's left</h2>"
        f"<div class='kv'><span>papers without usable full text</span><b>{len(stub_pmids)}</b></div>"
        f"<div class='mut' style='font-size:12px;margin-top:6px'>re-run <code>gvf gvf-run {esc(gene)} …</code> "
        f"(corpus cache skips what's done; a new publisher key re-fetches stubs)</div></div>"
        f"</div>"
        f"{_worklist_card(gene, s, db, stub_pmids, zero_var_extracted, usable_unextracted, cap_set, recall, per_pmid)}"
        f"<h2 style='margin-top:20px'>Papers</h2>"
        f"<div class='mut' style='font-size:12px;margin-bottom:6px'>"
        f"Variants = extracted; Gold/Matched/Δ vs the gold standard (where available); "
        f"Carriers/Aff/Unaff = cohort counts from penetrance_data; records = per-person rows. "
        f"<span class='flag'>⚠</span> = usable full text on disk but 0 variants extracted "
        f"(possible boilerplate/page-shell scrape or an extraction miss — open it to check)."
        + (
            " Check = per-paper final-check verdict (flag −N = checker thinks N carriers were missed)."
            if has_fc
            else ""
        )
        + "</div>"
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
        # trust + final-check + since-last-run badges (the no-gold progress read)
        trust = s.get("trust") or {}
        fcc = s.get("final_check_counts") or {}
        badges = []
        # only badge trust when the gate actually tiered facts into trusted/
        # quarantine; a present column with only untiered rows is not "healthy".
        if (
            trust.get("tiered")
            and ((trust.get("trusted", 0) or 0) + (trust.get("quarantine", 0) or 0)) > 0
        ):
            qr = trust.get("quarantine_rate") or 0
            tcls = "bad" if qr >= 0.35 else "warn" if qr >= 0.15 else "ok"
            badges.append(
                f"<span class='tag {tcls}'>trust {trust.get('trusted', 0)}✓ / "
                f"{trust.get('quarantine', 0)} quar</span>"
            )
        # flag rate over reviewable (ok+flag), not raw row count
        reviewable = fcc.get("ok", 0) + fcc.get("flag", 0)
        if reviewable:
            fcls = "bad" if fcc.get("flag", 0) > 0 else "ok"
            badges.append(
                f"<span class='tag {fcls}'>check {fcc.get('flag', 0)} flagged</span>"
            )
        badges_html = s.get("delta_compact") or ""
        badges_html += "".join(badges)
        badges_html = (
            f"<div style='margin-top:6px'>{badges_html}</div>" if badges_html else ""
        )
        cards.append(
            f"<div class='card'><h2><a href='{esc(gene)}/index.html'>{esc(gene)}</a></h2>"
            f"<div class='kv'><span>papers in corpus</span><b>{s['in_corpus']}</b></div>"
            f"{_barbar(s['usable'], s['in_corpus'])}"
            f"<div class='kv'><span>usable full text</span><b>{s['usable']} ({pct(s['usable'], s['in_corpus']):.0f}%)</b></div>"
            f"<div class='kv'><span>extracted</span><b>{s['extracted']}</b></div>"
            f"<div class='kv'><span>papers w/ variants</span><b>{s['with_variants']}</b></div>"
            f"<div class='kv'><span>unique variants</span><b>{s['unique']}</b></div>{gold}{badges_html}"
            f"<div class='mut' style='font-size:11px;margin-top:6px'>DB: <code>{esc(s.get('db_source') or 'none')}</code></div>"
            f"<div style='margin-top:8px'><a href='{esc(gene)}/index.html'>status</a> · "
            f"<a href='{esc(gene)}/variants.html'>variants</a> · "
            f"<a href='{esc(gene)}/process.html'>process ({s.get('artifact_count', 0)})</a>"
            "</div></div>"
        )
    tot = {
        k: sum(s[k] for s in summaries.values())
        for k in ("in_corpus", "usable", "extracted", "variant_rows")
    }
    score_note = (
        "gold scores computed during this build from the selected DBs"
        if any(s.get("recall") for s in summaries.values())
        else "gold scores not computed in this build"
    )
    body = (
        f"<header><h1>GeneVariantFetcher — status & provenance dashboard</h1>"
        f"<div class='sub'>generated {esc(generated)} · {esc(score_note)} · policy baseline: docs/RECALL_STATUS.md · "
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
    # dedup while preserving order — a duplicated --gene would otherwise diff a
    # run against the snapshot it just appended and report a spurious "no change".
    target_genes = list(dict.fromkeys(target_genes))
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
        # per-run delta: diff this run's snapshot against the last one recorded in
        # the (gitignored) corpus history, so the gene page can show what moved.
        hist_dir = corpus_dir / "dashboard_history"
        snap = _gene_snapshot(gene_stats(gene, corpus_rows, db), db, recall)
        prev = _load_prev_snapshot(hist_dir, gene)
        # carry the last SCORED recall forward through unscored (--no-score) builds
        # so an interleaved fast pass doesn't wipe the recall baseline.
        if (
            snap.get("recall_uv_pct") is None
            and prev
            and prev.get("recall_uv_pct") is not None
        ):
            snap["recall_uv_pct"] = prev["recall_uv_pct"]
        delta_html = _delta_html(prev, snap)
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
            delta_html=delta_html,
        )
        (gene_dir / "index.html").write_text(page, encoding="utf-8")
        artifacts = discover_process_artifacts(gene, db_path, corpus_dir)
        (gene_dir / "process.html").write_text(
            render_process_page(gene, artifacts, gene_dir, db_path), encoding="utf-8"
        )
        s["artifact_count"] = len(artifacts)
        stats["artifact_files"] = stats.get("artifact_files", 0) + len(artifacts)
        if _snapshot_changed(prev, snap):
            _append_snapshot(hist_dir, gene, snap, generated)
        s["recall"] = recall
        s["db_source"] = _display_db_source(db_path) if db_path else ""
        s["gold_cov"] = (
            gold_coverage_stats(set(per_pmid), corpus_rows) if per_pmid else None
        )
        # stash the no-gold progress signals for the overview cards
        s["trust"] = (db or {}).get("trust", {"tiered": False})
        s["final_check_counts"] = dict(
            (db or {}).get("final_check", {}).get("counts", {})
        )
        s["delta_compact"] = _delta_compact(prev, snap)
        summaries[gene] = s
        stats["genes"] += 1
        if db:
            papers_dir = gene_dir / "papers"
            papers_dir.mkdir(parents=True, exist_ok=True)
            paper_artifacts = build_paper_process_index(
                gene, db_path, corpus_dir, set(paper_pmids), artifacts
            )
            for pmid in paper_pmids:
                (papers_dir / f"{pmid}.html").write_text(
                    render_paper_page(gene, pmid, db, corpus_dir, papers_dir),
                    encoding="utf-8",
                )
                (papers_dir / f"{pmid}_process.html").write_text(
                    render_paper_process_page(
                        gene,
                        pmid,
                        db,
                        paper_artifacts.get(pmid, []),
                        papers_dir,
                        score=per_pmid.get(pmid) if per_pmid else None,
                    ),
                    encoding="utf-8",
                )
                stats["paper_pages"] += 1
                stats["paper_process_pages"] = stats.get("paper_process_pages", 0) + 1
            (gene_dir / "variants.html").write_text(
                render_variants_page(gene, db, set(paper_pmids)), encoding="utf-8"
            )
    (out_dir / "index.html").write_text(
        render_index(summaries, generated), encoding="utf-8"
    )
    stats["summaries"] = summaries
    return stats
