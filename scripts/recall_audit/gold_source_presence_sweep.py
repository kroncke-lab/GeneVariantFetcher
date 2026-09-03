#!/usr/bin/env python3
"""Where does each gold variant live in the source we actually hold on disk?

The mixed-gold tranches score the *reading protocol*, but part of every miss is
an *acquisition* ceiling: a gold row whose variant string appears nowhere in the
article text, any supplement, or any PDF we hold cannot be produced by any
prompt, router, or parser. This sweep classifies every gold row of the mixed
inventory against everything on disk for its paper, blind to predictions, so
the ceiling can be reported as a bounded stratum instead of a guess.

It deliberately searches more than ``tier_source_reachability.py`` does. That
probe reads the article body, plain-text supplements, and PDFs. Most supplement
bytes in the corpus are binary — ``.docx``, ``.xlsx``, ``.doc``, ``.pptx``,
``.zip`` bundles — and the fold step does not always land them in
``FULL_CONTEXT``. Declaring a row "absent from acquired source" without opening
those files would overstate the ceiling, so each binary supplement is converted
through the production ``FormatConverter`` (cached outside the corpus) and
searched too.

Classes, in evaluation order
----------------------------
``present_in_body``
    The variant is in ``<pmid>_FULL_CONTEXT.md`` or ``<pmid>_CLEANED.md``.
``present_in_supplement_only``
    Absent from the body, present in a supplement file after conversion. A
    reader that only sees the body cannot find it; the fold or the source
    selector owns the miss.
``present_in_pdf_only``
    Absent from body and supplements, present in a top-level article PDF.
``source_absent``
    Nothing searchable is on disk for this paper: no usable body, no
    convertible supplement, no PDF.
``text_absent_stub_body``
    Not in any searchable text, and everything we hold for the paper is
    shorter than an article (abstract page, landing page, caption stubs with
    no table bodies). The full text never landed. Acquisition.
``text_absent_garbled_body``
    Not in any searchable text, and the body is dominated by unmapped PDF
    glyph codes (``(cid:NN)``) or has almost no letters. The bytes landed but
    are unreadable as text. Acquisition (re-render/OCR), not extraction.
``text_absent_figures_present``
    Not in any searchable text, but figure images are on disk. A PNG cannot
    be string-searched; the source *may* hold the row. Bound, not proof.
``text_absent_notation_inconclusive``
    Not in any searchable text and the notation is an indel/frameshift/
    structural form the probe cannot reliably find. Unknown, not absent.
``text_absent_substitution``
    A plain substitution absent from every searchable representation with no
    figure images on disk. This is the hard acquisition ceiling.

Integrity properties: classification reads gold rows and on-disk source only,
never a prediction or a score. Nothing under ``corpus/`` is written; converted
text is cached under ``--cache-dir``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from benchmarks.codex_paper_eval.run_eval import extract_pdf_text, load_gold  # noqa: E402
from scripts.recall_audit.common import unresolved_variant_supplement_refs  # noqa: E402
from scripts.recall_audit.source_stratified_metrics import (  # noqa: E402
    ABSTRACT_ONLY_MARKER,
    MIN_FULL_TEXT_CHARS,
)
from scripts.recall_audit.tier_source_reachability import (  # noqa: E402
    SEARCHABLE_SUFFIXES,
    FormIndex,
    _squash,
    is_substitution,
)

MIXED_GOLD = REPO_ROOT / "benchmarks" / "evaluation_tiers" / "mixed_gold"
DEFAULT_INVENTORY = MIXED_GOLD / "inventory.tsv"
DEFAULT_REGISTRY = MIXED_GOLD / "registry.json"
DEFAULT_GOLD_ROOT = MIXED_GOLD / "answer_key"

IMAGE_SUFFIXES = frozenset({".png", ".jpg", ".jpeg", ".tif", ".tiff", ".gif", ".bmp"})
MEDIA_SUFFIXES = frozenset({".mp4", ".avi", ".mov", ".mpg", ".mpeg", ".wmv", ".mp3"})
CONVERTIBLE_SUFFIXES = frozenset(
    {".docx", ".doc", ".xlsx", ".xls", ".pptx", ".ppt", ".pdf", ".zip"}
)
import re  # noqa: E402

VARIANT_TOKEN_RE = re.compile(
    r"\b(?:p\.)?[A-Z][a-z]{2}\d{1,4}(?:[A-Z][a-z]{2}|Ter|X|\*|fs)"
    r"|\b[ACDEFGHIKLMNPQRSTVWY]\d{1,4}(?:[ACDEFGHIKLMNPQRSTVWY]|X|\*)\b"
    r"|\bc\.\d+[+-]?\d*[ACGT]>[ACGT]"
)
MAX_CONVERT_BYTES = 80 * 1024 * 1024
MAX_SEARCH_BYTES = 96 * 1024 * 1024

CLASSES = (
    "present_in_body",
    "present_in_supplement_only",
    "present_in_pdf_only",
    "source_absent",
    "text_absent_stub_body",
    "text_absent_garbled_body",
    "text_absent_figures_present",
    "text_absent_notation_inconclusive",
    "text_absent_substitution",
)
# The two classes a recall denominator may exclude under the existing
# reachability policy: the source could not have supported the row. Everything
# else stays penalized because a weak probe must never manufacture credit.
HARD_CEILING = frozenset(
    {
        "source_absent",
        "text_absent_stub_body",
        "text_absent_garbled_body",
        "text_absent_substitution",
    }
)
# Below this many searchable characters (body + supplements + PDFs) the paper
# is an abstract/landing page, not an article. Real bodies in this corpus are
# tens of KB; abstract captures are 2-5 KB.
STUB_CHARS = 6000
# pdftotext emits ``(cid:NN)`` for glyphs with no Unicode mapping. One code
# per ~40 characters means the running text is unreadable.
GARBLED_CID_RATIO = 1 / 40
# Upper bound on the ceiling: also count the rows the probe cannot decide.
WIDE_CEILING = HARD_CEILING | frozenset(
    {"text_absent_figures_present", "text_absent_notation_inconclusive"}
)


@dataclass
class PaperIndex:
    """Every searchable representation on disk for one paper."""

    body: str = ""
    body_chars: int = 0
    body_state: str = ""
    supplement_text: str = ""
    supplement_hits: dict = field(default_factory=dict)  # squashed text -> file
    supplement_text_files: int = 0
    supplement_converted_files: int = 0
    supplement_failed_files: int = 0
    supplement_unsearchable_files: int = 0
    pdf_text: str = ""
    article_pdf_files: int = 0
    figure_files: int = 0
    unfetched_links: int = 0
    advertised_unresolved_refs: bool = False
    paper_dir_exists: bool = False
    supplement_file_index: list = field(default_factory=list)
    searchable_chars: int = 0
    body_cid_codes: int = 0
    body_letter_ratio: float = 0.0
    body_variant_tokens: int = 0
    body_garbled: bool = False


def _cache_key(path: Path) -> str:
    st = path.stat()
    return hashlib.sha256(
        f"{path.resolve()}|{st.st_size}|{int(st.st_mtime)}".encode()
    ).hexdigest()


def convert_supplement(path: Path, cache_dir: Path, converter) -> str | None:
    """Convert a binary supplement to text, caching outside the corpus."""
    suffix = path.suffix.lower()
    try:
        if path.stat().st_size > MAX_CONVERT_BYTES:
            return None
    except OSError:
        return None
    cache_dir.mkdir(parents=True, exist_ok=True)
    cached = cache_dir / f"{_cache_key(path)}.txt"
    if cached.is_file():
        return cached.read_text(errors="replace")
    text = ""
    try:
        if suffix == ".pdf":
            text = extract_pdf_text([str(path.resolve())], MAX_SEARCH_BYTES) or ""
        elif suffix == ".docx":
            text = converter.docx_to_markdown(path) or ""
        elif suffix == ".doc":
            text = converter.doc_to_markdown(path) or ""
        elif suffix in {".xlsx", ".xls"}:
            text = converter.excel_to_markdown(path) or ""
        elif suffix in {".pptx", ".ppt"}:
            md = getattr(converter, "markitdown", None)
            text = md.convert(str(path)).text_content if md is not None else ""
        elif suffix == ".zip":
            dest = cache_dir / "zip" / _cache_key(path)
            _, text = converter.extract_zip_supplement(path, dest)
            text = text or ""
    except Exception as exc:  # noqa: BLE001 - one bad file must not stop a sweep
        text = ""
        (cache_dir / f"{_cache_key(path)}.err").write_text(
            f"{path}\n{type(exc).__name__}: {exc}\n"
        )
        return None
    cached.write_text(text)
    return text


def build_index(
    gene: str, pmid: str, corpus: Path, cache_dir: Path, converter
) -> PaperIndex:
    paper = corpus / gene / pmid
    out = PaperIndex(paper_dir_exists=paper.is_dir())
    if not out.paper_dir_exists:
        out.body_state = "no source on disk"
        return out

    body_parts: list[str] = []
    states: list[str] = []
    for name in (f"{pmid}_FULL_CONTEXT.md", f"{pmid}_CLEANED.md"):
        candidate = paper / name
        if not candidate.is_file():
            continue
        try:
            raw = candidate.read_text(errors="replace")
        except OSError as exc:
            states.append(f"unreadable: {exc.strerror or exc}")
            continue
        if raw.lstrip().startswith(ABSTRACT_ONLY_MARKER):
            states.append("abstract-only fallback")
            continue
        if len(raw) < MIN_FULL_TEXT_CHARS:
            states.append("source too short to be a body")
            # Still searchable: a stub can carry a variant in its abstract.
        out.body_chars = max(out.body_chars, len(raw))
        body_parts.append(raw)
    raw_body = "\n".join(body_parts)
    out.body = _squash(raw_body)
    if raw_body:
        out.body_cid_codes = raw_body.count("(cid:")
        letters = sum(ch.isalpha() for ch in raw_body)
        out.body_letter_ratio = round(letters / max(1, len(raw_body)), 3)
        out.body_variant_tokens = len(VARIANT_TOKEN_RE.findall(raw_body))
        out.body_garbled = (
            out.body_cid_codes / max(1, len(raw_body)) >= GARBLED_CID_RATIO
            or out.body_letter_ratio < 0.35
        )
    if not body_parts:
        out.body_state = states[0] if states else "no article body on disk"
    elif states:
        out.body_state = states[0]
    if raw_body:
        try:
            out.advertised_unresolved_refs = bool(
                unresolved_variant_supplement_refs(raw_body, gene)
            )
        except Exception:  # noqa: BLE001
            out.advertised_unresolved_refs = False

    artifacts = paper / f"{pmid}_artifacts.json"
    if artifacts.is_file():
        try:
            meta = json.loads(artifacts.read_text(errors="replace"))
            out.unfetched_links = len(meta.get("supplement_links_unfetched") or [])
        except Exception:  # noqa: BLE001
            pass

    figures = 0
    for f in paper.rglob("*"):
        if f.is_file() and f.suffix.lower() in IMAGE_SUFFIXES:
            figures += 1
    out.figure_files = figures

    supplements = paper / f"{pmid}_supplements"
    supp_parts: list[str] = []
    budget = MAX_SEARCH_BYTES
    if supplements.is_dir():
        for f in sorted(supplements.rglob("*")):
            if not f.is_file() or budget <= 0:
                continue
            suffix = f.suffix.lower()
            if f.name.startswith("."):
                continue
            text: str | None = None
            if suffix in SEARCHABLE_SUFFIXES:
                try:
                    if f.stat().st_size <= budget:
                        text = f.read_text(errors="replace")
                        out.supplement_text_files += 1
                except OSError:
                    text = None
            elif suffix in CONVERTIBLE_SUFFIXES:
                text = convert_supplement(f, cache_dir, converter)
                if text is None:
                    out.supplement_failed_files += 1
                elif text:
                    out.supplement_converted_files += 1
            elif suffix in IMAGE_SUFFIXES or suffix in MEDIA_SUFFIXES:
                continue
            else:
                out.supplement_unsearchable_files += 1
            if text:
                squashed = _squash(text)
                if squashed:
                    supp_parts.append(squashed)
                    out.supplement_file_index.append(
                        (str(f.relative_to(paper)), squashed)
                    )
                    budget -= len(text)
    out.supplement_text = "\n".join(supp_parts)

    pdfs = sorted(f for f in paper.glob("*.pdf") if f.is_file() and f.parent == paper)
    out.article_pdf_files = len(pdfs)
    if pdfs:
        parts = []
        for pdf in pdfs:
            text = convert_supplement(pdf, cache_dir, converter)
            if text:
                parts.append(_squash(text))
        out.pdf_text = "\n".join(parts)
    out.searchable_chars = len(out.body) + len(out.supplement_text) + len(out.pdf_text)
    return out


_NONSENSE_RE = re.compile(r"^(P\.)?([A-Z]{1,3})(\d+)(X|\*|TER|STOP)$")


def expand_forms(forms: frozenset[str]) -> frozenset[str]:
    """Add the nonsense spellings the scorer already bridges (X / * / Ter / Stop).

    ``VariantNormalizer.get_all_forms`` yields ``p.Arg222*`` for ``R222X`` but
    papers write ``Arg222Ter`` and ``Arg222Stop``; the scorer's ``matches``
    accepts those, so a presence probe that cannot see them would report a
    reachable row as absent.
    """
    out = set(forms)
    for form in forms:
        m = _NONSENSE_RE.match(form)
        if not m:
            continue
        prefix, ref, pos, _ = m.groups()
        for tail in ("X", "*", "TER", "STOP"):
            out.add(f"{ref}{pos}{tail}")
            out.add(f"P.{ref}{pos}{tail}")
    return frozenset(out)


def classify(variant: str, forms: frozenset[str], index: PaperIndex) -> tuple[str, str]:
    forms = expand_forms(forms)
    for form in forms:
        if form in index.body:
            return "present_in_body", "body"
    for form in forms:
        for rel, squashed in index.supplement_file_index:
            if form in squashed:
                return "present_in_supplement_only", f"supplement:{rel}"
    for form in forms:
        if form in index.pdf_text:
            return "present_in_pdf_only", "article_pdf"
    searchable = bool(index.body or index.supplement_text or index.pdf_text)
    if not searchable:
        return "source_absent", index.body_state or "nothing searchable on disk"
    if index.body_garbled and not index.supplement_text and not index.pdf_text:
        return (
            "text_absent_garbled_body",
            f"body has {index.body_cid_codes} unmapped glyph codes, letter ratio "
            f"{index.body_letter_ratio}",
        )
    if index.searchable_chars < STUB_CHARS:
        return (
            "text_absent_stub_body",
            f"only {index.searchable_chars} searchable chars on disk",
        )
    if index.figure_files:
        return (
            "text_absent_figures_present",
            f"{index.figure_files} figure image(s) on disk",
        )
    if not is_substitution(variant):
        return "text_absent_notation_inconclusive", "non-substitution notation"
    return "text_absent_substitution", "substitution absent from all searchable text"


def sweep_paper(task: dict) -> dict:
    """Worker: classify every gold row for one gene-paper attempt."""
    gene, pmid = task["gene"], task["pmid"]
    corpus = Path(task["corpus"])
    cache_dir = Path(task["cache_dir"])
    gold_root = Path(task["gold_root"])
    from harvesting.format_converters import FormatConverter

    converter = FormatConverter()
    started = time.time()
    rows_out: list[dict] = []
    try:
        gold_rows = load_gold(gold_root, gene, pmid)
    except Exception as exc:  # noqa: BLE001
        return {
            "gene": gene,
            "pmid": pmid,
            "error": f"gold load failed: {type(exc).__name__}: {exc}",
            "rows": [],
        }
    index = build_index(gene, pmid, corpus, cache_dir, converter)
    forms = FormIndex()
    for gold_row in gold_rows:
        variant = (gold_row.get("variant") or "").strip()
        if not variant:
            continue
        klass, where = classify(variant, forms.get(gene, variant), index)
        rows_out.append(
            {
                "gene": gene,
                "pmid": pmid,
                "variant": variant,
                "class": klass,
                "where": where,
                "is_substitution": is_substitution(variant),
                "gold_carriers": gold_row.get("carriers"),
                "gold_affected": gold_row.get("affected"),
                "gold_unaffected": gold_row.get("unaffected"),
            }
        )
    paper_meta = asdict(index)
    paper_meta.pop("body", None)
    paper_meta.pop("supplement_text", None)
    paper_meta.pop("pdf_text", None)
    paper_meta.pop("supplement_hits", None)
    paper_meta.pop("supplement_file_index", None)
    paper_meta["seconds"] = round(time.time() - started, 2)
    return {"gene": gene, "pmid": pmid, "rows": rows_out, "paper": paper_meta}


def read_inventory(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def tranche_membership(registry: Path) -> dict[tuple[str, str], str]:
    reg = json.loads(registry.read_text())
    out: dict[tuple[str, str], str] = {}
    for tier in reg.get("tiers", []):
        manifest = registry.parent / tier["manifest"]
        for line in manifest.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            gene, pmid = line.split("\t")[:2]
            out[(gene.strip(), pmid.strip())] = tier["id"]
    return out


def summarize(rows: list[dict], papers: dict) -> dict:
    total = len(rows)
    by_class = Counter(r["class"] for r in rows)
    by_gene: dict = defaultdict(Counter)
    by_tranche: dict = defaultdict(Counter)
    by_prov: dict = defaultdict(Counter)
    for r in rows:
        by_gene[r["gene"]][r["class"]] += 1
        by_tranche[r.get("tranche") or "unassigned"][r["class"]] += 1
        by_prov[r.get("gold_provenance") or "unknown"][r["class"]] += 1

    hard = sum(by_class[c] for c in HARD_CEILING)
    wide = sum(by_class[c] for c in WIDE_CEILING)

    absent_by_paper: Counter = Counter()
    for r in rows:
        if r["class"] in WIDE_CEILING:
            absent_by_paper[(r["gene"], r["pmid"], r["class"])] += 1
    top_papers = []
    paper_absent_total: Counter = Counter()
    for (gene, pmid, klass), n in absent_by_paper.items():
        paper_absent_total[(gene, pmid)] += n
    for (gene, pmid), n in paper_absent_total.most_common(40):
        meta = papers.get((gene, pmid), {})
        classes = {
            k: v for (g, p, k), v in absent_by_paper.items() if g == gene and p == pmid
        }
        top_papers.append(
            {
                "gene": gene,
                "pmid": pmid,
                "text_absent_rows": n,
                "gold_rows": sum(
                    1 for r in rows if r["gene"] == gene and r["pmid"] == pmid
                ),
                "classes": classes,
                "figure_files": meta.get("figure_files"),
                "unfetched_links": meta.get("unfetched_links"),
                "advertised_unresolved_refs": meta.get("advertised_unresolved_refs"),
                "supplement_converted_files": meta.get("supplement_converted_files"),
                "body_chars": meta.get("body_chars"),
            }
        )

    def pct(n: int, d: int) -> float:
        return round(100.0 * n / d, 3) if d else 0.0

    # Count-bearing gold: how much of the count signal sits behind the ceiling?
    def nonzero(r: dict, f: str) -> bool:
        v = r.get(f"gold_{f}")
        return isinstance(v, int) and v > 0

    carrier_rows = [r for r in rows if nonzero(r, "carriers")]
    carrier_behind_hard = sum(1 for r in carrier_rows if r["class"] in HARD_CEILING)
    carrier_behind_wide = sum(1 for r in carrier_rows if r["class"] in WIDE_CEILING)

    return {
        "gold_rows": total,
        "papers": len(papers),
        "by_class": dict(by_class),
        "by_class_pct": {k: pct(v, total) for k, v in by_class.items()},
        "hard_ceiling_rows": hard,
        "hard_ceiling_pct": pct(hard, total),
        "wide_ceiling_rows": wide,
        "wide_ceiling_pct": pct(wide, total),
        "max_reachable_recall_if_hard_excluded_pct": pct(total - hard, total),
        "max_reachable_recall_if_wide_excluded_pct": pct(total - wide, total),
        "nonzero_carrier_rows": len(carrier_rows),
        "nonzero_carrier_rows_behind_hard_ceiling": carrier_behind_hard,
        "nonzero_carrier_rows_behind_wide_ceiling": carrier_behind_wide,
        "by_gene": {g: dict(c) for g, c in sorted(by_gene.items())},
        "by_tranche": {t: dict(c) for t, c in sorted(by_tranche.items())},
        "by_gold_provenance": {p: dict(c) for p, c in sorted(by_prov.items())},
        "top_text_absent_papers": top_papers,
        "papers_with_conversion_failures": sum(
            1 for m in papers.values() if m.get("supplement_failed_files")
        ),
        "supplement_files_converted": sum(
            m.get("supplement_converted_files", 0) for m in papers.values()
        ),
        "supplement_files_failed": sum(
            m.get("supplement_failed_files", 0) for m in papers.values()
        ),
    }


def render_md(summary: dict, args_desc: str) -> str:
    lines = [
        "# Gold-row source presence sweep",
        "",
        f"Scope: {args_desc}",
        "",
        "Classification is blind to predictions and reads only the gold row and",
        "the source on disk (article body, converted supplements incl. .docx/.xlsx/",
        ".doc/.pptx/.zip, article PDFs). Nothing under `corpus/` was written.",
        "",
        f"- Gold rows: **{summary['gold_rows']}** across {summary['papers']} gene-paper attempts",
        f"- Hard ceiling (source absent or substitution absent from all text, no figures): "
        f"**{summary['hard_ceiling_rows']} rows ({summary['hard_ceiling_pct']}%)** — "
        f"max paper-derived recall if excluded: {summary['max_reachable_recall_if_hard_excluded_pct']}%",
        f"- Wide ceiling (also figures-present and non-searchable notation): "
        f"**{summary['wide_ceiling_rows']} rows ({summary['wide_ceiling_pct']}%)** — "
        f"max recall if excluded: {summary['max_reachable_recall_if_wide_excluded_pct']}%",
        f"- Non-zero carrier gold rows: {summary['nonzero_carrier_rows']}; behind hard ceiling "
        f"{summary['nonzero_carrier_rows_behind_hard_ceiling']}, behind wide ceiling "
        f"{summary['nonzero_carrier_rows_behind_wide_ceiling']}",
        f"- Supplement files converted: {summary['supplement_files_converted']}; failed: "
        f"{summary['supplement_files_failed']} (papers with failures: {summary['papers_with_conversion_failures']})",
        "",
        "## By class",
        "",
        "| class | rows | % |",
        "|---|---:|---:|",
    ]
    for klass in CLASSES:
        n = summary["by_class"].get(klass, 0)
        lines.append(f"| `{klass}` | {n} | {summary['by_class_pct'].get(klass, 0.0)} |")
    lines += [
        "",
        "## By gene",
        "",
        "| gene | rows | " + " | ".join(CLASSES) + " |",
        "|---|---:|" + "---:|" * len(CLASSES),
    ]
    for gene, counts in summary["by_gene"].items():
        total = sum(counts.values())
        lines.append(
            f"| {gene} | {total} | "
            + " | ".join(str(counts.get(c, 0)) for c in CLASSES)
            + " |"
        )
    lines += [
        "",
        "## By tranche",
        "",
        "| tranche | rows | hard ceiling | wide ceiling | present_in_supplement_only |",
        "|---|---:|---:|---:|---:|",
    ]
    for tranche, counts in summary["by_tranche"].items():
        total = sum(counts.values())
        hard = sum(counts.get(c, 0) for c in HARD_CEILING)
        wide = sum(counts.get(c, 0) for c in WIDE_CEILING)
        lines.append(
            f"| {tranche} | {total} | {hard} | {wide} | {counts.get('present_in_supplement_only', 0)} |"
        )
    lines += [
        "",
        "## By gold provenance",
        "",
        "| provenance | rows | hard ceiling | wide ceiling |",
        "|---|---:|---:|---:|",
    ]
    for prov, counts in summary["by_gold_provenance"].items():
        total = sum(counts.values())
        hard = sum(counts.get(c, 0) for c in HARD_CEILING)
        wide = sum(counts.get(c, 0) for c in WIDE_CEILING)
        lines.append(f"| {prov} | {total} | {hard} | {wide} |")
    lines += [
        "",
        "## Papers with the most text-absent gold rows (wide ceiling)",
        "",
        "| gene | PMID | text-absent rows | gold rows | classes | figures | unfetched links | advertised unresolved supp refs | converted supps | body chars |",
        "|---|---:|---:|---:|---|---:|---:|---|---:|---:|",
    ]
    for p in summary["top_text_absent_papers"]:
        classes = ", ".join(f"{k}={v}" for k, v in sorted(p["classes"].items()))
        lines.append(
            f"| {p['gene']} | {p['pmid']} | {p['text_absent_rows']} | {p['gold_rows']} | {classes} | "
            f"{p['figure_files']} | {p['unfetched_links']} | {p['advertised_unresolved_refs']} | "
            f"{p['supplement_converted_files']} | {p['body_chars']} |"
        )
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    ap.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    ap.add_argument("--gold-root", type=Path, default=DEFAULT_GOLD_ROOT)
    ap.add_argument("--corpus", type=Path, default=REPO_ROOT / "corpus")
    ap.add_argument(
        "--cache-dir",
        type=Path,
        required=True,
        help="conversion cache, outside corpus/",
    )
    ap.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="restrict to one gene<TAB>pmid manifest",
    )
    ap.add_argument(
        "--include-unavailable",
        action="store_true",
        help="also sweep inventory rows marked source_unavailable",
    )
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    if not (args.corpus.is_dir()):
        print(f"corpus not reachable: {args.corpus}", file=sys.stderr)
        return 2
    args.out_dir.mkdir(parents=True, exist_ok=True)

    inventory = read_inventory(args.inventory)
    tranches = tranche_membership(args.registry)
    wanted: set[tuple[str, str]] | None = None
    if args.manifest:
        wanted = set()
        for line in args.manifest.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            g, p = line.split("\t")[:2]
            wanted.add((g.strip(), p.strip()))

    attempts: list[dict] = []
    for row in inventory:
        key = (row["gene"], row["pmid"])
        if wanted is not None and key not in wanted:
            continue
        if row["status"] == "quarantined":
            continue
        if row["status"] == "source_unavailable" and not args.include_unavailable:
            continue
        attempts.append(row)
    if args.limit:
        attempts = attempts[: args.limit]

    tasks = [
        {
            "gene": a["gene"],
            "pmid": a["pmid"],
            "corpus": str(args.corpus),
            "cache_dir": str(args.cache_dir),
            "gold_root": str(args.gold_root),
        }
        for a in attempts
    ]
    meta_by_key = {(a["gene"], a["pmid"]): a for a in attempts}

    rows: list[dict] = []
    papers: dict = {}
    errors: list[dict] = []
    started = time.time()
    done = 0
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(sweep_paper, t) for t in tasks]
        for fut in as_completed(futures):
            result = fut.result()
            done += 1
            key = (result["gene"], result["pmid"])
            inv = meta_by_key.get(key, {})
            if result.get("error"):
                errors.append(
                    {"gene": key[0], "pmid": key[1], "error": result["error"]}
                )
                continue
            papers[key] = result["paper"]
            for r in result["rows"]:
                r["tranche"] = tranches.get(key, "")
                r["gold_provenance"] = inv.get("gold_provenance", "")
                r["inventory_status"] = inv.get("status", "")
                r["body_chars"] = result["paper"].get("body_chars")
                r["searchable_chars"] = result["paper"].get("searchable_chars")
                r["figure_files"] = result["paper"].get("figure_files")
                rows.append(r)
            if done % 50 == 0 or done == len(tasks):
                print(
                    f"[{done}/{len(tasks)}] {time.time() - started:.0f}s rows={len(rows)}",
                    file=sys.stderr,
                    flush=True,
                )

    rows.sort(key=lambda r: (r["gene"], r["pmid"], r["variant"]))
    summary = summarize(rows, papers)
    summary["errors"] = errors
    summary["elapsed_seconds"] = round(time.time() - started, 1)
    summary["arguments"] = {
        "inventory": str(args.inventory),
        "registry": str(args.registry),
        "gold_root": str(args.gold_root),
        "corpus": str(args.corpus.resolve()),
        "manifest": str(args.manifest) if args.manifest else None,
        "include_unavailable": args.include_unavailable,
        "attempts": len(attempts),
    }

    fieldnames = [
        "gene",
        "pmid",
        "tranche",
        "gold_provenance",
        "inventory_status",
        "variant",
        "class",
        "where",
        "is_substitution",
        "gold_carriers",
        "gold_affected",
        "gold_unaffected",
        "body_chars",
        "searchable_chars",
        "figure_files",
    ]
    with (args.out_dir / "gold_rows.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore"
        )
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with (args.out_dir / "papers.json").open("w") as fh:
        json.dump(
            {f"{g}:{p}": m for (g, p), m in sorted(papers.items())},
            fh,
            indent=1,
            sort_keys=True,
        )
    with (args.out_dir / "summary.json").open("w") as fh:
        json.dump(summary, fh, indent=1, sort_keys=True)
    desc = (
        f"{len(attempts)} gene-paper attempts from {args.inventory.name}"
        + (f" restricted to {args.manifest.name}" if args.manifest else "")
        + (" (including source_unavailable)" if args.include_unavailable else "")
    )
    (args.out_dir / "summary.md").write_text(render_md(summary, desc))
    print(
        json.dumps(
            {
                k: summary[k]
                for k in (
                    "gold_rows",
                    "by_class",
                    "hard_ceiling_pct",
                    "wide_ceiling_pct",
                    "elapsed_seconds",
                )
            },
            indent=1,
        )
    )
    return 0


if __name__ == "__main__":
    os.environ.setdefault("GVF_DISABLE_LOCAL_DATA", "1")
    raise SystemExit(main())
