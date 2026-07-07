#!/usr/bin/env python3
"""Rebuild the curated-extraction-eval fixture from the registry + gold + corpus.

The set is defined by **registry.tsv** (the editable source of truth: one row per
paper — gene, pmid, strategy, note). This script reads it and deterministically
regenerates every derived file, so the set can be rebuilt or audited from scratch
and is trivial to GROW: append a row to registry.tsv (or use add_paper.py) and
re-run.

It writes (all relative to this script's directory):
  - gold/normalized/<GENE>_recall_input.csv  frozen gold SUBSET (only set PMIDs).
        Same schema as the repo gold, so scripts/run_recall_suite.py scores it
        unmodified. Rows come from the repo gold standard when the PMID is there,
        else from gold_overrides/<GENE>_recall_input.csv (see below).
  - pmids/<GENE>.txt                          PMID-per-line lists for
        `gvf gvf-run <GENE> --pmid-file ...` (the re-extraction path).
  - manifest.csv / manifest.md                one row per paper (strategy, gold
        counts, where the gold came from, title, corpus path, PubMed URL).
  - sources/<GENE>/<PMID>/...                  LOCAL mini-corpus snapshot copied
        from corpus/ (gitignored: paper full text is never committed). This has
        the same layout as corpus/<GENE>/<PMID>/ so run_benchmark.py --mode
        extract can feed the normal gvf-run cache path. Skip with --no-sources.

ADDING PAPERS (the set is meant to grow as you find problem cases)
-----------------------------------------------------------------
1. Append a row to registry.tsv (gene<TAB>pmid<TAB>strategy<TAB>note), or run
   `python add_paper.py <GENE> <PMID> --strategy table --note "..."`.
2. The paper needs a gold answer to be scoreable:
   - If the gene already has a repo gold standard
     (gene_variant_fetcher_gold_standard/normalized/<GENE>_recall_input.csv) and
     the PMID has rows there, nothing else is needed.
   - Otherwise (a new gene, or a paper the curator never scored) supply the
     expected variants in gold_overrides/<GENE>_recall_input.csv (same schema:
     variant,pmid,carriers,affected,unaffected[,gold_v2_*]). This is how you add
     "problem" papers from genes without a full gold standard. See
     gold_overrides/README.md.
3. The source must be cached at corpus/<GENE>/<PMID>/ (true after any gvf-run that
   fetched it). build_fixture warns if it is missing.
4. Re-run `python build_fixture.py`. It reports any paper missing gold or source.

HOW THE ORIGINAL 24 WERE CHOSEN
-------------------------------
Goal: a small, strategy-diverse set of gold papers the current pipeline already
extracts WELL, so prompt/harness/guardrail changes can be measured on recall +
MAE cheaply. Candidate pool = every gold PMID scored against the canonical DBs; a
paper qualified when it had >=6 gold variants, >=85% row recall, and accurate
counts (summed |count error| <= 1.0 per matched row), and each pick was then
re-validated live against the canonical DBs (>=95% row recall). Tagged by
dominant extraction strategy (table/figure/text/mixed) and curated to balance
gene x strategy and span small-to-large with some count-disagreement cases so MAE
has signal. As the set grows with problem cases, those new papers will (by design)
score lower — that is the point of a regression set.
"""

from __future__ import annotations

import csv
import shutil
import sqlite3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
GOLD_SRC = REPO / "gene_variant_fetcher_gold_standard" / "normalized"
OVERRIDES_DIR = HERE / "gold_overrides"
CORPUS = REPO / "corpus"
REGISTRY = HERE / "registry.tsv"

VALID_STRATEGIES = {"table", "text", "figure", "mixed"}

# Most benchmark source is synced into corpus/<GENE>/<PMID>/, but a few no-gold
# new-gene runs only have their reusable cache under results/.../pmc_fulltext.
RUN_LOCAL_SOURCE_DIRS = {
    "BRCA1": [REPO / "results/BRCA1/20260616_132646/pmc_fulltext"],
    "BRCA2": [
        REPO
        / "results/BRCA2/20260606_134517_hereditary_breast_cancer_500/BRCA2/20260606_134519/pmc_fulltext"
    ],
    "MYBPC3": [REPO / "results/MYBPC3/20260616_132646/pmc_fulltext"],
    "APOE": [REPO / "results/APOE/20260616_132646/pmc_fulltext"],
}

# Canonical gold-CSV schema (matches the repo normalized recall inputs). The
# first five are required; gold_v2_* are optional adjudication-overlay columns.
CANON_HEADER = [
    "variant",
    "pmid",
    "carriers",
    "affected",
    "unaffected",
    "gold_v2_carriers",
    "gold_v2_affected",
    "gold_v2_unaffected",
    "gold_v2_status",
    "gold_v2_note",
    "gold_v2_source",
]

# Canonical DBs (used here only for titles; run_benchmark.py defaults to these
# for score mode). Kept in sync with docs/RECALL_STATUS.md.
CANONICAL_DBS = {
    "KCNH2": REPO / "results/KCNH2/e2e_working_20260529_full/02_strict/KCNH2.db",
    "KCNQ1": REPO
    / "validation_runs/20260517_203904/results/KCNQ1/20260517_204424/KCNQ1.db",
    "SCN5A": REPO
    / "validation_runs/turnkey_e2e_20260518_213934/results/SCN5A/20260518_213938/SCN5A.db",
    "RYR2": REPO
    / "validation_runs/turnkey_e2e_20260518_213934/results/RYR2/20260518_213938/RYR2.db",
}

PUBMED = "https://pubmed.ncbi.nlm.nih.gov/{pmid}/"


def load_registry() -> list[dict[str, str]]:
    """Read registry.tsv into a list of {gene,pmid,strategy,note}."""
    if not REGISTRY.exists():
        raise SystemExit(f"registry not found: {REGISTRY}")
    picks: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for lineno, raw in enumerate(REGISTRY.read_text().splitlines(), 1):
        line = raw.rstrip("\n")
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            raise SystemExit(
                f"registry.tsv line {lineno}: need at least gene<TAB>pmid<TAB>"
                f"strategy, got {line!r}"
            )
        gene = parts[0].strip().upper()
        pmid = parts[1].strip()
        strat = parts[2].strip().lower()
        note = parts[3].strip() if len(parts) > 3 else ""
        if strat not in VALID_STRATEGIES:
            print(
                f"  WARNING registry line {lineno}: unknown strategy {strat!r} "
                f"(expected {sorted(VALID_STRATEGIES)})"
            )
        key = (gene, pmid)
        if key in seen:
            print(
                f"  WARNING registry line {lineno}: duplicate {gene}/{pmid}, skipping"
            )
            continue
        seen.add(key)
        picks.append({"gene": gene, "pmid": pmid, "strategy": strat, "note": note})
    return picks


def _read_gold_csv(path: Path) -> dict[str, list[dict[str, str]]]:
    """pmid -> list of rows (each normalized to CANON_HEADER keys)."""
    out: dict[str, list[dict[str, str]]] = {}
    if not path.exists():
        return out
    with path.open(newline="") as f:
        for r in csv.DictReader(f):
            pmid = (r.get("pmid") or "").strip()
            if not pmid:
                continue
            row = {k: (r.get(k) or "") for k in CANON_HEADER}
            row["pmid"] = pmid
            out.setdefault(pmid, []).append(row)
    return out


def gold_for_gene(gene: str) -> tuple[dict[str, list[dict[str, str]]], dict[str, str]]:
    """Merge repo gold + local override gold for one gene.

    Returns (pmid -> rows, pmid -> source) where source is 'repo' or 'override'.
    Repo gold takes precedence for a PMID present in both.
    """
    repo = _read_gold_csv(GOLD_SRC / f"{gene}_recall_input.csv")
    override = _read_gold_csv(OVERRIDES_DIR / f"{gene}_recall_input.csv")
    merged: dict[str, list[dict[str, str]]] = dict(override)
    merged.update(repo)
    source = {p: "override" for p in override}
    source.update({p: "repo" for p in repo})
    return merged, source


def gold_stats(rows: list[dict[str, str]]) -> dict[str, int]:
    car = aff = unaff = 0
    for r in rows:
        for key, name in (
            ("carriers", "car"),
            ("affected", "aff"),
            ("unaffected", "unaff"),
        ):
            try:
                v = int(float(r.get(key) or 0))
            except ValueError:
                v = 0
            if name == "car":
                car += v
            elif name == "aff":
                aff += v
            else:
                unaff += v
    return {
        "gold_rows": len(rows),
        "carriers": car,
        "affected": aff,
        "unaffected": unaff,
    }


def title_for(gene: str, pmid: str) -> str:
    db = CANONICAL_DBS.get(gene)
    if db and db.exists():
        try:
            conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
            try:
                cols = [r[1] for r in conn.execute("PRAGMA table_info(papers)")]
                if "title" in cols:
                    row = conn.execute(
                        "SELECT title FROM papers WHERE pmid=?", (pmid,)
                    ).fetchone()
                    if row and (row[0] or "").strip():
                        t = row[0].strip()
                        if t.lower() not in {"unknown title", f"paper {pmid}"}:
                            return t
            finally:
                conn.close()
        except sqlite3.Error:
            pass
    fc = source_full_context(gene, pmid)
    if fc.exists():
        try:
            skip = {
                "abstract",
                "introduction",
                "methods",
                "results",
                "discussion",
                "main text",
                "conclusions",
                "references",
            }
            with fc.open() as f:
                for line in f:
                    if line.startswith("## ") and not line.startswith("### "):
                        t = line[3:].strip()
                        if t and t.lower() not in skip and len(t) > 12:
                            return t
        except OSError:
            pass
    return ""


def source_full_context(gene: str, pmid: str) -> Path:
    """Return the best cached FULL_CONTEXT path for a benchmark paper."""
    corpus_fc = CORPUS / gene / pmid / f"{pmid}_FULL_CONTEXT.md"
    if corpus_fc.exists():
        return corpus_fc
    for root in RUN_LOCAL_SOURCE_DIRS.get(gene, []):
        fc = root / f"{pmid}_FULL_CONTEXT.md"
        if fc.exists():
            return fc
    return corpus_fc


def copy_source_snapshot(corpus_dir: Path, dst_dir: Path, pmid: str) -> int:
    """Copy one corpus PMID directory into the benchmark's local mini corpus."""
    if dst_dir.exists():
        shutil.rmtree(dst_dir)
    dst_dir.mkdir(parents=True, exist_ok=True)

    copied_bytes = 0
    for name in (
        f"{pmid}_FULL_CONTEXT.md",
        f"{pmid}_CLEANED.md",
        f"{pmid}_artifacts.json",
    ):
        src = corpus_dir / name
        if src.is_file():
            dst = dst_dir / name
            shutil.copy2(src, dst)
            copied_bytes += dst.stat().st_size

    for suffix in ("_figures", "_supplements"):
        src = corpus_dir / f"{pmid}{suffix}"
        if src.is_dir():
            dst = dst_dir / src.name
            shutil.copytree(src, dst)
            copied_bytes += sum(f.stat().st_size for f in dst.rglob("*") if f.is_file())
    return copied_bytes


def main(argv: list[str]) -> int:
    materialize_sources = "--no-sources" not in argv
    picks = load_registry()

    by_gene: dict[str, list[dict[str, str]]] = {}
    for p in picks:
        by_gene.setdefault(p["gene"], []).append(p)

    (HERE / "gold" / "normalized").mkdir(parents=True, exist_ok=True)
    (HERE / "pmids").mkdir(parents=True, exist_ok=True)

    manifest_rows: list[dict[str, object]] = []
    src_bytes_total = 0
    missing_src: list[str] = []
    missing_gold: list[str] = []

    if materialize_sources:
        shutil.rmtree(HERE / "sources", ignore_errors=True)

    for gene, items in by_gene.items():
        gold_rows_by_pmid, gold_source = gold_for_gene(gene)

        # frozen gold subset CSV (canonical schema)
        out_csv = HERE / "gold" / "normalized" / f"{gene}_recall_input.csv"
        with out_csv.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=CANON_HEADER)
            w.writeheader()
            for it in items:
                for row in gold_rows_by_pmid.get(it["pmid"], []):
                    w.writerow(row)

        # pmid list
        (HERE / "pmids" / f"{gene}.txt").write_text(
            "\n".join(it["pmid"] for it in items) + "\n"
        )

        for it in items:
            pmid, strat, note = it["pmid"], it["strategy"], it["note"]
            rows = gold_rows_by_pmid.get(pmid, [])
            stats = gold_stats(rows)
            gsrc = gold_source.get(pmid, "MISSING")
            if gsrc == "MISSING" or stats["gold_rows"] == 0:
                missing_gold.append(f"{gene}/{pmid}")

            source_fc = source_full_context(gene, pmid)
            src_status = "ok" if source_fc.exists() else "MISSING_IN_CORPUS"
            if not source_fc.exists():
                missing_src.append(f"{gene}/{pmid}")

            if materialize_sources and source_fc.exists():
                src_bytes_total += copy_source_snapshot(
                    source_fc.parent,
                    HERE / "sources" / gene / pmid,
                    pmid,
                )

            manifest_rows.append(
                {
                    "gene": gene,
                    "pmid": pmid,
                    "strategy": strat,
                    "gold_variant_rows": stats["gold_rows"],
                    "gold_carriers": stats["carriers"],
                    "gold_affected": stats["affected"],
                    "gold_unaffected": stats["unaffected"],
                    "gold_source": gsrc,
                    "title": title_for(gene, pmid),
                    "why_selected": note,
                    "corpus_source": str(source_fc.relative_to(REPO))
                    if source_fc.exists()
                    else "",
                    "source_status": src_status,
                    "pubmed_url": PUBMED.format(pmid=pmid),
                }
            )

    manifest_rows.sort(key=lambda r: (str(r["gene"]), str(r["pmid"])))
    fields = [
        "gene",
        "pmid",
        "strategy",
        "gold_variant_rows",
        "gold_carriers",
        "gold_affected",
        "gold_unaffected",
        "gold_source",
        "title",
        "why_selected",
        "corpus_source",
        "source_status",
        "pubmed_url",
    ]
    with (HERE / "manifest.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(manifest_rows)

    # manifest.md (human-readable; never hand-edit)
    md = [
        "# Curated Extraction Eval — paper manifest",
        "",
        "Auto-generated by `build_fixture.py` from `registry.tsv`. Do not hand-edit.",
        "Machine-readable copy: `manifest.csv`.",
        "",
        "`strategy` = how the gold variants are primarily encoded, so the set",
        "exercises every extraction route: **table**, **text**, **figure** (read by",
        "the vision layer), **mixed** (table+figure routing). `gold` = where the",
        "scoring answer comes from: `repo` (gene_variant_fetcher_gold_standard) or",
        "`override` (gold_overrides/, curator-supplied for problem/no-gold papers).",
        "",
        "| gene | pmid | strategy | gold | gold rows | car/aff/unaf | title |",
        "| --- | --- | --- | --- | ---: | --- | --- |",
    ]
    for r in manifest_rows:
        counts = f"{r['gold_carriers']}/{r['gold_affected']}/{r['gold_unaffected']}"
        title = (str(r["title"]) or "—").replace("|", "/")
        md.append(
            f"| {r['gene']} | [{r['pmid']}]({r['pubmed_url']}) | {r['strategy']} | "
            f"{r['gold_source']} | {r['gold_variant_rows']} | {counts} | {title} |"
        )
    by_strat_md: dict[str, int] = {}
    for r in manifest_rows:
        by_strat_md[str(r["strategy"])] = by_strat_md.get(str(r["strategy"]), 0) + 1
    md += [
        "",
        f"**{len(manifest_rows)} papers** · "
        f"{sum(int(r['gold_variant_rows']) for r in manifest_rows)} gold variant-rows · "
        f"by strategy: {by_strat_md} · "
        f"by gene: {{ {', '.join(f'{g}: {len(v)}' for g, v in by_gene.items())} }}",
        "",
    ]
    (HERE / "manifest.md").write_text("\n".join(md))

    total_gold = sum(int(r["gold_variant_rows"]) for r in manifest_rows)
    by_strat: dict[str, int] = {}
    for r in manifest_rows:
        by_strat[str(r["strategy"])] = by_strat.get(str(r["strategy"]), 0) + 1
    print(f"Wrote {len(manifest_rows)} papers across {len(by_gene)} genes.")
    print(f"  total gold variant-rows: {total_gold}")
    print(f"  by strategy: {by_strat}")
    print(f"  by gene: {{ {', '.join(f'{g}: {len(v)}' for g, v in by_gene.items())} }}")
    if materialize_sources:
        print(f"  materialized sources: {src_bytes_total / 1024:.0f} KB (gitignored)")
    rc = 0
    if missing_gold:
        print(f"  WARNING no gold answer for: {missing_gold}")
        print(
            "    -> add their variants to gold_overrides/<GENE>_recall_input.csv "
            "(see gold_overrides/README.md), or remove them from registry.tsv."
        )
        rc = 1
    if missing_src:
        print(f"  WARNING missing corpus source for: {missing_src}")
        print(
            "    -> fetch them first (e.g. a gvf-run that includes the PMID), "
            "then re-run. score mode does not need source; extract mode does."
        )
        rc = 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
