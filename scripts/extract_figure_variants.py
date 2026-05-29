#!/usr/bin/env python3
"""CLI: extract gene variants from figure images using a vision LLM.

Walks ``{pmc_fulltext_dir}/{pmid}_figures/`` for each requested PMID, asks a
multimodal model to read every image for variants of the given gene, and
writes one JSON report per PMID. Optionally appends the discovered variants
to a SQLite DB so they show up in recall comparisons.

Typical use::

    # Read figures for one PMID, just print JSON
    python scripts/extract_figure_variants.py \\
        --gene KCNH2 --pmid 24667783 \\
        --pmc-dir results/KCNH2/20260517_074737/pmc_fulltext \\
        --out /tmp/figure_reads

    # Process the top blockers and inject into the recall DB
    python scripts/extract_figure_variants.py \\
        --gene KCNH2 \\
        --pmid 29650123 --pmid 24667783 --pmid 19038855 \\
        --pmc-dir results/KCNH2/20260517_074737/pmc_fulltext \\
        --out recall_metrics/figure_reads_$(date +%Y%m%d) \\
        --db results/KCNH2/20260517_074737/KCNH2.db
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sqlite3
import sys
from pathlib import Path
from typing import Iterable, List, Tuple

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from harvesting.figure_variant_reader import (  # noqa: E402
    PMIDFigureReport,
    find_pmid_figures,
    read_figures_for_pmid,
)
from utils.variant_normalizer import VariantNormalizer  # noqa: E402


logger = logging.getLogger("extract_figure_variants")


# ---------------------------------------------------------------------------
# Figure-variant precision gate
#
# Vision/figure reads are noisy: the model sees pedigree legends, comparator
# genes, axis labels, and OCR artifacts, and the reader returns *no* confidence
# score. Without a gate, anything with a non-empty cdna/protein string is
# written RAW into the DB, and because figures run on every normal ``gvf-run``
# (via the recovery layer) this contaminates every new DB with potential false
# variants.
#
# ``GVF_FIGURE_VARIANT_GATE`` controls how aggressively figure variants are
# filtered before insertion. It is read fresh on every ingest call so tests and
# callers can flip it per-run:
#
#   off         - today's raw behavior: accept any variant with a non-empty
#                 cdna OR protein string. No validation. Use to reproduce the
#                 pre-gate baseline.
#   validate    - (DEFAULT) cheap sanity gate. Reject variants that (a) carry no
#                 well-formed cdna or protein notation, (b) sit at a protein
#                 position outside the gene's known length, or (c) are known
#                 non-target / oncology hotspot artifacts. Pure-text, no network.
#   corroborate - validate + a documented hook for future text-corroboration
#                 (require the variant to also appear in the assembled full
#                 text). NOT YET IMPLEMENTED: currently behaves exactly like
#                 ``validate``. See ``_figure_variant_passes_gate``.
#
# Unknown / empty values fall back to ``validate`` (the safe default).
# ---------------------------------------------------------------------------
FIGURE_VARIANT_GATE_ENV = "GVF_FIGURE_VARIANT_GATE"
_FIGURE_GATE_MODES = {"off", "validate", "corroborate"}
_FIGURE_GATE_DEFAULT = "validate"


def _figure_variant_gate_mode() -> str:
    """Return the active gate mode, defaulting to ``validate``.

    Read live (not cached) so a single process can score multiple runs under
    different policies. Unknown values warn once and fall back to the default.
    """
    raw = (os.environ.get(FIGURE_VARIANT_GATE_ENV) or "").strip().lower()
    if not raw:
        return _FIGURE_GATE_DEFAULT
    if raw not in _FIGURE_GATE_MODES:
        logger.warning(
            "Unknown %s=%r; falling back to %r (valid: %s)",
            FIGURE_VARIANT_GATE_ENV,
            raw,
            _FIGURE_GATE_DEFAULT,
            ", ".join(sorted(_FIGURE_GATE_MODES)),
        )
        return _FIGURE_GATE_DEFAULT
    return raw


def _figure_variant_passes_gate(
    gene: str,
    cdna: str | None,
    protein: str | None,
    mode: str | None = None,
) -> Tuple[bool, str]:
    """Decide whether a single figure-reader variant may enter the DB.

    Returns ``(passed, reason)``. ``reason`` is a short, stable token used for
    the dropped-variant breakdown log; on pass it is ``"ok"`` (or ``"gate_off"``
    when the gate is disabled).

    Reuses the canonical validators in ``utils.variant_normalizer``:
      * ``normalize_cdna`` / ``normalize_to_single_letter`` -> well-formedness.
      * ``is_non_target_variant`` -> hotspot + out-of-range rejection.
      * ``validate_position`` (via ``PROTEIN_LENGTHS``) -> protein-position
        bounds. No-ops for genes of unknown length (cold-start genes), so the
        gate never over-rejects an unfamiliar gene.
    """
    if mode is None:
        mode = _figure_variant_gate_mode()

    # ``off`` reproduces the historical behavior exactly: a non-empty cdna OR
    # protein string is enough. The empty-check itself lives in the caller, so
    # by the time we are here at least one of the two is non-empty.
    if mode == "off":
        return True, "gate_off"

    normalizer = VariantNormalizer(gene)

    # 1) Well-formed notation. A non-empty string is not enough; it must parse
    #    as a c.-style cdna OR a protein substitution/indel/frameshift/stop.
    cdna_ok = bool(cdna) and normalizer.normalize_cdna(cdna) is not None
    protein_ok = bool(protein) and (
        normalizer.normalize_to_single_letter(protein) is not None
        or normalizer.normalize_protein(protein) is not None
    )
    if not (cdna_ok or protein_ok):
        return False, "malformed_notation"

    # 2) Non-target / out-of-range protein artifacts. ``is_non_target_variant``
    #    folds in both the oncology-hotspot list and the position>length check.
    if protein:
        is_non_target, _reason = normalizer.is_non_target_variant(protein)
        if is_non_target:
            return False, "non_target"

        # 3) Explicit protein-position bound (belt-and-suspenders; also catches
        #    in-range hotspots' siblings that parse but fall outside length).
        pos = normalizer.extract_position(protein)
        if pos is not None and not normalizer.validate_position(pos):
            return False, "position_out_of_range"

    if mode == "corroborate":
        # TODO(corroborate): require the variant to also appear in the assembled
        # full text for this PMID before accepting. Extension point: thread the
        # PMID's text here and confirm presence of the normalized cdna/protein
        # form. Until that lands, ``corroborate`` == ``validate`` (the checks
        # above), so enabling it is safe but not yet stricter.
        return True, "ok"

    return True, "ok"


def _discover_pmids_with_figures(pmc_dir: Path) -> List[str]:
    """Find every PMID under *pmc_dir* whose ``{PMID}_figures/`` has images."""
    from harvesting.figure_variant_reader import is_image_path  # noqa: E402

    pmids: List[str] = []
    if not pmc_dir.is_dir():
        return pmids
    for fig_dir in sorted(pmc_dir.glob("*_figures")):
        if not fig_dir.is_dir():
            continue
        if any(is_image_path(p) for p in fig_dir.iterdir() if p.is_file()):
            pmid = fig_dir.name.replace("_figures", "")
            if pmid.isdigit():
                pmids.append(pmid)
    return pmids


def _load_pmid_list(args: argparse.Namespace, pmc_dir: Path) -> List[str]:
    pmids: List[str] = list(args.pmid or [])
    if args.pmid_file:
        for line in Path(args.pmid_file).read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                pmids.append(line.split()[0])
    if args.missing_csv:
        import csv

        with Path(args.missing_csv).open() as f:
            seen = set()
            for row in csv.DictReader(f):
                pid = (row.get("pmid") or "").strip()
                if pid and pid not in seen:
                    seen.add(pid)
                    pmids.append(pid)
    if args.auto_pmids:
        pmids.extend(_discover_pmids_with_figures(pmc_dir))
    if not pmids:
        sys.exit(
            "No PMIDs provided. Use --pmid, --pmid-file, --missing-csv, "
            "or --auto-pmids to walk the pmc-dir for any PMID with figures."
        )
    seen = set()
    out: List[str] = []
    for p in pmids:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def _ensure_paper(con: sqlite3.Connection, pmid: str, gene: str) -> None:
    if con.execute("SELECT 1 FROM papers WHERE pmid=?", (pmid,)).fetchone():
        return
    con.execute(
        "INSERT OR IGNORE INTO papers (pmid, gene_symbol, extraction_summary) "
        "VALUES (?, ?, 'figure-reader stub')",
        (pmid, gene),
    )


def _ensure_variant(
    con: sqlite3.Connection, gene: str, cdna: str | None, protein: str | None
) -> int:
    row = con.execute(
        """SELECT variant_id FROM variants
           WHERE gene_symbol=? AND cdna_notation IS ?
             AND protein_notation IS ? AND genomic_position IS NULL""",
        (gene, cdna, protein),
    ).fetchone()
    if row:
        return row[0]
    cur = con.execute(
        "INSERT INTO variants (gene_symbol, cdna_notation, protein_notation) VALUES (?, ?, ?)",
        (gene, cdna, protein),
    )
    return int(cur.lastrowid)


def ingest_report(
    report: PMIDFigureReport, db_path: Path, source_tag: str = "figure-reader"
) -> int:
    """Write extracted variants to *db_path*. Returns number of new variant_papers rows."""
    if not report.distinct_variants:
        return 0
    return ingest_cached_variants(
        pmid=report.pmid,
        gene=report.gene,
        distinct=report.distinct_variants,
        db_path=db_path,
        source_tag=source_tag,
    )


def ingest_cached_variants(
    pmid: str,
    gene: str,
    distinct: list,
    db_path: Path,
    source_tag: str = "figure-reader (cached)",
) -> int:
    """Write a pre-deduped variant list into *db_path*.

    Each candidate passes through ``_figure_variant_passes_gate`` (controlled by
    ``GVF_FIGURE_VARIANT_GATE``, default ``validate``) *before* it is written.
    Dropped variants are counted by reason and logged — never silently
    truncated. With ``GVF_FIGURE_VARIANT_GATE=off`` the gate is a no-op and the
    historical raw behavior is reproduced.
    """
    if not distinct:
        return 0
    gate_mode = _figure_variant_gate_mode()
    con = sqlite3.connect(str(db_path))
    try:
        _ensure_paper(con, pmid, gene)
        added = 0
        dropped: dict[str, int] = {}
        for v in distinct:
            cdna = (v.get("cdna") or "").strip() or None
            protein = (v.get("protein") or "").strip() or None
            if not (cdna or protein):
                dropped["empty_notation"] = dropped.get("empty_notation", 0) + 1
                continue
            passed, reason = _figure_variant_passes_gate(
                gene, cdna, protein, mode=gate_mode
            )
            if not passed:
                dropped[reason] = dropped.get(reason, 0) + 1
                logger.debug(
                    "PMID %s: gate(%s) dropped figure variant cdna=%r protein=%r (%s)",
                    pmid,
                    gate_mode,
                    cdna,
                    protein,
                    reason,
                )
                continue
            vid = _ensure_variant(con, gene, cdna, protein)
            exists = con.execute(
                "SELECT 1 FROM variant_papers WHERE variant_id=? AND pmid=?",
                (vid, pmid),
            ).fetchone()
            if exists:
                continue
            note = json.dumps(
                {k: v.get(k) for k in ("carriers", "affected", "unaffected", "context")}
            )
            con.execute(
                """INSERT INTO variant_papers
                   (variant_id, pmid, source_location, additional_notes, key_quotes)
                   VALUES (?, ?, ?, ?, ?)""",
                (vid, pmid, source_tag, note, "[]"),
            )
            added += 1
        con.commit()
        if dropped:
            total_dropped = sum(dropped.values())
            breakdown = ", ".join(
                f"{reason}={count}" for reason, count in sorted(dropped.items())
            )
            logger.info(
                "PMID %s: figure gate(%s) accepted %d, dropped %d (%s)",
                pmid,
                gate_mode,
                added,
                total_dropped,
                breakdown,
            )
        return added
    finally:
        con.close()


def run(args: argparse.Namespace) -> int:
    pmc_dir = Path(args.pmc_dir).expanduser()
    pmids = _load_pmid_list(args, pmc_dir)
    out_dir = Path(args.out).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)

    summary = {
        "gene": args.gene,
        "pmc_dir": str(pmc_dir),
        "model": args.model,
        "pmids": [],
    }

    total_variants = 0
    total_db_added = 0
    cached_pmids = 0
    for pmid in pmids:
        images = find_pmid_figures(pmc_dir, pmid)
        if not images:
            logger.warning("PMID %s: no figures on disk under %s", pmid, pmc_dir)
            summary["pmids"].append(
                {
                    "pmid": pmid,
                    "figures": 0,
                    "variants": 0,
                    "db_added": 0,
                    "cached": False,
                }
            )
            continue

        if args.max_images is not None:
            images = images[: args.max_images]

        report_file = out_dir / f"{pmid}.json"
        if report_file.exists() and not args.force:
            try:
                cached = json.loads(report_file.read_text())
                logger.info(
                    "PMID %s: cached report at %s, skipping (use --force to re-run)",
                    pmid,
                    report_file,
                )
                cached_pmids += 1
                variants = cached.get("distinct_variants", [])
                total_variants += len(variants)
                db_added = 0
                if args.db and variants:
                    # Rebuild a minimal report object for ingestion
                    report = PMIDFigureReport(
                        pmid=cached.get("pmid", pmid),
                        gene=cached.get("gene", args.gene),
                    )
                    # We only need distinct_variants on the report for ingest
                    report.per_figure = []  # not needed for ingest
                    # Inject the deduped variant list directly
                    setattr(report, "_cached_distinct", variants)
                    db_added = ingest_cached_variants(
                        pmid=cached.get("pmid", pmid),
                        gene=cached.get("gene", args.gene),
                        distinct=variants,
                        db_path=Path(args.db),
                    )
                    total_db_added += db_added
                summary["pmids"].append(
                    {
                        "pmid": pmid,
                        "figures": len(images),
                        "variants": len(variants),
                        "db_added": db_added,
                        "cached": True,
                    }
                )
                continue
            except Exception as exc:
                logger.warning(
                    "PMID %s: cached read failed (%s); re-running", pmid, exc
                )

        logger.info("PMID %s: reading %d figure(s)", pmid, len(images))
        report = read_figures_for_pmid(
            pmid=pmid,
            gene=args.gene,
            pmc_fulltext_dir=pmc_dir,
            model=args.model,
            max_images=args.max_images,
        )

        report_file.write_text(json.dumps(report.to_dict(), indent=2))
        variants = report.distinct_variants
        total_variants += len(variants)
        db_added = 0
        if args.db and variants:
            db_added = ingest_report(report, Path(args.db))
            total_db_added += db_added
        logger.info(
            "PMID %s: %d distinct variant(s), %d added to DB",
            pmid,
            len(variants),
            db_added,
        )
        summary["pmids"].append(
            {
                "pmid": pmid,
                "figures": len(images),
                "variants": len(variants),
                "db_added": db_added,
                "cached": False,
            }
        )
    summary["cached_pmids"] = cached_pmids

    summary["total_variants"] = total_variants
    summary["total_db_added"] = total_db_added
    summary_file = out_dir / "summary.json"
    summary_file.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--gene", required=True, help="Gene symbol, e.g. KCNH2")
    p.add_argument(
        "--pmc-dir",
        required=True,
        help="Directory containing {PMID}_figures/ subdirectories",
    )
    p.add_argument(
        "--out", required=True, help="Output directory for per-PMID JSON reports"
    )
    p.add_argument(
        "--pmid",
        action="append",
        default=None,
        help="PMID to process (repeat for multiple)",
    )
    p.add_argument(
        "--pmid-file", help="File with one PMID per line ('#' comments allowed)"
    )
    p.add_argument(
        "--missing-csv",
        help="CSV with a 'pmid' column (e.g. recall missing_in_sqlite.csv)",
    )
    p.add_argument(
        "--auto-pmids",
        action="store_true",
        help="Auto-discover PMIDs by walking pmc-dir for any {PMID}_figures/ "
        "subdirectory with at least one image file. Composable with the other "
        "PMID sources.",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Re-run the vision LLM even when a cached per-PMID report exists "
        "at {out}/{PMID}.json. Default: skip cached PMIDs.",
    )
    p.add_argument(
        "--db",
        help="Optional SQLite DB to insert extracted variants into",
    )
    p.add_argument(
        "--model",
        default=None,
        help="Vision model override (defaults to config.settings.get_vision_model)",
    )
    p.add_argument(
        "--max-images",
        type=int,
        default=None,
        help="Per-PMID image cap (default: all)",
    )
    p.add_argument("--verbose", "-v", action="store_true")
    return p


def main() -> int:
    args = build_parser().parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
