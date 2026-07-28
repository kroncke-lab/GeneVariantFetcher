#!/usr/bin/env python
"""Fill missing per-variant counts in an existing run database.

Standalone entry point for :mod:`pipeline.count_recovery` — the additive
counterpart to the carrier guard / outlier guard / count classifier, which only
ever remove counts. Use this to backfill a database that has already been
extracted, without re-running extraction:

    scripts/recover_counts.py --db results/KCNH2/<run>/KCNH2.db --gene KCNH2 --dry-run
    scripts/recover_counts.py --db results/KCNH2/<run>/KCNH2.db --gene KCNH2

``--dry-run`` opens the database read-only and reports exactly what a real run
would write, including counts it would skip because the slot is already
populated. Always run it first.

Writes are confined to NULL count columns, land as ``trust_tier='quarantine'``
for Step 3.7 to promote, and are logged to ``count_recovery_log`` with the
grounding quote, the structured count role, and the evidence locator — so the
pass is auditable and reversible. The database is copied to
``<db>.before_count_recovery.db`` before the first write.

Per-call LLM traces are written to a run-scoped directory beside the database
(override with ``--trace-dir``); the manifest and a self-contained HTML report
are built at the end.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from config.settings import get_settings  # noqa: E402

# The reusable helpers live in the PACKAGED module: `scripts` is excluded from
# the wheel (pyproject [tool.setuptools.packages.find]), so anything cli/ or
# pipeline/ needs at runtime cannot live here. Re-exported for backward
# compatibility with existing callers and docs.
from pipeline.count_recovery import (  # noqa: E402
    SOURCE_BASENAMES,
    corpus_root,
    default_source_roots,
    make_source_resolver,
    recover_counts,
)

logger = logging.getLogger("recover_counts")

__all__ = [
    "SOURCE_BASENAMES",
    "corpus_root",
    "default_source_roots",
    "make_source_resolver",
    "main",
]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--db", type=Path, required=True)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--pmid-file", type=Path, help="Restrict to these PMIDs.")
    ap.add_argument(
        "--source-root",
        type=Path,
        action="append",
        default=None,
        help="Extra source directory; repeatable. Defaults to the run "
        "dir's pmc_fulltext then corpus/.",
    )
    ap.add_argument("--model", help="Override COUNT_RECOVERY_MODEL.")
    ap.add_argument(
        "--reasoning-effort", help="Override COUNT_RECOVERY_REASONING_EFFORT."
    )
    ap.add_argument("--fields", help="Comma list: carriers,affected,unaffected.")
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Report gaps and grounded answers without writing.",
    )
    ap.add_argument(
        "--include-linkage-rows",
        action="store_true",
        help="Inspection mode: also ask about ClinVar/PubTator linkage-only "
        "variants the paper may never mention. Recovery normally restricts "
        "itself to paper-derived rows.",
    )
    ap.add_argument(
        "--trace-dir",
        type=Path,
        help="Where to write per-call LLM traces. Defaults to "
        "<db parent>/llm_traces_count_recovery_<run id>. Set "
        "GVF_LLM_TRACE_ENABLED=false to disable retention entirely.",
    )
    ap.add_argument(
        "--no-backup",
        action="store_true",
        help="Skip the pre-mutation database copy (not recommended).",
    )
    ap.add_argument("--json-out", type=Path, help="Write the summary as JSON.")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )
    if not args.db.is_file():
        raise SystemExit(f"no such database: {args.db}")

    settings = get_settings()
    model = args.model or settings.count_recovery_model
    effort = args.reasoning_effort or settings.count_recovery_reasoning_effort
    fields = tuple(
        f.strip()
        for f in (args.fields or settings.count_recovery_fields).split(",")
        if f.strip()
    )
    pmids = None
    if args.pmid_file:
        pmids = [p.strip() for p in args.pmid_file.read_text().split() if p.strip()]

    roots: list[Path] = []
    for root in list(args.source_root or []) + default_source_roots(args.db):
        if root not in roots:
            roots.append(root)

    # Standalone recovery configures its OWN trace root so its calls are not
    # written into (and do not invalidate) the manifest of the extraction run
    # that produced this database.
    from utils.llm_trace import (
        TRACE_MANIFEST_NAME,
        build_trace_manifest,
        configure_llm_tracing,
        resolve_trace_location,
        safe_run_id,
        tracing_enabled_by_environment,
    )
    from utils.llm_trace_html import TRACE_REPORT_NAME, build_trace_html_report
    from utils.llm_utils import litellm_completion

    run_id = (
        f"recover-counts-{datetime.now(timezone.utc):%Y%m%dT%H%M%S%f}Z-{os.getpid()}"
    )
    trace_enabled = tracing_enabled_by_environment()
    # --trace-dir is this run's directory verbatim; GVF_LLM_TRACE_DIR is a
    # storage BASE and gets a per-run child, so two sequential recoveries against
    # the same database never share a trace tree.
    if args.trace_dir:
        trace_root, trace_run_id = Path(args.trace_dir), run_id
    else:
        location = resolve_trace_location(
            run_id, default_root=args.db.parent / f"llm_traces_{safe_run_id(run_id)}"
        )
        trace_root, trace_run_id = location.root, location.run_id
    configure_llm_tracing(trace_root, run_id=trace_run_id, enabled=trace_enabled)
    if trace_enabled:
        logger.info("LLM traces -> %s (run %s)", trace_root, trace_run_id)

    stats = recover_counts(
        args.db,
        args.gene,
        source_for_pmid=make_source_resolver(args.gene, roots),
        llm_caller=litellm_completion,
        model=model,
        reasoning_effort=effort,
        pmids=pmids,
        fields=fields,
        max_variants_per_call=settings.count_recovery_max_variants_per_call,
        max_source_chars=settings.count_recovery_max_source_chars,
        dry_run=args.dry_run,
        paper_derived_only=not args.include_linkage_rows,
        backup=not args.no_backup,
    )

    stats.pop("result_objects", None)
    papers = stats.pop("results", [])
    logger.info(
        "count recovery %s: %d gap(s) across %d paper(s) -> %d grounded, %d written "
        "(%d rejected, %d paper(s) without source, %d failed) via %s@%s",
        "DRY RUN" if args.dry_run else "applied",
        stats["gaps_found"],
        stats["papers_with_gaps"],
        stats["counts_accepted"],
        stats["counts_written"],
        stats["counts_rejected"],
        stats["papers_without_source"],
        stats["papers_failed"],
        model,
        effort,
    )
    if trace_enabled and Path(trace_root).is_dir():
        try:
            manifest = build_trace_manifest(
                trace_root,
                output_path=Path(trace_root) / TRACE_MANIFEST_NAME,
                run_id=trace_run_id,
            )
            report = (
                Path(trace_root).parent
                / f"{safe_run_id(trace_run_id)}_{TRACE_REPORT_NAME}"
            )
            build_trace_html_report(
                trace_root,
                output_path=report,
                title=f"{args.gene} · count recovery trace review",
                run_id=trace_run_id,
            )
            logger.info(
                "trace manifest: %d call(s), %d decision(s) -> %s",
                manifest["llm_call_count"],
                manifest["decision_event_count"],
                report,
            )
        except Exception as exc:  # noqa: BLE001 - tracing must not change the outcome
            logger.warning("could not finalize count-recovery trace artifacts: %s", exc)
    if args.json_out:
        payload = dict(stats)
        payload["papers"] = papers
        args.json_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        logger.info("summary -> %s", args.json_out)

    # A pass that attempted papers and grounded nothing while every batch failed
    # is a failure, not a no-op. Returning 0 there made a total model/parse
    # outage indistinguishable from "there was nothing to recover".
    if stats.get("all_batches_failed"):
        logger.error(
            "count recovery failed on every attempted paper (%d batch failure(s))",
            stats.get("batch_failures", 0),
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
