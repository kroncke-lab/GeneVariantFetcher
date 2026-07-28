#!/usr/bin/env python3
"""CLI wrapper for the packaged LLM-trace HTML viewer.

The implementation lives in :mod:`utils.llm_trace_html` because `cli/gvf_run.py`
and `cli/automated_workflow.py` build the report on every run, and `scripts` is
excluded from the wheel — importing the builder from here made the tracing
feature's headline artifact silently dead in an installed `gvf`.

Usage:
    .venv/bin/python scripts/build_llm_trace_html.py /path/to/run
    .venv/bin/python scripts/build_llm_trace_html.py /path/to/run/llm_traces \
        --output /path/to/trace_review.html
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from utils.llm_trace_html import (  # noqa: E402
    DEFAULT_MAX_FIELD_CHARS,
    DEFAULT_MAX_PAPERS_PER_FILE,
    DEFAULT_MAX_RECORDS_PER_PAPER,
    TRACE_REPORT_NAME,
    build_trace_html_report,
    collect_trace_report_data,
    render_trace_report,
)

__all__ = [
    "TRACE_REPORT_NAME",
    "build_trace_html_report",
    "collect_trace_report_data",
    "render_trace_report",
    "main",
]


def _resolve_paths(value: Path, output: Path | None) -> tuple[Path, Path, Path]:
    candidate = value.expanduser().resolve()
    if (candidate / "llm_traces").is_dir():
        run_dir = candidate
        trace_root = candidate / "llm_traces"
    else:
        trace_root = candidate
        run_dir = candidate.parent
    destination = (
        output.expanduser().resolve()
        if output is not None
        else run_dir / TRACE_REPORT_NAME
    )
    return run_dir, trace_root, destination


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("path", type=Path, help="Run directory or llm_traces directory")
    parser.add_argument("--output", type=Path, help="Output HTML file")
    parser.add_argument("--title", help="Reader-facing report title")
    parser.add_argument("--run-id", help="Run id to attribute this report to")
    parser.add_argument(
        "--max-field-chars",
        type=int,
        default=DEFAULT_MAX_FIELD_CHARS,
        help="Cap on each embedded prompt/response body (default %(default)s).",
    )
    parser.add_argument(
        "--max-papers-per-file",
        type=int,
        default=DEFAULT_MAX_PAPERS_PER_FILE,
        help="Above this many trace groups the report shards per paper "
        "(default %(default)s).",
    )
    parser.add_argument(
        "--max-records-per-paper",
        type=int,
        default=DEFAULT_MAX_RECORDS_PER_PAPER,
        help="Cap on records embedded per paper (default %(default)s).",
    )
    parser.add_argument(
        "--locked",
        action="store_true",
        help="Assert that a lock file covers this manifest, so the integrity "
        "banner may report the 'locked' level.",
    )
    args = parser.parse_args(argv)

    run_dir, trace_root, output_path = _resolve_paths(args.path, args.output)
    if not trace_root.is_dir():
        parser.error(f"trace directory not found: {trace_root}")
    data = build_trace_html_report(
        trace_root,
        output_path=output_path,
        run_dir=run_dir,
        title=args.title,
        run_id=args.run_id,
        max_field_chars=args.max_field_chars,
        max_papers_per_file=args.max_papers_per_file,
        max_records_per_paper=args.max_records_per_paper,
        locked=args.locked,
    )
    summary = data["summary"]
    print(f"LLM trace report: {output_path}")
    print(
        "  "
        f"{summary['paper_count']} papers, {summary['llm_call_count']} calls, "
        f"{summary['decision_event_count']} decisions"
    )
    integrity = data["integrity"]
    if integrity["valid"]:
        print(f"  trace manifest integrity: {integrity['level']}")
        print(f"    {integrity['level_description']}")
    else:
        print(f"  trace manifest integrity: {len(integrity['errors'])} issue(s)")
        for error in integrity["errors"][:10]:
            print(f"    - {error}")
    missing = data["coverage"]["missing_decision_links"]
    if missing:
        print(f"  route-coverage gaps: {len(missing)}")
        for gap in missing[:10]:
            print(f"    - {gap['stage']}: {gap['reason']}")
    if data.get("omissions"):
        print(f"  omitted from the page: {len(data['omissions'])} item(s)")
    if data.get("sharded"):
        print(f"  sharded per paper under: {data['shard_dir']}")
    print("  Open the HTML file in a web browser; no server is required.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
