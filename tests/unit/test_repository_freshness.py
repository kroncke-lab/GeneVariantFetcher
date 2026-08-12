"""Mechanical guards against stale repository references and status claims."""

from __future__ import annotations

import json
import re
from pathlib import Path
from urllib.parse import unquote


REPO = Path(__file__).resolve().parents[2]
MARKDOWN_LINK_RE = re.compile(r"(?<!!)\[[^\]]+\]\(([^)]+)\)")


def _documentation_files() -> list[Path]:
    roots = list(REPO.glob("*.md"))
    for directory in (
        "benchmarks",
        "data",
        "docs",
        "gene_variant_fetcher_gold_standard",
        "scripts",
    ):
        roots.extend((REPO / directory).rglob("*.md"))
    embedded_evidence_dirs = {"runs", "sources", "pmc_fulltext", "scout_output"}
    return sorted(
        path
        for path in roots
        if path.is_file()
        and not embedded_evidence_dirs.intersection(path.relative_to(REPO).parts)
    )


def test_local_markdown_links_resolve_without_editor_line_suffixes():
    broken: list[str] = []
    for document in _documentation_files():
        for raw_target in MARKDOWN_LINK_RE.findall(
            document.read_text(encoding="utf-8", errors="replace")
        ):
            target = raw_target.strip().strip("<>").split("#", 1)[0]
            if not target or re.match(r"^[a-z][a-z0-9+.-]*:", target, re.I):
                continue
            if re.search(r":\d+$", target):
                broken.append(f"{document.relative_to(REPO)} -> {raw_target}")
                continue
            resolved = (document.parent / unquote(target)).resolve()
            if not resolved.exists():
                broken.append(f"{document.relative_to(REPO)} -> {raw_target}")
    assert broken == [], "broken local Markdown links:\n" + "\n".join(broken)


def test_public_dashboard_is_explicitly_labeled_historical_until_rescore():
    status = (REPO / "docs" / "RECALL_STATUS.md").read_text(encoding="utf-8")
    dashboard = (REPO / "docs" / "dashboard" / "index.html").read_text(encoding="utf-8")
    strategy = json.loads(
        (REPO / "docs" / "dashboard" / "strategy.json").read_text(encoding="utf-8")
    )

    assert "Historical view" in status
    assert 'data-dashboard-status="historical"' in dashboard
    assert strategy["status"] == "historical pre-correction snapshot"


def test_runtime_sources_do_not_reference_retired_claude_35_haiku():
    offenders: list[str] = []
    for root_name in (
        "pipeline",
        "harvesting",
        "config",
        "utils",
        "gene_literature",
        "cli",
        "scripts",
    ):
        for path in (REPO / root_name).rglob("*.py"):
            if "claude-3-5-haiku-20241022" in path.read_text(
                encoding="utf-8", errors="replace"
            ):
                offenders.append(str(path.relative_to(REPO)))
    assert offenders == []


def test_tasks_is_the_single_forward_checklist_without_completed_history():
    tasks = (REPO / "TASKS.md").read_text(encoding="utf-8")
    headings = re.findall(r"^## (.+)$", tasks, flags=re.MULTILINE)

    assert headings == [
        "1. Re-establish the scientific baseline",
        "2. Recover missing source before adding inference",
        "3. Measure count semantics and recovery",
        "4. Finish the trust and fleet boundary",
        "5. Engineering handoff follow-ups",
        "Deliberate decisions and non-goals",
    ]
    assert "- [x]" not in tasks.lower()
    assert "START HERE" not in tasks
    assert "## Completed" not in tasks
    assert "## Active Tasks" not in tasks


def test_current_docs_do_not_restore_retired_cohort_or_output_labels():
    current_docs = [
        REPO / "README.md",
        REPO / "CLAUDE.md",
        REPO / "CONTRIBUTING.md",
        REPO / "TASKS.md",
        REPO / "docs" / "README.md",
        REPO / "docs" / "ARCHITECTURE.md",
        REPO / "docs" / "QUICKSTART.md",
        REPO / "docs" / "VARIANT_BROWSER_INTEGRATION.md",
        REPO / "benchmarks" / "curated_extraction_eval" / "README.md",
        REPO / "benchmarks" / "codex_paper_eval" / "README.md",
    ]
    retired = (
        "101-paper",
        "azure_first_101",
        "gvf_curated_101",
        "Current-production baselines",
        "./output/KCNH2/20260210_143022",
    )

    offenders: list[str] = []
    for document in current_docs:
        text = document.read_text(encoding="utf-8")
        for phrase in retired:
            if phrase in text:
                offenders.append(f"{document.relative_to(REPO)} -> {phrase}")
    assert offenders == []


def test_changelog_has_one_unreleased_section():
    changelog = (REPO / "CHANGELOG.md").read_text(encoding="utf-8")
    assert len(re.findall(r"^## \[Unreleased\]", changelog, flags=re.MULTILINE)) == 1
