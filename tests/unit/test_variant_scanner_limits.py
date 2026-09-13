"""Offline regression coverage for scanner resource limits and durable skips."""

import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
from pathlib import Path

import pytest
from pydantic import ValidationError

from config.settings import Settings, get_settings
from utils import variant_scanner as module
from utils.variant_scanner import ScanResult, VariantScanner


def test_size_check_precedes_normalization_and_worker(monkeypatch, tmp_path):
    from utils import source_text

    def forbidden(*args, **kwargs):
        pytest.fail("oversized input must not reach normalization or a worker")

    monkeypatch.setattr(source_text, "normalize_source_text", forbidden)
    monkeypatch.setattr(module, "_scan_in_subprocess", forbidden)
    result = VariantScanner().scan(
        "A" * 101, "PMID_123", max_chars=100, audit_dir=tmp_path
    )
    assert result.status == "skipped"
    assert result.reason == "input_exceeds_max_chars"
    assert result.variants == []
    assert result.unique_normalized == set()
    assert result.get_hints_for_prompt() == ""
    record = json.loads(next(tmp_path.glob("*.json")).read_text())
    assert record == result.to_metadata()
    assert record["input_chars"] == 101
    assert record["max_chars"] == 100
    assert record["source"] == "PMID_123"
    assert record["elapsed_seconds"] >= 0


def test_exact_size_limit_is_allowed(monkeypatch):
    seen = []

    def worker(text, gene, source, budget):
        seen.append(text)
        return ScanResult()

    monkeypatch.setattr(module, "_scan_in_subprocess", worker)
    scanner = VariantScanner()
    assert scanner.scan("abc", max_chars=3).status == "complete"
    assert scanner.scan("abcd", max_chars=3).status == "skipped"
    assert seen == ["abc"]


@pytest.mark.parametrize(
    "field,value",
    [
        ("scanner_max_chars", 0),
        ("scanner_max_chars", -1),
        ("scanner_max_chars", 1.5),
        ("scanner_budget_seconds", 0),
        ("scanner_budget_seconds", -1),
        ("scanner_budget_seconds", float("inf")),
        ("scanner_budget_seconds", float("nan")),
    ],
)
def test_settings_reject_invalid_limits(field, value):
    with pytest.raises(ValidationError):
        Settings(_env_file=None, **{field: value})


def test_limits_are_environment_backed(monkeypatch):
    monkeypatch.setenv("SCANNER_MAX_CHARS", "12345")
    monkeypatch.setenv("SCANNER_BUDGET_SECONDS", "2.5")
    get_settings.cache_clear()
    assert get_settings().scanner_max_chars == 12345
    assert get_settings().scanner_budget_seconds == 2.5


@pytest.mark.parametrize("value", [0, -1, True, 1.5])
def test_explicit_invalid_character_limit_is_rejected(value):
    with pytest.raises(ValueError):
        VariantScanner().scan("A", max_chars=value)


def test_real_stuck_regex_is_killed_from_thread_and_next_scan_works(
    monkeypatch, tmp_path
):
    processes = []
    real_popen = subprocess.Popen

    def tracked_popen(*args, **kwargs):
        process = real_popen(*args, **kwargs)
        processes.append(process)
        return process

    monkeypatch.setattr(module.subprocess, "Popen", tracked_popen)
    with monkeypatch.context() as patch:
        # Tiny input, deliberately exponential regex: a size cap cannot help.
        patch.setattr(
            module,
            "_SCAN_WORKER_CODE",
            'import re; re.search(r"(a+)+$", "a" * 128 + "!")',
        )
        with ThreadPoolExecutor(max_workers=1) as pool:
            result = pool.submit(
                VariantScanner().scan,
                "small input",
                "PMID_456",
                budget_seconds=0.2,
                audit_dir=tmp_path,
            ).result(timeout=10)
    assert result.status == "timed_out"
    assert result.reason == "wall_clock_budget_exceeded"
    assert result.variants == []
    assert processes[0].poll() is not None  # killed AND reaped, not a runaway thread
    assert (
        json.loads(next(tmp_path.glob("*.json")).read_text())["status"] == "timed_out"
    )
    following = VariantScanner().scan("KCNH2 p.Arg534Cys was found.")
    assert following.status == "complete"
    assert "R534C" in following.unique_normalized


def test_worker_failure_is_visible_and_does_not_raise(monkeypatch, tmp_path):
    monkeypatch.setattr(module, "_SCAN_WORKER_CODE", 'print("not valid JSON")')
    result = VariantScanner().scan("KCNH2 R534C", audit_dir=tmp_path)
    assert result.status == "failed"
    assert result.reason == "worker_failed"
    assert result.variants == []
    assert json.loads(next(tmp_path.glob("*.json")).read_text())["status"] == "failed"


def test_unwritable_audit_is_logged_without_aborting_scan(
    monkeypatch, tmp_path, caplog
):
    def unavailable(*args, **kwargs):
        raise PermissionError("audit unavailable")

    monkeypatch.setattr(Path, "mkdir", unavailable)
    result = VariantScanner().scan("XX", max_chars=1, audit_dir=tmp_path / "audit")
    assert result.status == "skipped"
    assert "Could not persist variant scan audit" in caplog.text
    assert "audit_write_error" in result.stats


@pytest.mark.parametrize(
    "gene,text",
    [
        (
            "KCNH2",
            "KCNH2 p.Arg534Cys and c.1600C>T were found in the patient. IVS9+1G>A.",
        ),
        ("BRCA2", "BRCA2 deletion of exons 3-5 was identified in the patient."),
        ("SCN5A", "SCN5A p.(Arg18Gln) and p.Arg1623fs*10 were identified. ΔKPQ."),
    ],
)
def test_worker_round_trip_preserves_variants_context_offsets_and_stats(gene, text):
    scanner = VariantScanner(gene)
    expected = scanner._scan_unbounded(text)
    actual = scanner.scan(text)
    assert actual.status == "complete"
    assert [asdict(v) for v in actual.variants] == [
        asdict(v) for v in expected.variants
    ]
    assert actual.unique_normalized == expected.unique_normalized
    for key, value in expected.stats.items():
        assert actual.stats[key] == value


def test_pathological_suffixes_finish_under_external_watchdog():
    # Exercise the actual active regex objects without a scanner timeout masking
    # a regression in the regex repair. This process is killed on test failure.
    code = module._SCAN_WORKER_CODE[
        : module._SCAN_WORKER_CODE.index("from utils.variant_scanner")
    ]
    code += """
from utils.variant_scanner import VariantScanner
for pattern, prefix, suffix in [
    (VariantScanner.PROTEIN_HGVS_FULL, "p.Arg123fs", "_"),
    (VariantScanner.PROTEIN_HGVS_PAREN, "p.(Arg123fs", "_"),
    (VariantScanner.PROTEIN_THREE_LETTER, "Arg123fs", "_"),
]:
    assert pattern.search(prefix + "1" * 100_000 + suffix) is None
    assert pattern.search(prefix + "*10" + (")" if "(" in prefix else "")) is not None
from utils.structural_alleles import DELTA_RE
assert DELTA_RE.search('delta' + ' ' * 100_000 + '1') is None
assert DELTA_RE.search('delta  -  KPQ').group(1) == 'KPQ'
"""
    subprocess.run(
        [sys.executable, "-c", code, str(Path(module.__file__).parent)],
        check=True,
        timeout=5,
        capture_output=True,
        text=True,
    )


@pytest.mark.parametrize("llm_fails", [False, True])
def test_extraction_skip_audit_survives_llm_failure(monkeypatch, tmp_path, llm_fails):
    from pipeline.extraction import ExpertExtractor
    from utils.models import Paper

    monkeypatch.setenv("SCANNER_MAX_CHARS", "20")
    get_settings.cache_clear()
    extractor = ExpertExtractor(
        models=["test-model"], tier_threshold=1, fulltext_dir=str(tmp_path)
    )
    extractor.enable_ensemble_qa = False
    paper = Paper(
        pmid="26669661",
        gene_symbol="KCNH2",
        title="Clinical observations",
        full_text="We studied patients with KCNH2 variants and recorded clinical findings. "
        * 100,
    )
    monkeypatch.setattr(extractor, "_try_table_router", lambda *_: None)
    monkeypatch.setattr(extractor, "_extract_variants_from_tables", lambda *_: [])

    def call(*args, **kwargs):
        if llm_fails:
            raise RuntimeError("synthetic LLM failure after scanner skip")
        return {"variants": [], "extraction_metadata": {}}, False, "{}"

    monkeypatch.setattr(extractor, "call_llm_json_with_status", call)
    # A later adjudicator may replace the payload; the scanner status must
    # still be stamped on the final retained extraction.
    monkeypatch.setattr(
        extractor,
        "_maybe_adjudicate_extraction",
        lambda **kwargs: {"variants": [], "extraction_metadata": {}},
    )
    result = extractor.extract(paper)
    records = [
        json.loads(path.read_text())
        for path in (tmp_path / "variant_scan_audit").glob("*.json")
    ]
    assert records
    assert all(record["status"] == "skipped" for record in records)
    assert all(record["source"] == "PMID_26669661" for record in records)
    if llm_fails:
        assert not result.success
    else:
        assert result.success
        assert (
            result.extracted_data["extraction_metadata"]["variant_scan"]["status"]
            == "skipped"
        )


def test_possessive_repair_preserves_match_spans_and_capture_groups():
    import itertools
    import re

    for pattern, prefix, ending in [
        (VariantScanner.PROTEIN_HGVS_FULL, "p.Arg123", ""),
        (VariantScanner.PROTEIN_HGVS_PAREN, "p.(Arg123", ")"),
        (VariantScanner.PROTEIN_THREE_LETTER, "Arg123", ""),
    ]:
        original = re.compile(pattern.pattern.replace(r"\d*+", r"\d*"), pattern.flags)
        for suffix, extension, boundary in itertools.product(
            ["Cys", "fs", "fs*", "fs10", "fs*10", "del", "dup", "ins", "Ter", "*", "="],
            ["", "1", "123", "*", "*10", "X10"],
            ["", " ", ".", "_", "A"],
        ):
            text = prefix + suffix + extension + ending + boundary
            before, after = original.search(text), pattern.search(text)
            assert bool(before) == bool(after), text
            if before:
                assert (before.span(), before.groups()) == (
                    after.span(),
                    after.groups(),
                ), text
