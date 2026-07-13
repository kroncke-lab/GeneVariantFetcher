import json
import os
import time
from pathlib import Path
from types import SimpleNamespace

import cli.gvf_run as gvf_run


def test_csv_has_nonheader_rows(tmp_path: Path):
    missing = tmp_path / "missing.csv"
    assert gvf_run._csv_has_nonheader_rows(missing) is False

    header_only = tmp_path / "header.csv"
    header_only.write_text("pmid,doi\n", encoding="utf-8")
    assert gvf_run._csv_has_nonheader_rows(header_only) is False

    populated = tmp_path / "populated.csv"
    populated.write_text("pmid,doi\n12345,10.1/example\n", encoding="utf-8")
    assert gvf_run._csv_has_nonheader_rows(populated) is True


def test_step_source_recovery_continues_after_partial_fetch_and_refreshes(
    tmp_path: Path, monkeypatch
):
    run_dir = tmp_path / "results" / "GENE" / "run"
    source_qc = run_dir / "source_qc"
    source_qc.mkdir(parents=True)
    (source_qc / "source_acquisition_worklist.csv").write_text(
        "gene,pmid,action,route\nGENE,12345,fetch,fetch_elsevier_insttoken\n",
        encoding="utf-8",
    )
    (source_qc / "fetch_input.csv").write_text(
        "pmid,doi,route\n12345,10.1016/example,fetch_elsevier_insttoken\n",
        encoding="utf-8",
    )
    (source_qc / "source_override.csv").write_text(
        "gene,pmid,action,route,available_context_path\n"
        "GENE,67890,refresh_replay,existing,/tmp/67890_FULL_CONTEXT.md\n",
        encoding="utf-8",
    )

    calls = []
    refreshed_db = run_dir / "GENE.refresh_test.db"

    def fake_run(cmd, capture_output, text):
        del capture_output, text
        calls.append(cmd)
        rendered = " ".join(str(part) for part in cmd)
        if "fetch_paywalled.py" in rendered:
            return SimpleNamespace(returncode=1, stdout="", stderr="partial fetch")
        if "summarize_acquisition_outcome.py" in rendered:
            out = Path(cmd[cmd.index("--out") + 1])
            out.write_text(
                json.dumps(
                    {
                        "pmid_coverage": {
                            "selected_for_fetch_download": {
                                "pmids": 1,
                                "total_pmids": 2,
                                "coverage": 0.5,
                            }
                        }
                    }
                ),
                encoding="utf-8",
            )
            override_out = Path(cmd[cmd.index("--source-override-out") + 1])
            override_out.write_text(
                "gene,pmid,action,route,available_context_path\n"
                "GENE,12345,refresh_replay,fetched,/tmp/12345_FULL_CONTEXT.md\n",
                encoding="utf-8",
            )
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        if "refresh_run_db.py" in rendered:
            refresh_dir = run_dir / "refresh_test"
            layers_dir = refresh_dir / "layers"
            layers_dir.mkdir(parents=True)
            (layers_dir / "progression.json").write_text(
                json.dumps({"scoring_enabled": False, "progression": []}),
                encoding="utf-8",
            )
            refreshed_db.write_bytes(b"sqlite")
            summary = refresh_dir / "refresh_summary.json"
            summary.write_text(
                json.dumps({"active_db": str(refreshed_db)}),
                encoding="utf-8",
            )
            future = time.time() + 1
            os.utime(summary, (future, future))
            return SimpleNamespace(returncode=0, stdout="{}", stderr="")
        raise AssertionError(f"unexpected command: {cmd}")

    monkeypatch.setattr(gvf_run.subprocess, "run", fake_run)

    result = gvf_run.step_source_recovery(
        gene="GENE",
        run_dir=run_dir,
        source_qc_dir=source_qc,
        gold=None,
        run_recovery_layers=True,
        timeout_s=12,
    )

    assert result is not None
    assert result.active_db == refreshed_db
    assert result.layer_outdir == run_dir / "refresh_test" / "layers"
    assert result.refresh_summary == run_dir / "refresh_test" / "refresh_summary.json"

    rendered_calls = [" ".join(str(part) for part in cmd) for cmd in calls]
    assert any("fetch_paywalled.py" in call for call in rendered_calls)
    assert (
        sum("summarize_acquisition_outcome.py" in call for call in rendered_calls) == 2
    )
    refresh_call = next(call for call in rendered_calls if "refresh_run_db.py" in call)
    assert "--stage-extractions" in refresh_call
    assert "--only-forced-pmids" in refresh_call
    assert refresh_call.count("--source-override-csv") == 2


def test_step_source_recovery_keeps_supplement_only_rows_out_of_paywall_fetch(
    tmp_path: Path, monkeypatch
):
    run_dir = tmp_path / "results" / "SCN5A" / "run"
    source_qc = run_dir / "source_qc"
    source_qc.mkdir(parents=True)
    (run_dir / "pmc_fulltext").mkdir()
    (source_qc / "source_acquisition_worklist.csv").write_text(
        "gene,pmid,action,route\n"
        "SCN5A,29325976,fetch_supplement_only,fetch_elsevier_supplements_only\n",
        encoding="utf-8",
    )
    (source_qc / "fetch_input.csv").write_text("PMID,DOI,route\n", encoding="utf-8")
    (source_qc / "supplement_input.csv").write_text(
        "PMID,DOI,route\n"
        "29325976,10.1016/j.hrthm.2018.01.014,fetch_elsevier_supplements_only\n",
        encoding="utf-8",
    )
    (source_qc / "source_override.csv").write_text(
        "gene,pmid,action,route,available_context_path\n", encoding="utf-8"
    )

    calls = []

    def fake_run(cmd, capture_output, text):
        del capture_output, text
        calls.append(cmd)
        rendered = " ".join(str(part) for part in cmd)
        if "fetch_elsevier_supplements.py" in rendered:
            Path(cmd[cmd.index("--source-override-out") + 1]).write_text(
                "gene,pmid,action,route,available_context_path\n"
                "SCN5A,29325976,refresh_replay,fetch_elsevier_supplements_only,"
                f"{run_dir / 'pmc_fulltext' / '29325976_FULL_CONTEXT.md'}\n",
                encoding="utf-8",
            )
            return SimpleNamespace(returncode=0, stdout="no new files", stderr="")
        if "summarize_acquisition_outcome.py" in rendered:
            Path(cmd[cmd.index("--out") + 1]).write_text("{}", encoding="utf-8")
            Path(cmd[cmd.index("--source-override-out") + 1]).write_text(
                "gene,pmid,action,route,available_context_path\n",
                encoding="utf-8",
            )
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        if "refresh_run_db.py" in rendered:
            refresh_dir = run_dir / "refresh_test"
            refresh_dir.mkdir()
            (refresh_dir / "refresh_summary.json").write_text(
                '{"active_db": "SCN5A.refreshed.db"}', encoding="utf-8"
            )
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        raise AssertionError(f"unexpected command: {cmd}")

    monkeypatch.setattr(gvf_run.subprocess, "run", fake_run)

    result = gvf_run.step_source_recovery(
        gene="SCN5A",
        run_dir=run_dir,
        source_qc_dir=source_qc,
        gold=None,
        run_recovery_layers=False,
        timeout_s=12,
    )

    assert result is not None
    rendered_calls = [" ".join(str(part) for part in cmd) for cmd in calls]
    assert any("fetch_elsevier_supplements.py" in call for call in rendered_calls)
    assert not any("fetch_paywalled.py" in call for call in rendered_calls)
    refresh_call = next(call for call in rendered_calls if "refresh_run_db.py" in call)
    assert "supplement_source_override.csv" in refresh_call
