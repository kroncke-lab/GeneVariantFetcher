"""Source-frozen, gold-free production comparison; CLI spend is excluded.

Run from the repository root. Credentials are obtained in memory from the
existing Azure login. No key is written to disk or passed in command arguments.
"""

import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from benchmarks.codex_paper_eval.setup_production_eval import runtime_fingerprint
from pipeline.source_snapshot import freeze_sources

HERE = Path(__file__).resolve().parent
PYTHON = ROOT / ".venv/bin/python"
HARNESS = ROOT / "benchmarks/codex_paper_eval"
GOLD = ROOT / "benchmarks/evaluation_tiers/mixed_gold_continuation_120/answer_key"
FROZEN = ROOT / "validation_runs/model_routing_20260906/frozen_sources"
ARMS = {
    "grok43": "azure_ai/grok-4.3",
    "grok46": "azure_ai/grok-4.6",
    "astra": "azure_ai/gpt-6-astra",
}
SHARED = {
    "MODEL_PROVIDER": "azure",
    "ANTHROPIC_API_KEY": "",
    "TIER2_MODEL": "azure_ai/gpt-5.6-luna",
    "TIER2_REASONING_EFFORT": "xhigh",
    "TABLE_ROUTER_MODEL": "azure_ai/Kimi-K2.6-1",
    "TIER3_ADJUDICATOR_MODELS": "azure_ai/gpt-5.6-sol",
    "TIER3_ADJUDICATOR_REASONING_EFFORT": "medium",
    "TIER3_VERIFIER_REASONING_EFFORT": "medium",
    "VISION_MODEL": "azure_ai/gpt-5.6-sol",
    "VISION_REASONING_EFFORT": "",
    "COUNT_RECOVERY_ENABLED": "0",
    "COUNT_RECOVERY_MODEL": "azure_ai/gpt-5.6-luna",
    "TIER3_MAX_TOKENS": "32000",
    "MAX_WORKERS": "1",
    "FILTER_MAX_WORKERS": "2",
    "AZURE_MAX_WORKERS": "2",
    "AZURE_RPM": "30",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_dir(arm):
    return HARNESS / "runs" / ("20260906_model12_" + arm)


def config(arm):
    return dict(
        SHARED,
        TIER3_MODELS=ARMS[arm],
        TIER3_REASONING_EFFORT="" if arm == "grok43" else "high",
    )


def run(command, logfile=None, env=None):
    if logfile:
        with Path(logfile).open("w") as log:
            subprocess.run(
                [str(x) for x in command],
                cwd=ROOT,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
    else:
        subprocess.run([str(x) for x in command], cwd=ROOT, env=env, check=True)


def prepare(arm):
    folder = run_dir(arm)
    assert not folder.exists(), folder
    run(
        [
            PYTHON,
            HARNESS / "run_eval.py",
            "prepare",
            "--run-id",
            folder.name,
            "--seed",
            "2026090601",
            "--paper-manifest",
            HERE / "paper_manifest.tsv",
            "--corpus-root",
            FROZEN,
            "--gold-root",
            GOLD,
            "--eligibility-mode",
            "variant",
            "--minimum-chars",
            "1000",
        ],
        HERE / (arm + "_prepare.log"),
    )
    selection = json.loads((folder / "selection.json").read_text())
    expected = {
        tuple(line.split())
        for line in (HERE / "paper_manifest.tsv").read_text().splitlines()
        if line.strip() and not line.startswith("#")
    }
    assert {(p["gene"], str(p["pmid"])) for p in selection["papers"]} == expected
    freeze_sources(selection["papers"], folder / "frozen_corpus")
    (folder / "pmids").mkdir()
    (folder / "production_runs").mkdir()
    (folder / "operator_logs").mkdir()
    for gene in sorted({g for g, p in expected}):
        (folder / "pmids" / f"{gene}.txt").write_text(
            "".join(p + "\n" for g, p in sorted(expected) if g == gene)
        )
    setup = {
        "classification": "Opened calibration; source-frozen primary model configuration comparison, not registered discovery or confirmation",
        "runtime": runtime_fingerprint(),
        "model_configuration": config(arm),
        "source_snapshot_sha256": sha(folder / "frozen_corpus/source_snapshot.json"),
        "campaign_hook_sha256": sha(HERE / "budget_guard.py"),
        "campaign_script_sha256": sha(Path(__file__)),
        "arm": arm,
        "attempts": len(expected),
        "budget_envelope_usd": 150,
        "gold_root": str(GOLD),
    }
    (folder / "analysis_setup.json").write_text(json.dumps(setup, indent=2) + "\n")
    (HERE / (arm + "_setup.json")).write_text(json.dumps(setup, indent=2) + "\n")
    print("PREPARED", arm, flush=True)


def environment(arm):
    from config.settings import get_settings

    s = get_settings()
    env = dict(os.environ, **config(arm))
    env["AZURE_AI_API_BASE"] = s.azure_ai_api_base
    env["AZURE_AI_API_KEY"] = s.azure_ai_api_key
    key = subprocess.check_output(
        [
            "az",
            "cognitiveservices",
            "account",
            "keys",
            "list",
            "--name",
            "magen-api-2-resource",
            "--resource-group",
            "MAGen",
            "--subscription",
            "3c2867cc-afb1-4ff8-96e8-4cee3e09d869",
            "--query",
            "key1",
            "-o",
            "tsv",
        ],
        text=True,
    ).strip()
    env["GVF_NEW_AZURE_API_KEY"] = key
    env["AZURE_AI_MODEL_ROUTES"] = json.dumps(
        {
            m: {
                "api_base": "https://magen-api-2-resource.services.ai.azure.com/openai/v1",
                "api_key_env": "GVF_NEW_AZURE_API_KEY",
            }
            for m in ("gpt-6-astra", "grok-4.6")
        }
    )
    env["GVF_EXPERIMENT_BUDGET"] = str(HERE / "budget.json")
    env["GVF_EXPERIMENT_ARM"] = arm
    env["GVF_FROZEN_SOURCE_MANIFEST"] = str(
        run_dir(arm) / "frozen_corpus/source_snapshot.json"
    )
    return env


def check(arm):
    folder = run_dir(arm)
    setup = json.loads((folder / "analysis_setup.json").read_text())
    assert runtime_fingerprint() == setup["runtime"], "Runtime drift"
    assert sha(HERE / "budget_guard.py") == setup["campaign_hook_sha256"], "Hook drift"
    assert (
        sha(folder / "frozen_corpus/source_snapshot.json")
        == setup["source_snapshot_sha256"]
    ), "Source manifest drift"
    return folder


def extract(arm):
    folder = check(arm)
    env = environment(arm)

    def gene_run(file):
        gene = file.stem
        print("START", arm, gene, flush=True)
        run(
            [
                PYTHON,
                HERE / "run_instrumented_cli.py",
                "gvf-run",
                gene,
                "--email",
                "brett.kroncke@gmail.com",
                "--output",
                folder / "production_runs",
                "--pmid-file",
                file,
                "--no-source-recovery",
                "--no-corpus-sync",
                "--no-publish-review",
                "--gold-free-run",
            ],
            folder / "operator_logs" / f"{gene}.log",
            env,
        )
        print("COMPLETED", arm, gene, flush=True)

    with concurrent.futures.ThreadPoolExecutor(3) as pool:
        list(pool.map(gene_run, sorted((folder / "pmids").glob("*.txt"))))
    check(arm)


def lock(arm):
    folder = check(arm)
    for path in sorted((folder / "production_runs").glob("*/*/RUN_STATUS.json")):
        status = json.loads(path.read_text())
        assert (
            status["status"] == "completed"
            and not status["stage_failures"]
            and status["gold_access"]["disabled"]
        ), path
    run(
        [
            PYTHON,
            HARNESS / "rebind_production_sources.py",
            "--run-dir",
            folder,
            "--production-root",
            folder / "production_runs",
        ],
        HERE / (arm + "_rebind.log"),
    )
    run(
        [
            PYTHON,
            HARNESS / "db_to_predictions.py",
            "--run-dir",
            folder,
            "--production-root",
            folder / "production_runs",
            "--trust-mode",
            "trusted",
            "--identity-mode",
            "trusted",
            "--paper-primary",
            "--out",
            folder / "predictions.json",
        ],
        HERE / (arm + "_projection.log"),
    )
    run(
        [PYTHON, HARNESS / "run_eval.py", "lock", "--run-dir", folder],
        HERE / (arm + "_lock.log"),
    )
    print("LOCKED", arm, flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("action", choices=["prepare", "extract", "lock", "score"])
    p.add_argument("arm", choices=ARMS)
    args = p.parse_args()
    if args.action == "score":
        run(
            [
                PYTHON,
                HARNESS / "run_eval.py",
                "score",
                "--run-dir",
                run_dir(args.arm),
                "--gold-root",
                GOLD,
            ],
            HERE / (args.arm + "_score.log"),
        )
    else:
        globals()[args.action](args.arm)
