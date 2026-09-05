"""Rebuild the compact campaign ledger from immutable, completed run artifacts."""

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
RUNS = ROOT / "benchmarks/codex_paper_eval/runs"
IDS = (
    "20260905_mechanism10_candidate",
    "20260905_mechanism4_final",
    "20260905_repository1_final",
)


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    pricing = read(ROOT / "benchmarks/evaluation_tiers/cost_calibration.json")
    rates = pricing["model_rates_per_million_tokens"]
    runs = []
    for rid in IDS:
        directory = RUNS / rid
        lock = read(directory / "LOCK.json")
        for name in ("selection", "predictions"):
            assert sha(directory / f"{name}.json") == lock[f"{name}_sha256"]
        report = read(directory / "report.json")
        predictions = read(directory / "predictions.json")
        cost = 0
        for model, usage in predictions["token_usage"]["models"].items():
            key = next(k for k in rates if k.lower() in model.lower())
            cost += (
                usage["input_tokens"] * rates[key]["input"]
                + usage["output_tokens"] * rates[key]["output_or_reasoning"]
            ) / 1_000_000
        runs.append(
            {
                "run_id": rid,
                "overall": report["overall"],
                "cost_proxy_usd": cost,
                "calls": predictions["token_usage"]["llm_calls"],
                "papers": [
                    {k: p[k] for k in ("gene", "pmid", "tp", "fp", "fn", "count")}
                    for p in report["papers"]
                ],
                "sha256": {
                    name: sha(directory / name)
                    for name in (
                        "LOCK.json",
                        "selection.json",
                        "predictions.json",
                        "report.json",
                        "analysis_setup.json",
                        "frozen_corpus/source_snapshot.json",
                    )
                },
            }
        )
    reference = {
        (p["gene"], p["pmid"]): p
        for p in read(RUNS / IDS[0] / "selection.json")["papers"]
    }
    pairs = []
    for p in read(RUNS / IDS[1] / "selection.json")["papers"]:
        old = reference[(p["gene"], p["pmid"])]
        pairs.append(
            {
                "gene": p["gene"],
                "pmid": p["pmid"],
                "actual_rendering_hash_equal": p["source_sha256"]
                == old["source_sha256"],
                "prototype_sha256": old["source_sha256"],
                "final_sha256": p["source_sha256"],
            }
        )
    prior = read(OUT.parent / "phenotype_failure_panel_20260905/manifest.json")
    for path, expected in prior["input_sha256"].items():
        assert sha(ROOT / path) == expected, path
    result = {
        "classification": "opened calibration; no registered acceptance decision",
        "runs": runs,
        "matched_actual_sources": pairs,
        "prior_locked_inputs_verified": len(prior["input_sha256"]),
        "pricing_as_of": pricing["pricing_as_of"],
        "test_cost_proxy_usd": sum(r["cost_proxy_usd"] for r in runs),
    }
    (OUT / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    budget = read(OUT / "budget.json")
    budget["test_spend_proxy_usd"] = result["test_cost_proxy_usd"]
    budget["completed_tests"] = [
        {k: r[k] for k in ("run_id", "cost_proxy_usd", "calls")} for r in runs
    ]
    budget["next_authorized_step"] = None
    budget["remaining_uncommitted_test_envelope_usd"] = (
        budget["ceiling_usd"] - result["test_cost_proxy_usd"]
    )
    (OUT / "budget.json").write_text(json.dumps(budget, indent=2) + "\n")
    print(json.dumps(budget, indent=2))


if __name__ == "__main__":
    main()
