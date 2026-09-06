"""Run the medium-effort arm one gene at a time within live reservation headroom."""

import json
import time
from campaign import HERE, PYTHON, check, environment, run

folder = check("astra_medium")
env = environment("astra_medium")
for gene in ("MYBPC3", "RYR2", "SCN5A"):
    assert not (folder / "production_runs" / gene).exists(), (
        f"Refusing duplicate gene run: {gene}"
    )
    while True:
        data = json.loads((HERE / "budget.json").read_text())
        committed = (
            sum(row.get("accounted_usd", row["reserved_usd"]) for row in data["calls"])
            + data["smoke_uncertainty_reserve_usd"]
        )
        if committed + 25 <= data["api_ceiling_usd"]:
            break
        time.sleep(5)
    check("astra_medium")
    print("START astra_medium", gene, flush=True)
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
            folder / "pmids" / f"{gene}.txt",
            "--no-source-recovery",
            "--no-corpus-sync",
            "--no-publish-review",
            "--gold-free-run",
        ],
        folder / "operator_logs" / f"{gene}.log",
        env,
    )
    print("COMPLETED astra_medium", gene, flush=True)
check("astra_medium")
