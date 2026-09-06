"""Conservatively correct Azure Grok's separate reasoning counter in old hook rows.
Runs alongside the frozen extraction hook; never changes scientific runtime.
"""

import sys
import os
import time
from pathlib import Path

sys.path[:0] = [
    str(Path.cwd()),
    str(Path.cwd() / "docs/evidence/model_routing_20260906"),
]
from budget_guard import update, RATES

os.environ["GVF_EXPERIMENT_BUDGET"] = str(
    Path.cwd() / "docs/evidence/model_routing_20260906/budget.json"
)


def reconcile(data):
    data["additional_accounting_uncertainty_reserve_note"] = (
        "Smoke reserve includes $20 extra headroom while normalizing Grok separate reasoning token counters. Corrected every second, with no scientific runtime changes."
    )
    data["smoke_uncertainty_reserve_usd"] = 25
    for row in data["calls"]:
        if "gpt-6-astra" in row["model"] and row["status"] != "returned":
            original = row.setdefault("original_reserved_usd", row["reserved_usd"])
            bound_in = (original * 1e6 / 6 - row["output_cap"] * 50) / 10
            if bound_in <= 272000:
                # A byte-bound below the long-context threshold cannot incur
                # that premium. Still reserve THREE complete SDK attempts.
                row["reserved_usd"] = original / 2
                row["reservation_basis"] = (
                    "three full attempts at standard rates; input byte bound below Astra long-context threshold"
                )
        if "grok-" not in row["model"] or row["status"] != "returned":
            continue
        usage = row["usage"]
        inp = usage.get("prompt_tokens", usage.get("input_tokens", 0))
        raw = usage.get("completion_tokens", usage.get("output_tokens", 0))
        out = max(raw, usage.get("total_tokens", inp + raw) - inp)
        rin, rout = next(v for k, v in RATES.items() if k in row["model"].lower())
        cost = (inp * rin + out * rout) / 1e6 * (2 if inp > 200000 else 1)
        row["reported_completion_tokens"] = raw
        row["billable_output_tokens_conservative"] = out
        row["api_proxy_usd"] = max(cost, row.get("api_proxy_usd", 0))
        row["accounted_usd"] = max(cost, row.get("accounted_usd", 0))


if __name__ == "__main__":
    while True:
        update(reconcile)
        time.sleep(1)
