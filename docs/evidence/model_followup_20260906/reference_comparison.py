"""Post-lock diagnostic to existing gold; no source-count acceptance by gold match."""

import json
import sys
import time
from pathlib import Path
from client import HERE, sha, write

ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
from benchmarks.codex_paper_eval.run_eval import load_gold, twin_identical
from score import norm

assert (HERE / "outputs_locked.json").exists() and (
    HERE / "grok_outputs_locked.json"
).exists()
scores = json.loads((HERE / "source_scores.json").read_text())
prior = (
    ROOT / "benchmarks/codex_paper_eval/runs/20260906_model12_grok43/predictions.json"
)
baseline = {(p["gene"], p["pmid"]): p for p in json.loads(prior.read_text())["papers"]}
goldroot = ROOT / "benchmarks/evaluation_tiers/mixed_gold_continuation_120/answer_key"
rows = []
# Source-exact numerical candidates are identical across the three arms, including
# the extra DOC representation. Use one copy; do not pool that repeated table.
for r in scores["results"]:
    if (
        r["arm"] != "astra_low"
        or not r["reference_people"]
        or r["packet"] == "myb204_structured"
    ):
        continue
    packet = json.loads((HERE / "packets" / (r["packet"] + ".json")).read_text())
    gene, pmid = packet["gene"], packet["pmid"]
    gold = load_gold(goldroot, gene, pmid)
    refs = json.loads((HERE / "source_reference.json").read_text())["rows"][r["packet"]]
    spelling = {norm(v): v for person in refs for v in person["variants"]}
    for v, counts in r["counts"].items():
        label = spelling.get(v, v)
        gs = [g for g in gold if twin_identical(label, g["variant"], gene)]
        bs = [
            b
            for b in baseline[gene, pmid]["variants"]
            if twin_identical(label, b["variant"], gene)
        ]
        row = dict(
            gene=gene,
            pmid=pmid,
            source_variant=label,
            counts=counts,
            gold_matches=gs,
            baseline_matches=bs,
            matching="unique strict equivalent allele only; ambiguous/unmapped retained",
            endpoint_comparable_for_AU=pmid != "30403697",
            source_qualified_variant=r[
                "source_adjudicated_complete_variant_counts"
            ].get(v)
            is not None,
        )
        comparisons = []
        if len(gs) == 1:
            for sourcefield, field in [
                ("carriers", "carriers"),
                ("positive", "affected"),
                ("negative", "unaffected"),
            ]:
                if field != "carriers" and pmid == "30403697":
                    continue
                g = gs[0][field]
                b = bs[0].get(field) if len(bs) == 1 else None
                comparisons.append(
                    dict(
                        field=field,
                        source_value=counts[sourcefield],
                        gold=g,
                        baseline=b,
                        gold_exact=counts[sourcefield] == g,
                        baseline_exact=b == g,
                        baseline_missing=b is None,
                        source_qualified=row["source_qualified_variant"],
                        interpretation="Packet/cohort/table endpoint may still differ from existing gold; inspect disagreements.",
                    )
                )
        row["comparisons"] = comparisons
        rows.append(row)
cs = [c for r in rows for c in r["comparisons"]]
summary = {
    field: dict(
        matched=sum(c["field"] == field for c in cs),
        source_exact=sum(c["field"] == field and c["gold_exact"] for c in cs),
        baseline_exact=sum(c["field"] == field and c["baseline_exact"] for c in cs),
        new_exact_fills=sum(
            c["field"] == field and c["gold_exact"] and c["baseline_missing"]
            for c in cs
        ),
        new_disagreeing_fills=sum(
            c["field"] == field and not c["gold_exact"] and c["baseline_missing"]
            for c in cs
        ),
    )
    for field in ["carriers", "affected", "unaffected"]
}
write(
    HERE / "existing_gold_diagnostic.json",
    dict(
        graded_unix=time.time(),
        source_score_sha256=sha(HERE / "source_scores.json"),
        baseline_sha256=sha(prior),
        classification="Post-lock source candidate overlap, strict matching; not whole-paper precision/recall, no production merge",
        rows=rows,
        summary=summary,
        unmapped_or_ambiguous=[
            (r["gene"], r["pmid"], r["source_variant"])
            for r in rows
            if len(r["gold_matches"]) != 1
        ],
    ),
)
print(json.dumps(summary, indent=2))
print(
    "unmapped",
    [(r["pmid"], r["source_variant"]) for r in rows if len(r["gold_matches"]) != 1],
)
print(
    "disagreements",
    [
        (r["pmid"], r["source_variant"], c)
        for r in rows
        for c in r["comparisons"]
        if not c["gold_exact"]
    ],
)
