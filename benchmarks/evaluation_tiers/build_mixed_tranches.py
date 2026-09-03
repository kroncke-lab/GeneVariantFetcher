#!/usr/bin/env python3
"""Build the complete, source-available gold inventory into mixed tranches.

The builder reads gold only to determine gene/PMID membership.  It never uses a
variant identity, count value, or gold row count to rank or balance a tranche.
Every source-available gene-paper attempt is assigned exactly once, and all
attempts for the same PMID stay in one tranche so a paper cannot leak across
protocol comparisons under another gene.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from utils.gold_standard import gold_row_excluded  # noqa: E402

DEFAULT_OUT = HERE / "mixed_gold"
DEFAULT_CORPUS = REPO / "corpus"
DEFAULT_COST = HERE / "cost_calibration.json"
REPO_GOLD = REPO / "gene_variant_fetcher_gold_standard" / "normalized"
OVERRIDE_GOLD = REPO / "benchmarks" / "curated_extraction_eval" / "gold_overrides"

QUARANTINE = {
    ("KCNH2", "10086972"): "not a genetics paper; gold row belongs to 10086971",
    (
        "KCNH2",
        "14642689",
    ): "angiotensin-II receptor expression study, not KCNH2 genetics",
    ("BRCA2", "19944633"): "non-human canine BRCA2 ortholog; negative control only",
}
LEAD_APPROVED_BRCA2 = {"26833046", "26848529"}


class BuildError(RuntimeError):
    """The mixed suite cannot be built reproducibly from the current inputs."""


def sha256_file(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def stable_rank(seed: int, *parts: object) -> str:
    material = "|".join([str(seed), *(str(part) for part in parts)])
    return hashlib.sha256(material.encode()).hexdigest()


def provenance_for(source_kind: str, gene: str, pmid: str) -> str:
    if source_kind == "repo_gold":
        return "human_curated_ryr2_pilot" if gene == "RYR2" else "human_curated_cardiac"
    if gene == "BRCA2" and pmid in LEAD_APPROVED_BRCA2:
        return "collaborator_approved_nonexhaustive"
    return "curated_override_mixed_provenance"


def gold_inputs() -> list[tuple[str, str, Path]]:
    inputs: list[tuple[str, str, Path]] = []
    for path in sorted(REPO_GOLD.glob("*_recall_input.csv")):
        inputs.append(("repo_gold", path.name.removesuffix("_recall_input.csv"), path))
    for path in sorted(OVERRIDE_GOLD.glob("*_recall_input.csv")):
        gene = path.name.removesuffix("_recall_input.csv")
        if gene.startswith("_"):
            continue
        inputs.append(("curated_override", gene, path))
    return inputs


def read_gold_inventory() -> tuple[
    dict[str, list[dict[str, str]]],
    dict[tuple[str, str], dict[str, str]],
    list[dict[str, Any]],
]:
    """Merge answer keys with repo gold taking precedence per gene/PMID."""

    rows_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    attempt_meta: dict[tuple[str, str], dict[str, str]] = {}
    input_meta: list[dict[str, Any]] = []
    for source_kind, gene, path in gold_inputs():
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            fieldnames = reader.fieldnames or []
            rows = list(reader)
        if not {"variant", "pmid"} <= set(fieldnames):
            raise BuildError(f"gold input lacks variant/pmid columns: {path}")
        source_attempts: set[tuple[str, str]] = set()
        source_name = str(path.relative_to(REPO))
        for row in rows:
            if gold_row_excluded(row):
                continue
            pmid = str(row.get("pmid") or "").strip()
            variant = str(row.get("variant") or "").strip()
            if not pmid.isdigit() or not variant:
                continue
            key = (gene, pmid)
            # The complete repo gold is authoritative when an override happens
            # to repeat one of its paper attempts.  Overrides add only papers a
            # gene's repo key does not already contain.
            if (
                key in attempt_meta
                and attempt_meta[key]["gold_source"] != source_name
                and source_kind == "curated_override"
            ):
                continue
            if key not in attempt_meta:
                attempt_meta[key] = {
                    "gene": gene,
                    "pmid": pmid,
                    "gold_provenance": provenance_for(source_kind, gene, pmid),
                    "gold_source": source_name,
                }
            source_attempts.add(key)
            rows_by_gene[gene].append(row)
        input_meta.append(
            {
                "kind": source_kind,
                "gene": gene,
                "path": str(path.relative_to(REPO)),
                "sha256": sha256_file(path),
                "eligible_attempts_contributed": len(source_attempts),
            }
        )
    if not attempt_meta:
        raise BuildError("no named gold variant attempts found")
    return rows_by_gene, attempt_meta, input_meta


def usable_source(
    corpus: Path, gene: str, pmid: str, minimum_chars: int
) -> Path | None:
    directory = corpus / gene / pmid
    candidates = [
        directory / f"{pmid}_FULL_CONTEXT.md",
        directory / f"{pmid}_CLEANED.md",
    ]
    usable = [
        path
        for path in candidates
        if path.is_file() and path.stat().st_size >= minimum_chars
    ]
    if not usable:
        return None
    # The run harness may choose a richer CLEANED rendering only when it Pareto
    # dominates FULL_CONTEXT.  The tranche inventory needs presence, not that
    # extraction-time choice; selection.json binds the exact chosen bytes later.
    return usable[0]


def model_rate(model: str, rates: dict[str, dict[str, float]]) -> dict[str, float]:
    matches = [rate for key, rate in rates.items() if key in model.lower()]
    if len(matches) != 1:
        raise BuildError(
            f"expected one price match for model {model!r}, got {len(matches)}"
        )
    return matches[0]


def cost_rates(path: Path) -> tuple[dict[str, float], dict[str, Any]]:
    profile = json.loads(path.read_text())
    rates = profile["model_rates_per_million_tokens"]
    unit_costs: dict[str, float] = {}
    for gene, calibration in profile["calibrations"].items():
        attempts = int(calibration["attempts"])
        if attempts <= 0:
            raise BuildError(f"{gene}: invalid cost calibration attempt count")
        cost = 0.0
        for model, usage in calibration["models"].items():
            rate = model_rate(model, rates)
            cost += (
                int(usage["input_tokens"]) * float(rate["input"])
                + int(usage["output_tokens"]) * float(rate["output_or_reasoning"])
            ) / 1_000_000
        unit_costs[gene] = cost / attempts
    return unit_costs, profile


def partition_attempts(
    attempts: list[dict[str, str]], seed: int, target_size: int
) -> list[list[dict[str, str]]]:
    """Article-atomic, gene/provenance-balanced deterministic bin packing."""

    by_pmid: dict[str, list[dict[str, str]]] = defaultdict(list)
    for attempt in attempts:
        by_pmid[attempt["pmid"]].append(attempt)
    groups = sorted(
        by_pmid.values(),
        key=lambda rows: (-len(rows), stable_rank(seed, rows[0]["pmid"])),
    )
    count = math.ceil(len(attempts) / target_size)
    tranches: list[list[dict[str, str]]] = [[] for _ in range(count)]
    genes = [Counter() for _ in range(count)]
    provenances = [Counter() for _ in range(count)]
    for group in groups:
        group_genes = Counter(row["gene"] for row in group)
        group_provenances = Counter(row["gold_provenance"] for row in group)
        candidates = [
            index
            for index, tranche in enumerate(tranches)
            if len(tranche) + len(group) <= target_size
        ]
        if not candidates:
            raise BuildError("article-atomic packing exceeded the tranche size")

        def placement_score(index: int) -> tuple[int, int, int, str]:
            gene_pressure = sum(genes[index][gene] for gene in group_genes)
            provenance_pressure = sum(
                provenances[index][name] for name in group_provenances
            )
            return (
                gene_pressure,
                provenance_pressure,
                len(tranches[index]),
                stable_rank(seed, "bin", index, group[0]["pmid"]),
            )

        destination = min(candidates, key=placement_score)
        tranches[destination].extend(group)
        genes[destination].update(group_genes)
        provenances[destination].update(group_provenances)

    for index, tranche in enumerate(tranches):
        tranche.sort(
            key=lambda row: stable_rank(seed, "order", index, row["gene"], row["pmid"])
        )
    return tranches


def write_answer_key(
    root: Path,
    rows_by_gene: dict[str, list[dict[str, str]]],
    attempts: dict[tuple[str, str], dict[str, str]],
) -> None:
    root.mkdir(parents=True, exist_ok=True)
    standard_fields = [
        "variant",
        "pmid",
        "carriers",
        "affected",
        "unaffected",
        "gold_v2_carriers",
        "gold_v2_affected",
        "gold_v2_unaffected",
        "gold_v2_status",
        "gold_v2_note",
        "gold_v2_source",
    ]
    for gene, rows in sorted(rows_by_gene.items()):
        rows.sort(
            key=lambda row: (str(row.get("pmid") or ""), str(row.get("variant") or ""))
        )
        target = root / f"{gene}_recall_input.csv"
        with target.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=standard_fields, lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(
                {key: row.get(key, "") for key in standard_fields} for row in rows
            )
    with (root / "provenance.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["gene", "pmid", "gold_provenance", "gold_source"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(attempts[key] for key in sorted(attempts))


def digest_answer_key(root: Path) -> str:
    value = hashlib.sha256()
    files = sorted(root.glob("*_recall_input.csv")) + [root / "provenance.tsv"]
    for path in files:
        value.update(path.name.encode())
        value.update(b"\0")
        value.update(sha256_file(path).encode())
        value.update(b"\n")
    return value.hexdigest()


def build(args: argparse.Namespace) -> dict[str, Any]:
    rows_by_gene, attempt_meta, input_meta = read_gold_inventory()
    args.out.mkdir(parents=True, exist_ok=True)
    answer_key = args.out / "answer_key"
    write_answer_key(answer_key, rows_by_gene, attempt_meta)

    included: list[dict[str, str]] = []
    inventory: list[dict[str, str]] = []
    for key in sorted(attempt_meta):
        meta = dict(attempt_meta[key])
        source = usable_source(args.corpus, *key, args.minimum_chars)
        if key in QUARANTINE:
            status = "quarantined"
            reason = QUARANTINE[key]
        elif source is None:
            status = "source_unavailable"
            reason = f"no FULL_CONTEXT/CLEANED source >= {args.minimum_chars} bytes"
        else:
            status = "included"
            reason = ""
            included.append(meta)
        inventory.append(
            {
                **meta,
                "status": status,
                "reason": reason,
                "source": (str(source.relative_to(REPO)) if source is not None else ""),
            }
        )

    tranches = partition_attempts(included, args.seed, args.tranche_size)
    unit_costs, cost_profile = cost_rates(args.cost_calibration)
    missing_cost = sorted({row["gene"] for row in included} - set(unit_costs))
    if missing_cost:
        raise BuildError(f"no observed cost calibration for genes {missing_cost}")

    for stale in args.out.glob("tranche_*.tsv"):
        stale.unlink()
    gold_root_digest = digest_answer_key(answer_key)
    tier_records = []
    for number, tranche in enumerate(tranches, 1):
        tier_id = f"mixed_gold_{number:02d}"
        manifest = args.out / f"tranche_{number:02d}.tsv"
        lines = [
            f"# {tier_id}: mixed all-gold protocol tranche; seed {args.seed}.",
            "# Article-disjoint; variant identities/counts never participate in assignment.",
            "# Primary score: paper-derived. ClinVar/PubTator: secondary linkage lane.",
            *(f"{row['gene']}\t{row['pmid']}" for row in tranche),
        ]
        manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")
        estimate = sum(unit_costs[row["gene"]] for row in tranche)
        tier_records.append(
            {
                "id": tier_id,
                "order": number,
                "manifest": manifest.name,
                "sha256": sha256_file(manifest),
                "attempt_count": len(tranche),
                "unique_pmid_count": len({row["pmid"] for row in tranche}),
                "gene_attempt_counts": dict(
                    sorted(Counter(row["gene"] for row in tranche).items())
                ),
                "gold_provenance_attempt_counts": dict(
                    sorted(Counter(row["gold_provenance"] for row in tranche).items())
                ),
                "selection_seed": args.seed,
                "eligibility_mode": "variant",
                "gold_root": "answer_key",
                "gold_root_sha256": gold_root_digest,
                "primary_score_lane": "paper_derived",
                "comparison_score_lanes": ["linkage_assisted"],
                "role": "mixed_gold_protocol_regression",
                "estimated_cost_usd": round(estimate, 4),
                "budget_with_headroom_usd": round(
                    estimate * float(cost_profile["budget_headroom_multiplier"]), 4
                ),
                "paired_estimated_cost_usd": round(estimate * 2, 4),
                "paired_budget_with_headroom_usd": round(
                    estimate * 2 * float(cost_profile["budget_headroom_multiplier"]),
                    4,
                ),
            }
        )

    inventory_path = args.out / "inventory.tsv"
    with inventory_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "gene",
                "pmid",
                "gold_provenance",
                "gold_source",
                "status",
                "reason",
                "source",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(inventory)

    estimated_total = sum(float(tier["estimated_cost_usd"]) for tier in tier_records)
    consumption_log = args.out / "consumption_log.jsonl"
    if not consumption_log.exists():
        consumption_log.write_text("", encoding="utf-8")
    registry = {
        "schema_version": 1,
        "as_of": "2026-09-03",
        "suite_id": "mixed_all_gold",
        "count_unit": "gene_paper_attempt",
        "selection_seed": args.seed,
        "target_tranche_size": args.tranche_size,
        "assignment": "article_atomic_gene_and_provenance_balanced_sha256_v1",
        "gold_value_blinded": True,
        "primary_score_lane": "paper_derived",
        "comparison_score_lanes": ["linkage_assisted"],
        "consumption_log": {
            "path": consumption_log.name,
            "schema_version": 1,
            "required_arms": ["baseline", "candidate"],
        },
        "evaluation_design": {
            "comparison": "paired_frozen_baseline_and_candidate_on_same_manifest",
            "primary_endpoint": "paper_derived_micro_variant_identity_recall",
            "guardrails": [
                "paper_derived_micro_variant_identity_precision_noninferiority",
                "per_gene_results_reported",
                "per_gold_provenance_results_reported_without_scientific_pooling",
            ],
            "consume_order": [tier["id"] for tier in tier_records],
            "burn_policy": (
                "Any score inspected before a protocol decision burns that tranche "
                "for confirmation; confirm on the next still-unopened tranche."
            ),
            "cache_policy": "frozen_source_bytes_only_no_prior_prediction_or_db_reuse",
            "cluster_unit": "PMID",
            "decision_rule": {
                "delta_definition": "candidate_minus_baseline_on_the_same_tranche",
                "primary": {
                    "metric": "paper_derived_micro_variant_identity_recall",
                    "minimum_observed_delta": 0.01,
                    "noninferiority_margin": -0.01,
                    "one_sided_confidence_level": 0.95,
                },
                "precision_guardrail": {
                    "metric": "paper_derived_micro_variant_identity_precision",
                    "noninferiority_margin": -0.02,
                    "one_sided_confidence_level": 0.95,
                },
                "confidence_interval": {
                    "method": "paired_cluster_bootstrap_nearest_rank",
                    "cluster_unit": "PMID",
                    "resamples": 10000,
                    "seed": args.seed,
                },
                "pass_condition": (
                    "Observed recall delta >= +0.01 and its one-sided 95% lower "
                    "bound >= -0.01, while the precision-delta one-sided 95% "
                    "lower bound is >= -0.02."
                ),
                "confirmation": (
                    "A discovery pass advances only to the next unopened tranche; "
                    "accept the candidate only if the same frozen rule passes there."
                ),
            },
        },
        "inventory": {
            "path": inventory_path.name,
            "sha256": sha256_file(inventory_path),
            "gold_attempts": len(attempt_meta),
            "source_available_attempts": len(included),
            "source_unavailable_attempts": sum(
                row["status"] == "source_unavailable" for row in inventory
            ),
            "quarantined_attempts": sum(
                row["status"] == "quarantined" for row in inventory
            ),
            "unique_source_available_pmids": len({row["pmid"] for row in included}),
        },
        "gold_inputs": input_meta,
        "answer_key": {
            "path": "answer_key",
            "sha256": gold_root_digest,
        },
        "cost_model": {
            "calibration": str(args.cost_calibration.relative_to(REPO)),
            "calibration_sha256": sha256_file(args.cost_calibration),
            "method": "sum of observed per-gene production cost per gene-paper attempt",
            "pricing_basis": cost_profile["pricing_basis"],
            "pricing_as_of": cost_profile["pricing_as_of"],
            "estimated_suite_cost_usd": round(estimated_total, 4),
            "budget_with_headroom_usd": round(
                estimated_total * float(cost_profile["budget_headroom_multiplier"]),
                4,
            ),
            "paired_estimated_suite_cost_usd": round(estimated_total * 2, 4),
            "paired_budget_with_headroom_usd": round(
                estimated_total * 2 * float(cost_profile["budget_headroom_multiplier"]),
                4,
            ),
            "unit_cost_usd_by_gene": {
                gene: round(cost, 6) for gene, cost in sorted(unit_costs.items())
            },
        },
        "tiers": tier_records,
    }
    (args.out / "registry.json").write_text(
        json.dumps(registry, indent=2) + "\n", encoding="utf-8"
    )
    readme_lines = [
        "# Mixed all-gold protocol tranches",
        "",
        (
            f"This seeded suite assigns all **{len(included):,} source-available** "
            f"gene-paper attempts from the repository's **{len(attempt_meta):,}** "
            "named-variant gold inventory to exactly one of "
            f"**{len(tier_records)}** mixed tranches."
        ),
        "",
        (
            "The primary score is `paper_derived`. Rows originating in ClinVar or "
            "PubTator remain locked as `external_linkage_variants` and contribute "
            "only to the secondary `linkage_assisted` diagnostic. A database "
            "citation is therefore never counted as evidence that the protocol "
            "found a variant in the paper. Ambiguous `mixed`, legacy, and unknown "
            "origins are locked in `unattributed_variants` but excluded from both "
            "scored lanes."
        ),
        "",
        (
            f"Assignment is deterministic from seed `{args.seed}`, article-atomic "
            "(the same PMID cannot appear in different tranches under different "
            "genes), and balanced by gene and gold provenance without reading "
            "variant identities, count values, or gold row counts. `inventory.tsv` "
            f"records the {registry['inventory']['source_unavailable_attempts']} "
            "attempts without usable local source and the "
            f"{registry['inventory']['quarantined_attempts']} quarantined "
            "wrong-paper attempt; they are not silently treated as extraction "
            "failures."
        ),
        "",
        (
            "Gold provenance remains a reporting stratum even though workloads are "
            "mixed. Do not pool `human_curated_cardiac`, the RYR2 spreadsheet "
            "pilot, lead-approved non-exhaustive BRCA2 records, and "
            "mixed-provenance curated overrides into one scientific headline. "
            "`run_eval.py` emits `by_gold_provenance` for that reason."
        ),
        "",
        "## Running one tranche",
        "",
        "```bash",
        ".venv/bin/python benchmarks/codex_paper_eval/setup_production_eval.py create \\",
        "  --tier-id mixed_gold_01 \\",
        "  --paper-manifest benchmarks/evaluation_tiers/mixed_gold/tranche_01.tsv \\",
        "  --registry benchmarks/evaluation_tiers/mixed_gold/registry.json \\",
        f"  --seed {args.seed} --comparison-arm baseline \\",
        "  --run-id YYYYMMDD_protocol_mixed01_baseline \\",
        "  --email brett.kroncke@gmail.com",
        "```",
        "",
        (
            "The generated extraction remains gold-free and nonpublishing. The "
            "registry resolves the composite answer key only for the post-lock "
            "score. Reuse frozen source bytes, never prior predictions, extraction "
            "JSON, or databases. Compare the frozen baseline and candidate on the "
            "same manifest; attempts are clustered by PMID. Once a tranche's score "
            "has been inspected, treat it as calibration and confirm a change on "
            "the next still-unopened tranche in registry order."
        ),
        "",
        (
            "For this paired suite, create and score `--comparison-arm baseline` "
            "first, then create `--comparison-arm candidate` after the protocol "
            "change. The scorer appends both openings to `consumption_log.jsonl`; "
            "setup refuses a repeated arm or an out-of-order tranche."
        ),
        "",
        "## Preregistered paired decision",
        "",
        (
            "Deltas are candidate minus baseline on the same tranche. A change "
            "passes discovery only when observed paper-derived recall improves by "
            "at least 1 percentage point, the one-sided 95% PMID-cluster-bootstrap "
            "lower bound for recall is at least -1 point, and the corresponding "
            "precision lower bound is at least -2 points. A discovery pass is not "
            "acceptance: the identical rule must pass on the next unopened tranche. "
            "These are preregistered engineering thresholds, not clinical-effect "
            "thresholds. Use `run_eval.py compare` to apply the registry rule."
        ),
        "",
        (
            "Gold-row denominators vary from 110 to 593 per tranche and rare "
            "non-cardiac provenance strata are sparse. Paired inference stays "
            "within one tranche; per-provenance values must show sample size and "
            "remain diagnostic. Absolute pooled precision is not an exhaustive "
            "estimate for non-exhaustive collaborator gold. Until an append-only "
            "abandonment event exists, do not open a baseline unless its candidate "
            "arm will be completed."
        ),
        "",
        "## Per-tranche cost estimate",
        "",
        (
            "These are planning proxies, not invoices. Each estimate sums observed "
            "per-gene production token cost per gene-paper attempt from "
            "`../cost_calibration.json` using the dated 2026-08-24 public price "
            "card. The primary protocol comparison is paired, so it costs two arms; "
            "headroom adds 25% for retry/variance. Source acquisition and human "
            "review are excluded."
        ),
        "",
        "| Tranche | Attempts | Genes | One arm | Paired A/B | Paired +25% |",
        "|---|---:|---|---:|---:|---:|",
    ]
    for tier in tier_records:
        genes = ", ".join(
            f"{gene} {count}" for gene, count in tier["gene_attempt_counts"].items()
        )
        readme_lines.append(
            f"| `{tier['id']}` | {tier['attempt_count']} | {genes} | "
            f"${tier['estimated_cost_usd']:.2f} | "
            f"${tier['paired_estimated_cost_usd']:.2f} | "
            f"${tier['paired_budget_with_headroom_usd']:.2f} |"
        )
    readme_lines.extend(
        [
            f"| **Suite** | **{len(included)}** |  | "
            f"**${estimated_total:.2f}** | "
            f"**${registry['cost_model']['paired_estimated_suite_cost_usd']:.2f}** | "
            f"**${registry['cost_model']['paired_budget_with_headroom_usd']:.2f}** |",
            "",
        ]
    )
    (args.out / "README.md").write_text("\n".join(readme_lines), encoding="utf-8")
    return registry


def parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=2026090301)
    ap.add_argument("--tranche-size", type=int, default=50)
    ap.add_argument("--minimum-chars", type=int, default=2000)
    ap.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS)
    ap.add_argument("--cost-calibration", type=Path, default=DEFAULT_COST)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    return ap


def main() -> int:
    args = parser().parse_args()
    try:
        registry = build(args)
    except BuildError as exc:
        print(f"error: {exc}")
        return 2
    inventory = registry["inventory"]
    print(
        f"{len(registry['tiers'])} tranches; "
        f"{inventory['source_available_attempts']}/{inventory['gold_attempts']} "
        "gold attempts source-available; estimated "
        f"${registry['cost_model']['estimated_suite_cost_usd']:.2f} "
        f"(${registry['cost_model']['budget_with_headroom_usd']:.2f} with headroom)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
