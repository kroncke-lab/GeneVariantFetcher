# Mixed all-gold protocol tranches

This seeded suite assigns all **1,422 source-available** gene-paper attempts from the repository's **1,534** named-variant gold inventory to exactly one of **29** mixed tranches.

The primary score is `paper_derived`. Rows originating in ClinVar or PubTator remain locked as `external_linkage_variants` and contribute only to the secondary `linkage_assisted` diagnostic. A database citation is therefore never counted as evidence that the protocol found a variant in the paper. Ambiguous `mixed`, legacy, and unknown origins are locked in `unattributed_variants` but excluded from both scored lanes.

Assignment is deterministic from seed `2026090301`, article-atomic (the same PMID cannot appear in different tranches under different genes), and balanced by gene and gold provenance without reading variant identities, count values, or gold row counts. `inventory.tsv` records the 111 attempts without usable local source and the 1 quarantined wrong-paper attempt; they are not silently treated as extraction failures.

Gold provenance remains a reporting stratum even though workloads are mixed. Do not pool `human_curated_cardiac`, the RYR2 spreadsheet pilot, lead-approved non-exhaustive BRCA2 records, and mixed-provenance curated overrides into one scientific headline. `run_eval.py` emits `by_gold_provenance` for that reason.

## Running one tranche

```bash
.venv/bin/python benchmarks/codex_paper_eval/setup_production_eval.py create \
  --tier-id mixed_gold_01 \
  --paper-manifest benchmarks/evaluation_tiers/mixed_gold/tranche_01.tsv \
  --registry benchmarks/evaluation_tiers/mixed_gold/registry.json \
  --seed 2026090301 --comparison-arm baseline \
  --run-id YYYYMMDD_protocol_mixed01_baseline \
  --email brett.kroncke@gmail.com
```

The generated extraction remains gold-free and nonpublishing. The registry resolves the composite answer key only for the post-lock score. Reuse frozen source bytes, never prior predictions, extraction JSON, or databases. Compare the frozen baseline and candidate on the same manifest; attempts are clustered by PMID. Once a tranche's score has been inspected, treat it as calibration and confirm a change on the next still-unopened tranche in registry order.

For this paired suite, create and score `--comparison-arm baseline` first, then create `--comparison-arm candidate` after the protocol change. The scorer appends both openings to `consumption_log.jsonl`; setup refuses a repeated arm or an out-of-order tranche.

## Preregistered paired decision

Deltas are candidate minus baseline on the same tranche. A change passes discovery only when observed paper-derived recall improves by at least 1 percentage point, the one-sided 95% PMID-cluster-bootstrap lower bound for recall is at least -1 point, and the corresponding precision lower bound is at least -2 points. A discovery pass is not acceptance: the identical rule must pass on the next unopened tranche. These are preregistered engineering thresholds, not clinical-effect thresholds. Use `run_eval.py compare` to apply the registry rule.

Gold-row denominators vary from 110 to 593 per tranche and rare non-cardiac provenance strata are sparse. Paired inference stays within one tranche; per-provenance values must show sample size and remain diagnostic. Absolute pooled precision is not an exhaustive estimate for non-exhaustive collaborator gold. Until an append-only abandonment event exists, do not open a baseline unless its candidate arm will be completed.

## Per-tranche cost estimate

These are planning proxies, not invoices. Each estimate sums observed per-gene production token cost per gene-paper attempt from `../cost_calibration.json` using the dated 2026-08-24 public price card. The primary protocol comparison is paired, so it costs two arms; headroom adds 25% for retry/variance. Source acquisition and human review are excluded.

| Tranche | Attempts | Genes | One arm | Paired A/B | Paired +25% |
|---|---:|---|---:|---:|---:|
| `mixed_gold_01` | 49 | KCNH2 9, KCNQ1 10, MYBPC3 1, RYR2 6, SCN5A 23 | $5.13 | $10.25 | $12.81 |
| `mixed_gold_02` | 49 | BRCA1 1, KCNH2 8, KCNQ1 10, RYR2 6, SCN5A 24 | $5.21 | $10.41 | $13.01 |
| `mixed_gold_03` | 49 | KCNH2 8, KCNQ1 10, MYBPC3 1, RYR2 6, SCN5A 24 | $5.15 | $10.29 | $12.87 |
| `mixed_gold_04` | 49 | APOE 1, KCNH2 9, KCNQ1 10, RYR2 6, SCN5A 23 | $4.99 | $9.97 | $12.47 |
| `mixed_gold_05` | 49 | KCNH2 9, KCNQ1 10, MYBPC3 1, RYR2 6, SCN5A 23 | $5.13 | $10.25 | $12.81 |
| `mixed_gold_06` | 49 | KCNH2 8, KCNQ1 10, MYBPC3 1, RYR2 6, SCN5A 24 | $5.15 | $10.29 | $12.87 |
| `mixed_gold_07` | 49 | APOE 1, BRCA2 1, KCNH2 8, KCNQ1 11, RYR2 5, SCN5A 23 | $5.11 | $10.22 | $12.77 |
| `mixed_gold_08` | 49 | BRCA2 1, KCNH2 8, KCNQ1 10, RYR2 6, SCN5A 24 | $5.14 | $10.27 | $12.84 |
| `mixed_gold_09` | 49 | BRCA2 1, KCNH2 9, KCNQ1 10, RYR2 6, SCN5A 23 | $5.12 | $10.23 | $12.79 |
| `mixed_gold_10` | 49 | BRCA2 1, KCNH2 9, KCNQ1 10, RYR2 6, SCN5A 23 | $5.12 | $10.23 | $12.79 |
| `mixed_gold_11` | 49 | APOE 1, KCNH2 9, KCNQ1 10, RYR2 6, SCN5A 23 | $4.99 | $9.97 | $12.47 |
| `mixed_gold_12` | 49 | APOE 1, BRCA1 1, KCNH2 9, KCNQ1 10, RYR2 5, SCN5A 23 | $5.18 | $10.35 | $12.94 |
| `mixed_gold_13` | 49 | BRCA1 1, KCNH2 9, KCNQ1 10, RYR2 6, SCN5A 23 | $5.18 | $10.37 | $12.96 |
| `mixed_gold_14` | 50 | BRCA2 1, KCNH2 9, KCNQ1 10, RYR2 6, SCN5A 24 | $5.22 | $10.45 | $13.06 |
| `mixed_gold_15` | 49 | BRCA1 1, BRCA2 1, KCNH2 9, KCNQ1 10, RYR2 5, SCN5A 23 | $5.30 | $10.61 | $13.26 |
| `mixed_gold_16` | 49 | KCNH2 9, KCNQ1 10, MYBPC3 1, RYR2 6, SCN5A 23 | $5.13 | $10.25 | $12.81 |
| `mixed_gold_17` | 49 | BRCA1 1, KCNH2 9, KCNQ1 11, RYR2 5, SCN5A 23 | $5.16 | $10.32 | $12.90 |
| `mixed_gold_18` | 49 | BRCA1 1, KCNH2 9, KCNQ1 10, RYR2 6, SCN5A 23 | $5.18 | $10.37 | $12.96 |
| `mixed_gold_19` | 49 | BRCA1 1, BRCA2 1, KCNH2 9, KCNQ1 10, RYR2 5, SCN5A 23 | $5.30 | $10.61 | $13.26 |
| `mixed_gold_20` | 49 | BRCA1 1, KCNH2 9, KCNQ1 10, RYR2 6, SCN5A 23 | $5.18 | $10.37 | $12.96 |
| `mixed_gold_21` | 49 | APOE 1, KCNH2 9, KCNQ1 10, RYR2 5, SCN5A 24 | $4.98 | $9.96 | $12.45 |
| `mixed_gold_22` | 49 | BRCA1 1, KCNH2 9, KCNQ1 10, RYR2 6, SCN5A 23 | $5.18 | $10.37 | $12.96 |
| `mixed_gold_23` | 49 | KCNH2 9, KCNQ1 10, MYBPC3 1, RYR2 6, SCN5A 23 | $5.13 | $10.25 | $12.81 |
| `mixed_gold_24` | 49 | BRCA2 1, KCNH2 9, KCNQ1 10, RYR2 6, SCN5A 23 | $5.12 | $10.23 | $12.79 |
| `mixed_gold_25` | 49 | KCNH2 8, KCNQ1 11, MYBPC3 1, RYR2 6, SCN5A 23 | $5.13 | $10.26 | $12.82 |
| `mixed_gold_26` | 49 | KCNH2 9, KCNQ1 10, RYR2 6, SCN5A 24 | $4.99 | $9.98 | $12.47 |
| `mixed_gold_27` | 49 | APOE 1, KCNH2 9, KCNQ1 10, RYR2 6, SCN5A 23 | $4.99 | $9.97 | $12.47 |
| `mixed_gold_28` | 49 | KCNH2 9, KCNQ1 10, MYBPC3 1, RYR2 5, SCN5A 24 | $5.12 | $10.24 | $12.80 |
| `mixed_gold_29` | 49 | BRCA1 1, KCNH2 9, KCNQ1 11, RYR2 5, SCN5A 23 | $5.16 | $10.32 | $12.90 |
| **Suite** | **1422** |  | **$148.83** | **$297.66** | **$372.08** |
