# Mixed all-gold protocol continuation

This seeded suite assigns all **1,324 source-available** gene-paper attempts from the repository's **1,534** named-variant gold inventory to exactly one of **11** mixed tranches.

This continuation excludes **98** gene-paper attempts across **90** PMIDs listed in the consumed manifests bound in `registry.json`. The prior locked runs remain immutable calibration evidence.

The primary score is `paper_derived`. Rows originating in ClinVar or PubTator remain locked as `external_linkage_variants` and contribute only to the secondary `linkage_assisted` diagnostic. A database citation is therefore never counted as evidence that the protocol found a variant in the paper. Ambiguous `mixed`, legacy, and unknown origins are locked in `unattributed_variants` but excluded from both scored lanes.

Assignment is deterministic from seed `2026090302`, article-atomic (the same PMID cannot appear in different tranches under different genes), and balanced by gene and gold provenance without reading variant identities, count values, or gold row counts. `inventory.tsv` records the 111 attempts without usable local source and the 1 quarantined wrong-paper attempt; they are not silently treated as extraction failures.

Gold provenance remains a reporting stratum even though workloads are mixed. Do not pool `human_curated_cardiac`, the RYR2 spreadsheet pilot, lead-approved non-exhaustive BRCA2 records, and mixed-provenance curated overrides into one scientific headline. `run_eval.py` emits `by_gold_provenance` for that reason.

## Running one tranche

```bash
.venv/bin/python benchmarks/codex_paper_eval/setup_production_eval.py create \
  --tier-id mixed_gold_cont120_01 \
  --paper-manifest benchmarks/evaluation_tiers/mixed_gold_continuation_120/tranche_01.tsv \
  --registry benchmarks/evaluation_tiers/mixed_gold_continuation_120/registry.json \
  --seed 2026090302 --comparison-arm baseline \
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
| `mixed_gold_cont120_01` | 120 | APOE 1, BRCA1 1, KCNH2 21, KCNQ1 25, MYBPC3 1, RYR2 14, SCN5A 57 | $12.56 | $25.11 | $31.39 |
| `mixed_gold_cont120_02` | 120 | APOE 1, BRCA2 1, KCNH2 21, KCNQ1 25, MYBPC3 1, RYR2 14, SCN5A 57 | $12.49 | $24.97 | $31.22 |
| `mixed_gold_cont120_03` | 120 | BRCA1 1, BRCA2 1, KCNH2 21, KCNQ1 25, MYBPC3 1, RYR2 14, SCN5A 57 | $12.68 | $25.37 | $31.71 |
| `mixed_gold_cont120_04` | 120 | APOE 1, BRCA1 1, KCNH2 22, KCNQ1 25, RYR2 14, SCN5A 57 | $12.40 | $24.79 | $30.99 |
| `mixed_gold_cont120_05` | 121 | BRCA1 1, BRCA2 1, KCNH2 22, KCNQ1 25, MYBPC3 1, RYR2 14, SCN5A 57 | $12.77 | $25.54 | $31.93 |
| `mixed_gold_cont120_06` | 121 | APOE 1, BRCA1 1, KCNH2 22, KCNQ1 25, MYBPC3 1, RYR2 14, SCN5A 57 | $12.64 | $25.29 | $31.61 |
| `mixed_gold_cont120_07` | 120 | BRCA1 1, BRCA2 1, KCNH2 21, KCNQ1 25, MYBPC3 1, RYR2 14, SCN5A 57 | $12.68 | $25.37 | $31.71 |
| `mixed_gold_cont120_08` | 120 | BRCA1 1, BRCA2 1, KCNH2 22, KCNQ1 25, RYR2 14, SCN5A 57 | $12.53 | $25.05 | $31.31 |
| `mixed_gold_cont120_09` | 120 | BRCA1 1, BRCA2 1, KCNH2 22, KCNQ1 25, RYR2 14, SCN5A 57 | $12.53 | $25.05 | $31.31 |
| `mixed_gold_cont120_10` | 121 | APOE 1, BRCA2 1, KCNH2 22, KCNQ1 25, MYBPC3 1, RYR2 14, SCN5A 57 | $12.57 | $25.15 | $31.44 |
| `mixed_gold_cont120_11` | 121 | APOE 1, BRCA1 1, BRCA2 1, KCNH2 22, KCNQ1 24, RYR2 14, SCN5A 58 | $12.65 | $25.30 | $31.63 |
| **Suite** | **1324** |  | **$138.50** | **$277.00** | **$346.25** |
