# Canonical evaluation tiers

These are the only active paper cohorts for protocol rollout. Dated benchmark
runs and older manifests remain reproducibility artifacts; they are not extra
gates and must not be described as current.

| Order | Tier | Gene–paper attempts | Unique PMIDs | Use |
|---:|---|---:|---:|---|
| 1 | `gold_50` | 50 | 50 | Scored gate: fixed cardiac 48 plus Nate's two lead-approved BRCA2 papers |
| 2 | `cardiac_120` | 120 | 98 | Cardiac-only expansion: the first 30 ranked review papers in each cardiac workspace |
| 3 | `reviewer_396` | 396 | 357 | Full established reviewer backlog: seven 50-paper queues plus the 46-paper BRCA2 queue |

Counts are **gene–paper extraction/review attempts**, not necessarily distinct
articles. A paper can appear in more than one gene–disease workspace. This is
why the full backlog has 396 attempts but 357 unique PMIDs, close to the
remembered “about 340 papers.”

## Gate definitions

1. `tier1_gold_50.tsv` is the fixed gold-scored strategy gate. Cardiac and
   BRCA2 metrics remain separate because their gold provenance differs.
2. `tier2_cardiac_120.tsv` is the operational cardiac expansion. It takes the
   first 30 entries from each ranked 50-paper cardiac review list. It is a
   reviewer/rollout cohort, not a new answer key.
3. `tier3_reviewer_396.tsv` is the exact live paper membership of the eight
   established `review-50-all-genes-20260713` workspaces on 2026-08-11:
   APOE, BRCA1, BRCA2, KCNH2, KCNQ1, MYBPC3, RYR2, and SCN5A. BRCA2 has 46
   papers after the provenance exclusion.

Tier 2 is a strict subset of Tier 3. Tier 1 is deliberately a separate scored
benchmark rather than a subset: 27 of its high-carrier cardiac papers are not in
the operational reviewer queues.

The private Variant Browser dashboard also contains separate BMPR2 cold-start
(50 papers), LMNA/TTN temporal experiments (99 each), and empty legacy
workspaces. They are deliberately excluded here: they were not part of the
established 50-per-pair cohort and would turn this into an unplanned 644-attempt
rollout.

`registry.json` is the machine-readable index and pins each manifest SHA-256.
The unit test
`tests/unit/test_evaluation_tiers.py` pins counts, membership derivation, PMID
validity, and the live SCN5A reconciliation so the active tiers cannot drift
silently.

## Historical experiments

Keep dated contents under `benchmarks/*/runs/` and the append-only protocol
history. They support scientific reconstruction of what worked and failed.
Do not add another active cohort directory for an experiment; put experimental
outputs in a dated run directory and promote them here only after a deliberate
cohort decision.
