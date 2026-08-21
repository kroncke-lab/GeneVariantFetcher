# Canonical evaluation tiers

These are the only active paper cohorts for protocol rollout. Dated benchmark
runs and older manifests remain reproducibility artifacts; they are not extra
gates and must not be described as current.

| Order | Tier | Gene–paper attempts | Unique PMIDs | Use |
|---:|---|---:|---:|---|
| 1 | `gold_50` | 50 | 50 | Scored gate: fixed cardiac 48 plus Nate's two lead-approved BRCA2 papers |
| 2 | `gold_120` | 118 | 114 | Scored cardiac expansion: originally 30/gene; KCNH2 is 28 after quarantining two invalid PMIDs |
| 3 | `reviewer_545` | 545 | 506 | Full reviewer backlog: ten 50-paper queues plus the 45-paper BRCA2 queue |

Counts are **gene–paper extraction/review attempts**, not necessarily distinct
articles. A paper can appear in more than one gene–disease workspace. This is
why the full backlog has 545 attempts but 506 unique PMIDs.

## Gate definitions

1. `tier1_gold_50.tsv` is the fixed gold-scored strategy gate. Cardiac and
   BRCA2 metrics remain separate because their gold provenance differs.
2. `tier2_gold_120.tsv` is the scored cardiac expansion. It is a seeded
   (`2026081301`), gold-value-blinded sample of source-available papers per
   cardiac gene, restricted to PMIDs with explicit gold assertions for carriers,
   affected, and unaffected. The original draw was 30 per gene. On 2026-08-15
   KCNH2 10086972 was removed because that PMID is not a genetics paper (the
   gold row belonged to 10086971, already in the sample). On 2026-08-21 KCNH2
   14642689 was also quarantined: PubMed identifies it as an angiotensin-II
   receptor expression study, not a KCNH2 genetics paper, and the intended
   source for its A561T gold row is unknown. Live membership is 118 attempts /
   114 unique PMIDs (KCNH2 28). Selection used only PMID
   eligibility and assertion presence; extraction does not see gold values or
   row counts.
3. `tier3_reviewer_545.tsv` is the full private-review population on
   2026-08-11: APOE, BMPR2, BRCA1, BRCA2, KCNH2, KCNQ1, LMNA, MYBPC3, RYR2,
   SCN5A, and TTN. Every workspace has 50 papers except BRCA2, which has 45
   after the provenance exclusion and the 2026-08-20 removal of the explicitly
   canine BRCA2 ortholog paper PMID 19944633. BMPR2 preserves its existing 50-paper
   cohort; LMNA and TTN are narrowed from 99 to the ranked 50-paper manifests
   in `reviewer_pmids_50_20260811/`.

Tiers 1 and 2 are scored benchmarks; Tier 3 is the independent operational
review backlog. The `gold_120` sample overlaps the cardiac arm of `gold_50` by
nine attempts and adds 110 attempts, avoiding a second high-carrier-only sample.

The empty ALPL, N4BP2L1, and TTR legacy workspaces are not paper cohorts.

`registry.json` is the machine-readable index and pins each manifest SHA-256.
The unit test
`tests/unit/test_evaluation_tiers.py` pins counts, gold eligibility, PMID
validity, and the live reviewer reconciliation so the active tiers cannot drift
silently.

## Historical experiments

Keep dated contents under `benchmarks/*/runs/` and the append-only protocol
history. They support scientific reconstruction of what worked and failed.
Do not add another active cohort directory for an experiment; put experimental
outputs in a dated run directory and promote them here only after a deliberate
cohort decision.
