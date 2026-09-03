# Canonical evaluation tiers

The four registry entries below preserve the staged rollout cohorts and their
historical gates. For ongoing protocol-change testing, the canonical broad
suite is now [`mixed_gold/`](mixed_gold/README.md): it inventories every named
variant gold source, assigns every source-available gene-paper attempt exactly
once across 29 heterogeneous tranches, makes paper-derived identity the primary
score, and records an observed-usage cost estimate for every tranche.
The estimate is recorded both per arm and for the paired frozen-baseline versus
candidate comparison required for protocol changes.

Dated benchmark runs and older manifests remain reproducibility artifacts; they
are not extra gates and must not be described as current.

| Order | Tier | Gene–paper attempts | Unique PMIDs | Use |
|---:|---|---:|---:|---|
| 1 | `gold_50` | 50 | 50 | Scored gate: fixed cardiac 48 plus Nate's two lead-approved BRCA2 papers |
| 2 | `gold_120` | 118 | 114 | Scored cardiac expansion: originally 30/gene; KCNH2 is 28 after quarantining two invalid PMIDs |
| 3 | `reviewer_545` | 545 | 506 | Full reviewer backlog: ten 50-paper queues plus the 45-paper BRCA2 queue |
| 4 | `gold_120b` | 125 | 124 | Scored replication tranche: 30/gene cardiac plus the 5 remaining curated BRCA2 papers, sharing no article with tiers 1-3 or the gold-150 rosters |

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

4. `tier4_gold_120b.tsv` is the second scored tranche, drawn on 2026-08-25 with
   seed `2026082501` by `select_tranche.py`. It answers a different question
   from `gold_120` — does a protocol tuned on tier 2 hold on gold papers it has
   never been scored against? — so it is **article-level** disjoint from tiers
   1-3 and from the preregistered gold-150 calibration/holdout rosters. Same
   article under a different gene still counts as used: a multi-gene paper
   already scored under KCNQ1 has had its tables optimised against, and
   BRCA1/BRCA2 output differing only in the gene column is a recorded failure
   here. Eligibility calls the same `run_eval.gold_count_eligible_pmids` helper
   tier 2 used, so the rule cannot drift between tranches, and selection reads
   gold presence only — never values or row counts. BRCA2 is 5 because only 8
   BRCA2 papers carry curated gold at all and three are already spent in tier 1
   and gold-150. Provenance, per-paper source digests, and the rejection ledger
   are in `tier4_gold_120b_selection.json`; the frozen per-gene answer key is in
   `gold_120b_answer_key/`.

   Its BRCA2 gold comes from the 2026-07-06 `gold_overrides` reconciliation, not
   from the lead-approved Variant Browser snapshot — both lead-approved BRCA2
   papers are already in tier 1. Report the BRCA2 arm separately and name that
   provenance; do not fold it into a cardiac headline.

   Statistical note: the cardiac arm carries 286 gold rows over 120 attempts
   against tier 2's 632 over 118. Both draws are blind and their count-bearing
   share matches (79% vs 81%); gold rows per paper are long-tailed with a median
   of 1, and tier 2 simply caught more tail papers. Expect wider recall
   intervals here, and do not re-draw to chase row count — selecting on row
   count is selecting on the answer key.

Tiers 1, 2, and 4 are scored benchmarks; Tier 3 is the independent operational
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
