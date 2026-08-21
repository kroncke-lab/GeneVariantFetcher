# PMID 19944633 paper-scope failure and repair

Date: 2026-08-20. Publication status: **HOLD** pending a source-grounded
precision review and trusted-field-only import validation.

## What actually happened

PMID 19944633 is titled “Single nucleotide variation in exon 11 of canine
BRCA2 in healthy and cancerous mammary tissue.” It should never be a paper in a
human BRCA2 clinical review queue.

Two distinct runs exposed two harness generations:

- The 2026-08-13 experimental run admitted the explicit PMID and sent it to an
  LLM, which returned 19 explicitly canine variants. That was a harness failure,
  not a hard biomedical judgment.
- In the 2026-08-20 no-publish run, the extractor's deterministic ortholog gate
  worked: its decision trace records zero provider attempts,
  `model_used=skipped-nonhuman-ortholog`, and zero accepted variants. The DB
  persisted that marker. PubTator recovery then ignored it and attached four
  `variant_papers` links, which is why the final structural audit failed.

The bad paper also remained in the explicit 46-paper BRCA2 allowlist, so the
current system was spending source/recovery effort on an input that should have
been rejected before extraction.

## Repair

- Remove 19944633 from the active BRCA2 allowlist; retain it only in the
  wrong-species negative controls and explicit review-scope exclusions.
- Run a deterministic paper-scope gate even for explicit PMID manifests.
- Persist one paper-level exclusion reason from both the early pre-model skip
  and the late post-extraction gene filter.
- Make source replay refuse excluded papers even when forced.
- Make ClinVar, PubTator, and figure recovery honor the DB exclusion and purge
  legacy PMID-linked evidence while retaining the paper/extraction audit rows.
- Add PubMed/PubTator title preflights as defense in depth for legacy DBs.

The active canonical Tier 3 is now `reviewer_545`: 545 gene-paper attempts / 506
unique PMIDs, including 45 BRCA2 papers.

## Verification

- Synthetic coverage pins explicit admission, early/late exclusion persistence,
  migration, forced replay refusal, ClinVar, PubTator, figures, source audit,
  active cohort membership, and legacy evidence purge.
- A copy of the real Aug-20 BRCA2 DB started with four PMID 19944633
  `variant_papers` links; the shared repair removed four, left zero, and retained
  the `skipped-nonhuman-ortholog` audit metadata.
- Focused regression matrix: 83 passed. Full offline suite: 2029 passed / 1
  skipped; changed files pass Ruff, format, and `git diff --check`.
- Grok 4.6 `xhigh` adversarial review first returned NO-GO for a stale cohort
  test, an early-marker-only durable key, replay of late drops, and the missing
  ClinVar title belt. After those were fixed and tested, the same review returned
  **GO** with no remaining P0/P1 bypass.

## Corrected rerun

The active 45-paper manifest was rerun from scratch at
`results/vb_reextract_20260820_scopefix/BRCA2/20260820_194856`. It completed
45/45 extractions and the final validation run exited 0 with no stage failures
or warnings. The final DB contains 1,917 variants, complete VariantFeatures
coverage for every live variant, and 137 quarantined false positives. PMID
19944633 is absent from the staged manifest, extraction files, papers,
extraction metadata, variant links, fact provenance, penetrance, individual,
functional, and phenotype tables. The regenerated three-gene collaborator
audit passes BRCA1 50/50, BRCA2 45/45, and BMPR2 50/50 with zero live
species-scope links. The final offline unit suite is 2,041 passed / 1 skipped.

This makes the corrected dataset structurally eligible for source adjudication;
it does not make either the rejected Aug-20 output or the new output
publication-ready. A source-grounded precision sample and trusted-field-only
Variant Browser import validation remain required.
