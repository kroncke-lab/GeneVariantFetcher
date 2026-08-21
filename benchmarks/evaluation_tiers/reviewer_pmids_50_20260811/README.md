# Additional 50-paper reviewer cohorts (2026-08-11)

These files add BMPR2, LMNA, and TTN to the canonical full reviewer tier.

- BMPR2 preserves the exact 50-paper cohort already live in Variant Browser,
  including its existing review order. Its 2026-08-08 source audit found 42
  full-text papers and eight abstract-only papers; the latter remain reviewer
  and source-recovery work, not evidence that the extractor saw the complete
  article.
- LMNA and TTN were legacy 99-paper temporal-mix snapshots. Each is narrowed to
  50 by descending count-bearing VariantEvidence rows, then descending total
  VariantEvidence rows, then ascending PMID. This prioritizes papers with
  reviewable carrier/count evidence while filling TTN's remaining slots with
  its richest variant-evidence papers.
- At selection time none of the three workspaces had reviewer adjudications or
  current gold records.

Live reconciliation was verified on 2026-08-11. BMPR2 already matched its
manifest and retained its original order. LMNA and TTN now match their manifests
exactly, have `review_order` 1–50, and record the manifest path, checksum,
selection policy, and retained PMIDs in Variant Browser source runs labeled
`review_scope_50_20260811`. All 150 staged papers expose a non-empty
paper-specific reviewer summary.

The files preserve ranking order. The combined full-tier manifest is
`../tier3_reviewer_545.tsv`.
