# BRCA2 collaborator-gold scope correction — methods record

## Reason for the correction

The 2026-08-10 56-paper study treated all eight BRCA2 override papers as a
single provisional stratum. A provenance audit on 2026-08-11 found that only
PMIDs 26833046 and 26848529 were present in the lead-approved Variant Browser
gold snapshot with reviewer metadata. The other six originated in an internal
candidate re-derivation and review pass and therefore did not meet the intended
collaborator-gold inclusion rule.

## Inclusion rule

A BRCA2 paper is eligible for the active scored comparison cohort only when
its current reference records are traceable to a Variant Browser collaborator
review and lead approval. Applying this rule retained:

- PMID 26833046 — source reviewer `nate`, approving lead `kronckbm`.
- PMID 26848529 — source reviewer `nate`, approving lead `kronckbm`.

It excluded PMIDs 10398279, 15365993, 18489799, 21356067, 22655046, and
25802882 from active review and new strategy comparisons. In the live BRCA2
review queue, four were present and removed; 10398279 and 21356067 were already
absent. The queue changed from 50 to 46 papers, retained both Nate-reviewed
papers and all 87 associated current gold records, and removed no adjudications
or gold. Historical benchmark records were not deleted or rewritten.

## Metric projection

The original cardiac and BRCA2 prediction files remained locked. No model call,
extraction, or adjudication was rerun. Per-paper contributions from the six
excluded BRCA2 papers were removed from the historical aggregate, leaving the
cardiac 48 and the two eligible BRCA2 papers. Conditional MAE and count recall
retain the definitions in `METHODS_20260810.md`.

| Field | Before absolute error / MAE | After absolute error / MAE | Count recall after |
| --- | ---: | ---: | ---: |
| Carriers | 298 / 0.9058 | **20 / 0.0608** | 329/1007 (32.67%) |
| Affected | **193 / 0.7148** | 194 / 0.7185 | 270/1006 (26.84%) |
| Unaffected | 39 / 0.1566 | **25 / 0.1004** | 249/1005 (24.78%) |

The two-paper BRCA2 carrier stratum has only three supplied carrier counts and
MAE 1.3333 (4/3). The combined 0.0608 value is dominated by the cardiac stratum
and must always be reported with this decomposition.

## Reproducibility

- Active manifests:
  `benchmarks/codex_paper_eval/highcarrier48_plus_brca2_collaborator2_20260811.tsv`
  and `brca2_2_collaborator_reviewed_20260811.tsv`.
- Machine-readable projection:
  `runs/20260811_collaborator_gold_50/metrics.json`.
- Historical source:
  `runs/20260810_luna_xhigh_56/metrics.json` and its locked prediction digests.
- Review-publishing exclusions:
  `benchmarks/curated_extraction_eval/review_scope_exclusions.tsv`.
- Live Variant Browser audit record: SourceRun
  `brca2_review_scope_exclusion_20260811`, allowlist SHA-256
  `32d9dc822106da9a78aeef469f7bbc38428754e35bbf8da3362e11e5e41f5cd0`.

## Interpretation

This correction strengthens provenance but reduces BRCA2 sample size. It does
not validate BRCA2 generalization. Future BRCA2 papers enter the active cohort
only after collaborator review and lead approval are synchronized from Variant
Browser.
