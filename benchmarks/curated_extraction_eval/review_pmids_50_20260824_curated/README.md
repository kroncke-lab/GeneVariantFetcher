# Source-screened BRCA1/BRCA2/BMPR2 collaborator cohort

This is the presentation-ready private-review cohort: exactly 50 papers for
each of BRCA1, BRCA2, and BMPR2. It supersedes the earlier three-gene staging
membership for collaborator review, but it does not rewrite the frozen
historical manifests or authorize public annotation publication.

## Paper-level inclusion contract

- BRCA1/BRCA2: a primary human study with current-study germline carrier,
  family, case/control, or cohort observations for the target gene.
- BMPR2: a primary human PAH or closely related pulmonary-vascular study with
  current-study BMPR2 carrier observations.
- Reviews, literature/database compilations, method-only validation sets,
  wrong-gene interaction papers, somatic-only studies, explicitly negative
  target-gene cohorts, off-disease association studies, and cell-model-only
  papers are out of scope.
- Mixed studies remain only when the target-gene clinical observations are a
  direct part of the current cohort and can be separated during review.

All 29 replacement papers were checked against their PubMed title/abstract and
an identity-matched local source before selection. All replacement sources are
full text rather than abstract-only fallbacks. The final three manifests are
also 150/150 full text. `membership_changes.tsv` records the 29 removals and 29
additions with paper-specific reasons.

The manifests are inputs to a fresh gold-disabled extraction. Publication must
use the resulting completed active DBs through Variant Browser's trusted,
login-gated staging importer. Public `publish_annotations` remains out of
scope.
