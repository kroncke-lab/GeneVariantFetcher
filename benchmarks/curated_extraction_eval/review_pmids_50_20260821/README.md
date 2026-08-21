# Corrected 50-paper BRCA2 collaborator cohort

This versioned manifest restores a 50-paper BRCA2 review queue without
reintroducing the five papers removed from the original July cohort. It keeps
the exact active 45-paper order from
`../review_pmids_20260811_brca2_provenance/BRCA2.txt` and appends five primary
human clinical studies selected by a deterministic walk down the canonical
carrier-bearing ranking.

The four collaborator-provenance exclusions (`15365993`, `18489799`,
`22655046`, and `25802882`) remain excluded. Canine BRCA2 ortholog paper
`19944633` remains a negative control and must never enter the clinical review
database.

`selection_walk.tsv` records every decision from canonical rank 51 through the
fifth eligible replacement. Eligibility requires primary human clinical BRCA2
carrier observations from this study. Functional-only, review, computational,
somatic-only, wrong-gene-ranked, methods-validation, and other-study-only papers
are rejected. The ranking determines inspection order, not automatic inclusion.

The two papers with known richer same-PMID source copies under the BRCA1 corpus
(`26183948` and `25884701`) require supplement recovery or an identity-verified
source override before the BRCA2 run can pass readiness. `20380699` and
`26843898` are explicit attribution regression targets: historical rows contain
wrong-gene or other-study contamination that current extraction must not repeat.

This is a collaborator-review cohort, not a gold standard and not a headline
metric set. It does not replace the frozen historical manifests. Public
publication remains opt-in and requires the complete source, trust, identity,
and Variant Browser import gates.
