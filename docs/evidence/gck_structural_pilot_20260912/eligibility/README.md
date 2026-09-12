# GCK pilot identity, eligibility, and endpoint audit

Date: 2026-09-12. This is a read-only audit of the frozen count inputs and their
archived feature/source lineage. No counts, empirical hyperparameters, source
identities, or predictor scores were corrected or refitted here.

## Eligibility used by the exploratory pilot

All **398 count-bearing keys** remain eligible for the common empirical prior,
including the seven synonymous keys. They contain 932 affected observations,
82 literature unaffected observations, and 3,620 gnomAD observations treated as
unaffected under the user's assumption. The restored synonymous keys account for
2,553 gnomAD observations omitted from the earlier 391-row histogram subset.
Affected plus literature unaffected exactly reproduces the archived literature
denominator for every key. The frozen empirical prior remains the one computed
from all 398 keys, before the geometry exclusions below.

For the local structural kernel, require the archived consequence to be
`missense`, a canonical protein position, and an exact reference-amino-acid match
against **P35557-1**. Initiator methionine changes are ineligible; there are no
position-1 rows in this input. This yields **249 of 268 missense keys** spanning
182 canonical residues. Nineteen missense keys have a reference mismatch and are
explicitly ineligible for geometry. The other 130 rows have non-missense
consequences and remain outside this kernel. Availability and confidence of an
actual structural context are additional checks performed by the geometry
runner; `geometry_eligible_pre_structure` is not a claim that a residue is
experimentally resolved.

The 19 mismatched keys are W289C, C253Y, G73E, I160M, V61G, A209V, E275V,
G262R, G73R, M35I, N257T, R187W, R225H, R272C, R423P, S442P, T210M, T260M,
and T432K. See [the exact reference comparisons](wt_mismatches.csv). These rows
are retained in the prior and output universe, with a geometry exclusion reason.
They have not been renamed, shifted, or merged on the basis of a plausible
correction. Source/transcript review is needed to resolve them.

All archived feature rows name `ENST00000403799.8`. That field records the
feature-join transcript, rather than proving that every original paper used that
transcript. The canonical sequence is frozen by the geometry audit in
`../geometry/P35557-1.fasta`; this audit checked the sequence, not the labels on
experimental author numbering.

## Identity and variant-only exclusion

There is one row per archived `key`; all 398 keys are unique. All 268 missense
`(reference, position, alternate)` tuples are also unique. The pilot identifier
is `GCK:<key>`. This is a **protein-key aggregate**, not a newly adjudicated
genomic-allele identity. Excluding this identifier removes the whole archived
aggregate. It does not claim allele-level leave-one-out when a key combines or
ambiguously maps different DNA alleles.

Eight missense keys map to more than one warehouse allele in the frozen feature
join: W99R, D278E, D274E, H317Q, V91L, W99C, M197I, and M251I. Ten missense
keys have more than one cDNA notation in retained source observations: A259T,
A449T, D278E, F152L, F171L, G261R, G318R, L146P, V455L, and W99R. Those source
notations can represent distinct alleles or extraction/notation errors; their
truth has not been adjudicated. They are preserved in
[identity_alias_audit.csv](identity_alias_audit.csv). Do not split their counts
or treat those aliases as separate spatial donors.

Two warehouse IDs occur under two different missense keys in the archived join:
A209V/T209M and W289C/Y289C. In both pairs, the first key fails canonical WT
validation and is excluded from geometry; the valid members do not share an
allele ID with another retained geometry-eligible key. The exact records are in
[shared_allele_ids.csv](shared_allele_ids.csv). A shared feature ID is evidence of
an identity conflict, not sufficient evidence to merge the source counts.

The upstream count builder already normalizes protein notation and folds
cDNA-only aliases into a protein key when linked by a source row. Within a
key–PMID group it keeps the variant ID with the largest summed denominator;
586 of the 616 retained observations record merged duplicate source IDs. This
audit preserves that frozen aggregation. It does not independently establish
person/family uniqueness across papers, or recover allele-specific counts lost
by aggregation.

Among the 249 WT-valid missense keys, **56 residues have multiple substitutions,
containing 123 keys**. These are retained as distinct donors when another
substitution at the same residue is the target, following the user's
variant-only exclusion rule. See
[same_residue_alternatives.csv](same_residue_alternatives.csv).

## Endpoint evidence and limits

The structural pilot must be labeled **exploratory pooled GCK clinical endpoint**.
The archived cohort count table has affected/unaffected counts and trust fields,
but no per-count disease or direction-of-effect field. The per-variant feature
table also lacks an adjudicated phenotype label. A trusted extraction flag does
not establish that all rows measure MODY or the same clinical endpoint.

The source `individual_records` table contains `phenotype_details` and
`evidence_sentence`. Inspection restricted to the exact retained database
variant ID and PMID identifies clear hypoglycemia/hyperinsulinism examples:

- W99R, PMID 34680961: the retained individual-record observation contains four
  affected relatives, with extracted descriptions of hyperinsulinemic
  hypoglycemia from early life.
- V389L, PMID 24890200: the retained cohort observation is 5 affected / 0
  unaffected; associated extracted individual evidence describes a previously
  characterized activating GCK mutation and hyperinsulinemic hypoglycemia.

These records demonstrate that the frozen affected counts are not uniformly
GCK-MODY outcomes. This audit does not determine corrected counts or remove
variants automatically. A disease-specific comparison requires source-level
endpoint, inheritance, and effect-direction curation. In particular, some
keywords occur in medication-related hypoglycemia, negated statements, or
hyperglycemia with insulin resistance; a keyword is not an exclusion rule.

[phenotype_keyword_review.csv](phenotype_keyword_review.csv) records a bounded
screen of associated extracted individual evidence for retained observations:
29 key–PMID observations across 22 keys contain `hypoglyc`, `hyperinsulin`,
`activating`, or `gain.of.function`. The screen is incomplete for cohort rows
without individual evidence. Its excerpts are source-extraction material, not
independently adjudicated diagnoses. No subjects' individual identifiers are
included in this report.

V455E is retained with its frozen PMID 36208030 observation, 3 affected and 25
literature unaffected; that cohort row provides no phenotype-specific field.
Its endpoint or activating status cannot be inferred from that row, the
AlphaMissense score, or memory. It is not automatically excluded.

## Archived AlphaMissense comparator

[archived_alphamissense.csv](archived_alphamissense.csv) contains exactly one
row for every one of the 398 frozen keys, joined one-to-one. It copies the
archived `GCK_protocol/variants_features.csv` score without a warehouse refresh,
imputation, re-selection, or averaging in this audit. Coverage is 255/398 across
all consequence classes, 252/268 for missense, and **246/249** for canonical
WT-valid missense. Missing scores remain empty with `am_missing=True`.

The feature-join method, transcript, mapped allele count, and archived
`alphamissense_allele_range` are retained. For the eight multi-allele missense
keys, that archived range is zero. Identical scores do not resolve their DNA
identity ambiguity. Six WT-invalid missense keys also have archived scores;
score presence must not override the canonical geometry exclusion.

## Files and checks

- `variant_eligibility.csv`: all 398 keys; provisional geometry eligibility,
  canonical WT, explicit exclusion reason, identity grain, ambiguity flags,
  and endpoint label. Every row remains `prior_eligible=True`.
- `archived_alphamissense.csv`: one archived comparator score or explicit
  missing value per key, with source hash and join provenance.
- `wt_mismatches.csv`, `same_residue_alternatives.csv`,
  `identity_alias_audit.csv`, `shared_allele_ids.csv`: inspectable identity
  checks without rewriting the original count aggregates.
- `phenotype_keyword_review.csv`: retained-observation source evidence for
  endpoint review, with database and PMID lineage.
- `source_provenance.csv`: source paths and SHA-256 hashes. The six read-only
  GVF databases still match the hashes in the original observation summary;
  all listed file hashes were verified during this audit.

Checks passed: unique key universes across source counts, empirical posteriors,
eligibility and archived scores; exact archived affected/unaffected/literature
denominator agreement; count partitions; canonical sequence bounds and WT
comparison; one-to-one score join; and original source-file hashes. Geometry,
clinical endpoint adjudication, and independence of observations remain separate
from these identity and arithmetic checks.
