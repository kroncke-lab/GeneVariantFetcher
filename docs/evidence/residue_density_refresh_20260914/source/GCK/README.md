# GCK source and endpoint correction — 2026-09-14

The refreshed primary is a **GCK-MODY / mild-hyperglycemia analysis with a
bounded source correction**, not a fully adjudicated penetrance cohort. It
excludes verified hyperinsulinemic-hypoglycemia observations and incompatible
genotype observations. The phenotype review also corrects two families whose
heterozygous MODY relatives had been mixed with neonatal-diabetes probands.
No unknown phenotype is converted to unaffected, and no gnomAD count is changed.

## Inputs for the new population union

| File | Meaning | Canonical missense clinical A | Canonical missense clinical U |
|---|---|---:|---:|
| Frozen input | Before this review | 670 | 72 |
| [clinical_input.csv.gz](clinical_input.csv.gz) | Primary: source/zygosity corrections; diagnosis-only proxy quarantined | 603 | 32 |
| [clinical_input_diabetes_proxy.csv.gz](clinical_input_diabetes_proxy.csv.gz) | Broader sensitivity: same corrections, retaining the flagged diagnosis-only proxy | 617 | 66 |

Both files preserve all **398 original GCK clinical keys** and every original
identity column. `affected_literature` and `unaffected_literature` are the
corrected input counts. Old counts, removed counts, restored counts, review
status and strict-endpoint columns preserve the decisions. **Filter to
`n_literature > 0` before rebuilding the union.** There are 28 keys without
remaining clinical evidence in the primary and 21 in the broader sensitivity.
Their population alleles must remain in the observed population universe and
be reassigned by the normal union rules; do not delete population alleles or
retain an unsupported clinical protein aggregate.

These are clinical-key totals before the population union's identity rules.
The complete primary ledger across all variant types has A=860 and U=42;
the broader ledger has A=879 and U=76. Missense and nonsense priors must still
be fit separately, with canonical-WT eligibility applied after the union.

## Verified corrections

[GCK_observation_corrections.csv](GCK_observation_corrections.csv) records 35
exact frozen variant-by-paper rows with source coordinates and SHA-256 pins.
The affected/unaffected replacement accounting removes the original A=59/U=6
and restores verified A=6/U=0. This is **net A−53/U−6**, not 59 discarded
unique affected people. Repeated family members across papers remain a separate
ownership question.

- **Hypoglycemia and hyperinsulinism:** the review extends the earlier W99R/V389L
  audit to the exact observed rows for A456V, Y214C, V452L, M197I, W99L, T103S,
  I211F, V455M, M197T, W99C, T65I, E67V, S64P, V91L, Y215C, M197V, R447L and
  V455L. The ledger lists the specific paper for each removal. Variants are not
  excluded simply because an in-vitro experiment is activating.
- **A456V, PMID 11916951:** the mother had no hypoglycemic symptoms but did have
  measured low fasting glucose. The original A1/U1 concerns two hypoglycemia
  carriers; neither becomes MODY-unaffected. See the
  [primary abstract](https://pubmed.ncbi.nlm.nih.gov/11916951/).
- **V389L, PMID 24890200:** the original A5 includes four biochemical
  hypoglycemia observations and one unassessed relative. All five leave the
  MODY analysis; the unassessed person is not a negative. This follows the
  [previous exact patient review](../../../structural_sanity_20260913/gck/endpoint_review.csv).
- **V455L:** one hypoglycemic child from
  [PMID 32928245](https://pmc.ncbi.nlm.nih.gov/articles/PMC7490857/) and four
  HI-ascertained family members from
  [PMID 34532767](https://pmc.ncbi.nlm.nih.gov/articles/PMC8563668/) leave the
  MODY endpoint. One adult in the latter family subsequently developed type 2
  diabetes; that does not establish four MODY observations.
- **R397L, PMID 15644838:** replace A1/U2 with **A2/U0**. The infant is homozygous
  with neonatal diabetes; both heterozygous parents have measured mild
  hyperglycemia, explicitly interpreted as GCK-MODY by the paper. The parents'
  fasting values are 5.8 and 6.1 mmol/L. See the
  [primary report](https://pubmed.ncbi.nlm.nih.gov/15644838/).
- **N254H, PMID 26587058:** replace A4/U0 with **A2/U0**. Two infants have
  compound heterozygous neonatal diabetes, while their two mothers carry N254H
  heterozygously and have mild fasting hyperglycemia. Retain the mothers' MODY
  observations. See [Results and Figure 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC4652399/).
- **T209M, PMID 25555642:** retain **A2/U0** for two constitutionally heterozygous
  siblings with MODY. Remove the original U1 for their father, whose blood
  variant is low-level mosaic; his normal measured glucose is preserved in the
  ledger, but he does not represent a constitutional heterozygous denominator.
  See the [primary report](https://pubmed.ncbi.nlm.nih.gov/25555642/).
- **H424Y, PMID 39610513:** quarantine A1/U2. The proband is homozygous with
  neonatal hyperglycemia; the available source does not resolve each
  heterozygous parent's mild-hyperglycemia status. The report is retained for
  source follow-up, not used to assert two unaffected parents. See the
  [primary abstract](https://pubmed.ncbi.nlm.nih.gov/39610513/).

Source-specific decisions preserve counterexamples. **V62M remains eligible**:
the original families co-segregate with measured hyperglycemia despite increased
activity in some in-vitro assays. **N180D, PMID 34496959 remains eligible** because
the study's clinical table explicitly identifies the observed child as GCK-MODY.
An activator-drug trial and a study mentioning iatrogenic hyperinsulinemia also
retain their GCK-MODY observations. This avoids treating a keyword or mechanism
as the patient's endpoint.

## Diabetes diagnosis is a separate endpoint

[GCK_diabetes_endpoint_pending.csv](GCK_diabetes_endpoint_pending.csv) preserves
the 15 observations from PMID 36208030: A19/U34 across all types, including
missense A14/U34. The source separates participants by diabetes diagnosis codes
and considers a secondary HbA1c threshold of 6.5%. It does not provide the
variant-specific split at the lower mild-hyperglycemia threshold appropriate to
GCK-MODY. These rows are therefore quarantined in the primary and included only
in the broader sensitivity. In particular, V455E's A3/U25 does not establish
25 people without mild GCK hyperglycemia. See the
[paper's methods and supplementary table](https://pmc.ncbi.nlm.nih.gov/articles/PMC9659663/).

PMID 36257325 defines GCK penetrance using fasting glucose ≥5.6 mmol/L and/or
HbA1c ≥39 mmol/mol. It contributes **zero usable observations to the frozen
GCK input**; its unsplit variant totals are not added or assigned affected
status here. This paper's high penetrance for curated pathogenic variants does
not measure every missense variant. See
[phenotype definitions and Figure 4](https://pmc.ncbi.nlm.nih.gov/articles/PMC9674944/).

Removing incorrect unaffected observations can increase a posterior or a prior.
Neither source correction nor this endpoint choice is justified by whether it
makes the residue plot lower. The broader sensitivity makes that distinction
inspectable.

## Remaining limits and reproduction

The [19 canonical-WT mismatches](GCK_canonical_identity_queue.csv) remain
excluded from the missense fit. This pass makes **no identity repairs**. Some
papers repeat inconsistent residue/cDNA notation; several patterns suggest
isoform numbering, but that is insufficient to authorize an offset. Exact
source/transcript alignment and nucleotide ownership are still required.
W99R's cDNA conflict also remains unresolved. A WT match alone cannot rule out
every isoform error, because neighboring residues can have the same amino acid.

[missense_source_inventory.csv](missense_source_inventory.csv) is the initial
paper-level triage inventory from the frozen observations and source databases.
Placeholder titles were checked against cached full texts where relevant; it
is not a record that all remaining sources or family ownership have been
adjudicated. Mixed diabetes etiologies, age, treatment, co-occurring variants,
case ascertainment and repeated people remain limitations of the pooled model.

Run [build_source_ledger.py](build_source_ledger.py) with the sibling
BayesianPenetranceEstimator Python environment. It checks that original clinical
key counts exactly reproduce the frozen variant-by-paper observations, asserts
unique observation identities and nonnegative corrections, pins reviewed
sources, and writes deterministic GZip/CSV with LF newlines. All generated
files are smaller than 1.2 MB. [source_checks.json](source_checks.json) records
the counts, source hashes and script hash. Source DBs, earlier frozen evidence,
and population counts are unchanged.
