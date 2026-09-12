# Population-inclusive union: independent identity and count audit

Date: 2026-09-12. This audit preserves the frozen literature data and supplies
identity flags for the new population-inclusive analysis. It does not rewrite
paper counts, split aggregate counts across alleles, or change the historical
results. gnomAD carriers are treated as unaffected, as requested.

## Source and transcript alignment

The five archived feature transcripts are HNF1A `ENST00000257555.11`, GCK
`ENST00000403799.8`, LDLR `ENST00000558518.6`, BRCA2 `ENST00000380152.8`, and
KCNQ1 `ENST00000155840.12`. The source collector independently confirms these
versions against the current gene API annotations. The warehouse's
`transcripts` registry contains MANE/canonical rows for HNF1A, GCK and LDLR,
but no BRCA2 or KCNQ1 rows. Its `genes.canonical_transcript` value for KCNQ1 is
the RefSeq label `NM_000218.3`, not the archived Ensembl identifier. Therefore
neither a first transcript row nor an unqualified cDNA alias is a safe join.
See [canonical_transcript_audit.csv](canonical_transcript_audit.csv).

Canonical WT checks use the five UniProt FASTA snapshots collected under
`../population/`, including their file hashes: HNF1A P20823, 631 residues; GCK
P35557, 465; LDLR P01130, 860; BRCA2 P51587, 3,418; and KCNQ1 P51787, 676.
Gene/canonical-transcript consequences must agree with those reference residues
before a protein-point-substitution join is accepted. Position-only agreement
is insufficient. A transcript version change needs sequence/coordinate
verification rather than silent version stripping.

## Frozen literature identity problems

There are 8,317 unique gene–key rows in the archived count universe. The table
below reports warnings about their identity, not evidence that their papers
or biological variants are false.

| Gene | Literature keys | Keys mapping multiple warehouse alleles | Shared warehouse IDs / literature keys | Recommended identity quarantine |
|---|---:|---:|---:|---:|
| HNF1A | 234 | 6 | 0 / 0 | 9 |
| GCK | 398 | 10 | 2 / 4 | 26 |
| LDLR | 880 | 10 | 43 / 88 | 209 |
| BRCA2 | 6,106 | 545 | 18 / 33 | 118 |
| KCNQ1 | 699 | 7 | 5 / 10 | 69 |

The quarantine recommendation is an explicit canonical WT mismatch or an
archived shared-allele join that points to a different canonical substitution.
It is broader than the earlier GCK missense-only audit because it includes all
consequences. No rows have an out-of-range position in the provided reference
sequences; cDNA-only records are separately flagged as having no protein
position. Such records require a validated canonical-transcript cDNA match.

Archived `vf_variant_ids` are **candidate feature joins**, not independent
genomic truth. The original join allowed a cDNA match when no compatible
protein match existed. Examples include GCK A209V/T209M and W289C/Y289C, and
BRCA2 frameshift G2313fs/missense R2336L. LDLR has many 21-residue-offset pairs,
such as V408M/V429M and W66G/W87G. These are compatible with a mature-protein
versus precursor-numbering issue, but this audit does not apply a universal
offset or redistribute their counts.

Four shared-ID/PMID groups contain multiple retained literature keys: LDLR
W23X/W44X and C201X/C222X in PMID 31102204, LDLR P518=/P539= in PMID 25606447,
and KCNQ1 R147H/R174H in PMID 30758498. Those observations must not simply be
summed after an alias merge. They are potential duplicate evidence, not
adjudicated duplicates; see
[shared_identity_pmid_overlap.csv](shared_identity_pmid_overlap.csv).

The non-destructive primary rule is to keep these source rows in a separate
ledger with the reasons disclosed, and exclude incompatible identities from
automatic count assignment. Canonical-compatible population variants remain
in the population inventory. Resolving mature numbering, contradictory
notations, or overlapping paper counts is a subsequent source-curation task.

## Variant universe and count ownership

Retain every actually observed, QC-eligible population genomic allele in an
inspectable inventory, including population-only alleles with zero assigned
literature affected observations. Do not restrict that inventory to variants
already found in literature, AlphaMissense-scored variants, missense variants,
or warehouse-enumerated SNVs. Include observed indels. Enumerated but
unobserved possible variants, and positive-AN/zero-AC sites, are not additional
observed variant rows. Missing AC is unknown, not evidence of absence.

Keep canonical coding/splice eligibility and consequence flags so two prior
universes can be reported transparently: all observed gene-associated alleles,
and coding/splice alleles plus compatible literature-only variants. The
all-gene distribution answers the broad observed-variation question. The
coding/splice distribution answers a different question and is a useful
comparison for a local protein-density feature. Intronic/UTR rows must not
silently enter or disappear from a purported coding-variant prior. The final
primary scope follows the user's choice; both inventories remain available.

For each population allele, use its normalized GRCh38 identity exactly once.
Attach a literature observation to an allele only when canonical genomic,
protein-point-substitution, or transcript-specific cDNA identity is compatible.
Do not reuse incompatible archived warehouse joins. Synonymous alleles are
particularly important: several distinct DNA changes can have the same
unchanged amino acid at one residue, so a protein key such as `P518=` is not a
unique genomic allele.

Where literature counts belong to a protein-level aggregate with several
compatible genomic alleles, retain an explicit aggregate identity and list its
member alleles. Attach the literature affected/unaffected counts once; attach
each member's population carrier count once; do not also count those members
as independent prior rows carrying the same literature evidence. Never copy
the full literature count onto each possible allele. If an aggregate cannot
be resolved without contradictory identities, quarantine it rather than
inventing an allocation. Distinct frameshifts must retain their full genomic
or HGVS identity; equal starting positions alone do not establish equivalence.

The resulting ledger can contain genomic alleles and explicitly identified
aggregate units. Record `prior_unit_grain` and member IDs, and report how many
population alleles each prior unit represents. This is preferable to silently
claiming allele-level estimates for protein-level evidence. A harmonized
protein-point-substitution aggregation sensitivity can quantify the effect of
that unit choice. Without individual genotypes, sums over different alleles
remain variant-carrier observations, not a proven count of unique people who
carry any allele in the aggregate.

## Carrier counts and the empirical prior

These five genes are autosomal. For a diploid allele in one coherent release,

```text
AC = heterozygote carriers + 2 * homozygote-alt carriers
population carriers = AC - homozygote-alt count
```

Use joint AC and joint homozygote counts from the same source release and
filter policy. Do not subtract a homozygote count from another assay or
release. Require finite integer counts, `AC >= 2 * n_homozygotes >= 0`, and
`AN >= AC`; retain AC, AN, homozygotes, filters and source metadata beside the
derived carrier count. Missing homozygote counts make exact carrier conversion
unknown; an AC-proxy result must be labeled as such. Report AC as the
prespecified historical sensitivity, rather than silently retaining the old
carrier approximation.

gnomAD's joint resource combines exome and genome information. Use joint
counts once; never add joint counts to exomes/genomes, add ancestry rows to
their `all` total, or add versions of the same dataset. Keep the source
collector's resolved joint/assay QC policy explicit. The joint callable-site
AN also differs from a sum restricted to datasets where a variant happened
to be observed. [gnomAD v4.1 release notes](https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1)

The old `gnomad_added` column must be **replaced**, not added again, when the
new release counts are assigned. Preserve it only as a historical comparison.
For one resolved prior unit:

```text
A = assigned literature affected
U = assigned literature unaffected + population carriers
n = A + U > 0
empirical posterior alpha = empirical prior alpha + A
empirical posterior beta  = empirical prior beta  + U
```

Population-only units have A=0 because no affected literature observation was
assigned; their population carriers are unaffected under the adopted
assumption. They receive a proper empirical posterior using positive common
prior parameters, not an improper Beta(0,U).

Refit the requested saturating-weight mean/weighted-MSE empirical Beta moments
on the expanded prior units, using the historical formulas from the structural
plan. Do not substitute `sum(A)/sum(n)` for the empirical-prior mean: that is a
pooled carrier fraction, a different estimand. Adding many observed
population-only variants changes the variant distribution as well as the
pooled counts; that is the correction being tested. Check the moment
conditions and fail explicitly if invalid. Retain the original literature-only
fit as historical evidence rather than overwriting it.

## Machine-readable artifacts

- `literature_identity_flags.csv.gz`: 8,317 unique gene–key rows, with canonical
  WT status, candidate join flags, quarantine recommendations and unchanged
  literature A/U/n. `protein_point_join_allowed` requires a compatible
  missense/nonsense/synonymous point key and WT identity.
- `archived_shared_allele_conflicts.csv`: each shared warehouse ID, its
  canonical enumerated consequence, source key and compatibility result.
- `archived_multi_allele_joins.csv`: all archived keys mapping multiple alleles.
- `shared_identity_pmid_overlap.csv`: possible same-allele/same-paper overlaps.
- `canonical_transcript_audit.csv`, `literature_join_summary.csv`, and
  `literature_identity_summary.csv`: compact source and count summaries.
- `population_count_invariants.csv`: independent checks of the collector's
  30,101 gene–genomic-allele rows. Identities are unique; joint AC/AN/homozygote
  counts are integer-valued and satisfy the count inequalities; all derived
  carrier counts equal AC minus homozygotes; and eligibility equals observed
  positive AC plus the resolved QC flag. QC-eligible counts are HNF1A 2,389,
  GCK 2,292, LDLR 4,218, BRCA2 9,037, and KCNQ1 3,912. Their canonical
  coding/splice subsets contain 1,302, 956, 2,216, 6,811, and 1,716 alleles,
  respectively. No annotated canonical WT mismatches were found in those
  eligible population rows. These are allele counts, not unique people.

All SQLite access was read-only. Canonical sequence hashes and archived
feature hashes are recorded in the tables. No original counts or files were
mutated by this audit.
