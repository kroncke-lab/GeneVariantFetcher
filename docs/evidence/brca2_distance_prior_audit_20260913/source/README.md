# BRCA2 source-count audit, 13 September 2026

**A real source error inflated the BRCA2 missense counts.** An exported ClinVar
variant catalog in PMID 40664060 was treated as a clinical patient table. Each
catalog row became one affected carrier even though it contained no human count.
This folder freezes the proof, exact observation quarantine, genuine clinical
restoration and clinical-key correction needed to rebuild the population union.
It does not establish that every other source is correctly ascertained or that a
structural neighborhood score is clinically calibrated disease penetrance.

## What the source actually contains

[Shi et al., PMID 40664060](https://pmc.ncbi.nlm.nih.gov/articles/PMC12281534/)
studied 100 subjects and reported 22 variants, including 15 BRCA2 variants. The
paper also downloaded a separate ClinVar catalog for classification comparisons
(Methods, “Genetic variants classification”; Supplementary Table 4). That catalog
is a list of variants, not an enumeration of affected people. Primary Table 2
(`tbl0002`) contains the actual clinical sample rows.

The exact frozen source example is BRCA2 A1043V, local DB variant 12496, PMID
40664060. Its source header contains variant identifiers, genomic coordinates,
classifications and review status; its row includes VCV001002452 and
`NC_000013.11:32337482:C:T`. There is no patient identifier or count. Nevertheless:

- `count_provenance.carriers_column_label` is
  `implicit one carrier per clinical row`.
- The upstream database stores total carriers 1, affected 1, unaffected NULL.
- The affected fact is explicitly recorded as `synthesized`, source row 6075.
- The downstream observation builder derives unaffected 0 from total minus
  affected. Its arithmetic is consistent with the incorrect upstream input.

The source DB SHA-256 matches the original frozen protocol summary exactly.
`source_provenance.json` gives the paths, hashes and public XML identity;
`paper_source_db_records.csv.gz` preserves the original records. Independent
current-parser reproduction is recorded in the adjacent `../parser/` folder.
The extraction label “Table 2” on a catalog record must not be confused with the
article's actual clinical Table 2.

## Exact correction and preservation of real patients

| Frozen observation category | Observations / A removed |
| --- | ---: |
| Confirmed catalog missense | 2,979 |
| Confirmed catalog nonsense | 196 |
| Confirmed catalog frameshift | 836 |
| Confirmed catalog in-frame indel | 59 |
| **All confirmed catalog observations** | **4,070** |

The paper supplied 4,071 selected observations in total. One A938fs observation
has a truncated source header and remains explicitly pending. It is outside the
missense/nonsense fits; it is not removed merely because it shares the PMID.
Of the 2,979 removed missense observations, 2,976 belonged to the frozen canonical
missense analysis. The other three were already excluded by identity eligibility.

Five real breast-cancer sample rows in primary Table 2 are restored separately:

| Sample | Protein variant | cDNA | Restored A | Restored U |
| --- | --- | --- | ---: | ---: |
| 6 | R118H | c.353G>A | 1 | 0 |
| 27 | P2276T | c.6826C>A | 1 | 0 |
| 58 | V2076I | c.6226G>A | 1 | 0 |
| 59 | S28N | c.83G>A | 1 | 0 |
| 63 | W2619C | c.7857G>C | 1 | 0 |

These counts come from explicit clinical sample rows with a breast-cancer
endpoint, not from a pathogenicity label. Canonical reference residues match
P51587. The frozen population matches for R118H, S28N and V2076I are consumed by
the corrected clinical union rather than duplicated. P2276T and W2619C have no
observed matching allele in the frozen population inventory. Absence does not
create an unaffected person. W2619C's GRCh38 identity, 13-32362574-G-C, is also
verified by the [NCBI ClinVar allele record](https://www.ncbi.nlm.nih.gov/clinvar/variation/438744/);
the classification itself is not used as a carrier count.

**Authoritative rebuilding contract:** apply the corrected clinical-key table
before rebuilding the clinical/population union. Drop clinical keys whose
corrected A + literature U is zero. Their observed gnomAD alleles remain in the
population inventory and become genomic population units where appropriate.
Then rejoin the five clinical observations using their actual identities.
Removal matches `(key, PMID, original observation variant ID)`; never remove a
variant ID across all papers or delete every observation from the paper.

| File | Role |
| --- | --- |
| `BRCA2_catalog_quarantine_observations.csv` | 4,070 exact original identities and source evidence |
| `BRCA2_clinical_restoration_observations.csv` | Five genuine clinical observations |
| `BRCA2_literature_identity_corrected.csv.gz` | Authoritative complete corrected clinical-key input: 6,110 keys, including 3,644 zero-clinical-evidence keys to drop before the union |
| `BRCA2_clinical_key_corrections.csv.gz` | Compact clinical-key delta view |
| `BRCA2_count_corrections.csv` | **Provisional old-unit accounting only; not a corrected union or fitting input** |

The old-unit accounting gives 4,574 remaining missense units and 457 nonsense
units with evidence. Those are not the final rebuilt counts. Releasing former
protein aggregates yields **4,584 missense units**, A **2,085**, with all
**1,516,774** missense gnomAD carriers retained. The final nonsense count is
**453**, with gnomAD U **13,375**: four released genomic alleles are canonical
frameshifts rather than nonsense alleles under the old clinical X labels. Their
511 population carriers remain in the complete union; they are excluded only
from the nonsense class. The adjacent rebuilt-analysis outputs are authoritative
for these final units and fitted parameters.

## Why the old prior stayed high despite many population carriers

There is no success/failure inversion or missing population denominator in the
frozen missense fit. It had 6,656 units, A 5,056, literature U 1,482 and gnomAD U
1,516,774. Posterior alpha adds A; posterior beta adds both U components.

The accepted historical fit uses `w = 1 - 1/(n + 0.01)` per variant and fits the
mean of `A/n` with those weights. Its variance divides the weighted squared
errors by the number of variants. `audit_source.py` independently reproduces
the old mean **0.14859817** and strength **6.10677443**, without modifying this
formula. A common allele receives approximately the same fitting weight as
another sufficiently observed allele; the weights are not carrier counts.

- The ten largest population denominators contain **93.29%** of gnomAD carriers
  but receive **0.483%** of the prior-fit weight. Their own posteriors are already
  close to zero where their observed fractions support that result.
- There are **2,296 A=1/U=0 units**; **2,084 (90.77%)** are entirely explained by
  the confirmed catalog error. All A=1/U=0 units together contribute 1.097
  percentage points to the old mean under the historical weights.
- The **1,020 A=1/U>0 units** contribute another **9.969 percentage points**.
  Catalog-derived A also contaminated such population-matched variants. The
  source problem therefore extends beyond visually obvious A=1/U=0 singletons.
- **3,041 A=0 gnomAD-observed units** already exist, including 2,971 strictly
  population-only units. The large population counts are present; they do not
  dominate the prior in proportion to their carrier totals.

`missense_count_weight_groups.csv`, `largest_population_denominators.csv` and
`zero_affected_population_units.csv.gz` expose these components. This diagnosis
does not justify changing the fit to force a desired 0.1% mean.

## Bounded residual source screen

`cross_gene_scan.py` checks the selected frozen observation source records against
the exact catalog-header plus implicit-carrier signature, verifying all original
DB hashes. PMID 40664060 occurs in none of the other four genes, and that exact
signature was not found in them.

| Gene | Frozen observations | Same-PMID rows | Exact catalog signature | Implicit clinical-row carrier inference |
| --- | ---: | ---: | ---: | ---: |
| HNF1A | 419 | 0 | 0 | 26 |
| GCK | 616 | 0 | 0 | 60 |
| LDLR | 1,554 | 0 | 0 | 102 |
| KCNQ1 | 1,909 | 0 | 0 | 0 |
| BRCA2 | 8,182 | 4,071 | 4,070 | 5,471 |

A negative signature is not validation. Source quotes can be empty/truncated;
an implicit clinical row can also be a legitimate patient. This screen makes no
additional changes to counts.

Two concrete remaining BRCA2 queues are recorded in
`remaining_BRCA2_source_queue.csv.gz`. They retain **227 canonical missense
observations / A** after the current correction:

- [PMID 36385461](https://pmc.ncbi.nlm.nih.gov/articles/PMC10098510/): 22 canonical
  missense observations, A 22/U 0; 762 observations across all types. The paper
  compiles variants from literature and genomic databases. Its Supplementary
  Table S3C lists variants by country with annotation/database fields and no
  person-count columns. A representative selected A1393fs record used implicit
  clinical-row inference, while that inventory row has no patient/count data.
  D2723G and E3002K also map to inventory entries, but the 22 selected missense
  rows have not been individually adjudicated. Review the exact source route,
  real cohort counts elsewhere and overlap with original publications before
  removing or restoring counts.
- [PMID 33054725](https://pmc.ncbi.nlm.nih.gov/articles/PMC7556962/): 205 canonical
  missense observations, A 205/U 0; 236 across all types. The source contains
  actual tumor sample IDs, including A1439T in TCGA-BR-4184-01 with stomach
  adenocarcinoma. Other supplementary rows explicitly identify somatic variants.
  These are not established invented people; germline status, cancer endpoint
  and overlap with other cohorts need row-level adjudication before they can be
  used as hereditary breast/ovarian disease carrier counts.

The current correction adjudicates one proven source block. It does not certify
the remaining BRCA2 sources or the other four genes. The next source task is
row-level review of these queues, including germline status, endpoint and cohort
overlap; no automatic blanket subtraction of the remaining 227 A is authorized
by this screen.

## Reproduce

Run `audit_source.py` with pandas/numpy to recheck the frozen observations,
quarantine, clinical restoration, corrected key counts, orientation and old fit
arithmetic. Its `checks.json` includes input hashes. Run `cross_gene_scan.py` to
repeat the bounded source-signature screen when the original read-only frozen
DBs are available. Both scripts write only this evidence folder. They perform no
network calls, DB writes, core changes, model refits or git operations.
