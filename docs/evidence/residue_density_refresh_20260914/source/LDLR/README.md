# Bounded LDLR source adjudication for the FH residue-density refresh

The three largest affected-count sources **among frozen canonical missense keys** cover 71 observations, 197 affected counts and 19 unaffected counts. Two have proven endpoint/partition errors; the third supports its retained clinical-FH counts. This is not a complete audit of the other 242 papers contributing to these canonical missense keys.

| PMID | Canonical missense rows | Original A/U | Adjudication |
| --- | ---: | ---: | --- |
| [29802317](https://pubmed.ncbi.nlm.nih.gov/29802317/) | 55 | 73/0 | MI case/control study; FH status unestablished, quarantine the observation contribution |
| [25463123](https://pubmed.ncbi.nlm.nih.gov/25463123/) | 13 | 62/19 | Detected-variant table partitions index patients/affected relatives; correct to 81/0 |
| [31491741](https://pubmed.ncbi.nlm.nih.gov/31491741/) | 3 | 62/0 | Retain explicit counts among clinically diagnosed FH patients |

The complete clinical input changes from **4,595 A / 216 U to 4,535 A / 193 U** across all types. Within the original canonical missense keys, totals change from **1,400 A / 114 U to 1,346 A / 95 U**, before the population union is rebuilt. Counts here are source observations, not a claim of unique individuals across publications.

## Authoritative correction contract

`clinical_input.csv.gz` is the complete 880-key LDLR table in the original literature-identity-flags schema, with additional baseline/delta audit columns. It retains 835 keys with clinical evidence. Its SHA-256 is:

`34807f665c9307f1cf725ddc7f2e77b6d6a799bbc884f10b880a375ec9d3eaa3`

The root rerun must drop zero-clinical-evidence keys before rebuilding the clinical/population union. All observed population alleles remain; they may change from clinical protein aggregates into genomic population units. Neither unknown FH status nor a missing clinical observation is a new unaffected person.

`observation_decisions.csv.gz` records exact original key + PMID + variant_id, original A/U, removal and restoration amounts, nullable adjudicated counts, source row/line, genomic allele or cDNA, endpoint, rationale and original extraction provenance. `corrected_observations.csv.gz` independently reconstructs the full key totals. `clinical_key_deltas.csv` is a convenience view of affected keys. Five unchanged Greek rows are explicit retain decisions, so source counts already containing affected relatives are not increased twice.

## MI endpoint proof: PMID 29802317

The full article describes 9,956 myocardial-infarction cases and 8,373 controls. Its Methods define MI using clinical records and cardiac tests. The Discussion explicitly says FH was not checked/excluded. The supplement's Cases and Controls therefore refer to MI status. An FH database identifier in another column establishes prior annotation, not the phenotype of the sampled carriers. [Primary article record](https://pubmed.ncbi.nlm.nih.gov/29802317/)

All **64 selected rows (59 missense, 5 frameshift; 83 A / 0 U)** match the exact frozen target row, including chromosome 19, position, REF, ALT, substitution and both count columns. The frozen key quotes name Table S8. The cleaned source contains the complete table, around lines 4298–4440. Exact matches prevent accidental joining to PCSK9 substitutions with the same protein label. Their contribution to the FH primary input is quarantined; the observed MI counts are preserved in the ledger. No control-only rows are newly extracted in this bounded repair.

## Greek index/relative partition proof: PMID 25463123

The Subjects section defines the cohort as clinically diagnosed heFH: 262 pediatric index patients and 299 relatives. Table 2 is explicitly the list of detected LDLR variants; its columns give the nucleotide change, precursor and mature protein changes, index count, and relative count. The **26 source rows sum to 140 index patients and 141 relatives**, exactly matching the paper's conclusion describing molecularly identified heFH patients. Thus these per-variant relatives are mutation-positive clinically affected carriers, not an unaffected denominator. [Primary article record](https://pubmed.ncbi.nlm.nih.gov/25463123/)

The exact cDNA and precursor protein reference/position are checked for each selected observation. Mature protein numbering is retained as source context and is never used as an assumed offset. `Greek_Table2_exact_rows.csv` freezes this proof. The **16 selected rows (13 missense, 2 nonsense, 1 frameshift)** change from **67 A / 23 U to 90 A / 0 U**. Eleven partitions change and five rows already contain the correct total. All observed carrier totals are preserved. This is a diagnosis-based cohort count; it does not estimate lifetime untreated biochemical penetrance or establish independence of relatives.

## Retained source: PMID 31491741

The Methods describe 801 clinically diagnosed HeFH patients, including 650 unrelated patients, with explicit diagnostic criteria. The Results give C338S=23, D433H=20 and L568V=19 among the unrelated FH cohort. These three exact missense observations remain 62 A / 0 U. `retained_source_review.csv` provides the original identities and source lines (59 and duplicated Results at 77). CAD analyses elsewhere in this paper are not substituted for FH case status. [Primary article record](https://pubmed.ncbi.nlm.nih.gov/31491741/)

## Scope, reproduction and remaining limits

The first metadata screen selected top sources across all variant types. For LDLR, two large contributors mainly use legacy protein labels already excluded by canonical WT checks. The relevant final screen instead selected the top three sources **after restricting to the frozen canonical missense key set**; `canonical_missense_top3_coverage.csv` and its original-observation companion document that denominator. No legacy identities are repaired here. Exact two-paper corrections also cover the selected other-type/noncanonical rows supported by the same tables; this does not place those rows into a missense prior.

Remaining clinical data are still ascertained from heterogeneous studies; other source endpoint errors and duplicate people may remain. This repair does not validate the entire literature, relax identity exclusions, modify original observations/DBs, refit a prior, or change the geometry/kernel. The resulting neighborhood density remains a model feature rather than clinically calibrated variant penetrance.

Reproduce from this repository:

```sh
/Users/kronckbm/GitRepos/BayesianPenetranceEstimator/.venv/bin/python docs/evidence/residue_density_refresh_20260914/source/LDLR/adjudicate.py
```

`checks.json` records baseline observations/identity flags, source hashes, exact count assertions, and the deterministic output hash. The compact data use LF and deterministic gzip; original full texts/DBs remain in ignored source caches.
