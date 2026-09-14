# Bounded source-context screen of HNF1A, LDLR and KCNQ1

This is the initial warning screen used to choose source adjudications, not a validation stamp for three unchanged genes. The selected observation inventories and original DBs were opened read-only and checked against their earlier provenance hashes. The screen combines the three largest sources by affected count across all types with explicit title/quote warnings for catalog, somatic, liver-neoplasm and opposite cardiac endpoints. A keyword hit is a review trigger, not an automatic exclusion.

| Gene | Original observations A/U | Initial screen rows A/U | Initial screen restricted to canonical missense rows A/U |
| --- | --- | --- | --- |
| HNF1A | 419; 762/115 | 18; 155/49 | 9; 131/35 |
| LDLR | 1,554; 4,595/216 | 16; 739/1 | 5; 64/1 |
| KCNQ1 | 1,909; 8,392/2,264 | 374; 2,311/200 | 244; 1,477/177 |

These are the **original frozen** counts, and row counts are observations rather than people. Coverage is available in `coverage.csv`; source titles/flags in `papers_in_scope.csv`; every selected original observation's source metadata in `source_metadata_screen.csv.gz`; warnings in `keyword_flagged_rows.csv.gz`. The snapshot contains all 3,882 original rows, while `reviewed_source_scope.csv.gz` contains only the initial selected subset. Some metadata fields are sparse and can miss errors.

## Adjudication handoffs

- **HNF1A:** the initial screen found type-2-diabetes counts labeled for MODY, somatic tumor observations, and liver-adenoma/diabetes phenotype partitions. The exact source adjudication and complete corrected clinical ledger are in `../HNF1A/`. It changes all-type totals from 762/115 to 633/83, while retaining the two established type-2-study counts only in a separately labeled broader-diabetes sensitivity input. Some flagged frameshift/splice records remain outside the current missense/nonsense fits; see that source report.
- **LDLR:** the initial raw ranking was not adequate for the plotted subset because large historical contributors failed canonical WT identity checks. A second bounded screen selected the top three sources *within canonical missense keys*: PMIDs 29802317, 25463123 and 31491741, covering **71 rows, 197 A and 19 U**. Exact source review proved MI/FH endpoint confusion and index/affected-relative partition errors. `../LDLR/` contains the authoritative corrections and retained-source proof. Large excluded legacy identities were not repaired.
- **KCNQ1:** the screen found explicit AF/short-QT endpoints and inconsistent use of all carriers, clinical diagnoses and symptomatic events. In particular, the G589D cohort in PMID 21244686 is not evidence that all 492 carriers had symptomatic events; other papers explicitly partition symptomatic/asymptomatic carriers. These warning rows were passed to a separate KCNQ1 adjudicator. The authoritative KCNQ1 decision ledger and final scope belong to its source folder, not this unchanged metadata snapshot.

This initial script intentionally does not mutate counts. Later corrections live in the gene-specific folders. `screen_checks.json` describes the screen at the time it ran; its statement that counts are unchanged refers to this scanner only, not the subsequent refresh.

Reproduce the scanner with:

```sh
/Users/kronckbm/GitRepos/BayesianPenetranceEstimator/.venv/bin/python docs/evidence/residue_density_refresh_20260914/source/other_genes/scan_source.py
```

No source re-extraction, source DB changes, external LLM calls or new variant identity repair were performed by this screen. Neither the initial coverage table nor the bounded repairs establish that all remaining counts across these genes share a validated clinical endpoint.
