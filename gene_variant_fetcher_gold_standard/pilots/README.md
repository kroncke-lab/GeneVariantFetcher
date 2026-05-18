# Gold-Standard Pilot Set

Deterministic pilot subsets generated from `gene_variant_fetcher_gold_standard`.
The default clean pilot genes are `KCNH2, KCNQ1, SCN5A` with up to `10`
PMIDs per gene.

## What To Compare Against

- Use `normalized/<GENE>_recall_input.csv` here for scoring pilots. These files
  have the same reduced schema as the full recall inputs: `variant`, `pmid`,
  `carriers`, `affected`, `unaffected`.
- Use the parent package's `normalized/<GENE>_clinical_counts.csv` for audit and
  provenance. Those files retain `source_row_id`, `source_row_ordinal`,
  `source_type`, raw/normalized variants, disease-specific affected columns,
  unaffected counts, and QC flags.
- `pilot_pmids.csv` records why each PMID was selected and gives aggregate
  counts for quick API smoke-test triage.

## Full vs Partial

The parent KCNH2, KCNQ1, and SCN5A `clinical_counts` files are full flat
snapshots of their exported source tables. The recall inputs are intentionally
partial: PubMed-indexed scoring rows only, collapsed to the counts GVF compares.
RYR2 is excluded by default because it is xlsx-derived and mixes literature rows
with population-only gnomAD rows; add it with `--include-ryr2` when testing that
edge case specifically.

## Regenerate

```bash
python scripts/build_gold_standard_pilots.py
```

Run live API smoke checks only when the local Python environment has GVF
dependencies installed:

```bash
python scripts/run_gold_standard_api_pilots.py --live
```
