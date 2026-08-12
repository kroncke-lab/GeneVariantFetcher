# GeneVariantFetcher Gold Standard

Committed baseline snapshot of the literature-curated variant counts used by
GeneVariantFetcher recall analyses. The three SQL-derived gene exports preserve
the legacy Variant_Browser source described below; RYR2 is spreadsheet-derived.
For current lead-approved review corrections, scoring syncs the Variant_Browser
machine API into `adjudications/review_gold.sqlite3` as documented in
[`../docs/VARIANT_BROWSER_INTEGRATION.md`](../docs/VARIANT_BROWSER_INTEGRATION.md).

## What's here

```
gene_variant_fetcher_gold_standard/
├── schema/clinical_counts.schema.json   Column contract for the normalized long format
├── raw/                                 Verbatim per-table dumps from Variant_Browser (audit trail)
├── normalized/                          Tolerant-parsed, source-classified, variant-normalized output
├── adjudications/                       Review correction overlay from the Variant_Browser round-trip
├── qc/                                  Row-count summaries, anomaly reports
├── manifest.json                        Generation metadata (commit, timestamp, row counts, filters)
└── README.md                            (this file)
```

### Genes covered

| Gene  | Source table                          | Per-PMID? |
|-------|---------------------------------------|-----------|
| KCNH2 | `[dbo].[KCNH2_clinical]`              | yes       |
| KCNQ1 | `[dbo].[KCNQ1_clinical_v10]`          | yes       |
| SCN5A | `[SCN5A].[SCN5A_papers]`              | yes       |
| RYR2  | `RYR2_20241129.xlsx` hand-curated literature subset | yes, xlsx-derived |

RYR2's per-paper carrier counts are **not** stored in Variant_Browser. Only an
aggregate variant-level browser table exists there. The shipped RYR2 gold
standard is therefore xlsx-derived: the literature rows produce
`normalized/RYR2_clinical_counts.csv` and `normalized/RYR2_recall_input.csv`,
while population-only gnomAD rows are preserved separately in
`raw/RYR2_gnomad_rows.csv`. Treat RYR2 as an edge-case pilot rather than the
clean Variant_Browser SQL baseline used for KCNH2/KCNQ1/SCN5A.

## Baseline-snapshot regeneration

This command regenerates the committed baseline exports from their legacy source
tables; it does not replace the live review-gold sync. Run it only when those
tables and credentials are intentionally in scope. The exporter uses
Variant_Browser's venv (ODBC + Django + credentials):

```bash
cd /Users/kronckbm/GitRepos/Variant_Browser
set -a; source .env; set +a
venv/bin/python /Users/kronckbm/GitRepos/GeneVariantFetcher/scripts/build_gold_standard_from_varbrowser.py \
  --out /Users/kronckbm/GitRepos/GeneVariantFetcher/gene_variant_fetcher_gold_standard
```

It is idempotent — re-running overwrites `raw/`, `normalized/`, `qc/`, and
`manifest.json` in place.

## Baseline snapshot provenance

- **Source DB:** `sandboxtest.database.windows.net::sandboxtest` (Azure SQL)
- **Source tables** (with schema brackets — important because the SCN5A table
  lives in a non-default schema):
  - `[dbo].[KCNH2_clinical]`
  - `[dbo].[KCNQ1_clinical_v10]`
  - `[SCN5A].[SCN5A_papers]`
- **KCNH2 cross-walk filter:** `[dbo].[all_vars_annotated]` is restricted to
  `isoform='A'` (the canonical NM_000238 transcript) when matching clinical
  variants to browser variants.
- **Generator:** `scripts/build_gold_standard_from_varbrowser.py` — exact git
  commit recorded in `manifest.json`.

## Column meanings (normalized long format)

See `schema/clinical_counts.schema.json` for the authoritative contract. Quick
reference:

| Column                | Meaning                                                                  |
|-----------------------|--------------------------------------------------------------------------|
| `gene`                | Gene symbol (KCNH2 / KCNQ1 / SCN5A).                                     |
| `source_row_ordinal`  | 1-based row position in the matching `raw/<table>.csv` (same ORDER BY). Guaranteed-unique audit key. |
| `source_row_id`       | Natural primary-key value from the source. **Not unique for KCNQ1** — use `source_row_ordinal` for traceability. |
| `pmid`                | Stringified PMID. `null` for non-pubmed rows.                            |
| `source_type`         | `pubmed`, `cohort`, `personal_communication`, or `blank`.                |
| `variant_raw`         | Variant string as written in the clinical table.                         |
| `variant_normalized`  | HGVS-style normalized name. Falls back to raw if normalization fails.    |
| `year`                | Publication / cohort year.                                               |
| `affected_*`          | Per-disease affected counts (e.g. `affected_lqt`, `affected_brs1`).      |
| `unaffected`          | Asymptomatic carriers.                                                   |
| `ambiguous`           | Ambiguous-phenotype carriers (KCNH2/KCNQ1 only).                         |
| `unknown_phenotype`   | Unknown-phenotype carriers (KCNH2/KCNQ1 only).                           |
| `homozygous`          | Homozygous carriers (KCNH2/KCNQ1 only).                                  |
| `total_stated`        | Total carriers as stated by source table (when present).                 |
| `total_computed`      | Sum of all component-count columns.                                      |
| `count_balance_delta` | `total_stated - total_computed`. Non-zero → flagged in `qc/`.            |
| `quality_flags`       | Pipe-delimited tags, e.g. `parse_int_failure:UNA`.                       |

## Recall input (the file GVF actually consumes)

`normalized/<GENE>_recall_input.csv` is the only file GeneVariantFetcher needs
for recall scoring. It carries one row per (variant, PMID) pair restricted to
`source_type='pubmed'`. Columns:

| Column     | Meaning                                                         |
|------------|-----------------------------------------------------------------|
| `variant`  | Normalized variant name (matches GVF's expected shape).         |
| `pmid`     | PMID as a string of digits.                                     |
| `carriers` | Carrier total used for matching. Sum of all affected/unaffected/ambiguous/unknown counts. |
| `affected` | Affected carriers only (sum across all disease-affected columns). |
| `unaffected` | Unaffected/asymptomatic carriers from the source clinical table. |
| `gold_v2_carriers` | Optional manually adjudicated carrier total. Original `carriers` is preserved. |
| `gold_v2_affected` | Optional manually adjudicated affected-carrier count. |
| `gold_v2_unaffected` | Optional manually adjudicated unaffected-carrier count; blank means explicit null when `gold_v2_status` is populated. |
| `gold_v2_status` | Populated only for rows that received manual adjudication or confirmation. |
| `gold_v2_note` | Short rationale for the v2 value. |
| `gold_v2_source` | Local adjudication artifact path. |

Cohort, personal-communication, and blank-PMID rows are dropped here on purpose
— GVF's literature-extraction recall is measured only against PubMed-indexed
sources.

Claim-verification audits default to the v2 value set when `gold_v2_status` is
populated, and otherwise fall back to the original gold columns. Use
`--gold-value-set original` in `scripts/recall_audit/run_claim_verification_pilot.py`
to reproduce original-gold scoring.

Current adjudication semantics:

- `carriers` means variant-positive people, not everyone enrolled, sampled, sequenced, or screened.
- `affected` means variant carriers meeting the paper's disease/phenotype definition, including ECG/QTc-defined affection when the paper defines it that way.
- `unaffected` means explicit unaffected/asymptomatic carriers. When v2 has a populated status and blank `gold_v2_unaffected`, the adjudicated value is intentionally null rather than the original unaffected count.

## Current review-gold overlay (`adjudications/`)

`adjudications/review_gold.sqlite3` holds the versioned current view, immutable
snapshots, record changes, tiers, and reversible exclusions pulled from the
Variant_Browser machine API. Release scoring uses
`scripts/run_recall_suite.py --review-gold-sync required`; raw, disputed,
withheld, stale, or checksum-invalid records fail closed.

CSV export ingest remains a compatibility path for an explicit offline export:

```bash
cd ~/GitRepos/Variant_Browser && set -a && source .env && set +a
python manage.py export_adjudications --out adjudications.csv
cd ~/GitRepos/GeneVariantFetcher
python scripts/ingest_review_adjudications.py --export-csv ~/GitRepos/Variant_Browser/adjudications.csv
```

This writes (idempotently), per gene present in the export:

- `adjudications/<GENE>_review_adjudications.csv` — one row per adjudication, keyed
  by `(gene, source_notation, pmid)` where `source_notation` is the GVF
  `protein_notation`. Keeps **both** the extracted counts and the adjudicated
  `corrected_*` values so metrics can be recomputed without mutating raw rows.
  Verdicts map to actions: `confirm`→`gold_confirmed`,
  `correct_counts`→`count_override`, `wrong_variant`→`false_positive`,
  `wrong_paper`→`excluded`, `missing`→`followup_missing`, `other`→`followup_other`.
- `adjudications/review_followup_queue.csv` — missing/other verdicts plus any
  unresolved matches, for human follow-up.
- `adjudications/review_adjudications_summary.json` — per-gene counts + net count
  deltas.

The scorer consumes the SQLite overlay by tier. The legacy per-gene CSV overlay
is opt-in compatibility input. See
[`../docs/VARIANT_BROWSER_INTEGRATION.md`](../docs/VARIANT_BROWSER_INTEGRATION.md)
for the authoritative sync, schema, and scoring contract.

## Audit / QC outputs

`qc/summary.json` is the headline. The other CSVs surface specific anomalies:

- `qc/duplicate_variant_source_pairs.csv` — the same (gene, variant, PMID) appearing in multiple rows. Some are legitimate (re-extraction, different cohorts in one paper); the report exists so you can decide.
- `qc/unmatched_variants.csv` — clinical rows whose variant didn't crosswalk to the browser table. Useful for finding name conflicts (HGVS vs legacy single-letter, isoform B residues, frameshift notation).
- `qc/count_balance_flags.csv` — rows where component counts (affected + unaffected + ambiguous + unknown + homozygous) don't sum to the stated total. These usually flag a parsing failure or a partial-count row.
- `qc/browser_total_reconciliation.csv` — per-variant: sum of literature carriers vs. the browser's adjusted `total_carriers`. The browser column is curator-adjusted so mismatches are expected; this file makes them inspectable.

## What this is **not**

- **Not** a database. It's a flat-file snapshot. Re-run the exporter for fresh data.
- **Not** the only source of literature signal for GVF. Per-paper FULL_CONTEXT
  files (under `results/<GENE>/<run>/`) remain the input the extractor reads.
- **Not** a CSV the recall analyzer parses directly — it parses the
  `<GENE>_recall_input.csv` form, which is a strict subset.
