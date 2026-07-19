# Variant_Browser Review Integration

GeneVariantFetcher round-trips with the sibling **Variant_Browser** curation app
(Azure SQL `vb-curation` on `variantbrowser-prod-sql.database.windows.net`). The two halves:

1. **Publish** — push a finished GVF run into the review database so collaborators
   can adjudicate it.
2. **Ingest adjudications** — pull collaborator verdicts back into GVF's gold
   standard as a durable correction overlay.

Both endpoints are **owned by Variant_Browser** (`scripts/gvf_publish.sh` and
the bearer-authenticated `/review/api/gold-standard/` read API). GVF calls the
contract rather than querying Azure SQL tables directly, so the browser remains
responsible for excluding raw, stale, disputed, withheld, and archived calls.

```
  GVF run (<GENE>.db)                              Variant_Browser (Azure vb-curation)
        │                                                       │
        │  gvf-run --publish-review                             │
        │  → gvf_publish.sh GENE DB [DISEASE]                   │
        ├──────────────────────────────────────────────────────▶  papers + carrier
        │                                                       │  evidence + individuals
        │                                                       │        │
        │                                                       │   collaborators adjudicate
        │                                                       │        │
        │  ingest_review_adjudications.py                       │
        │  ◀── authenticated JSON from live Azure DB ───────────┤
        ▼                                                       │
  review_gold.sqlite3   (approved gold + correction overlay + sync audit)
```

---

## 1. Publish a run (GVF → review DB)

Opt-in final step of a run. OFF by default.

```bash
gvf gvf-run KCNH2 --email you@example.com --output results/ \
  --disease "Long QT type 2" --publish-review
```

What it does:

- After the scored `<GENE>.db` is final (including any source-recovery DB), GVF
  shells out to `Variant_Browser/scripts/gvf_publish.sh <GENE> <db> [DISEASE]`.
- The browser script pushes the gene's run into `vb-curation`: papers, per-
  (variant, paper) carrier evidence (counts + source location + verbatim quotes),
  and per-patient phenotype/individual records. It matches GVF `protein_notation`
  to browser variant identities; non-missense (frameshift/nonsense/indel) become
  carrier-only entries; cross-gene and unparseable rows are dropped.
- `DISEASE` is optional (resolved from the gene's review pair if omitted). **The
  gene-disease pair must already exist** in the browser — it is created upstream
  from variantFeatures; `gvf_publish.sh` errors clearly if not.
- Idempotent on the browser side: re-running replaces that gene's carrier data on
  the current snapshot.

**Best-effort:** a missing repo, a non-zero exit, or a timeout is logged and
warned — it never fails the GVF run. So a run that scored fine still produces its
report even if the review DB is unreachable.

### Resolving the Variant_Browser checkout

`--publish-review` needs to find `gvf_publish.sh`. Resolution order:

1. `--review-repo /path/to/Variant_Browser`
2. `GVF_REVIEW_REPO` or `VARIANT_BROWSER_DIR` environment variable
3. The conventional sibling checkout `../Variant_Browser`

If none contains `scripts/gvf_publish.sh`, the step warns and skips.

GVF needs **no** Azure credentials for this — `gvf_publish.sh` loads them from
`Variant_Browser/.env` itself.

Wired in `cli/gvf_run.py` (`step_publish_review`, `_find_review_repo`).

### Curated 101-paper staging publish

For the fixed curated benchmark, publish the benchmark DBs only to the
Variant_Browser staging/review surface, not the public site:

```bash
python scripts/publish_curated_review_set.py \
  --run-root validation_runs/azure_first_101_YYYYMMDD \
  --dataset-label gvf_curated_101_YYYYMMDD \
  --create-pairs
```

Use `--dry-run` first to confirm the per-gene DBs and PMID manifests:

```bash
python scripts/publish_curated_review_set.py \
  --run-root validation_runs/azure_first_101_YYYYMMDD \
  --dataset-label gvf_curated_101_YYYYMMDD \
  --dry-run
```

The helper accepts repeatable `--db GENE=/path/to.db` overrides, or
`--extract-root` when the extraction directory is not under `RUN_ROOT/extract`.
It still calls `Variant_Browser/scripts/gvf_publish.sh`; staging dataset labels,
PMID restrictions, and optional pair creation are passed through environment
variables consumed by the browser repo.

---

## 2. Sync adjudications (live Azure review DB → GVF gold DB)

The normal path pulls the complete current lead-approved dataset directly from
the live Azure-backed review database. There is no browser download and no CSV
transport:

```bash
export GVF_REVIEW_GOLD_TOKEN='<machine credential>'
python scripts/ingest_review_adjudications.py \
  --source-url https://variantbrowser.org/review/api/gold-standard/
```

The token lives in the Variant Browser App Service setting and the GVF GitHub
Actions secret named `GVF_REVIEW_GOLD_TOKEN`; never put it in a URL, log, or
tracked file. The endpoint is HTTPS-only, read-only, `Cache-Control: no-store`,
and uses a SHA-256 dataset revision that GVF independently recomputes before
committing the cache. Redirects are rejected so the bearer token cannot be
forwarded to another host.

For offline recovery only, `--export-csv` still accepts a lead-approved
`export_gold_standard` file. Raw multi-reviewer `export_adjudications` files and
disputed/withheld rows fail closed on both paths.

### The round-trip key

`(gene, source_notation, pmid)`. `source_notation` **is** the GVF
`variants.protein_notation` (e.g. `p.Ser818Leu`). The ingest matches each export
row against a run DB's extracted `(variant, paper)` rows using the *same* canonical
notation bridging the recall scorer uses (`cli/compare_variants.to_canonical_form`),
so `p.Ser818Leu` and `S818L` resolve to the same bucket — and a match here means a
match in scoring.

Run DBs are taken from `--db GENE=path` (repeatable), else auto-discovered as the
newest `<GENE>.db` under `--results-dir` (default `results/`). With `--no-db`, or
when no DB is found for a gene, the overlay still records the adjudication with
`match_status=no_db` (the export is a valid correction record on its own).

### Verdict → action mapping

| Export `verdict`  | Overlay `action`   | Meaning                                            |
|-------------------|--------------------|----------------------------------------------------|
| `confirm`         | `gold_confirmed`   | Extracted (variant, paper) record is correct.      |
| `correct_counts`  | `count_override`   | Use the `corrected_affected/unaffected/total`.     |
| `wrong_variant`   | `false_positive`   | Extracted variant identity is wrong.               |
| `wrong_paper`     | `excluded`         | Variant/paper association is wrong; exclude it.    |
| `missing`         | `followup_missing` | GVF missed it — queue for a human.                 |
| `other`           | `followup_other`   | Free-text `comment` — queue for a human.           |

### Outputs

The live path atomically replaces
`gene_variant_fetcher_gold_standard/adjudications/review_gold.sqlite3`:

- `review_gold_records` — exact accepted Azure records, including browser
  `record_key`, revision/status, source-reviewer account, approving-lead account,
  and the full contract payload.
- `review_gold_overlays` — the resolved correction overlay keyed by
  `(gene, source_notation, pmid)`, with extracted and corrected counts side by
  side for recall, precision, MAE, and RMSE/RMSD.
- `review_gold_sync_state` — source URL, schema version, dataset checksum,
  record count, Azure generation time, and local sync time.
- `review_followup_queue.csv` — everything needing a human: `missing`/`other`,
  plus any actionable verdict whose round-trip key didn't resolve to an extracted
  row (`unmatched`) or that had no DB to verify against (`no_db`).
- `review_adjudications_summary.json` — per-gene counts by action and match status,
  plus net adjudicated count deltas (`net_affected_delta`, `net_unaffected_delta`)
  for matched `count_override` rows.

### Scoring and freshness

Approved gold calls are a **correction overlay**: the extracted value and the adjudicated
value are kept side by side. Nothing rewrites the raw extraction JSON or the run DB.

`scripts/run_recall_suite.py` defaults to `--review-gold-sync auto`: when
`GVF_REVIEW_GOLD_TOKEN` is configured it pulls a fresh snapshot before scoring
and reads the SQLite overlay directly. For release/headline measurements use:

```bash
python scripts/run_recall_suite.py --review-gold-sync required --score
```

`required` fails before scoring if authentication, network access, schema,
status gating, record count, or checksum validation fails. It never falls back
to a stale cache. `off` gives an explicit unadjudicated baseline.

The scheduled **Azure review gold sync** GitHub Action exercises the live
contract daily and on demand. Its ephemeral smoke cache is not a replacement
for the local canonical run DBs; the required sync immediately before a metric
run is the freshness guarantee.

The ingester verifies that the file contains gold approval metadata and only
accepted `gold_standard`/legacy `adjudicated` statuses. `disputed` and `withheld`
audit rows fail closed so they cannot silently alter evaluation metrics.

Approved gold remains a correction overlay: raw extraction JSON and run DB rows
are not rewritten. This keeps model output and human truth separately auditable.

---

## Files

| File | Role |
|------|------|
| `cli/gvf_run.py` (`step_publish_review`) | The opt-in publish step. |
| `cli/__init__.py` (`gvf-run --publish-review`) | CLI flags. |
| `scripts/ingest_review_adjudications.py` | Authenticated pull, validation, matching, and atomic cache sync. |
| `scripts/run_recall_suite.py` | Fresh-sync gate and live overlay consumer. |
| `gene_variant_fetcher_gold_standard/adjudications/review_gold.sqlite3` | Local gold DB + correction overlay + sync audit. |
| `.github/workflows/review-gold-sync.yml` | Daily/on-demand live contract smoke. |
| `Variant_Browser/scripts/gvf_publish.sh` | Publish endpoint (owned by the browser). |
| `Variant_Browser/review/api.py` (`live_gold_standard`) | Machine-authenticated Azure read contract. |

## Coverage

Only KCNH2, KCNQ1, SCN5A, RYR2, and BRCA2 have carrier data in the review DB today.
TTR / ALPL / N4BP2L1 have no GVF extraction yet.
