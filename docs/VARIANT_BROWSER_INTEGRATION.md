# Variant_Browser Review Integration

GeneVariantFetcher round-trips with the sibling **Variant_Browser** curation app
(Azure SQL `vb-curation` on `variantbrowser-prod-sql.database.windows.net`). The two halves:

1. **Publish** — push a finished GVF run into the review database so collaborators
   can adjudicate it.
2. **Ingest adjudications** — pull collaborator verdicts back into GVF's gold
   standard as a durable correction overlay.

Both endpoints are **owned by Variant_Browser** (`scripts/gvf_publish.sh` and
`manage.py export_adjudications`). GVF only calls them — it does not duplicate the
DB schema, Azure credentials, or the GVF→browser variant matching. Keep that
contract: changes to the publish/match logic belong in Variant_Browser.

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
        │  ◀── export_adjudications --out adjudications.csv ────┤
        ▼                                                       │
  gene_variant_fetcher_gold_standard/adjudications/   (correction overlay + follow-up queue)
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

## 2. Ingest adjudications (review → gold standard)

Export the collaborator verdicts from Variant_Browser, then fold them into GVF's
gold standard as a correction overlay.

```bash
# 1. Export from the browser (writes the round-trip CSV).
cd ~/GitRepos/Variant_Browser
set -a && source .env && set +a
python manage.py export_adjudications [--gene KCNH2] --out adjudications.csv

# 2. Ingest into GVF's gold standard.
cd ~/GitRepos/GeneVariantFetcher
python scripts/ingest_review_adjudications.py --export-csv ~/GitRepos/Variant_Browser/adjudications.csv
```

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

Written under `gene_variant_fetcher_gold_standard/adjudications/` (override with
`--out-dir`). Idempotent — re-running replaces every gene present in the export:

- `<GENE>_review_adjudications.csv` — the **correction overlay**. One row per
  adjudication, keyed by `(gene, source_notation, pmid)`. It keeps **both** the
  pipeline's extracted counts (`extracted_carriers/affected/unaffected`, captured
  from the run DB) **and** the adjudicated corrections (`corrected_*`), so recall /
  precision / MAE can be recomputed without mutating raw extracted rows.
- `review_followup_queue.csv` — everything needing a human: `missing`/`other`,
  plus any actionable verdict whose round-trip key didn't resolve to an extracted
  row (`unmatched`) or that had no DB to verify against (`no_db`).
- `review_adjudications_summary.json` — per-gene counts by action and match status,
  plus net adjudicated count deltas (`net_affected_delta`, `net_unaffected_delta`)
  for matched `count_override` rows.

### Correction overlay, not mutation

Adjudications are a **correction overlay**: the extracted value and the adjudicated
value are kept side by side. Nothing rewrites the raw extraction JSON or the run DB.

> **Note (current limitation):** the recall scorer (`cli/compare_variants.py`,
> `scripts/run_recall_suite.py`) does **not yet consume** this overlay — nor the
> older, hand-curated `gold_v2_*` columns in `normalized/<GENE>_recall_input.csv`.
> The overlay preserves everything needed to recompute metrics, but joining it into
> the live scorer is a deliberate next step, not automatic today.

This is the live, export-driven cousin of the older hand-curated `gold_v2_*`
overlay columns in `normalized/<GENE>_recall_input.csv` (which were populated by a
one-off manual probe, since retired).

---

## Files

| File | Role |
|------|------|
| `cli/gvf_run.py` (`step_publish_review`) | The opt-in publish step. |
| `cli/__init__.py` (`gvf-run --publish-review`) | CLI flags. |
| `scripts/ingest_review_adjudications.py` | The adjudication ingest. |
| `gene_variant_fetcher_gold_standard/adjudications/` | Correction overlay + queue + summary. |
| `Variant_Browser/scripts/gvf_publish.sh` | Publish endpoint (owned by the browser). |
| `Variant_Browser/.../export_adjudications.py` | Export endpoint (owned by the browser). |

## Coverage

Only KCNH2, KCNQ1, SCN5A, RYR2, and BRCA2 have carrier data in the review DB today.
TTR / ALPL / N4BP2L1 have no GVF extraction yet.
