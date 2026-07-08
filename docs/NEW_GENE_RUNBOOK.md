# New Gene Turn-Key Runbook

This is the cold-start flow for a gene with no manually curated gold standard.
The turnkey entry point is:

```bash
gvf gvf-run <GENE> --email you@example.com --output ./results [--disease "<phenotype>"]
```

Source recovery (paywall + supplement acquisition) runs by default in
`gvf-run`. Add `--no-source-recovery` for a fast PMC/free-text-only pass. The
optional `--disease` flag scopes discovery and Tier 2 filtering to a
gene-disease pair. The steps below detail what `gvf-run` runs.

Scope: use this for a new gene-disease pair where recall cannot be scored yet.
Use `RECALL_REFRESH_RUNBOOK.md` for an existing scored run and
`END_TO_END_RECALL_RUN.md` when recreating the recall workflow on another
machine.

## 1. Prepare

Use the project virtual environment. Local setup lives in `docs/QUICKSTART.md`;
credential details live in `docs/API_KEYS.md`.

```bash
source .venv/bin/activate
```

Recommended credentials:

- `NCBI_EMAIL` and `NCBI_API_KEY`: higher rate limits for PubMed, ClinVar, and metadata calls.
- One LLM provider key (`ANTHROPIC_API_KEY`, `OPENAI_API_KEY`, or Azure AI/OpenAI credentials):
  table/figure extraction and Tier 2 triage.
- `ELSEVIER_INSTTOKEN`: major gain for Elsevier-hosted full text and figures when institutional access is allowed.
- `WILEY_API_KEY` and `SPRINGER_API_KEY`: improves publisher full-text and supplement recovery.

Do not use SSO cookie scraping or institutional paywall credentials unless the run policy explicitly allows it.

## 2. Discover PMIDs From Scratch

Run without `--pmid-file`; explicit PMID lists are for validation against known gold standards and bypass the discovery/filtering behavior that matters for a new gene.

Prefer the `gvf <command>` console script. When it is not on `PATH`, use
`.venv/bin/python -m cli <command>` as the fallback.
Use the lower-level `extract` command here only when debugging individual stages;
the preferred production path is the `gvf-run` command at the top of this runbook.

```bash
gvf extract GENE --email "$NCBI_EMAIL" --output results/ --scout-first
```

The workflow should discover PMIDs through PubMind and PubMed (Europe PMC is added only when `USE_EUROPEPMC=1`), then write `results/GENE/<timestamp>/GENE_pmids.txt`.

## 3. Triage

Use Tier 1 keyword filtering and Tier 2 clinical LLM triage with the default thresholds first. Do not tune thresholds until you inspect false positives and false negatives from the new gene's first run.

Key outputs to inspect:

- `dropped_pmids.csv`: top reasons papers were removed.
- Tier 2 logs: malformed or low-confidence decisions should fail open or be reviewed.
- `GENE_pmids.txt` versus filtered PMID count: a very small retained fraction is a warning sign.

## 4. Harvest and Reuse Content

If resuming, point `GVF_RESUME_DIR` at the existing timestamped run directory.

```bash
GVF_RESUME_DIR=results/GENE/<timestamp> \
  gvf extract GENE --email "$NCBI_EMAIL" --output results/ --scout-first
```

The harvester should reuse non-empty `pmc_fulltext/*_FULL_CONTEXT.md` files and retry only empty, abstract-only, thin, or publisher-shell artifacts. Avoid re-downloading already valid content.

## 5. Extract and Migrate

After extraction JSONs are written, migrate to SQLite.

```bash
python harvesting/migrate_to_sqlite.py \
  --data-dir results/GENE/<timestamp> \
  --db results/GENE/<timestamp>/GENE.db
```

Check that the DB has non-zero `papers`, `variants`, `variant_papers`, and, for clinical papers, `individual_records` or `penetrance_data`.

If source recovery lands after the initial run, do not patch SQLite rows
directly. Refresh from source artifacts so the run remains reproducible:

```bash
python scripts/refresh_run_db.py \
  --gene GENE \
  --run-dir results/GENE/<timestamp> \
  --replace-db
```

This selects stale or under-counted PMIDs from reusable source files, rewrites
canonical extraction JSONs, rebuilds a fresh DB from the complete extraction
directory, then runs DB-observed recovery layers. A gold standard is optional;
without one, recall scoring is skipped but the same ClinVar, PubTator, and
figure recovery layers still run.

`gvf-run` also writes a no-gold source QC bundle under
`results/GENE/<timestamp>/source_qc/`:

- `source_acquisition_worklist.csv`: every run PMID classified as `fetch`,
  `refresh_replay`, `manual_or_blocked`, `inspect`, or `none`.
- `fetch_input.csv`: PMIDs without usable run-local full text; pass this to
  `scripts/fetch_paywalled.py`. Known blocked publisher cases are kept out of
  this file and remain in the worklist as `manual_or_blocked`. This also
  includes papers whose main text explicitly points target-gene variants to a
  missing supplemental table/body (`missing_variant_supplement` in the worklist).
- `source_override.csv`: PMIDs with usable source text that should be replayed
  with `refresh_run_db.py --source-override-csv`.
- `source_acquisition_summary.json`: coverage for usable full text, selected
  fetch PMIDs, selected source-refresh PMIDs, manual/blocked PMIDs, and
  zero-variant full-text PMIDs. It also counts
  `pmid_coverage.missing_variant_supplement`, which is a source-acquisition
  bucket rather than a generic extractor miss.

The default source QC denominator is the harvest/extraction surface, not raw
PubMed discovery. Pass `--include-discovery-pmids` only for diagnostic audits
that intentionally include every discovered PMID.

The no-gold acquisition and staged refresh loop runs by default as part of
`gvf-run`:

```bash
gvf gvf-run GENE --email you@example.com --output results/
```

This runs `source-qc`, fetches `source_qc/fetch_input.csv` with
`scripts/fetch_paywalled.py`, writes
`source_qc/acquisition_outcome_summary.json`, emits
`source_qc/fetched_source_override.csv`, and calls
`scripts/refresh_run_db.py --stage-extractions --only-forced-pmids` with both
existing-source and fetched-source override CSVs. Recovery layers run against the
refreshed DB unless `--skip layers` is also supplied.

For a fast PMC/free-text-only pass, skip source recovery:

```bash
gvf gvf-run GENE --email you@example.com --output results/ --no-source-recovery
```

For a manual no-gold audit on an existing run:

```bash
python scripts/recall_audit/source_acquisition_audit.py \
  --gene GENE \
  --run-dir results/GENE/<timestamp> \
  --out results/GENE/<timestamp>/source_qc/source_acquisition_worklist.csv \
  --summary-out results/GENE/<timestamp>/source_qc/source_acquisition_summary.json \
  --fetch-input-out results/GENE/<timestamp>/source_qc/fetch_input.csv \
  --source-override-out results/GENE/<timestamp>/source_qc/source_override.csv
```

After running `fetch_paywalled.py`, summarize selected-vs-successful
acquisition. With no gold, this reports worklist coverage; with `--gold` or
`--report`, it also reports PMID recall:

```bash
python scripts/recall_audit/summarize_acquisition_outcome.py \
  --gene GENE \
  --worklist results/GENE/<timestamp>/source_qc/source_acquisition_worklist.csv \
  --fetch-summary results/GENE/<timestamp>/source_qc/fetch/summary.json \
  --out results/GENE/<timestamp>/source_qc/acquisition_outcome_summary.json \
  --source-override-out results/GENE/<timestamp>/source_qc/fetched_source_override.csv
```

Replay successful acquisitions and existing refresh candidates without mutating
the canonical extraction directory:

```bash
python scripts/refresh_run_db.py \
  --gene GENE \
  --run-dir results/GENE/<timestamp> \
  --source-override-csv results/GENE/<timestamp>/source_qc/source_override.csv \
  --source-override-csv results/GENE/<timestamp>/source_qc/fetched_source_override.csv \
  --only-forced-pmids \
  --stage-extractions
```

After replay, rerun the outcome summary with the refresh summary attached to
record the PMIDs that were actually accepted into refreshed extraction JSON:

```bash
python scripts/recall_audit/summarize_acquisition_outcome.py \
  --gene GENE \
  --worklist results/GENE/<timestamp>/source_qc/source_acquisition_worklist.csv \
  --fetch-summary results/GENE/<timestamp>/source_qc/fetch/summary.json \
  --refresh-summary results/GENE/<timestamp>/refresh_<timestamp>/refresh_summary.json \
  --out results/GENE/<timestamp>/source_qc/acquisition_outcome_summary.json \
  --source-override-out results/GENE/<timestamp>/source_qc/fetched_source_override.csv
```

## 6. Cold-Start Recovery Layers

Run only layers that are gene-agnostic. `gvf-run` and `refresh_run_db.py` call
the combined layer driver automatically; for manual runs use:

```bash
python scripts/recall_recovery/run_all_layers.py \
  --gene GENE \
  --db results/GENE/<timestamp>/GENE.db \
  --pmc-dir results/GENE/<timestamp>/pmc_fulltext \
  --backup
```

The layer driver uses DB-observed PMIDs by default:

1. ClinVar PMID-citation recovery with `--pmid-source db`.
2. PubTator recovery with `--pmid-source db`; override `--gene-id` only if auto-resolution is wrong.
3. Figure reader on every PMID with a `*_figures/` directory.

Do not use `--pmid-source gold` for cold-start claims. That mode is diagnostic
only because it uses the target answer set to choose enrichment PMIDs.

Do not apply KCNH2's v12 manual recovery DB to a new gene. That layer is a curated KCNH2 artifact, not turn-key pipeline behavior.

## 7. QC Without a Gold Standard

Use internal signals instead of recall:

- Article coverage: downloaded plus reused full-context PMIDs divided by filtered PMIDs.
- Abstract-only fallback rate: should be reviewed when high, especially for high-priority papers.
- Zero-variant papers: inspect papers with good full text but no extracted variants.
- Supplement coverage: count expected supplement links, downloaded supplements, converted tables, and failed conversions.
- Figure coverage: count PMIDs with figure directories and distinct variants found by figure reader.
- Extraction density: variants per full-text paper and clinical count rows per variant-bearing paper.
- Source mix: distinguish PMC, publisher-free, API, abstract-only, and paywalled records.
- DB integrity: no duplicate canonical variants for the same gene/PMID; `variant_papers` rows should connect to valid `papers` and `variants`.

When a gold standard exists, distinguish PMID recall for PMIDs selected for
acquisition from PMID recall for PMIDs that actually landed usable full text.
`summarize_acquisition_outcome.py` reports both under
`pmid_recall.selected_for_fetch_download` and
`pmid_recall.usable_fulltext_downloaded`. After replay, pass
`--refresh-summary` as well; this adds `pmid_recall.source_refresh_attempted`
and `pmid_recall.source_refresh_successful`, the strict count of PMIDs accepted
into refreshed extraction JSON.

## 8. Expected First-Pass Behavior

For a genuinely new gene with no manual recovery DB, expect article/source
coverage to be better than variant-row and count recall. The current validation
status is summarized in `docs/RECALL_STATUS.md`, but do not
use those scored genes as hard expectations for a new disease area.

For no-gold runs, judge progress by internal QC: source completeness, supplement
coverage, zero-variant full-text papers, variant density, count plausibility,
and whether high-priority paywalled papers are represented by real full text
rather than abstract/paywall stubs.

## 9. Publish for Review and Ingest Adjudications (optional round-trip)

Once a run is scored, you can push it into the Variant_Browser review database so
collaborators can adjudicate it, then pull their verdicts back into the gold
standard. Both directions are opt-in.

Publish on the run (best-effort; a publish failure warns but never fails the run):

```bash
gvf gvf-run GENE --email you@example.com --output results/ \
  --disease "<phenotype>" --publish-review
```

Ingest the adjudications afterward (writes a correction overlay under
`gene_variant_fetcher_gold_standard/adjudications/` that keeps both the extracted
and the adjudicated counts):

```bash
cd ~/GitRepos/Variant_Browser && set -a && source .env && set +a
python manage.py export_adjudications [--gene GENE] --out adjudications.csv
cd ~/GitRepos/GeneVariantFetcher
python scripts/ingest_review_adjudications.py --export-csv ~/GitRepos/Variant_Browser/adjudications.csv
```

The gene-disease pair must already exist in the browser (created upstream from
variantFeatures). See `docs/VARIANT_BROWSER_INTEGRATION.md` for the full contract,
the round-trip key, and the verdict→action mapping.
