# New Gene Turn-Key Runbook

This is the proposed cold-start flow for a gene with no manually curated gold standard.

## 1. Prepare

Use the project virtual environment.

```bash
source .venv/bin/activate
```

Recommended credentials:

- `NCBI_EMAIL` and `NCBI_API_KEY`: higher rate limits for PubMed, ClinVar, and metadata calls.
- `OPENAI_API_KEY` or Anthropic vision-model credentials: table/figure extraction and Tier 2 triage.
- `ELSEVIER_INSTTOKEN`: major gain for Elsevier-hosted full text and figures when institutional access is allowed.
- `WILEY_API_KEY` and `SPRINGER_API_KEY`: improves publisher full-text and supplement recovery.

Do not use SSO cookie scraping or institutional paywall credentials unless the run policy explicitly allows it.

## 2. Discover PMIDs From Scratch

Run without `--pmid-file`; explicit PMID lists are for validation against known gold standards and bypass the discovery/filtering behavior that matters for a new gene.

```bash
python -m cli extract GENE --scout-first
```

The workflow should discover PMIDs through PubMind, PubMed, and Europe PMC, then write `results/GENE/<timestamp>/GENE_pmids.txt`.

## 3. Triage

Use Tier 1 keyword filtering and Tier 2 clinical LLM triage with the default thresholds first. Do not tune thresholds until you inspect false positives and false negatives from the new gene's first run.

Key outputs to inspect:

- `dropped_pmids.csv`: top reasons papers were removed.
- Tier 2 logs: malformed or low-confidence decisions should fail open or be reviewed.
- `GENE_pmids.txt` versus filtered PMID count: a very small retained fraction is a warning sign.

## 4. Harvest and Reuse Content

If resuming, point `GVF_RESUME_DIR` at the existing timestamped run directory.

```bash
GVF_RESUME_DIR=results/GENE/<timestamp> python -m cli extract GENE --scout-first
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

## 6. Cold-Start Recovery Layers

Run only layers that are gene-agnostic:

1. ClinVar PMID-citation recovery, after the script accepts `--gene`, `--db`, and auto-resolves the NCBI gene id.
2. PubTator recovery, after the script accepts `--gene`, `--db`, and auto-resolves `CorrespondingGene`.
3. Figure reader on every PMID with a `*_figures/` directory, after the CLI gains `--auto-pmids`.

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

## 8. Expected First-Pass Recall

For a genuinely new cardiac channel gene with no manual recovery DB, expect article coverage to be better than variant-row recall. Based on the current multi-gene validation, first-pass variant-row recall is likely around 30-45% without the recovery stack and could plausibly reach 40-60% after gene-agnostic ClinVar, PubTator, and figure layers. KCNH2's 82%+ result should not be treated as cold-start performance because it includes KCNH2-specific manual recovery.
