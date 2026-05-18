# GVF Turn-Key Audit

Audit timestamp: `20260518_115629`

## Executive Summary

The main discovery/harvest/extract workflow is mostly gene-parameterized, but the recall-recovery stack is not turn-key for a new gene. The current high KCNH2 recall depends on KCNH2-only recovery assets and scripts. A cold-start gene can run through discovery and extraction, but it will not get the same recall lift unless the ClinVar, PubTator, figure-reader, and recovery-driver layers are made gene-agnostic.

SCN5A did not have a usable extraction DB in the existing run artifacts, so this audit could not score a full SCN5A end-to-end DB. The prior SCN5A junk-supplement failure mode is covered by the current denylist test, which passes.

## Reused Artifacts and Scores

| Gene | Reused run | Audit DB | PMIDs | variant_rows | unique_variants | patients | affected | unaffected |
|---|---|---|---:|---:|---:|---:|---:|---:|
| KCNH2 | `results/KCNH2/20260517_074737` | `dbs/KCNH2_base.db` | 69.5% | 35.9% | 44.0% | 53.7% | 43.5% | 67.8% |
| KCNQ1 | `validation_runs/20260517_203904/results/KCNQ1/20260517_204424` | `dbs/KCNQ1_base.db` | 80.0% | 27.7% | 37.7% | 45.2% | 30.2% | 68.4% |
| SCN5A | `validation_runs/20260517_203904/results/SCN5A/20260518_075754` | none | n/a | n/a | n/a | n/a | n/a | n/a |
| RYR2 | `validation_runs/20260517_203904/results/RYR2/20260517_204708` | `dbs/RYR2_base.db` | 57.3% | 30.9% | 34.5% | 47.3% | 40.0% | 79.5% |

Reference prior recovered baselines:

- KCNH2 live recovered DB: 90.8% PMIDs, 82.8% variant_rows, 83.2% unique_variants, 86.0% patients, 87.5% affected, 81.8% unaffected. This matches `docs/NEXT_STEPS.md:7-12` and `recall_metrics/kcnh2_FINAL_BEST_20260517_230019/KCNH2/summary.json`.
- KCNQ1 patched reference: 80.3% PMIDs, 41.5% variant_rows, 59.7% unique_variants, 54.4% patients, 46.8% affected, 68.4% unaffected.
- RYR2 reference is unchanged from audit base.

## A. KCNH2-Specific or Gold-Standard-Specific Dependencies

1. `scripts/recall_recovery/ingest_clinvar.py:13-15` hard-codes the KCNH2 DB, KCNH2 gold CSV, and a Vanderbilt email. It also queries `KCNH2[gene]` at `scripts/recall_recovery/ingest_clinvar.py:36-47` and inserts/selects only `gene_symbol='KCNH2'` at `scripts/recall_recovery/ingest_clinvar.py:108-127`. This is not a multi-gene recovery layer.

2. `scripts/recall_recovery/ingest_pubtator.py:9-14` loads only `KCNH2_recall_input.csv`, filters only `CorrespondingGene:3757` at `scripts/recall_recovery/ingest_pubtator.py:40-47`, and writes to the hard-coded KCNH2 DB at `scripts/recall_recovery/ingest_pubtator.py:58-91`.

3. `scripts/recall_recovery/merge_v12_db.py:1-7` is explicitly a KCNH2 v12 manual-recovery merge. The code is generic once given DBs, but the required input artifact is manually curated KCNH2 state. This layer should never be counted as cold-start capability.

4. `utils/variant_normalizer.py:132-149` contains hard-coded KCNH2 aliases, `utils/variant_normalizer.py:157-172` loads only `data/kcnh2_variant_aliases.json`, and `_lookup_alias` returns no aliases unless `gene.upper() == "KCNH2"` at `utils/variant_normalizer.py:175-181`. The code degrades gracefully for other genes, but other genes lose alias bridging. The public function default is also `gene_symbol="KCNH2"` at `utils/variant_normalizer.py:1320-1342`.

5. KCNH2 IVS splice mapping is special-cased in `utils/variant_normalizer.py:196-205` and applied only when `gene_symbol.upper() == "KCNH2"` at `utils/variant_normalizer.py:1411-1423`. Other genes keep raw IVS notation rather than cDNA-equivalent notation.

6. `cli/compare_variants.py` is a scorer-side matcher, not the extraction path, but it uses a simplified local normalizer (`cli/compare_variants.py:698-711`) rather than `utils.variant_normalizer`. Fuzzy matching is guarded by positional digits (`cli/compare_variants.py:1107-1130`) and defaults to fuzzy mode (`cli/compare_variants.py:2156-2165`). This is defensible for cardiac missense-heavy genes, but it should be explicitly documented as scoring policy, not evidence that the DB is clean.

7. Validation runs with `--pmid-file` exercise harvest/extract but intentionally bypass discovery and Tier 1/Tier 2 filtering: `cli/automated_workflow.py:219-228` skips discovery and `cli/automated_workflow.py:317-324` skips filtering. This is correct for gold validation but does not test the future cold-start path.

8. Reuse behavior exists, but thresholds are inconsistent. `pipeline/steps.py:606-688` copies prior `*_FULL_CONTEXT.md` artifacts larger than 500 bytes, while `harvesting/orchestrator.py:165-205` treats files above `THIN_FULL_CONTEXT_BYTES` as usable and `harvesting/orchestrator.py:79` sets that to 6000 bytes. The user-requested 5 KB reuse rule is not exposed as a flag.

9. The figure reader is gene-parameterized for DB writes (`scripts/extract_figure_variants.py:77-90`) but operationally requires hand-fed PMIDs: `_load_pmid_list` exits without `--pmid`, `--pmid-file`, or `--missing-csv` at `scripts/extract_figure_variants.py:49-67`. There is no `--auto-pmids` and no skip-existing JSON behavior; each run rereads images at `scripts/extract_figure_variants.py:158-185`.

10. Pilot CSVs under `gene_variant_fetcher_gold_standard/pilots/` do not appear to leak into the main path. References are limited to pilot-building and pilot-running scripts unless explicitly passed as `--gold-dir`.

## B. Cold-Start Minimum Flow

See `docs/NEW_GENE_RUNBOOK.md`. The short version:

1. Activate `.venv`.
2. Run without a PMID file so discovery, abstracts, and Tier 1/Tier 2 triage are exercised.
3. Use `GVF_RESUME_DIR` only to resume a timestamped run and reuse existing full context artifacts.
4. Migrate extraction JSONs into `results/GENE/<timestamp>/GENE.db`.
5. Run only gene-agnostic recovery layers. KCNH2 v12 merge is not cold-start.
6. Without gold, judge quality using article coverage, abstract-only fallback rate, zero-variant paper rate, supplement/figure coverage, extraction density, and DB integrity.

## C. Proposed Minimal Patches

1. Parameterize ClinVar recovery.
   - Patch sketch: add `argparse` options `--gene`, `--db`, `--gold`, `--email`, `--api-key`; replace `KCNH2[gene]` with `{gene}[gene]`; pass the gene into all SQL inserts/selects.
   - Rationale: makes the public ClinVar layer usable for any gene.

2. Parameterize PubTator recovery and auto-resolve NCBI gene id.
   - Patch sketch: add `--gene`, `--db`, `--gold`; call NCBI gene search or maintain a small cache mapping symbol to GeneID; filter `CorrespondingGene:<id>` dynamically.
   - Rationale: removes the current `CorrespondingGene:3757` KCNH2 lock.

3. Add a recovery driver.
   - Patch sketch: `scripts/recall_recovery/run_all_layers.py --gene G --db DB --pmc-dir DIR --score-out OUT`, running ClinVar, PubTator, figure reader, and scoring after each. Gate `merge_v12_db.py` behind explicit `--v12-db`.
   - Rationale: a cold-start user should not need to know the historical KCNH2 recovery order.

4. Add `--auto-pmids` and cache behavior to `scripts/extract_figure_variants.py`.
   - Patch sketch: when `--auto-pmids` is set, scan `pmc_fulltext/*_figures/`; skip `out/<pmid>.json` unless `--force`.
   - Rationale: figure recovery should not require a gold missing list.

5. Make variant aliases per-gene.
   - Patch sketch: load `data/{gene_lower}_variant_aliases.json` when present; keep KCNH2 aliases as one file, but do not default all top-level helpers to KCNH2 silently.
   - Rationale: preserves graceful fallback while letting KCNQ1, SCN5A, RYR2, and future genes add curated alias bridges.

6. Expose full-context reuse policy.
   - Patch sketch: replace fixed `THIN_FULL_CONTEXT_BYTES = 6_000` with a setting/CLI option and align prior-run consolidation with `full_context_needs_retry`.
   - Rationale: avoids silently reusing or redownloading thin artifacts differently across layers.

7. Add a first-run QC report.
   - Patch sketch: `scripts/report_cold_start_qc.py --run-dir results/GENE/<timestamp>` summarizing full-context status, abstract-only fallback rate, zero-variant papers, supplement conversions, figure directories, extraction JSON count, and DB row counts.
   - Rationale: gives a no-gold standard acceptance test.

## SCN5A Note

The current denylist includes `Media-Pack-2024.pdf` and `JournalCatalog2026.pdf` patterns in `harvesting/supplement_scraper.py:1398-1402`, checked by `tests/unit/test_generic_scraper.py:267-280`. The targeted test passes. The existing SCN5A artifacts still do not contain a completed extraction DB, so full SCN5A recall remains unscored in this audit.
