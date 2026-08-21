# Scripts

Prefer the `gvf` CLI for normal work. These scripts are support tools for
refreshing runs, scoring recall, recovering source, and preparing review data.

## Routine Recall And Refresh

Use these when updating an existing scored run or measuring recall.

| Script | Purpose |
| --- | --- |
| `refresh_run_db.py` | Rebuild a run DB from selected source artifacts and recovery layers. |
| `refresh_recall.py` | Refresh and score the tracked recall genes. |
| `run_recall_suite.py` | Score one or more gene DBs against normalized gold inputs. |
| `manage_review_gold.py` | Inspect immutable Azure gold snapshots, changes, tiers, and reversible exclusions. |
| `recall_mae.py` | Compare rows-mode MAE and run-to-run score changes. |
| `trust_report.py` | Inspect the trust-gate two-tier DB (trusted/quarantine counts, rate, by reason). |
| `precision_sample.py` | Draw + score an adjudicated precision-calibration sample (Wilson CI) for gate thresholds. |
| `phenotype_value_precision.py` | Measure exact-match precision of supplied carrier/affected/unaffected values on an existing locked prediction set. |
| `build_status_dashboard.py` | Render the self-contained historical strategy/status dashboard from scored artifacts. |

## Source And Corpus

Use these for full-text, supplement, figure, and corpus maintenance.

| Script | Purpose |
| --- | --- |
| `fetch_paywalled.py` | Authenticated paywall/full-text recovery. |
| `fetch_elsevier_supplements.py` | Elsevier supplement recovery. |
| `fetch_linked_supplements.py` | Recover supplement links already recorded in paper markup but missed by publisher/API acquisition. |
| `fold_supplements.py` | Fold downloaded supplements into full-context source files. |
| `fulltext_acquisition_pass.py` | Source-acquisition-only pass without extraction. |
| `build_source_corpus.py` | Build or update the local ignored source corpus. |
| `corpus_to_harvest.py` | Copy corpus source back into harvest-shaped run folders. |
| `stage_corpus_figures.py` | Stage corpus figure artifacts for extraction/review. |
| `extract_figure_variants.py` | Extract variant evidence from figures. |
| `fetch_gold_figures.py` | Fetch figures for gold-standard audit papers. |
| `fetch_reference_sequences.py` | Fetch reference protein sequences used by validation code. |
| `gold_source_worklist.py` | Build source recovery worklists from gold inputs. |
| `ezproxy_relogin.py` | Bootstrap or refresh the persistent EZproxy session used by authenticated recovery. |

## Recovery Layers

These mutate or enrich an existing SQLite DB. Back up the DB first unless the
calling workflow already does it.

| Script | Purpose |
| --- | --- |
| `recall_recovery/run_all_layers.py` | Run DB-observed recovery layers. |
| `recall_recovery/ingest_clinvar.py` | Ingest ClinVar-observed variants. |
| `recall_recovery/ingest_pubtator.py` | Ingest PubTator-observed variants. |
| `apply_count_classifier.py` | Apply count-role classifier output. |
| `apply_count_outlier_guard.py` | Quarantine implausible carrier counts. |
| `apply_somatic_germline_qc.py` | Flag somatic/germline ambiguity. |
| `backfill_paper_metadata.py` | Backfill paper metadata columns. |
| `backfill_source_layers.py` | Backfill source-layer provenance. |
| `dedup_db.py` | Deduplicate DB rows. |
| `quarantine_fp.py` | Quarantine known false-positive rows. |

## Review And Curation

Use these when preparing human review packets or Variant Browser round trips.

| Script | Purpose |
| --- | --- |
| `ingest_review_adjudications.py` | Import Variant Browser adjudication exports. |
| `build_curation_packet.py` | Build cold-start manual-curation packets. |
| `score_curation_packet.py` | Convert reviewed packets into recall inputs and score them. |
| `build_adjudication_sheet.py` | Build reviewer spreadsheet inputs. |
| `build_adjudication_html.py` | Build static reviewer HTML from adjudication data. |
| `build_evidence_cards.py` | Build evidence-card rows for review. |
| `build_llm_trace_manifest.py` | Hash/index per-call LLM and decision traces for adjudication (thin CLI over `utils.llm_trace`). |
| `build_llm_trace_html.py` | Build a self-contained per-paper browser view of LLM traces (thin CLI over `utils.llm_trace_html`; bounded bodies, shards large runs). |
| `recover_counts.py` | Standalone additive fill of NULL per-variant counts (`pipeline.count_recovery`). `--dry-run` first; refuses nothing but writes only NULL slots, lands them as `quarantine`, and backs up the DB. |
| `publish_curated_review_set.py` | Publish the curated extraction-eval subset to Variant Browser staging after scope exclusions; not the full-gene publish path. |

## Recall Audit

The `recall_audit/` subdirectory contains report builders and targeted audit
pilots. Start with `recall_audit/README.md`.

| Script | Purpose |
| --- | --- |
| `collaborator_readiness_audit.py` | Audit fixed collaborator cohorts for structural/source/trace/trust readiness while leaving publication on hold. |
| `final_check_triage_report.py` | Produce a no-LLM shadow report for the parked final-check triage predicate. |

## Gold Input Builders

These regenerate committed gold-standard inputs from external curation sources.
They are not part of normal extraction.

| Script | Purpose |
| --- | --- |
| `build_gold_standard_from_varbrowser.py` | Export gold counts from Variant Browser SQL. |
| `build_ryr2_gold_standard_from_xlsx.py` | Build RYR2 gold counts from the source spreadsheet. |
| `build_gold_standard_pilots.py` | Build pilot gold-standard packages. |
| `run_gold_standard_api_pilots.py` | Run API pilot checks for gold-standard packages. |

## Dev Or Historical Tools

These are retained for diagnosis or historical comparison. Do not use them for
new recall work unless a runbook points here.

| Script | Status |
| --- | --- |
| `test_insttoken_unlock.py` | Manual credential probe. |
| `check_ezproxy.py` | Manual institutional-access probe. |
| `check_wiley_api.py` | Manual Wiley TDM credential/connectivity probe; requires network. |
| `discover_recall.py` | Historical recall-discovery helper. |
| `enrich_from_variantfeatures.py` | Manual enrichment helper for legacy VariantFeatures data. |
| `run_priority_extraction.py` | Historical priority-walk driver support; prefer `gvf gvf-run` unless reproducing that experiment. |
| `targeted_land.py` | Historical acceptance-gated DB promotion helper. |
| `replay_cap_trip_extractions.py` | Deprecated 2026-08-20; retained only to reproduce existing dense-table/scanner cap-trip run artifacts until `gvf-run` owns equivalent targeted replay. |
| `retry_failed_extractions.py` | Deprecated 2026-08-20; retained for targeted repair of an existing run's failure CSV. New runs use the retry path in `pipeline.steps.extract_variants`. |
| `smoke_azure_models.py` | Opt-in live connectivity/JSON smoke test for configured Azure LLM deployments. |

## Cleanup Rule

New scripts should be added to exactly one section above and linked from a
runbook or test. One-off local drivers belong outside Git.
