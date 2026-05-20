# Recall Recovery Scripts

Scripts used to enrich a GVF extraction DB with source-independent recovery
signals and score each layer against a gold standard.

## Recovery layers (run in this order — each is additive, none destructive)

After the main extraction pipeline produces a base SQLite DB, layer in:

1. **`recover_paywall_oa.py`** — Explicit-batch source acquisition helper for
   Unpaywall + Europe PMC fallback. It writes recovered `*_FULL_CONTEXT.md`
   files but is not a default recall layer.

2. **`ingest_clinvar.py`** — Query ClinVar via E-Utils for variants citing
   PMIDs already present in the extraction DB. Use `--pmid-source gold` only
   for explicitly labeled diagnostics.

3. **`ingest_pubtator.py`** — NCBI PubTator3 variant-annotation API for
   text-mined mutations from DB-observed PMIDs.

4. **`../extract_figure_variants.py`** — Vision-LLM figure reader. Pulls
   variant tables and mutation maps out of image-only figures the HTML
   stripper dropped. Wires through `harvesting/figure_variant_reader.py`.

5. **`merge_v12_db.py`** — KCNH2-only historical rescue from an older,
   hand-curated DB. This is opt-in only and should not be counted as
   cold-start capability.

## Current ceiling

Do not maintain live recall ceilings in this README. Current missing-PMID and
gap rankings live in `docs/CURRENT_RECALL_STATUS_2026-05-20.md`. In general,
the ceiling is set by source acquisition first: missing Elsevier/ScienceDirect,
Wiley, Karger, Sage/Liebert, or supplement artifacts must become real
`FULL_CONTEXT.md` content before extraction fixes can help.

## Order matters

Migration (`harvesting/migrate_to_sqlite.py`) is destructive: re-inserting
a paper cascade-deletes its variant_papers. If you re-migrate, rerun the
recovery layers afterward. `run_all_layers.py --backup` snapshots the DB before
mutation and `gvf gvf-run` passes that flag by default.
