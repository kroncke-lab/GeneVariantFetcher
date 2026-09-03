# Current-code mixed_gold_cont120_01 production evaluation

This scaffold is pinned to `tranche_01.tsv` at
`8eb431aa8e2e4ea5cc1723b8d98ba959a3fca8ba5f25d28051a82f3424f33119`: **120 gene-paper attempts** /
**110 unique PMIDs** (APOE 1, BRCA1 1, KCNH2 21, KCNQ1 25, MYBPC3 1, RYR2 14, SCN5A 57).

## 1. Extract without opening gold values

```bash
/Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/runs/20260903_protocol_cont120_01_baseline/run_extraction.sh
```

The 7 commands use production `gvf-run` with exact PMID files. Source
recovery and corpus sync are disabled so this is a calibrated comparison over
the frozen source-available cohort; publication is explicitly disabled. The
gold-free wrapper also disables file-based alias maps whose provenance includes
benchmark gold. The gene processes run concurrently; each writes a separate
`operator_logs/<GENE>.log`, and the launcher fails if any
gene process fails.

## 2. Inspect production completion

Require one successful `RUN_STATUS.json`, final database, and finalized
write-time-verified trace manifest for every gene before proceeding.

## 3. Rebind exact inputs, project, lock, then score

```bash
/Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/runs/20260903_protocol_cont120_01_baseline/lock_and_score.sh
```

This second script binds the exact run-local source material and production
trace manifests into `predictions.json`, applies the collaborator-facing trusted
count and identity projection, and makes **paper-derived rows the primary
score**. ClinVar/PubTator citation-linkage rows stay in the same locked artifact
as an `external_linkage_variants` audit lane and a secondary
`linkage_assisted` comparison. Both views are locked before any gold value is
read. A raw (`--trust-mode all --identity-mode all`) projection may be generated
later only as a clearly labeled diagnostic; it must not replace the locked
primary.
