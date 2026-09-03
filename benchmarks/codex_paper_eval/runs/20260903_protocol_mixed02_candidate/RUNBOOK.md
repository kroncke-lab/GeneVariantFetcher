# Current-code mixed_gold_02 production evaluation

This scaffold is pinned to `tranche_02.tsv` at
`26f2368b19684ef79c38cb65de8ad243648b3a0b291f81f2de77a1b9f81f32ba`: **49 gene-paper attempts** /
**45 unique PMIDs** (BRCA1 1, KCNH2 8, KCNQ1 10, RYR2 6, SCN5A 24).

## 1. Extract without opening gold values

```bash
/Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/runs/20260903_protocol_mixed02_candidate/run_extraction.sh
```

The 5 commands use production `gvf-run` with exact PMID files. Source
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
/Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/runs/20260903_protocol_mixed02_candidate/lock_and_score.sh
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
