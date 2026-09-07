# Current-code mixed_gold_cont120_04 production evaluation

This scaffold is pinned to `tranche_04.tsv` at
`d01a79df0fed29ad294b64593c10adc2d79b70c8956cdec5e9abc38a38baf3ba`: **120 gene-paper attempts** /
**111 unique PMIDs** (APOE 1, BRCA1 1, KCNH2 22, KCNQ1 25, RYR2 14, SCN5A 57).

## 1. Extract without opening gold values

```bash
/Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/runs/20260907_protocol_cont120_04_baseline/run_extraction.sh
```

The 6 commands use production `gvf-run` with exact PMID files. Source
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
/Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/runs/20260907_protocol_cont120_04_baseline/lock_and_score.sh
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
