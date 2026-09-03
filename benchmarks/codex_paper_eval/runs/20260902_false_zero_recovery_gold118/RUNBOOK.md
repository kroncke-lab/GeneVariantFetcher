# Current-code gold_120 production evaluation

This scaffold is pinned to `tier2_gold_120.tsv` at
`a70389566304ad6573d666721139aaf79c8185631b95c7e1ea86231fad53097e`: **118 gene-paper attempts** /
**114 unique PMIDs** (KCNH2 28, KCNQ1 30, RYR2 30, SCN5A 30). The name remains
`gold_120` for continuity; KCNH2 PMID 10086972 is intentionally absent after the
erratum identified it as the wrong paper.

## 1. Extract without opening gold values

```bash
/Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/runs/20260902_false_zero_recovery_gold118/run_extraction.sh
```

The four commands use production `gvf-run` with exact PMID files. Source
recovery and corpus sync are disabled so this is a calibrated comparison over
the frozen source-available cohort; publication is explicitly disabled. The
four gene processes run concurrently, matching the accepted 2026-08-20 test;
each writes a separate `operator_logs/<GENE>.log`, and the launcher fails if any
gene process fails.

## 2. Inspect production completion

Require one successful `RUN_STATUS.json`, final database, and finalized
write-time-verified trace manifest for every gene before proceeding.

## 3. Rebind exact inputs, project, lock, then score

```bash
/Users/kronckbm/GitRepos/GeneVariantFetcher/benchmarks/codex_paper_eval/runs/20260902_false_zero_recovery_gold118/lock_and_score.sh
```

This second script binds the exact run-local source material and production
trace manifests into `predictions.json`, applies the collaborator-facing trusted
count and identity projection, locks predictions before any gold value is read,
and only then invokes the scorer. A raw (`--trust-mode all --identity-mode all`)
projection may be generated later only as a clearly labeled diagnostic; it must
not replace the locked primary.
