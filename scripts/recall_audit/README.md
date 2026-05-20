# Recall Audit Scripts

Small tools for inspecting scored recall failures without re-running the full
pipeline. By default, scripts read the scored baseline from `docs/RECALL_STATUS.md`.
Override with `--run-dir` or `GVF_RECALL_RUN_DIR` when auditing a different run.

The run directory must contain `dbs/{GENE}.db` and `recall_score/`.

## Examples

```bash
.venv/bin/python scripts/recall_audit/pmid_side_by_side.py --gene KCNH2 --pmid 19160088
.venv/bin/python scripts/recall_audit/failure_taxonomy_report.py --gene KCNH2 --out /tmp/kcnh2_taxonomy.csv
.venv/bin/python scripts/recall_audit/gene_filter_audit.py --gene RYR2 --out /tmp/ryr2_gene_filter.csv
.venv/bin/python scripts/recall_audit/study_wide_n_detector.py --gene RYR2
.venv/bin/python scripts/recall_audit/multicohort_collapse_detector.py --gene KCNH2 --pmid 19160088
.venv/bin/python scripts/recall_audit/acquisition_status_per_pmid.py --gene SCN5A --summary
```

`end_to_end_pmid_replay.py` makes a real LLM call and should be used only for a
targeted PMID after a prompt/parser change:

```bash
.venv/bin/python scripts/recall_audit/end_to_end_pmid_replay.py \
  --gene KCNH2 --pmid 24973560 --out /tmp/replay_24973560
```
