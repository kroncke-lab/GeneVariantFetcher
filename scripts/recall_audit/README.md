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

Claim-verification pilots compare against adjudicated `gold_v2_*` columns when
present and otherwise fall back to the original gold values. Use
`--gold-value-set original` to reproduce the original-gold comparison. The
escalation-queue builder can also consume direct verification records so
high-risk trusted-consensus cases are queued even when debate agrees:

```bash
.venv/bin/python scripts/recall_audit/run_claim_verification_pilot.py \
  --cases-csv /path/to/paper_disagreement_report.csv \
  --run-dir /path/to/scored_run \
  --out-dir /tmp/claim_verify_azure \
  --model azure_ai/gpt-5.4
.venv/bin/python scripts/recall_audit/rescore_claim_verification_records.py \
  --records-jsonl /path/to/claim_verification_records.jsonl \
  --out-dir /tmp/claim_verify_rescore_v2 \
  --gold-value-set v2
.venv/bin/python scripts/recall_audit/run_claim_debate_pilot.py \
  --baseline-records /tmp/claim_verify_rescore_v2/claim_verification_records.jsonl \
  --out-dir /tmp/claim_debate_azure
.venv/bin/python scripts/recall_audit/rescore_claim_debate_records.py \
  --records-jsonl /path/to/claim_debate_records.jsonl \
  --out-dir /tmp/claim_debate_rescore_v2 \
  --gold-value-set v2
.venv/bin/python scripts/recall_audit/build_claim_debate_escalation_queue.py \
  --debate-records /path/to/claim_debate_records.jsonl \
  --verification-records /path/to/claim_verification_records.jsonl \
  --out-csv /tmp/escalation_queue.csv
.venv/bin/python scripts/recall_audit/run_claim_debate_pilot.py \
  --baseline-records /tmp/claim_verify_rescore_v2/claim_verification_records.jsonl \
  --queue-csv /tmp/escalation_queue.csv \
  --final-adjudicator \
  --out-dir /tmp/final_adjudication_sonnet
.venv/bin/python scripts/recall_audit/run_claim_debate_pilot.py \
  --baseline-records /tmp/claim_verify_rescore_v2/claim_verification_records.jsonl \
  --queue-csv /tmp/escalation_queue.csv \
  --final-arbiter \
  --out-dir /tmp/final_arbiter_opus
```

The intended routing is Azure-first for routine work, GPT-5.6 for the canonical
final sniff test, and Anthropic only for optional exception queues:

- Routine triage: `azure_ai/gpt-5.4` (`azure_ai/gpt-5.4-nano` only if deployed on the same endpoint)
- Table routing: `azure_ai/Kimi-K2.6-1`
- Main extraction: `azure_ai/grok-4.3`
- Internal claim verification / debate: `azure_ai/gpt-5.4`,
  `azure_ai/DeepSeek-V4-Pro`, and `azure_ai/Kimi-K2.6-1`
- Canonical final per-paper sniff test (Step 3.8):
  `azure_ai/gpt-5.6-sol` at `xhigh`, default-on, with soft persisted verdicts
- Optional exception adjudication: `FINAL_ADJUDICATOR_MODELS`, defaulting to
  `anthropic/claude-sonnet-5`
- Optional hard-case escalation: `FINAL_ARBITER_MODEL`, defaulting to
  `anthropic/claude-opus-4-8`

`end_to_end_pmid_replay.py` makes a real LLM call and should be used only for a
targeted PMID after a prompt/parser change:

```bash
.venv/bin/python scripts/recall_audit/end_to_end_pmid_replay.py \
  --gene KCNH2 --pmid 24973560 --out /tmp/replay_24973560
```
