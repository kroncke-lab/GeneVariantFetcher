# Recall Audit Scripts

Small tools for inspecting scored recall failures without re-running the full
pipeline. By default, DB-backed scripts use the canonical database paths listed
in `docs/RECALL_STATUS.md`, while score-backed scripts use its `Current local
scored artifact` directory. Override an individual input with `--db` or
`--recall-score`; use `--run-dir` or `GVF_RECALL_RUN_DIR` when auditing a
different self-contained run.

An explicit run directory must contain `dbs/{GENE}.db` and `recall_score/`.

## Examples

```bash
.venv/bin/python scripts/recall_audit/paper_disagreement_report.py --run-dir /path/to/scored_run
.venv/bin/python scripts/recall_audit/source_stratified_metrics.py --genes RYR2
.venv/bin/python scripts/recall_audit/build_acquisition_worklist.py --gene SCN5A
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
.venv/bin/python scripts/recall_audit/run_claim_debate_pilot.py \
  --baseline-records /tmp/claim_verify_rescore_v2/claim_verification_records.jsonl \
  --out-dir /tmp/claim_debate_azure
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

The shipped provider default is Anthropic. The measured Azure audit route below
is selected with `MODEL_PROVIDER=azure`; explicit per-tier model variables still
win. GPT-5.6 is used only when the parked final sniff-test experiment is enabled,
and Anthropic exception queues are separately configurable:

- Routine triage: `azure_ai/gpt-5.4` (`azure_ai/gpt-5.4-nano` only if deployed on the same endpoint)
- Table routing: `azure_ai/Kimi-K2.6-1`
- Main extraction: `azure_ai/grok-4.3`
- Internal claim verification / debate: `azure_ai/gpt-5.4`,
  `azure_ai/DeepSeek-V4-Pro`, and `azure_ai/Kimi-K2.6-1`
- Parked final per-paper sniff test (Step 3.8):
  `azure_ai/gpt-5.6-sol` at `xhigh`, default-off since 2026-07-26; re-enable
  together with the companion gate only for a measured experiment
- Optional exception adjudication: `FINAL_ADJUDICATOR_MODELS`, defaulting to
  `anthropic/claude-sonnet-5`
- Optional hard-case escalation: `FINAL_ARBITER_MODEL`, defaulting to
  `anthropic/claude-opus-4-8`

For a targeted single-PMID replay after a prompt/parser change, use the
production path rather than a bespoke script — it is the one that stays in
sync with the pipeline, and it re-uses already-acquired source:

```bash
.venv/bin/python scripts/refresh_run_db.py --gene KCNH2 --pmids 24973560 --dry-run
```

## `gold_source_presence_sweep.py` — acquisition ceiling, blind to predictions

Classifies every gold row of the mixed-gold inventory against everything on
disk for its paper: the article body, every supplement after conversion
(`.docx`/`.xlsx`/`.doc`/`.pptx`/`.zip` through the production
`FormatConverter`, cached outside `corpus/`), and article PDFs. Classes:
`present_in_body`, `present_in_supplement_only`, `present_in_pdf_only`,
`source_absent`, `text_absent_stub_body`, `text_absent_garbled_body`,
`text_absent_figures_present`, `text_absent_notation_inconclusive`,
`text_absent_substitution`. The first four "absent" classes are the hard
acquisition ceiling a recall denominator may exclude; the two unknown classes
must stay penalized. Reads gold rows and source bytes only; never a
prediction or score.

```bash
.venv/bin/python scripts/recall_audit/gold_source_presence_sweep.py \
  --include-unavailable --cache-dir /path/outside/corpus --out-dir <dir>
```

Latest result and interpretation:
`docs/evidence/gold_source_presence_sweep_20260903.md`.

## `fn_root_cause.py` — where did each paper-derived miss disappear?

For a scored `codex_paper_eval` production run, walks every gold row in
`report.json` `missed_gold` through the stages the run left on disk: the
sweep class (corpus), the run-local text, the LLM request (what the model
saw), the LLM response, the extraction JSON, the DB row with its source
layers, and the paper-derived / linkage lanes of `predictions.json`. The leaf
names the owner: `acquisition`, `unknown_notation`, `source_selection`,
`condensing`, `model_missed`, `parser_dropped`, `postprocess_dropped`,
`paper_row_lost_to_linkage_origin`, `projection_dropped`, `matcher`.

```bash
.venv/bin/python scripts/recall_audit/fn_root_cause.py \
  --run-dir benchmarks/codex_paper_eval/runs/<run_id>
```

Writes `diagnostics/fn_root_cause/{summary.md,summary.json,fn_root_cause.tsv}`
inside the run directory. Run it only after the run is locked and scored; it
reads gold via the report.
