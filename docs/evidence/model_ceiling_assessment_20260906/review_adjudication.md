# Adjudication of the three CLI reviews

Dated 2026-09-06. This is a review of advice against saved evidence and current
code, not a vote or a fresh model extraction benchmark. The prompts included
the prior source audits and model configuration. They did not include the
subsequent raw-response probe or fully describe the common patient-row hook;
those findings qualify the reviews below.

| Reviewer | Adopt or investigate | Reject or qualify |
| --- | --- | --- |
| Claude (`sonnet` resolved to `claude-sonnet-5`; receipt also lists a small Haiku call) | Check whether RYR2 18929323's numbers were omitted by the model or cleared by code; do not treat acquisition labels as validated causality; test model and workflow separately. | Flat results with the same model do not argue against model capability as a bottleneck. The 76.8% decomposition is **carrier error**, not A/U error. Count recovery is default off. A selected small experiment is not guaranteed to fit the remaining budget without usage estimates. |
| Grok (`grok-4.6-build`) | Separate source readiness, conditional reading and gold incompleteness; reject a corpus/human ceiling from the selected 76% scenario; preserve uncertainty about transfer. | The 9 and 16 model-missed rows are from different cohorts, not a worsening trend. “Models last” is stronger than the evidence supports; cross model and workflow factors. Its suggested one-point stopping rule is not a registered gate and is unsuitable as a corpus decision from a small targeted panel. “No Anthropic” is stronger than the user's preference to use it sparingly. |
| Agy (requested `gemini-3.1-pro-high`; response has no resolved-model receipt) | Make cohort semantics and damaged table structure explicit before counting; challenge the assumption that complete files alone solve reading. | Diagnostics do not prove models are irrelevant or that recall will stall at 79%. The 74 BRCA2 extra rows are not all validated. The 15% error reduction is measured but concentrated, not an illusion. The hypothesized rule clearing the RYR2 carrier counts is contradicted by the raw response. New clinical-reader calls are not zero-cost merely because source/output caches exist. Its proposed 145 cases mix row-level diagnostics from different cohorts and cannot be treated as 145 distinct paper attempts. |

All three reviews require an architecture qualification: deterministic success
already reaches the specialized patient-row phenotype derivation hook. It is a
narrow closed-table audit, not a general clinical missing-fact reader. The claim
that there is no post-identity phenotype processing is incorrect.

The direct [RYR2 probe](primary_model_abstention_probe.json) resolves the most
useful disagreement. Grok 4.3's raw structured output already omits the carrier
counts, despite retaining n=13 and n=6 in its notes. The final extraction does
not clear them; verification is skipped below the risk threshold. The diagnosis
therefore includes primary-model over-abstention and incomplete triggering.
Whether stronger reasoning, clearer field rules, or both improve this on fresh
runs remains untested.

## Invocation and receipts

Each CLI was called once with a bounded prompt and a 270-second subprocess
timeout. No CLI requested subagents or ran a production extraction. The Python
launcher is local scratch under `validation_runs/model_ceiling_assessment_20260906/`;
the durable artifacts below contain the review inputs and outputs.
The repository hook normalized the prompt files' final blank line; their
substantive content is unchanged.

- Claude: `--print --model sonnet --effort high --max-budget-usd 1.00 --tools '' --strict-mcp-config --safe-mode --no-session-persistence --output-format json`.
  [Prompt](claude_prompt.txt), [response/usage](claude_raw.json),
  [stderr](claude_stderr.txt). Successful one-turn response;
  `total_cost_usd=0.12369680000000001`.
- Grok: `--verbatim --output-format json --max-turns 1 --reasoning-effort low --no-subagents --disable-web-search --tools '' --permission-mode plan` with a prompt file and scratch working directory.
  [Prompt](grok_prompt.txt), [response/usage](grok_raw.json),
  [stderr](grok_stderr.txt). Successful one-turn response;
  `total_cost_usd=0.00881076`.
- Agy: `--print <prompt> --model gemini-3.1-pro-high --effort high --mode plan --sandbox --disable-slash-commands --output-format json --print-timeout 4m`.
  [Prompt](agy_prompt.txt), [response/usage](agy_raw.json),
  [stderr](agy_stderr.txt). `SUCCESS`, one turn, 18,563 reported total tokens;
  dollar cost unavailable.

The CLI estimates are not reconciled API invoices. The active extraction
campaign ledger is unchanged; no unavailable cost has been recorded as zero.
