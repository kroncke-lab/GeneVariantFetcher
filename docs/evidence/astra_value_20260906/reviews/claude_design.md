## 1. Main design flaw / repair

The proposal is **budget-blind on output caps**. Astra output is $50/M; a repeat of the earlier 32K-cap failure mode on even 2 tasks would exhaust the entire $2.174 envelope before Grok is touched. Conversely, if you shrink caps to fit the budget, any Astra failure on a "harder" task becomes indistinguishable from truncation-starvation rather than a genuine reasoning failure — this is the confound that would invalidate the whole comparison.

**Repair:** compute the cap arithmetically before running anything. Using the 7-paper diagnostic as an empirical anchor (168 fields, patient-row JSON was compact), a per-task output cap of ~3–4K tokens is enough for the described task types (no whole-roster restatement, no 32K regime). Freeze: `Σ(calls × (input_est × $rate_in + cap × $rate_out)) ≤ $2.174` *before* any paid call, and treat that arithmetic as part of the pre-registration, not a runtime adjustment.

A second, subtler risk: because all tasks come from already-opened papers, the human author of the task specs already knows which fields Astra previously got right (78/93 vs Grok's 71/93 source-bound records). That knowledge can unconsciously shape predicate wording toward Astra's known success style. **Repair:** write task specs generically (as already intended) and have a second pass check specs *before* seeing which model is "arm A/B" — i.e., author blind to model assignment, not just blind to output.

## 2. Minimal paired arm/task plan

- **Arm A:** Astra low effort, no retry (medium showed no numeric gain over low — don't pay for it).
- **Arm B:** Grok 4.6 low effort, no retry.
- **Arm C:** Grok 4.6 low effort + one cheap retry, *same* task/source, same schema.

Retry trigger must be **self-referential**, not comparative: fire only when the model's own output contains a null/error on a required field or a transport/cap failure — never triggered by disagreement with Astra or with any gold value. That keeps the retry fair and leak-free. Randomize A/B/C call order per task; lock the task JSON (source excerpt + schema + instructions) before any call is made, identical text served to all three arms.

Do not add a fourth "Astra + retry" arm — the budget doesn't support it, and Astra's retry value is not the open question here (cost is).

## 3. Gate, penalty, stopping rule

- **Utility score per task:** `+1` per unique correct accepted field, `−3` per wrong/unsupported accepted field, `−1` per null miss, `−5` per patient/variant identity error. Freeze these weights now, before results exist.
- **Inclusion gate:** a task counts only if (a) it's genuinely new content not previously scored in any prior calibration (no duplicate packets), (b) it strict-schema-parses without manual repair, (c) the source section it depends on was fully present in the packet (no silently dropped section).
- **Stopping rule:** run all 6–8 predeclared tasks as one fixed batch, no adaptive/early stopping after partial results (avoids optional-stopping bias). If the budget runs out mid-batch, stop and report **inconclusive**, not extrapolated.
- This is an operational go/no-go signal from a small, non-random, opened sample — explicitly not a population estimate.

## 4. Decision mapping and remaining unknowns

- **Default Astra (every paper):** justified only if Astra uniquely resolves a majority of tasks (e.g., ≥5/8) that Grok+retry fails, with zero `−5` identity errors, and the margin plausibly exceeds its ~5–8× cost multiple.
- **Selective Astra (predeclared-unresolved only):** justified if the advantage concentrates in specific task types (e.g., cross-table relative linkage, ambiguous timepoint/endpoint tasks) while Grok+retry matches elsewhere.
- **No Astra:** justified if Grok+retry matches Astra's utility score within noise, given the cost asymmetry.
- **Remains unknown regardless of outcome:** true population-level precision/recall, performance on genuinely unopened papers, any human-reader ceiling comparison, and whether an observed Astra edge on these 6–8 opened, failure-enriched tasks generalizes at all. This experiment can license an operational policy choice under the current budget constraint; it cannot certify accuracy or scale claims beyond that.
