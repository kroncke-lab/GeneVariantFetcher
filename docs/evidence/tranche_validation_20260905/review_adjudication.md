# Independent CLI reviews and corrections

Grok (`grok-4.6-build`) and Agy (`gemini-3.1-pro-high`) each completed a single,
facts-only review after pair 02, while pair 03 was running. No tools or code
edits were requested. Input is `review_prompt_02.txt`; raw responses and usage
are retained in `grok_review_02_raw.json` and `agy_review_02_raw.json`.
Grok reported $0.00476986 as CLI cost; Agy did not report dollars. These are
recorded separately from Azure trace-derived API estimates because a CLI's
notional usage price does not establish an additional API charge. No Anthropic
model was used.

Both reviews support keeping acquisition upside separate from frozen-source
reading results, avoiding a causal addition to the old 67.97% cohort score,
and checking whether gains transfer beyond previously scored papers. Their
recommendations are advisory, not evidence or new acceptance criteria.

Corrections to the reviewers' responses:

- Grok incorrectly says the -0.17 pp recall / -1.29 pp precision lower bounds
  fail non-inferiority limits of -1 / -2 pp. Both bounds pass. The identity
  endpoint fails because observed recall gain is below +1 pp.
- Agy incorrectly describes omissions as “zero error.” They are evaluated as
  a zero **prediction**, so a positive reference incurs its full count as
  error. Zero-valued references and abstaining from a wrong supplied value can
  lower error; supply and exact supplied counts are reported alongside MAE.
  This is not a scoring loophole unique to either arm.
- Agy's claim that one +0.70 pp comparison “proves” a plateau is unsupported.
  The modest observed effect and uncertainty warrant a restrained forecast,
  not proof that reading improvements are exhausted.
- Neither reviewer may invent an acceptance rule for pooled tranches or the
  exploratory A/U endpoint. The registered per-tranche rules remain unchanged.
- The prompt abbreviated the BRCA2 conflict as main-table count 2 versus
  supplement count 3. Grok misread those values as table numbers. Actual
  locations are main Table 3 and Supplementary Table 4.
- The prompt did not state that the RYR2 97/62/26 patient-row calculation is
  already implemented. Subsequent code/source cross-check confirmed it is;
  this is an existing endpoint disagreement, not future parser upside. The
  revised forecast must not count matching that reference as a promised gain.

The final estimate will use both locked pairs and source-supported remaining
work, not the reviewers' numerical opinions.

## Final review after both pairs

Both CLIs completed another facts-only review with the same respective models
and no tools/code edits: `review_prompt_final.txt`, `grok_review_final_raw.json`
and `agy_review_final_raw.json`. Grok reports 24,409 tokens and $0.00995214;
Agy reports 14,673 tokens and no dollar amount. Across these four new reviews,
Grok's two notional charges sum to $0.014722; Agy billing is unavailable. This
is not added to the Azure API ledger without evidence of an additional invoice.
No Anthropic model was used for either review round.

Accepted: neither paired result supports a general recall increase; distinguish
selected source substitutions from transfer; report the H558R influence and
unresolved person/allele semantics; lower the original cohort's A/U forecast.
The final forecast keeps 76% only as a rounded selected-recovery scenario with
net-neutral remaining identity work, moves its A/U central value to the already
observed 12.15%, and uses a conditional 10–17% range. For broader existing-source
reading, 0–5% count improvement is only a planning allowance and no reliable
recall lift is budgeted. Unimplemented source work has no corpus-wide yield
estimate. The initially proposed 15% A/U central / 10–20% range was not retained.

Corrections and judgments:

- Agy calls the 25.2% FN reduction an artifact of 25163546 alone. It is arithmetic
  from **two** selected runs: +11 TP from 20031634 and +20 from 25163546. It is
  useful for the named old cohort but cannot establish automated generalization.
- Agy says 11.55 percentage points of the pooled improvement come from H558R.
  That subtracts percentages with different denominators. The row contributes
  408/3,412 = 11.96 points on the original baseline denominator; removing the
  row entirely yields a different baseline, 3,004, and a 3.43% residual gain.
- Agy generalizes tranche-03 supply decline to overall coverage regression.
  Tranche 03 does decline, but pooled A/U supply rises 414 → 427 and exact
  supplied values rise 349 → 366. Omissions still incur full positive-reference
  error; withdrawing a previously wrong count is one possible improvement.
- “Fails to preserve baseline precision” is too broad as a conclusion from a
  failed non-inferiority test. Pair 03 fails that bound, while pair 02 passes;
  pooled precision is nearly unchanged with uncertainty. Neither demonstrates
  a reliable precision benefit or guarantees preservation.
- Agy's proposed return to 68–70% mixes the old candidate with the new paired
  baseline and ignores source recoveries for that specific old cohort. Its
  proposed 2–8% general A/U range is also uncalibrated. Neither is adopted as a
  measured estimate. Grok's narrower distinction between selected-cohort
  bookkeeping and remaining-work benefit is the appropriate interpretation.
- Grok's final response correctly distinguishes pair-02 non-inferiority passes
  from its failed +1-point gain threshold. Gate decisions are unchanged.

Reviewer responses are advisory. All final arithmetic is reproduced from
locked artifacts rather than accepted on model authority.
