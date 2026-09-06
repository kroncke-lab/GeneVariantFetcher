# Design consultation notes

All three source-contained CLI consultations completed successfully. Their exact
common brief is `design_prompt.txt`; readable responses and raw receipts sit
beside this file. The Claude consultation used Sonnet with a $3 CLI limit; Grok
and Agy used their CLI paths. No extraction API calls or credential reads were
performed by the consultation launcher. CLI cost is outside the user's API
envelope. These are advisory opinions, not source adjudications or results.

## Shared useful recommendations

- Compare Astra low with both a first Grok pass and a predeclared cheap second
  pass. The earlier different-example schema probes cannot establish model
  superiority.
- Freeze sources, generic prompt, schema, adequate output caps, cost bounds and
  scoring before paid calls. Reserve the full cap and retain failed-call cost
  uncertainty. A cap failure counts as an operational failure; it does not prove
  the model intrinsically cannot reason through the source.
- Separate unique correct source-supported candidates from values the unchanged
  production validator accepts. Only the latter directly establish value in the
  current protocol. A proposed derived-count path is a separate integration
  question.
- A small opened, selected difficult sample cannot justify Astra on every paper,
  population precision/recall, or a human-performance ceiling. Selective utility
  needs an explicit trigger that can run without gold.

## Recommendations requiring correction or qualification

- Agy and Grok propose stopping after the first three tasks match reference.
  This requires looking at reference scores before the fixed batch is locked;
  do not use it in the planned all-output-lock design. Claude's fixed-batch
  recommendation avoids that conflict.
- A null is an allowed, often correct result. Claude and Agy's language about
  nulls on mandatory endpoints must not convert missing evidence into a defect.
  Missing evidence is a retrieval problem; a retry should not be instructed to
  invent an implicit aggregate. Grok explicitly warns against retrying correct
  missing-source abstention, although deciding that correctness itself requires
  care without reference access.
- A parser-only retry predicate cannot detect confident semantic mistakes.
  Thus Grok's proposed Astra-only-after-parser-failure routing may miss the very
  cohort and timepoint errors being tested. A fixed independent second pass or a
  source-based ambiguity trigger can compare semantic recovery without gold
  leakage; the chosen policy must be frozen by the experiment owner.
- The suggested utility weights (usually +1 correct and -3 harmful accepted
  values, sometimes -5 identity errors) are proposed policy preferences, not
  empirically established clinical utilities. Agy's '+20% absolute margin'
  lacks a defined denominator. Record raw gains, regressions and costs alongside
  any explicitly chosen utility sensitivity.
- Claude suggests default Astra if it uniquely solves a majority of the eight
  selected tasks. That result alone would still not justify a default on the
  ordinary paper stream. Grok correctly rejects a universal default inference
  from this sample.
- Reviewer assertions that six to eight calls cannot fit are conditional on
  source length and output cap. Actual serialized input bounds and cap reserves
  decide feasibility. Do not accept a blanket budget claim or starve the output
  merely to fit a desired sample size.
- Grok's one-gain selective gate would show an example of incremental value,
  not reliable ongoing selective benefit. Report its fragility and avoid
  treating a single success as an established general routing category.

This note records consultation disposition at the design stage. It does not
replace the root experiment's frozen execution contract or its final decision.
