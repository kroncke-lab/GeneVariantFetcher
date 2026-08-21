# Diagnostic lock — not the acceptance run

This gold-value-blind production run is hash-locked and scored, but it is not
the acceptance run for the table-attribution change.

While the run was still blind, PMID 22677073 exposed a general false-negative
path: markdown table candidates retain `Header`, `Previous row`, `Target row`,
and `Next row` in one audit packet, while the current-study candidate filter
searched that entire packet for gene symbols. A valid target-gene row could
therefore be rejected as `gene_ambiguous` solely because an adjacent row named
another gene. The main LLM recovered the observed variants, but deterministic
short-circuit, hint, and additive-recovery paths cannot rely on that recovery.

The locked diagnostic result is 546 TP / 288 FP / 87 FN (65.47% precision,
86.26% recall, 74.44% F1). It used 2,606,143 tokens. Counted-extra precision is
97.85%; count-bearing-only precision is 94.47%.

The lock also confirms a pre-existing benchmark defect that must not be solved
in runtime extraction code: KCNH2 PMID 14642689 is an unrelated angiotensin-II
receptor / atrial-fibrillation paper, contains no KCNH2 variant study, and was
correctly rejected by the circuit breaker, yet the locked gold assigns it one
KCNH2 false negative. Preserve this run for historical comparison; adjudicate
the benchmark record separately with source provenance.

Required next step: implement a gene-agnostic, candidate-local table scope fix,
cover neighboring rows, same-row modifier genes, and blank gene-cell
continuations with synthetic tests, then run a fresh blind lock before any
BRCA1/BRCA2/BMPR2 collaborator extraction.
