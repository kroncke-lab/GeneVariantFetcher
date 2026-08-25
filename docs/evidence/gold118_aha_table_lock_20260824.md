# AHA source-bound gold-118 promotion — 2026-08-24

## Disposition

PROMOTE `benchmarks/codex_paper_eval/runs/20260824_aha_table_sourcebound_gold118`
as the authoritative 118-attempt cardiac lock. Predictions were generated with
gold access disabled, rebound to the exact production sources, hash-locked, and
only then scored. The candidate clears every preregistered non-regression gate
against `20260824_postfix_gold118`.

## What changed

The only new protocol candidate was the scoped AHA collapsed-table repair in
commit `e8e443008af2f7c955c82a2d890ea72bf4b3afd0`. The official AHA DOM for RYR2
PMID 15466642 contains a patient table beneath `figure.table`; the old renderer
kept only its caption. The repaired corpus source is 22,440 bytes, includes all
hidden patient rows, and is bound in the evaluation selection as SHA-256
`59482dba028834c285aa77aa6459998604cd9bb2216862ffeaca67f09b8a4367`.

The first generated scaffold, `20260824_aha_table_gold118`, was never run
because its selection still referenced the stale corpus source. The replacement
source-bound scaffold passed `setup_production_eval.py check` before spend.

## Locked result and promotion gates

| Metric | Accepted postfix lock | Source-bound AHA lock | Gate result |
| --- | ---: | ---: | --- |
| TP / FP / FN | 546 / 284 / 86 | **554 / 283 / 78** | pass |
| Variant recall | 86.3924% | **87.6582%** | pass (+1.2658 pp) |
| Raw precision | 65.7831% | **66.1888%** | pass (+0.4056 pp) |
| F1 | 74.6922% | **75.4255%** | pass (+0.7333 pp) |
| Counted-extra precision | 97.5000% | **97.8799%** | pass |
| Count-bearing-only precision | 93.6937% | **94.5946%** | pass |
| Count-bearing matched rows | 208 | **210** | pass |
| Carrier supplied / MAE | 206 / 0.3301 | **207 / 0.1932** | pass |
| Affected supplied / MAE | 49 / 0.5510 | **56 / 0.3214** | pass |
| Unaffected supplied / MAE | 18 / 0.5000 | **28 / 0.3214** | pass |

The target paper itself scored 8 TP / 0 FP / 0 FN. Unlike the earlier
held-fixed diagnostic, the fresh production run also recovered count fields for
that paper; they are included in the aggregate coverage and MAE above. The
largest remaining identity blocker is unchanged: KCNH2 PMID 29650123 supplies
2/22 identities because the publisher supplement contains no mutation roster.

## Lock and trace integrity

- 118 papers: KCNH2 28, KCNQ1 30, RYR2 30, SCN5A 30;
- all four production manifests completed and are write-time verified;
- 579 calls, 569 successful;
- 2,730,460 provider tokens: 1,951,910 input and 778,550
  total-minus-input-inclusive output/reasoning;
- prediction SHA-256:
  `101076cee9cc5c97900180f352b5c21a1422a2397372c721101f49cba81f7e15`;
- selection SHA-256:
  `cacc9de133408a6d44d87dd6b3a9faa6fd302386bd5e801918983c0f2a613e86`;
- lock SHA-256:
  `97796d25d7e436f1a5f4968e37e8880b8cc26b3ba8fe153226338e182a50e6d6`;
- report SHA-256:
  `7475561d8785d1df8f56d701d317a9a7c842d8db17edc192e5c5d14b17b28b17`.

The four gene jobs ran concurrently. The slowest job, SCN5A, completed in 43.4
minutes. Gold remained disabled through extraction, recovery, trust, source QC,
and report generation; `lock_and_score.sh` opened gold only after exporting and
locking the final predictions.

## Cost and budget

Trace-derived public-price proxy, using the repository's established rate
method:

| Model | Calls (successful) | Input | Visible output | Proxy |
| --- | ---: | ---: | ---: | ---: |
| Kimi K2.6 | 25 (25) | 42,969 | 48,896 | $0.23640455 |
| Grok 4.3 | 122 (115) | 1,219,508 | 308,581 | $2.29583750 |
| GPT-5.6 Sol | 432 (429) | 689,433 | 175,370 | $8.70826500 |
| **Total** | **579 (569)** | **1,951,910** | **532,847** | **$11.24050705** |

The conservative alternate that charges Grok's provider-reported
total-minus-input as output/reasoning is $11.85476455. The ledger uses the same
visible-output method as prior locks for comparability.

Prior charged work was $33.27855555, including the $0.20972125 production micro
replay. Adding this sealed run yields **$44.51906260 used / $55.48093740
remaining** from the user's $100 ceiling. Grok/AGY consultation CLI billing is
not exposed and remains excluded, not asserted to be zero.

## Verification

- focused AHA/browser tests: 30 passed;
- complete unit suite: 2,478 passed;
- Ruff check and format check: passed;
- AGY / Gemini 3.1 Pro maximum supported `high`: SHIP after AHA scoping;
- Grok 4.6 `xhigh`: SHIP after nested-table ownership filtering;
- sealed production wrapper: exit 0;
- lock and score wrapper: exit 0.

Exact BRCA1/BRCA2/BMPR2 150-paper P/R/F1/MAE remain undefined until the
preregistered exhaustive human answer keys are completed. This cardiac
promotion does not alter or release that holdout.
