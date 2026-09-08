# Phenotype safety and fast-path coverage audit — 2026-09-07

Opened-source calibration and regression work, starting from main `60952867`.
User authorized approximately $50 of API costs and requested adversarial Grok
and AGY CLI assistance. No new holdout or paid extraction arm was opened.

## Outcome

The starting handoff's “identity-only fixed-width fast path” diagnosis for
SCN5A PMID 30059973 is false. A fresh call through `_attempt_extraction`, with
both model-call methods patched to fail, returns **185 variant rows, all with
carrier counts, and no affected/unaffected values**. Table 5's 445 variant
occurrences represent 442 people; 185 is the number of distinct variants.

The source has useful unconsumed phenotype evidence, but it does not define a
single affected/unaffected answer. For E1784K, Table 14 gives 69 carriers:
29 ECG-negative, 13 isolated LQT3, 0 isolated BrS, 17 isolated PCCD, 0 SSS,
0 DCM and 10 overlap. Table 11 gives 60 asymptomatic and 9 syncope at
presentation, but 10 major cardiac events during follow-up. These counts
answer different questions. Copying 69 into affected would fit the reference
and misrepresent these source endpoints. No phenotype synthesis was added.
Fresh success-boundary integration: [source_integration.json](source_integration.json).
Source digest and bounded coverage output: [30059973_coverage.json](30059973_coverage.json).

## Changes

- Selected count headers now pass negative-status, clinical-endpoint,
  timepoint, population and assay exclusions. Previously “No. of unaffected
  patients”, “No. of patients without Brugada syndrome” and “Number of control
  cells” could acquire code-owned phenotype stamps.
- Patient/proband count columns require disease context in the table caption
  or count header. An explicit `Cases (n=...)` column retains its existing
  meaning. A patient noun alone cannot establish the phenotype of a mutation
  catalogue. This also prevents certifying unsupported parser copies of N.
- Caption resolution refuses contradictory reuse of a table anchor. Empty
  anchors do not compete with descriptive captions; repeated identical
  captions and linearized header fragments are handled separately. Resolution
  is cached once per distinct label, instead of rescanning the source per row.
- Deterministic and router short-circuits record pre-guard carrier/A/U supply,
  a source hash, and up to 20 labelled clinical-table candidates with LF-based
  source coordinates. Form feeds do not shift the line numbers. The audit
  makes no completeness or eligibility claim and triggers no model call.
  [`scripts/audit_table_phenotype_coverage.py`](../../../scripts/audit_table_phenotype_coverage.py)
  replays it from explicit source/extraction paths without modifying either.

Literal existing phenotype partitions, carriers and identities are preserved.
Changes apply to new extraction/refresh operations; historical locked runs and
production databases were not rewritten.

## Adversarial review

[Grok's source/code review](reviews/grok_review.md) independently confirmed the
carrier-versus-phenotype distinction and recommended coverage metadata and a
mutation-catalogue refusal before any endpoint-aware enrichment.

AGY's first headless run could not obtain read permission and returned no
substantive review. A [self-contained retry](reviews/agy_selfcontained.md)
completed without tool access. It supported keeping the endpoints separate,
but its response conflated **185 variants with 185 carriers**, invented local
file links, and suggested relaxing intentional mixed-cohort/clinical-header
refusals. Those claims and suggestions were rejected. Its warning that missing
unaffected values need not indicate an error is reflected in the audit's
explicit diagnostic-only interpretation. Integral floats are valid coverage
observations; rejecting `1.0` was not adopted.

[Grok's final review](reviews/grok_final_review.md) agreed with the cohort
changes and raised caption-prefix and audit-bound questions. Empty anchors and
linearized headers now have explicit regression checks; the expanded archive
check confirms those examples resolve. Arbitrarily choosing the longest
nonempty caption was **not** adopted: a longer caption can name a disease
subset of the other table's cohort. A regression test pins that counterexample.
Audit candidates and signals are bounded, with truncation flags. Original
unaugmented source coordinates are intentional and source-hashed, including
CRLF/form-feed tests; synthetic reconstructed-table line numbers would not
identify the cached source file. The audit's incompleteness is explicit. Raw
label whitespace is already normalized by `_table_label` before cache lookup.

Reviewer outputs are advice, not validation. [Receipts](reviews/receipts.json)
retain usage, denied actions, output hashes and reported costs.

## Measurement and limits

The zero-LLM replay uses the two already opened continuation-02/03 candidate
locks, with exactly the same archived predictions, source and scorer. See
[replay/replay_summary.json](replay/replay_summary.json). The cardiac subset is
234 of the 240 gene-paper attempts, and contains 1,118 positive affected
reference rows. These overlapping, previously examined locks are calibration,
not independent confirmation.

The stricter implementation preserves **567/1,118** exact affected values with
the projection, versus **203/1,118** when disabled: **364 new exact values and
46 new wrong values**, plus 13 new values on unmatched/reference-null rows.
Identity, carriers, unaffected and counted extras remain unchanged in this
replay. These are the previous projection's gains, **not a gain from this
iteration**. The attachment's 366-new-exact attribution was already corrected
by the repository's preceding review.

The preregistered 90% precision gate is **not passed** by these opened results:
364/(364+46) = **88.78%** for newly supplied reference-matched values; all supplied
cardiac affected values are 568/658 = **86.32%** exact. Restricting to positive
reference affected values gives 567/631 = **89.86%**. Excluding a difficult
paper after inspecting errors does not convert those into a prospective pass.
The 20129283 source/reference granularity disputes remain curator questions.
No headline, default-off stage or larger-cohort promotion follows.

A [gold-free archive sweep](transfer_smoke.json) scans **2,093 extraction
artifacts** across 11 gene filters, deduplicated by the existing smoke script's
gene/PMID/file-size key. These are not 2,093 distinct papers; archived versions
and duplicate representations remain. The earlier loose patient-header rule
would stamp 263 RYR2 variant rows in PMID 40875405. Its cached source explicitly
says it is a review compiling 221 publications and 964 patients, and its T4
has only `Variant | Number of Patients`. The stricter rule refuses it.
This is a verified unsafe projection prevented, not a gold recall gain.

The sweep also refuses patient-only mutation/annotation catalogues outside the
cardiac genes. The final sweep records 13 classified table entries and 675
projection/stamping decisions. These totals include repeated versions; they
are not independent variants or patients. The withheld projection opportunities are not all established
errors: requiring table-local disease context trades potential supply for
source specificity. No non-cardiac accuracy claim is made.

The canonical [phenotype figure](../../figures/evaluated_phenotype_counts/phenotype_count_recovery_stratified.png)
and its [run membership](../../figures/evaluated_phenotype_counts/phenotype_count_recovery_stratified.json)
were regenerated and visually inspected from unchanged locked inputs. They
represent historical cohorts, not this free replay. No new scored run requires
a companion difference figure.

## Remaining work

A phenotype reader needs a source-bound contract for variant ownership,
cohort, endpoint, timepoint, overlapping categories and disease specificity.
Current disease-word recognition remains broad and is not an ontology. The
coverage heuristic can miss unlabelled tables, wrapped headers beyond its
12-line window, and gene-specific relevance; it cannot certify completeness.
Neither raw metadata nor replay scores are clinical validation. The active
forward checklist remains [TASKS.md](../../../TASKS.md).

## Reproduction

```bash
.venv/bin/python -m pytest tests/unit -q
.venv/bin/python scripts/replay_table_cohort_phenotype.py \
  --run-dir benchmarks/codex_paper_eval/runs/20260905_protocol_cont120_02_candidate \
  --run-dir benchmarks/codex_paper_eval/runs/20260905_protocol_cont120_03_candidate \
  --out-dir tmp/fastpath_replay
.venv/bin/python scripts/table_cohort_smoke.py --root results \
  --genes KCNH2 KCNQ1 SCN5A RYR2 BMPR2 BRCA1 BRCA2 MYBPC3 APOE TTN LMNA \
  --out tmp/fastpath_smoke.json
```

**Validation:** 3,013 offline unit tests and three curated negative cases pass;
all configured formatting/lint hooks and `git diff --check` pass. Test results
and fresh API accounting are recorded in `validation.json` and
`budget.json` beside this report. They are separate from earlier $100/$150
campaign ledgers.
