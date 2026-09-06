# Two further tranches: recall did not improve reliably; count gains are concentrated

Completed 2026-09-06. Two paired 120-attempt tranches cover 240 gene-paper
attempts / 220 distinct articles, each read once per protocol. Both fail the
registered identity and carrier-count acceptance rules. Pooled recall is almost
unchanged; the apparent 15% A/U error improvement is dominated by one disputed
408-count row. A separate publisher-supplement recovery finds 20 previously
missed identities but supplies no counts. The accepted headline and default-off
recovery stages remain unchanged.

## Measured effects

Every cell below compares the frozen historical nine-file reader with the same
current candidate. These are fresh extractions from fixed source availability.
A/U error includes omitted identities/counts as zero predictions; NULLs remain
NULL in stored output. Relative error changes use identical asserted gold rows
within each comparison.

| Cohort | TP / FP / FN, baseline → candidate | Recall | Precision | Combined A/U absolute error | Relative A/U error reduction |
| --- | --- | --- | --- | --- | ---: |
| Tranche 02: 120 attempts, 571 gold identities | 441 / 182 / 130 → 445 / 180 / 126 | 77.23% → 77.93% | 70.79% → 71.20% | 1,531 → 1,487 | 2.87% |
| Tranche 03: 120 attempts, 959 gold identities | 773 / 133 / 186 → 767 / 135 / 192 | 80.60% → 79.98% | 85.32% → 85.03% | 1,881 → 1,414 | 24.83% |
| Pooled, descriptive: 240 attempts, 1,530 identities | 1,214 / 315 / 316 → 1,212 / 315 / 318 | 79.35% → 79.22% | 79.40% → 79.37% | 3,412 → 2,901 | 14.98% |

Pooled recall change is **−0.13 percentage points**, with a descriptive
PMID-cluster 95% interval of **−0.98 to +0.61 points**. Precision change is
−0.03 points (interval −1.44 to +1.12). The first tranche's +0.70-point recall
change misses the +1-point discovery threshold despite passing both
non-inferiority bounds. The second tranche's −0.63-point change also fails
both non-inferiority bounds. Both carrier-MAE upper bounds cross zero.
Tranche 03 is additional discovery/calibration, not confirmation of a success.

| Pooled count field | Absolute error, baseline → candidate | Supplied values | Exact supplied values |
| --- | ---: | ---: | ---: |
| Carriers | 2,401 → 1,859 | 861 → 875 | 770 → 785 |
| Affected | 2,616 → 2,526 | 241 → 249 | 191 → 205 |
| Unaffected | 796 → 375 | 173 → 178 | 158 → 161 |
| Affected + unaffected | 3,412 → 2,901 | 414 → 427 | 349 → 366 |

Affected error fell only **3.44%**, versus 52.89% for unaffected. Tranche 03
supplied fewer A/U values (273 → 264) and fewer exact ones (240 → 237) despite
lower error; pooled supply rose. Positive affected omissions fell 961 → 948,
but erroneous supplied affected zeros rose 0 → 1. Positive unaffected omissions
fell 119 → 113, while erroneous supplied unaffected zeros rose 2 → 4. These
are separate error types, not evidence that all missing-count problems improved.

The predeclared previously-unscored subset contains 196 attempts / 182 PMIDs;
the full tranches include prior exposure. Its pooled recall rises 78.61% →
78.91% (+0.30 points; interval −0.37 to +1.48), precision falls 84.58% → 84.45%,
and A/U error falls 1,713 → 1,228 (28.31%). This does not establish transfer:
the same dominant 408 row is previously unscored. Subsets and pooling are
exploratory and cannot replace the registered gate.

## Why the large count improvement needs caution

SCN5A PMID 20129283 H558R changes from missing carrier/unaffected values to
408/408. The table cell is real and matches gold, but Table 2 is framed around
2,600 reference alleles and the methods describe 1,300 healthy volunteers from
previously reported controls. The unit and study attribution need adjudication
before calling 408 a validated count of people. A structural table locator
alone does not resolve those semantics.

This row contributes **408 of the 511 pooled A/U error units recovered (79.8%)**.
A post-hoc sensitivity removing its count fields, while leaving all official
scores and decisions intact, gives:

| Sensitivity | A/U absolute error | Relative reduction |
| --- | ---: | ---: |
| Tranche 03, excluding H558R | 1,473 → 1,414 | 4.01% |
| Pooled, excluding H558R | 3,004 → 2,901 | 3.43% |
| Previously-unscored pooled subset, excluding H558R | 1,305 → 1,228 | 5.90% |

Carrier gains are similarly concentrated: 536 of 542 recovered error units
come from this paper and PMID 27566755. Excluding both entire papers leaves
carrier error 1,611 → 1,605, a 0.37% reduction. These influence checks are not
new acceptance tests. [Source audit](source_review_03.md),
[reproducible sensitivity](outlier_sensitivity.py),
[complete sensitivity output](outlier_sensitivity.json).

## Revised forecast

The previous 72–76% recall and 15–25% A/U error-reduction forecast concerned
the **original hard continuation-01 cohort**, whose candidate scored 67.97%
recall. It did not concern these new tranches or the whole corpus. That old
candidate also differs from the historical baseline used here, so the new
paired deltas cannot simply be added to 67.97%.

| Scope | Revised planning estimate | What supports it |
| --- | --- | --- |
| Broader reading on available sources | Budget for **no established recall lift** and only **0–5% A/U error reduction**; regression remains possible. | Paired recall is flat and residual count improvement is a few percent. The 0–5% allowance is a judgment, not a confidence interval or proven bound. No reliable FP reduction is forecast. |
| Original hard cohort, retaining the two selected source recoveries | Approximately **76% recall**, with a **74–78% planning range**; roughly 19–31% fewer identity misses than its old candidate. | Replacing only PMIDs 20031634 and 25163546 yields 292 / 143 / 92 TP/FP/FN, recall 76.04%, precision 67.13%. This is bookkeeping over selected reruns, not a fresh 120-paper result. The central estimate assumes remaining identity work is net neutral. |
| A/U error on that original hard cohort | Approximately **12% lower**, with a **10–17% planning range**, reduced from the earlier 15–25%. | The older 20031634 check already recovered 70/576 error units (12.15%); the new roster adds zero count gain. PMID 21302287 contains at most another 29 units (5.03%) of opportunity, conditional on correct source-resolvable joins. None of that unimplemented upside is in the central estimate. |

These ranges are engineering scenarios, not statistically calibrated intervals.
A corpus-wide yield from the remaining downloader and clinical-table changes
cannot yet be estimated reliably: the source successes were deliberately
selected problem papers. Grok and Agy both challenged broad extrapolation;
the final estimate separates selected recovery from transferable reading
improvement. [Machine-readable forecast](updated_forecast.json) and
[review dispositions](review_adjudication.md).

## Acquisition and reading changes worth doing next

1. **Validate the actual article component and let source quality govern cache updates.**
   The [publisher article](https://academic.oup.com/eurheartj/article/36/18/1123/2293182)
   links PMID 25163546's real Supplementary Data ZIP. Two unrelated cached files
   had obscured the gap. Its 53-page roster supplement supplies all 20 SCN5A
   identities in one fresh run: **20 TP / 0 FP / 0 FN**, with every count NULL.
   It has no per-variant count column. The verified archive and members are now
   in the external corpus; normal reuse preserves the 20 source cDNA strings
   and folds the two components exactly once. This separate opened-paper check
   is outside both tranches.
2. **Continue reading clinical tables after identities have been extracted.**
   Route unresolved count fields to already-acquired clinical components, with
   complete headers, footnotes, units, endpoints and cohort ownership. PMID
   30059973 Tables 11/14 are a concrete development case. A successful identity
   table shortcut must not imply that phenotype extraction is finished.
3. **Repair table grids and join explicit people before counting.**
   PMID 20433692's merged DOC cells lose family/mutation alignment during text
   conversion. Preserve row spans and deduplicate family + person + genotype;
   retain uncertain phenotype and affected non-carriers. Test explicit-ID joins
   on 21302287 and 30403697 without mixing their evaluation denominators.
4. **Abstain separately for carrier, affected and unaffected fields.**
   PMID 18929323 quotes 13 and 6 carriers but clears structured counts when the
   A/U split is unavailable. Preserve a supported carrier number while leaving
   A/U unknown. Do not derive unaffected by subtraction without a complete
   partition, or treat alleles, families or aggregate occurrences as people.

Endpoint disagreements such as RYR2 25814417's already implemented 97/62/26
calculation require adjudication, not an attempt to force agreement with a
different gold definition. See the [ranked recommendations](recommendations.md)
and [first-tranche source audit](source_review_02.md).

## Integrity, cost and artifacts

Candidate runtime stayed fixed at `e9fef8e…` over 251 files; the historical
baseline is nine files restored from `506a949c` with common d299 infrastructure
and a shared source-freezing overlay. This is not a whole-checkout reproduction.
All initial frozen snapshots match. Actual primary text hashes match on
**239/240** attempts; RYR2 31970460 uses abstract versus staged full text because
of protocol validation behavior. Its scored output is identical. One baseline
BRCA1 operational failure was preserved and retried once before gold access;
treating it as empty output leaves the final comparison numerically unchanged.
Gold-derived aliases were disabled, every successful gene job completed before
lock, and all 34 original audit input hashes remain intact. See
[METHODS](METHODS.md), [source receipt](source_review_03_receipt.json),
[retry receipt](operational_retry_03_baseline.json) and [verification](verification.json).

The four arms cost **$39.89468**, the failed attempt **$0.09305**, and the
supplement check **$0.05516**: **$40.04289 in new API proxy cost**. Including the
prior targeted work, the active campaign is **$44.73454 used / $55.26546 left**
of $100. These are dated list-price estimates, not invoices. All new extraction
used Azure; no Anthropic was used for these tests or the new reviews. CLI
notional prices are recorded separately in the review dispositions.

[Full results](results.json) include registered decisions, all/previously-unscored
summaries, count-error types and PMID-bootstrap intervals. Per-paper changes:
[02](paper_count_changes_02.json), [03](paper_count_changes_03.json).
Reproduce the descriptive calculations with:

```bash
.venv/bin/python docs/evidence/tranche_validation_20260905/summarize.py 02 03
.venv/bin/python docs/evidence/tranche_validation_20260905/outlier_sensitivity.py
```

Per-run paired figures:
[02](../../../benchmarks/codex_paper_eval/runs/20260905_protocol_cont120_02_candidate/figures/gold_difference.png),
[03](../../../benchmarks/codex_paper_eval/runs/20260905_protocol_cont120_03_candidate/figures/gold_difference.png).
The [canonical stratified figure](../../figures/evaluated_phenotype_counts/phenotype_count_recovery_stratified.png)
and its [run-membership JSON](../../figures/evaluated_phenotype_counts/phenotype_count_recovery_stratified.json)
retain legacy and opened-candidate views separately. Their 458 opened attempts
span earlier revisions and are not the pooled same-runtime estimate above.
Twelve PNGs were visually checked; [QA receipt](figure_qa.json).
