# New-model routing test — 6 September 2026

**Keep the current primary reader. Test Astra on bounded, source-linked patient
rosters and difficult evidence joins. The tested full-paper Astra configuration
is not ready for promotion, and the literal-only Astra overlay adds almost no
accepted counts.** The raw reader proposals show useful potential that the
current evidence contract cannot safely accept automatically.

This is a completed, failure-enriched **12-paper opened calibration** with the
same frozen sources across arms. It is not a random corpus sample, a blind
holdout, a model-weights-only comparison, or a comparison with a fresh human
reread. No default, headline metric or unopened tranche changes. All three final
arms locked before any new scores; abandoned pilot arms remain unscored.
[Methods](METHODS.md), [design and amendments](PLAN.md), [lock receipt](all_arms_locked_before_score.json).

## Measured results

All identity figures count variant-paper rows. “Exact/supplied” counts only
emitted values on matched reference rows; an abstention is not a correct supplied
zero. Error sums include every asserted reference row, evaluating missing
identities and abstentions as zero without changing stored nulls.

| Measure | Grok 4.3 control | Astra medium primary | Grok 4.3 + Astra clinical reader |
| --- | ---: | ---: | ---: |
| Identity TP / FP / FN | 699 / 36 / 89 | 608 / 21 / 180 | 699 / 36 / 89 |
| Identity recall | 88.71% | 77.16% | 88.71% |
| Identity precision | 95.10% | 96.66% | 95.10% |
| Carrier values, exact / supplied | 504 / 553 | 511 / 559 | 505 / 554 |
| Affected values, exact / supplied | 22 / 24 | 37 / 37 | 22 / 24 |
| Unaffected values, exact / supplied | 67 / 69 | 65 / 65 | 67 / 69 |
| Carrier absolute error | 337 | 512 | 335 |
| Combined affected/unaffected absolute error | 1,132 | 1,200 | 1,132 |
| Count-bearing identity extras | 13 | 15 | 13 |
| Returned API price proxy | $1.81 | $16.58 | $6.41 incremental |
| Unknown usage/retry reserve | $0 | $9.82 | $8.27 incremental |

The overlay reuses the control's identities. Its incremental $6.41 returned
proxy plus $8.27 unknown reserve yields **one accepted carrier fact and no A/U
facts**; it is not an independent identity observation. Its full-stack known
proxy is $8.22 including the $1.81 baseline.

The reference asserts 788 carrier fields, 788 affected fields and 787 unaffected
fields. Positive reference fields number 705, 650 and 86 respectively. Coverage
of those positive fields is **553/705, 24/650, 69/86** for control;
**559/705, 37/650, 65/86** for Astra; and **554/705, 24/650, 69/86** for the
accepted overlay. Thus high precision among supplied values coexists with very
large affected-count omissions. Count fields on identity extras, kept outside
that matched-value denominator, are carriers/affected/unaffected **13/3/1**,
**15/1/0**, and **13/3/1**. They are reference-membership extras, not automatically
proven source hallucinations. [Complete results](results.json).

Astra adds 15 exact affected values but loses two exact unaffected values and a net
91 matched identities. Carrier error rises 51.9%; combined A/U error rises 6.0%.
Its 100% exactness on supplied A/U values does not establish complete recovery.
The empty SCN5A 32533946 result also removes 20 control identity extras, helping
the apparent precision increase; empty output is not evidence of safer reading.

Two large SCN5A papers, 20129283 and 30059973, use deterministic extraction
shortcuts: together they account for 602/788 reference identities and 521 control
TPs. Their shared result does not demonstrate model equivalence. The explicitly
post-hoc slice excluding both is **178/26/8 → 87/11/99 TP/FP/FN**, recall
**95.70% → 46.77%**, for control → Astra. This is an additional descriptive
slice, not a replacement gate.

The two preflagged source/reference sensitivity exclusions, 20129283 and
25814417, answer a different question. Excluding both gives **359/29/11 →
269/14/101**, while combined A/U error improves **628 → 585**. The sign of the
count-error result therefore depends on disputed/influential rows. Preserve the
full result and the sensitivity together; neither provides a population forecast.

## What happened on the papers

| Gene / PMID | Control TP/FP/FN | Astra TP/FP/FN | Observation |
| --- | --- | --- | --- |
| RYR2 18929323 | 2/0/0 | 2/0/0 | Both capture explicit 13/6 carrier totals; variant-specific phenotype table is absent. |
| RYR2 19398417 | 1/0/0 | 1/0/0 | Frozen source is an abstract JSON rendering, despite a misleading runtime filename label. |
| RYR2 25435091 | 1/0/0 | 1/1/0 | Astra adds one identity extra. |
| RYR2 25814417 | 1/0/0 | 0/0/1 | Repeated 32k reasoning-only outputs; influential disputed phenotype counts. |
| RYR2 30403697 | 21/0/0 | 13/0/8 | Empty output, then a transport failure with SDK retries; later recovery retains some identities. |
| SCN5A 20031634 | 11/1/2 | 10/3/3 | Long Astra response succeeds, but final identity matching is worse. |
| SCN5A 20129283 | 339/7/78 | 339/7/78 | Shared deterministic path; study/unit issue remains. |
| SCN5A 25163546 | 20/0/0 | 20/0/0 | Both capture the recovered supplement identities, retaining null counts. |
| SCN5A 30059973 | 182/3/3 | 182/3/3 | Shared deterministic path. |
| SCN5A 32533946 | 83/20/0 | 0/0/83 | Astra exhausts 32k with no visible output; subsequent call denied by budget guard. |
| MYBPC3 20433692 | 13/0/3 | 13/0/3 | Full-paper Astra fails; downstream recovery still retains 13 identities. Compact roster reading succeeds separately. |
| MYBPC3 21302287 | 25/5/3 | 27/7/1 | Some identity/count gains, alongside two additional identity extras. |

The Astra primary stage made 15 calls including retries/JSON repair: six returned
no visible output at the length limit, eight ended with `stop`, and one failed
without usage. Grok's corresponding ten calls had nine `stop`, one `length`,
and no empty visible responses. These counts concern the extraction stage;
shared verification/routing calls are separate. The longest failed Astra request
included two SDK retries and took roughly 43 minutes. A returned API response,
valid extraction, trusted identity and accepted count are separate outcomes.
[Call outcomes](operational_summary.json), [failure mechanisms](failure_mechanisms.md).

The independent clinical reader triggered on nine papers, dispatched eight,
completed six without error, timed out on one and exhausted its output budget
on one. The final paper was refused **before dispatch** by the budget guard.
All 12 stay in the comparison, with baseline predictions retained on failure or
non-trigger. Its sole accepted fill is **RYR2 G1886S, PMID 30403697: 2 carriers**,
explicit in a current-study sentence and matching the reference. There are
**zero additional accepted affected/unaffected fields**.

The reader's **unvalidated raw additive lane** supplies 75 new fields: 71 agree
numerically with the reference and four do not. This includes **44 additional
exact A/U fields and three additional A/U disagreements**, plus 27 additional
exact carrier fields and one carrier disagreement. Combined A/U error would
fall **1,132 → 1,065** (5.9%) in that diagnostic lane. These are not accepted
improvements: most require person-level derivation or additional source binding.
The 54 “no unique retained identity” rejections count variant proposals; the
other rejection categories count fields and must not be added into a single
field-error denominator. The overlay cannot add missed identities or replace
wrong nonnull baseline counts. [Raw-lane analysis](raw_lane_results.json).

## Best use of the new models

1. **Use Astra for bounded evidence reading.** In a separate source-only probe,
   Astra low transcribed all 48 people from one compact MYBPC3 patient table in
   49.8 seconds for a $0.10346 API proxy. Original-DOC layout and converted-text
   audits agree in six fields on all people, with explicit qualifications for
   two blank family cells whose membership comes from prose. This changes task,
   schema, source scope and effort together; it does not isolate low versus
   medium effort or measure accepted clinical-count accuracy.
2. **Build a validated patient-record aggregation path.** Preserve people,
   genotype, variant, family/cohort, endpoint/timepoint, footnotes and source
   coordinates; deduplicate repeated tables/index cases; aggregate in code.
   Derived counts need an explicit evidence contract. The table's No Dx
   footnote says “Unaffected or healthy,” while prose distinguishes healthy
   carriers from suggestive ECG findings. Preserve those endpoint distinctions
   and contradictions. Do not weaken the literal validator to accept model sums.
3. **Recover missing evidence before escalating the reader.** RYR2 18929323's
   frozen body lacks Table 1; publisher/journal fetches returned 403 during a
   bounded follow-up. Another full-body model call does not supply the missing
   variant-specific patient rows. Keep table/figure/supplement completeness
   explicit and send caption, headers, footnotes and nearby scope text together.
4. **Gate Grok 4.6 on deployment health first.** Corrected Azure probes did not
   return usable answers within the tested bounds, so no corrected paper arm was
   dispatched and no accuracy ranking is available. Its functioning CLI used a
   different service path. Earlier API pilots silently lost the requested effort
   setting and are abandoned unscored. The SDK routing fix now preserves effort
   and is tested against actual serialized request bodies.

For planning, assign **no demonstrated overall A/U uplift to the accepted
literal-only overlay** and do not forecast a gain from replacing the primary
with the tested Astra configuration. The 5.9% raw A/U error reduction is an
opened-panel diagnostic opportunity, not an overall forecast. The bounded roster
result justifies the [next controlled test](roster_followup_design.md); it does
not justify a new corpus percentage. This experiment establishes **no 75%-of-human
ceiling**. Missing source, unfinished responses, endpoint ambiguity, identity
retention and the acceptance contract all limit this workflow before such a
ceiling could be measured.

## Budget, implementation and verification

All API pilots, probes and final arms are included: **$38.57768 returned-usage
proxy + $81.45984 unknown-usage/retry reserves + $25 uncertainty margin =
$145.03751 accounted envelope** against the approximately $150 authorization.
No calls remain in flight. CLI costs are excluded and no Anthropic extraction
was used. Unknown charges were not released as zero. Azure invoice reconciliation,
cache-price adjustments and exact billing for unreturned SDK attempts remain
unavailable; this envelope is not an exact invoice or guaranteed billing bound.
[Ledger](budget.json), [pricing basis](pricing_basis.md).

Implemented per-model Azure routes, actual transmission of new-model reasoning
effort, correct Astra completion caps/sampling parameters, and compatible vision
entry points. Current production defaults remain unchanged. Native and production
exports now retain authentic failed-call usage as unknown, with known subsets and
trace references. Those accounting amendments occurred after extraction and
before score: exact fingerprint reconstruction isolates the changes, the original
setups remain intact, scientific prediction arrays are unchanged, and the old
and new non-usage aggregate scores agree exactly. See [methods](METHODS.md),
[amendment receipt](production_usage_amendment.json) and
[score-equivalence audit](score_equivalence_audit.json).

Validation: **2,942 offline unit tests passed**, plus 18 focused projection/native
failure tests, 23 budget non-dispatch/native-failure tests, and eight existing
figure-layout tests. Compilation and Ruff pass. The current runtime was also
built and exercised as an isolated wheel; live vision accuracy was not tested.
A small post-score rendering fix keeps long source/status annotations inside
companion panels; it changes neither predictions nor score arithmetic.
[Verification](verification.json). Claude, Grok and Agy reviews and their
resolution are retained under [reviews](reviews/results_disposition.md).

The canonical [stratified figure](../../figures/evaluated_phenotype_counts/phenotype_count_recovery_stratified.png)
and [cohort manifest](../../figures/evaluated_phenotype_counts/phenotype_count_recovery_stratified.json)
were rebuilt with their existing sampling frames; this selected panel is not
pooled into them. Paired companions show the same reference rows for
[Astra primary](../../../benchmarks/codex_paper_eval/runs/20260906_model12_astra_medium_verified/figures/gold_difference.png)
and the [clinical overlay](../../../benchmarks/codex_paper_eval/runs/20260906_model12_grok43_astra_clinical_verified/figures/gold_difference.png).

Figure “exact” rates include evaluation-imputed zeros; the result table above
counts only explicitly supplied exact values. Do not interchange those denominators.
