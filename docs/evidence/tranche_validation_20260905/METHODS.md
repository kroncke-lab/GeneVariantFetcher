# Two additional registered paired tranches

This evaluation follows the prospective `PLAN.md` and the existing
`mixed_gold_continuation_120` registry. Each tranche contains 120 gene-paper
attempts and 110 PMIDs. Tranches 02 and 03 have disjoint article membership.
The intended experiment is current implemented reading versus the historical
nine-file reading protocol on identical frozen source assets.

## Runtime and source control

Candidate code is commit `d2996c74f993f9d66558d564c8e515272c679995`, fingerprint
`e9fef8ebe7e79457098e7ad6debf2c9cc6e9f601ae34cde659682d168cbab02b` over 251 files.
Baseline fingerprint is
`e620f9927cffc965293b1fa2feb235c412f67309c3a3492dca33dcdf31bd5d9d`, also 251 files.
It restores nine files from `506a949c`, with the same frozen-source staging branch
in `pipeline/steps.py` and common infrastructure from d2996c74. It is **not** a
full-checkout reproduction of the old commit. `baseline_runtime.json` identifies
the nine files and hashes; `baseline_frozen_source_overlay.patch` reproduces
the shared staging branch with `git apply --unidiff-zero`. The patch was applied
to an isolated copy and verified byte-for-byte against the archived baseline.

All work stayed in the canonical main checkout. Complete runtime archives were
made before switching the nine files. Switches refuse to run during extraction,
and launch/lock checks verify the full runtime fingerprint. The candidate is
unchanged across the two new tranches; no new reader behavior was fitted to
tranche 02's opened gold.

Baseline preparation copies source text and local PDF/figure/supplement assets
from the mounted corpus. Candidate preparation uses the baseline snapshot. The
initial manifest must match exactly; `summarize.py` additionally compares the
actual text input hashes rebound at lock for every attempt. PubMed and normal
production external metadata may still be loaded. No prior predictions or
SQLite databases are reused, gold-derived aliases are disabled, and corpus
sync and review publishing are disabled. These arms measure **reading with
fixed source availability**, not live acquisition success.

## Endpoint definitions

The registered identity endpoint requires observed recall gain ≥1 percentage
point, its one-sided 95% PMID-cluster bootstrap lower bound ≥−1 point, and the
precision lower bound ≥−2 points. Bootstrap: 10,000 resamples, seed 2026090302.
An identical candidate must pass the next unopened tranche to confirm a
successful discovery. Because tranche 02 fails discovery, tranche 03 is further
discovery/calibration, not confirmation of a pass.

The previously registered secondary carrier endpoint requires end-to-end MAE
delta ≤−0.05 and upper 95% bound <0, nondecreasing coverage on identity-matched
rows with its non-inferiority bound, and the identity safeguards. A/U endpoints
remain descriptive. Neither pooled results nor the previously-unscored subset
introduce a replacement acceptance threshold.

End-to-end count error uses every asserted reference value. A missing identity
or abstained count is evaluated as a zero prediction, so its absolute error is
the reference count. This never converts stored NULL into a literal zero.
Conditional MAE uses only supplied counts. We report both, plus supplied and
exact supplied values; a reduction in error can include withdrawal of an
incorrect count. Per-status transitions in `results.json` separate these cases.
The same result also separates omitted positive reference fields, supplied
zero against positive reference, supplied positive against reference zero,
and overcount/undercount units. These are reference-relative error categories,
not independent adjudications of the paper's clinical truth.

Raw identity extras and count-bearing extras retain the registry's accounting,
including known non-exhaustive BRCA gold. A/U-bearing extras are reported
separately from carrier-bearing extras. A count absent from a non-exhaustive
reference is not thereby independently proven scientifically false.

## Prior exposure and descriptive analysis

`prior_exposure.json` was produced from previously locked selection membership
before scoring the new pairs. Tranche 02 includes 18 previously scored PMIDs,
leaving 97 previously-unscored attempts / 92 PMIDs. Tranche 03 includes 20,
leaving 99 attempts / 90 PMIDs. The entire tranches are therefore not described
as wholly untouched holdouts. Registered full-cohort denominators remain intact;
the predeclared previously-unscored subset is a separate descriptive result.

`summarize.py 02 03` verifies prediction/selection locks, pairs identical rows,
reconciles extras with the official report, and checks the 34 original audit
input hashes. It also writes per-paper count changes. The pooled analysis uses
only these two disjoint tranches of the same runtime pair. It does not pool
older candidate revisions, the legacy linkage-assisted cohort, or targeted
calibration papers. Additional two-sided 95% PMID-cluster percentile intervals
use 10,000 resamples and seed 2026090502, and are labelled descriptive.

## Operational failure policy and costs

Before tranche 03 gold access, one baseline BRCA1 process stopped when a
string-valued `extraction_metadata.total_variants_found` could not be summed
with an integer. The failed source/output/trace directory was preserved outside
the scoring input with a full file-hash manifest. We allowed exactly one fresh
retry of that failed gene using the same runtime, models and frozen sources,
with the first operationally complete result selected. No quality score was
used to select it and no paper was dropped. See
`operational_retry_03_baseline.json`. The original driver still exits nonzero;
a resumed supervisor requires all seven successful completion records before
lock and scoring. This is reported as an operational retry, not an uninterrupted
120-attempt arm. Before the tranche 03 lock we also defined a descriptive
sensitivity that treats this one baseline paper as empty output, keeping the
candidate fixed. It does not replace the main comparison or its decision.
The metadata-type robustness issue remains future work.

All extraction launches explicitly select Azure and suppress the Anthropic API
credential. The fixed model configuration is in `model_configuration.json`.
`budget.json` includes completed-arm trace costs and separately retains the
failed operational attempt and the additional opened supplement check. Rates are the repository's dated list-price proxy,
not invoices. CLI-reported notional prices are listed with the reviews rather
than automatically counted as additional API charges. The existing campaign
ceiling remains $100, with $65 reserved for these two pairs.

## Forecast interpretation

The old 72–76% recall forecast concerned the hard, opened continuation-01
candidate (67.97% recall) after further acquisition and count-reading work.
That reference candidate is **not** the historical baseline in these new pairs.
Consequently current-minus-historical deltas cannot simply be added to 67.97%
as an incremental causal estimate.

The previously observed recovered-paper substitution (+11 identities out of 384;
70 fewer combined A/U error units out of 576) was included in the old forecast.
It is not a new full-cohort score and must not be counted twice. The proposed
clinical-table continuation and additional patient-ID joins remain unimplemented
and untested. During the new runs, a separate opened-source check recovered
the missing PMID 25163546 roster and all 20 of its gold identities, with no
extra variants and no supplied counts. Updated future-performance
ranges are engineering scenarios, not bootstrap confidence intervals or a
corpus-wide benchmark result.

Remaining-opportunity denominators must also stay aligned. MYBPC3 21302287
belongs to the old continuation-01 cohort and has 29 affected error units
(23 on matched identities, six on missed identities), at most 5.03% of its
576 combined A/U error units if every one were resolved correctly. This is an
error-budget upper bound, not a prediction that all rows are source-resolvable.
The larger SCN5A 30059973 clinical-table opportunity belongs to the earlier
mixed-02 cohort, and RYR2 30403697 belongs to gold-118. Their gains cannot be
added to continuation-01's numerator. The 22-paper panel's original
`panel.csv` retains all three cohort assignments and per-status errors.

## Separate publisher-supplement check

The ordinary article page for DOI 10.1093/eurheartj/ehu301 advertised its actual
Supplementary Data ZIP. Browser discovery supplied this article-bound link;
the earlier cached supplement filenames had instead represented unrelated
journal material. The archive contains supplementary methods and a 53-page
table PDF. Visual inspection of pages 6, 10 and 43 verified the table header,
20 SCN5A variant rows and a gene-level aggregate diagram. Table 6 has no
per-variant count column; the aggregate diagram cannot allocate individual
variant counts.

After acquiring and hashing the source, a separate $2 API reserve was recorded
before a single fresh source-only run. The same candidate runtime read the
repository body plus the complete supplement text. It used no prior database,
predictions or gold-derived aliases, was locked before scoring, and recovered
20 TP / 0 FP / 0 FN while leaving every carrier/A/U field NULL. Its measured
API proxy was $0.05516. This already-opened paper is outside both new registered
tranches and is never pooled into their effects. Source and result receipts are
`supplement_recheck_25163546.json` and `supplement_check_results.json`.
The final external-corpus layout retains the publisher ZIP, its two members
under the archive's stem directory, and a source-only acquisition receipt.
Standard fold sentinels replace the manual component headings, with the same
source text. A zero-API normal cache-reuse check preserves all 20 cDNA strings
read from the original PDF's SCN5A rows, folds two components exactly once and
is byte-identical on a second fold (`cache_reuse_check.json`). This local layout
check is not a second extraction score. The locked calibration source and
both registered tranche snapshots remain unchanged.

Substituting this run and the earlier 20031634 run into otherwise unchanged
continuation-01 predictions gives 292 TP / 143 FP / 92 FN: recall 76.04%,
precision 67.13%, 25.20% fewer identity misses and 12.15% less combined A/U
absolute error. These are selected-paper substitutions, not a fresh 120-attempt
score. The 20 newly recovered identities add no measured count improvement.

## Artifact preservation

The five new runs’ CSV artifacts explicitly disable Git line-ending conversion so
their byte hashes survive checkout. Tranche-03 prediction JSONs exceed the
normal 1,200-KB commit-hook limit; only those two exact run paths are exempted
from the size hook. All other hooks remain active, and locked predictions are
not shortened or rewritten. These repository artifact rules are outside the
251-file measured runtime; its fingerprint is unchanged.

## Final source parity and influence checks (post-lock)

Both pairs have identical initial source snapshots. Actual primary-text hashes
match on 120/120 attempts in 02 and 119/120 in 03. RYR2 31970460 differs because
the historical validator rejects the staged body and falls back to abstract
JSON, while candidate accepts the staged body. Both score 0/1/1 TP/FP/FN with
no counts; no new acquisition occurred. `source_review_03_receipt.json` binds
the actual inputs and driver logs. This fixed-availability comparison must not
be described as 240 identical actual primary inputs.

After inspecting the largest count change, `outlier_sensitivity.py` excludes
only SCN5A 20129283 H558R's three count fields from descriptive count summaries.
It neither removes the identity from recall nor changes the official gold or
registered tests. Table 2's 408 cell is reference-matching, but its person/allele
unit and prior-control-cohort attribution remain unresolved. The residual
pooled A/U reduction is 3.43% versus 14.98% with the row. A separate whole-paper
carrier sensitivity excludes 20129283 and 27566755. Neither sensitivity was a
prospective acceptance endpoint; both diagnose concentration of influence.
The operational BRCA1 empty-failure sensitivity produces the same final metrics
as the first successful retry, which also has no TP, FP or supplied count.
