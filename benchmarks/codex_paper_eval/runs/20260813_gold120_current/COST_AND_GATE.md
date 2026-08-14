# Gate 2 cost and acceptance record

This is the measured current-protocol run over `tier2_gold_120.tsv`: 120
gene-paper attempts / 116 unique PMIDs, 30 attempts per cardiac gene. The
locked primary is `predictions.json` + `LOCK.json`; `report.json` is the first
gold-reading score. `diagnostics/trusted_projection/` is a post-lock,
field-masked diagnostic and is not a second blinded primary.

## Result

**Gate 2 passed.** The repository's established counted-extra precision is
534/(534+58) = 90.20% raw and 534/(534+54) = 90.82% trusted, above the current
77.3% floor. The initially reported 73.15% is a stricter diagnostic with a
different numerator—count-bearing matched rows only—and was incorrectly
compared to that floor. Variant recall is 534/635 (84.09%); carrier MAE is 0.266
raw and 0.243 trusted.

| Projection | Variant precision | Variant recall | Precision vs counted extras | Count-bearing-only precision | Carrier supplied | Carrier MAE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Raw locked primary | 534/1440 (37.08%) | 534/635 (84.09%) | 534/(534+58) (90.20%) | 158/(158+58) (73.15%) | 154/635 (24.25%) | 0.266 |
| Trusted diagnostic | unchanged | unchanged | 534/(534+54) (90.82%) | 144/(144+54) (72.73%) | 140/635 (22.05%) | 0.243 |
| Vertical identity-only counterfactual | unchanged | unchanged | 534/(534+26) (95.36%) | 148/(148+26) (85.06%) | 144/635 (22.68%) | 0.278 |
| Vertical fix + trusted counterfactual | unchanged | unchanged | 534/(534+22) (96.04%) | 134/(134+22) (85.90%) | 130/635 (20.47%) | 0.254 |

The paper-only diagnostic (`diagnostics/paperonly_score.json`) removes pure
ClinVar/PubTator rows and scores 475 TP / 599 FP / 160 FN: 44.23% precision,
74.80% recall, 55.59% F1. External database layers therefore are not the sole
source of raw variant extras. The counted-extra acceptance metric remains above
its floor.

## Measured burn

| Gene | Active wall time | Calls (successful) | Tokens | Public-price proxy |
| --- | ---: | ---: | ---: | ---: |
| KCNH2 | 29.83 min | 112 (110) | 520,937 | $2.386 |
| KCNQ1 | 26.23 min | 167 (166) | 594,336 | $2.242 |
| RYR2 | 50.37 min | 132 (132) | 602,611 | $2.611 |
| SCN5A | 32.92 min | 121 (121) | 645,708 | $2.590 |
| **Total** | **139.34 min** | **532 (529)** | **2,363,592** | **$9.829** |

| Model role | Calls | Tokens | Summed provider time | Public-price proxy |
| --- | ---: | ---: | ---: | ---: |
| Kimi K2.6 table routing | 19 | 69,073 | 189.6 s | $0.171 |
| Grok 4.3 primary extraction | 114 | 1,615,045 | 4,404.9 s | $2.581 |
| GPT-5.6 Sol verification and figure vision | 399 | 679,474 | 2,965.4 s | $7.077 |

Per-gene model mix (calls / cost proxy):

| Gene | Kimi K2.6 | Grok 4.3 | GPT-5.6 Sol |
| --- | ---: | ---: | ---: |
| KCNH2 | 3 / $0.016 | 27 / $0.580 | 82 / $1.791 |
| KCNQ1 | 2 / $0.018 | 29 / $0.647 | 136 / $1.577 |
| RYR2 | 6 / $0.028 | 28 / $0.639 | 98 / $1.944 |
| SCN5A | 8 / $0.109 | 30 / $0.716 | 83 / $1.766 |

The proxy charges all traced input and output/reasoning tokens at public list
prices and ignores cached-input discounts: Sol $5/M input + $30/M output, Grok
$1.25/M + $2.50/M, Kimi $0.95/M + $4/M. Grok/Kimi are public provider proxies,
not the deployment-specific Azure invoice.

At the observed $0.0819 and 69.7 active seconds per attempt, BMPR2 50 + BRCA1
50 + BRCA2 46 would be about $11.96 and 2.83 active sequential hours. The exact
546-attempt reviewer tier would be about $44.72 and 10.57 active sequential
hours. A practical allowance is $12–$30 / 3–6 hours and $45–$100 / 11–20 hours,
respectively, plus source acquisition and human QC/curation.

## Extraction sequence measured here

1. Fix the PMID manifest and bind eligibility/source provenance without reading
   gold values.
2. Reuse full text, supplements, artifacts, PDFs, and figures; in a full uncached
   run, acquire missing material and fall back to an abstract when necessary.
3. Run deterministic text/table/variant scouting over untruncated material.
4. Use Kimi K2.6 to route usable variant/count tables.
5. Use Grok 4.3 for primary structured variant, count, and evidence extraction.
6. Use GPT-5.6 Sol for risk-ranked count-claim verification and figure vision.
7. Normalize notations, reject wrong-gene/junk rows, merge duplicates, and
   migrate raw observations to SQLite.
8. Add deterministic/lookup recovery layers (including ClinVar/PubTator and
   figure-derived observations), then adopt only grounded count repairs.
9. Apply the field-level trust gate; retain raw values for audit and expose a
   trusted projection separately.
10. Run source QC, finalize exact LLM traces, hash source material, and lock
    predictions before scoring.
11. Score variant precision/recall and count coverage/MAE against manual gold,
    with RYR2 and the trusted projection visible separately.
12. Advance to experimental extraction, QC, and Variant Browser publication
    only after the gold gate passes.

Gate 2 used the already staged corpus with source recovery and corpus sync
disabled, so its wall time does not measure fresh publisher acquisition.

## Excluded scope-correction burn

Before the user's cohort clarification, two launches were stopped and excluded
from scoring/publishing: a widened BMPR2 probe after 14 PMIDs (60 calls, 386,569
tokens, 924.0 summed provider-seconds, $1.357 proxy) and one paper from the old
non-gold-complete cardiac manifest (3 calls, 25,920 tokens, 90.0 seconds, $0.077
proxy). Together they used 412,489 tokens and about $1.43. They are not included
in the Gate 2 totals above.
