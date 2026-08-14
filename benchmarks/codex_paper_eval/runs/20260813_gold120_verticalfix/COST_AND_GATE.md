# Patched gold-120 gate and cost record

This is the fresh current-system revalidation of `tier2_gold_120.tsv`: 120
gene-paper attempts / 116 unique PMIDs, 30 attempts per cardiac gene. Exact
run-local source digests and all four write-time-verified production trace
manifests are bound into `predictions.json` and `LOCK.json` before the scorer
reads gold. `diagnostics/trusted_projection/` is a field-masked diagnostic of
the same raw DBs and is not a second blinded primary.

## Acceptance result

**Passed.** Counted-extra precision is 534/(534+24) = **95.70% raw** and
534/(534+23) = **95.87% trusted**, above the 77.3% floor. The separate
count-bearing-only diagnostic is 148/(148+24) = **86.05% raw** and
139/(139+23) = **85.80% trusted**. Variant recall remains 534/635 (**84.09%**).

Carrier MAE is **0.308 raw / 0.299 trusted**. This is about half the 0.614
canonical all-paper baseline but regresses from the preceding stochastic
gold-120 run's 0.266 / 0.243. Raw carrier absolute error is 45/146 supplied
matched counts versus 41/154 previously; trusted is 41/137. The advance is
therefore a measured precision improvement with acceptable low absolute MAE,
not a claim that every conditional count metric improved run-over-run.

The regression target is directly confirmed: PMID 26746457 retains all 42
eTable-6 KCNH2 variant identities and all 42 have null carrier counts because
the table classifies variants by laboratories rather than reporting one patient
per row.

## Measured LLM burn

| Gene | Summed provider time | Calls (successful) | Tokens | Cost proxy |
| --- | ---: | ---: | ---: | ---: |
| KCNH2 | 33.21 min | 108 (106) | 515,098 | $2.365 |
| KCNQ1 | 33.02 min | 166 (162) | 587,438 | $2.255 |
| RYR2 | 30.25 min | 132 (132) | 604,475 | $2.671 |
| SCN5A | 35.24 min | 121 (121) | 644,236 | $2.482 |
| **Total** | **131.73 min** | **527 (521)** | **2,351,247** | **$9.774** |

| Model role | Calls (successful) | Tokens | Provider time | Cost proxy |
| --- | ---: | ---: | ---: | ---: |
| Kimi K2.6 table routing | 19 (19) | 68,470 | 199.2 s | $0.168 |
| Grok 4.3 primary extraction | 115 (114) | 1,617,490 | 4,571.4 s | $2.587 |
| GPT-5.6 Sol verification and figure vision | 393 (388) | 665,287 | 3,133.0 s | $7.019 |

The four jobs ran concurrently and finished in 34.4–38.9 minutes each. The
cost proxy charges traced input/output at $5/M + $30/M for Sol, $1.25/M +
$2.50/M for Grok, and $0.95/M + $4/M for Kimi, ignoring cached-input discounts.
It is a public-list-price proxy, not the deployment-specific Azure invoice.

At the observed $0.0814 and 65.9 summed provider-seconds per attempt, the fixed
BMPR2 50 + BRCA1 50 + BRCA2 46 run is about **$11.89 / 2.67 active hours**.
The exact 546-attempt reviewer tier is about **$44.47 / 9.99 active hours**.
Practical allowances remain $12–$30 / 3–6 hours and $45–$100 / 11–20 hours,
respectively, plus source acquisition and human QC/curation.
