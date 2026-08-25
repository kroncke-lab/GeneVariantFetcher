# Gold-118 source-recovery validation and $100 ledger — 2026-08-24

## Disposition

The source-backed implementation is merged, but neither new blind lock replaces
the accepted `20260824_postfix_gold118` headline. The grouped/structural run
improved precision and count error while losing four identity matches. The
source-recovery run restored accepted recall and improved every conditional MAE,
but lost count coverage and missed the count-bearing precision floor. Selecting
the best metric from either stochastic run would be best-replicate selection,
so both locks remain failed validation evidence.

No accepted lock, gold file, public annotation, or collaborator database was
rewritten. The exact implementation commits are `6ce1349`, `c10cf94`, and
`ba1e6ee`; the failed source-recovery lock itself is preserved by `5e351de`.

## Locked comparison

All three rows use the same 118 gene-paper attempts / 114 unique PMIDs and the
same human cardiac gold denominator of 632 variant rows. Predictions and
write-time trace manifests were locked before gold values were read.

| Metric | Accepted postfix | Grouped/structural | Source recovery |
| --- | ---: | ---: | ---: |
| Run | `20260824_postfix_gold118` | `20260824_grouped_structural_gold118` | `20260824_source_recovery_gold118` |
| TP / FP / FN | 546 / 284 / 86 | 542 / 278 / 90 | 546 / 285 / 86 |
| Recall | **86.392%** | 85.759% | **86.392%** |
| Raw precision | **65.783%** | 66.098% | 65.704% |
| F1 | **74.692%** | 74.656% | 74.641% |
| Counted-extra precision | **97.500%** | 98.188% | **97.500%** |
| Count-bearing-only precision | **93.694%** | 95.413% | 93.396% |
| Carrier supplied / MAE / RMSE | **206 / 0.330 / 1.207** | 207 / 0.222 / 0.774 | 197 / **0.198 / 0.730** |
| Affected supplied / MAE | **49 / 0.551** | 56 / **0.339** | 56 / 0.429 |
| Unaffected supplied / MAE | **18 / 0.500** | 31 / **0.194** | 32 / 0.375 |

The preregistered promotion floors were the accepted count-bearing-only
precision (93.694%), accepted carrier coverage (206 rows), non-regressed
recall/precision, and lower carrier MAE. The source-recovery lock failed two:
93.396% count-bearing-only precision and 197 carrier rows. Its lower MAE is real
conditional improvement, but it partly reflects nine additional abstentions and
cannot be promoted alone.

## Source-recovery per-gene result

| Gene | TP / FP / FN | Raw P | Recall | Carrier supplied / MAE | Affected supplied / MAE | Unaffected supplied / MAE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| KCNH2 | 66 / 50 / 22 | 56.90% | 75.00% | 25 / 0.200 | 14 / 1.071 | 3 / 1.667 |
| KCNQ1 | 146 / 100 / 21 | 59.35% | 87.43% | 67 / 0.134 | 5 / 0.800 | 2 / 0.500 |
| RYR2 | 151 / 23 / 27 | 86.78% | 84.83% | 65 / 0.215 | 18 / 0.167 | 12 / 0.250 |
| SCN5A | 183 / 112 / 16 | 62.03% | 91.96% | 40 / 0.275 | 19 / 0.105 | 15 / 0.200 |

The next general levers are therefore acquisition for KCNH2 recall, false-row
control for KCNH2/KCNQ1/SCN5A raw precision, and source-explicit phenotype
counts rather than arithmetic completion. RYR2 already has the strongest raw
precision; broadening it is not the priority.

## What changed

The merged source-backed rules are general and fail closed:

- strip trailing footnote markers without treating them as identity digits;
- retain positive target-gene historical labels such as `2768Cdel` when the
  source establishes their coding role;
- reconstruct a flattened RYR2 genotype row only when its source columns and
  target-gene context are concrete;
- allow a residue-offset suspicion only when a same-row cDNA assertion provides
  a concrete check;
- match reference-less `p.` ranges only under the strict existing range rules.

No per-variant alias, codon-only identity inference, family-to-person
conversion, phenotype inference, or subtraction rule was added.

## Dominant source blocker: PMID 29650123

KCNH2 PMID 29650123 contributes 22 gold variants; the source-recovery lock finds
2 and misses 20, so one source accounts for 23.3% of all 86 remaining FNs. The
publisher page for DOI `10.1016/j.jacc.2018.01.078` exposes one supplement,
`mmc1.docx`, which is already acquired. Direct inspection shows Online Table 1
(cohort characteristics by diagnosis decade), Online Table 2 (beta-blocker
dosage), and three clinical/outcome figures. It contains no per-variant KCNH2
roster. The publisher DOM exposed no second mutation-table asset.

This is a provenance/source blocker, not permission to reconstruct 20 singleton
calls from the gold key or from genotype totals. The theoretical upside is
20 / 632 = 3.16 recall points, but no source-backed implementation can claim it from
the currently published assets.

Authoritative public records:

- <https://pubmed.ncbi.nlm.nih.gov/29650123/>
- <https://www.jacc.org/doi/10.1016/j.jacc.2018.01.078>
- <https://www.jacc.org/doi/suppl/10.1016/j.jacc.2018.01.078/suppl_file/mmc1.docx>

## Evaluated extraction strategy

The calibrated route uses exact frozen PMID files, gold-free production, and a
trusted count-and-identity projection. Deterministic paper census, scanner, and
table parsing run first. Kimi K2.6 routes tables; Grok 4.3 is the primary
per-paper extractor; GPT-5.6 Sol performs source-quoted claim verification and
caption-triaged figure reading. VariantFeatures and trust gates quarantine
unsupported identity/count fields. The route does not infer phenotype
partitions. Count recovery and the per-paper final-check/composer pair remain
default off. For a like-for-like evaluation, source recovery/corpus sync and
publication were disabled inside the sealed run.

The strategy remains source-first: acquire a missing mutation table, then make
one deterministic or source-quoted correction, freeze, and validate once. It
does not restore SCN5A PMID 22685113 carrier=1 values from an abstract that only
implies them, and it does not parse adjacent KCNQ1 PMID 25471708 `n` values when
compound-heterozygous configurations make allele attribution ambiguous.

## Independent max-effort review

Grok 4.6 ran at `xhigh` on the actual locks and implementation. Its strategy
verdict was `BLOCK` promotion and `SHIP` only the narrow source-first follow-up.
It specifically rejected best-replicate selection, the SCN5A singleton
inference, and a naive KCNQ1 adjacent-count parser. It recommended no immediate
third gold-118 run and required an exhaustive human answer key plus a hidden
holdout for the exact 150.

AGY was invoked with Gemini 3.1 Pro at its maximum supported `high` effort for
the same strategy review. It timed out without a substantive answer, so it is
not counted as approval. The earlier AGY `high` implementation review of the
source-backed code completed and returned `SHIP`.

## Cost and the active $100 ceiling

The ledger charges attributable work conservatively and excludes the accepted
baseline, which predates the improvement envelope.

| Item | Calls (success) | Tokens | Public-price proxy |
| --- | ---: | ---: | ---: |
| Targeted exploratory replays / missing telemetry reserve | — | — | $10.0000 |
| Grouped/structural blind lock | 580 (557) | 2,692,803 | $11.86999755 |
| Source-recovery blind lock | 547 (537) | 2,626,745 | $11.19883675 |
| **Charged** |  |  | **$33.06883430** |
| **Remaining** |  |  | **$66.93116570** |

The source-recovery model split was Kimi 26/26 calls ($0.26288925), GPT-5.6
Sol 400/396 ($8.12796000), and Grok 4.3 121/115 ($2.80798750). Deterministic
scoring and packet construction cost $0 in model calls. Grok/AGY consultation
CLI surfaces did not expose billing telemetry; their subscription/account cost
is therefore excluded, not asserted to be zero.

Public rates used for these trace-derived proxies:

- GPT-5.6 Sol: $5/M input, $30/M output —
  <https://azure.microsoft.com/en-us/blog/gpt-5-6-now-available-in-microsoft-foundry/>
- Grok 4.3 short context: $1.25/M input, $2.50/M output —
  <https://docs.x.ai/developers/models/grok-4.3>
- Kimi K2.6: $0.95/M input, $4/M output —
  <https://docs.fireworks.ai/serverless/pricing>

## Exact-150 status

The exact BRCA1 50 + BRCA2 50 + BMPR2 50 staged extraction is source-screened
and complete as a candidate set, but only 3/150 papers overlap any approved
curated fixture and BMPR2 has zero. Its precision, recall, F-score, and MAE are
undefined, not zero.

Commit `85f78a2` preregisters 30 calibration and 20 holdout papers per gene by
stable SHA-256 rank. It binds every source, keeps the holdout packet separate,
requires an exhaustive human `complete` or explicit `NONE` decision per paper,
retains somatic/NONE papers in the raw precision denominator, and excludes
family/unknown counts from person-level MAE. The scorer refuses blank variants,
mixed or count-bearing `NONE`, duplicate variant rows, incomplete status,
missing evidence, and PMID roster drift. No model may construct, approve, or
edit the key.

The next authorized step is human curation of the 90 calibration papers. Only
after calibration changes are frozen should the 60-paper holdout be released
and scored once. No additional paid exact-150 extraction is justified before
that answer key exists.
