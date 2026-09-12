# Grant end-to-end run — evidence record (started 2026-09-09)

Status: **in progress**. This README is filled in as stages complete; the
frozen protocol and launch contract are in [`FREEZE.md`](FREEZE.md).

## What this directory will hold

- `FREEZE.md` — frozen code state, model routing, disease clauses, launch table.
- `run_summaries/<GENE>.json` — copied `RUN_STATUS.json` + workflow statistics
  per gene once each `gvf-run` completes.
- `analysis/<GENE>/REPORT.md` + `figures/` — the per-gene Bayesian penetrance
  report produced by `BayesianPenetranceEstimator/iterations/grant_e2e_20260909/`.
- `analysis/SUMMARY.md` — cross-gene table.
- `cost.json` — dated-rate API cost proxy per gene (`summarize_pilot_cost.py`).

## Results (HNF1A, GCK, LDLR, BRCA2, KCNQ1 control) — updated 2026-09-12 11:30

Headline arm = **the lab's standing protocol** (`<GENE>_protocol`, `scripts/fit_protocol.py`): a
feature-conditioned prior fitted across variants, `logit(prior_v) = α + x_v·β`, updated by each variant's own
carriers as `p_v ~ Beta(prior_v·10 + affected_v, (1 − prior_v)·10 + unaffected_v)`. This is the Variant Browser
`p_mean_w` construction from `iterations/kcnq1_mave_20260813`. The prior enters as **10 pseudo-observations**, so
a variant known only from a handful of affected probands is tempered toward its feature-implied rate while a
variant with hundreds of carriers overwhelms the prior. Carriers include gnomAD carriers counted as unaffected.
Common variants are retained everywhere and reported as the known-outcome check. Earlier arms
(`_pooled`, `_gnomad_nomix`, `_gnomad`, `_sens`, and the literature-only `<GENE>`) are retained as sensitivity
analyses in `analysis/<label>/REPORT.md`; the reasoning that led to them and away from them is in FREEZE
deviations 3–8.

### All five genes under the protocol

| Gene | Papers | With counts | Variants | Literature carriers | gnomAD added | Eval set (≥5) | Common alleles | Common posterior median | Held-out prior AUC | ClinVar AUC | Discordant P/LP low | Discordant non-P/LP high |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| HNF1A | 1,000 | 173 | 233 | 877 | 179,610 | 65 | 3 | 0.0008 | 0.923 | 0.970 | 0 | 1 |
| GCK | 1,050 | 194 | 391 | 1,014 | 1,067 | 59 | 1 | 0.0038 | 0.784 | 0.932 | 1 | 1 |
| LDLR | 2,431 | 433 | 871 | 4,811 | 17,079 | 231 | 5 | 0.0053 | 0.881 | 0.751 | 1 | 1 |
| BRCA2 | 3,759 | 575 | 5,962 | 15,689 | 2,006,331 | 552 | 31 | 0.0020 | 0.885 | 0.896 | 4 | 0 |
| KCNQ1 (control) | 2,383 | — | 691 | 10,656 | 7,225 | 259 | 0.0164 | 0.803 | 0.814 | 4 | 9 |

The known-outcome check passes in every gene: **46 common alleles across the five, every one with its 97.5%
bound below 0.10**, against a literature fraction whose median is 1.00. The held-out feature prior — which never
sees a variant's own counts — ranks high-penetrance variants at AUC 0.78–0.92, beating ClinVar class on LDLR
(0.88 vs 0.75) and matching it on BRCA2 and KCNQ1. Machine-readable copy: `analysis/PROTOCOL_SUMMARY.csv`.

### KCNQ1 control (existing 2,383-paper extraction, curated Variant Browser penetrance as reference)

| Arm | Spearman vs curated | Pearson | median abs. difference | common-variant check (6 alleles, literature fraction 0.51) |
| --- | ---: | ---: | ---: | --- |
| **protocol (headline)** | 0.734 | 0.797 | **0.104** | posterior median 0.016, all 6 with 97.5% bound < 0.10 |
| pooled counts | 0.746 | 0.799 | 0.103 | 0.002, all 6 pass |
| anchored, no mixture | 0.816 | 0.832 | 0.118 | 0.002, all 6 pass |
| anchored, ascertainment mixture | 0.725 | 0.695 | 0.289 | 0.001, all 6 pass |
| literature only, mixture | 0.374 | 0.438 | 0.206 | 0.45 — **fails**, 0 of 6 below 0.10 |

134 variants with ≥5 carriers in both. Two things this table settles. First, the literature-only model fails the
known-outcome check because common alleles are reported almost only in affected probands: the population
denominator, not any ascertainment machinery, is what recovers their true outcome. Second, the comparison is
count-limited rather than model-limited — our extraction supplies 14,182 matched carriers against the
reference's 2,492, and the raw observed fractions themselves agree only at Spearman 0.771, so every anchored arm
is scoring near the ceiling that count difference imposes. The protocol arm's advantage is that it is the same
construction as the reference, so its numbers are directly interpretable against the published ones.

The protocol also shows its intended behaviour on this gene: V417M, ClinVar B/LB with 7 of 7 affected probands
and no gnomAD record, has a feature prior of 0.041 and a posterior of 0.44 [0.22, 0.67], where plain pooling
gives 0.84. A341V, with 1,438 carriers, moves only from a prior of 0.79 to 0.81 — the data overwhelm the prior,
as they should. G269S stays at 0.41 on 52 of 136 affected, a genuine reduced-penetrance finding.

### HNF1A (diabetes) — `analysis/HNF1A_protocol/REPORT.md`

- Collection → extraction: 1,416 candidates, 1,023 passed relevance filtering, 1,000 papers extracted across
  three shards and both passes, 173 with usable carrier counts, 419 variant × paper observations, **233
  variants, 877 literature carriers** (87% affected; 69% of observations are single affected probands),
  179,610 gnomAD carriers counted as unaffected. 547 count rows had no usable phenotype split, 218 of them the
  UK Biobank / Geisinger totals table (PMID 36257325).
- Known-outcome check: I27L (gnomAD AF 0.35, 70/93 affected in the literature) → 0.0008; S487N (AF 0.34) →
  0.0000; A98V (AF 0.029, 44/47) → 0.0063. All three 97.5% bounds below 0.10.
- Distribution: recurrent MODY alleles sit high — G292fs 0.96 [0.89, 1.00], R272H 0.97 — while the prior pulls
  thin frameshift series down, G288fs (7/26 affected) to 0.42 and P291fs (18/32) to 0.63. E508K, ClinVar B/LB
  with 56/68 affected in the literature but 122 gnomAD carriers, lands at 0.29 [0.23, 0.35] — a GWAS-scale
  rather than Mendelian effect, which the literature alone would have called penetrant.
- Ranking: held-out feature prior AUC **0.923**, ClinVar class 0.970 (65 variants with ≥5 carriers).
- Discordance: one variant, I213T, a ClinVar VUS at 6/6 affected, estimated 0.93 [0.76, 1.00].

### GCK (diabetes / GCK-MODY) — `analysis/GCK_protocol/REPORT.md`

- 2,609 candidates, 1,060 passed filtering, 1,050 papers extracted, 194 with usable counts, 616 observations,
  **391 variants, 1,014 literature carriers** (92% affected; 66% single probands), 1,067 gnomAD carriers.
  20,658 count rows lacked a phenotype split, 18,020 of them from the yeast deep-mutational-scanning paper
  (PMID 37101203) that the extractor read as one carrier per assay variant.
- Known-outcome check: the single common allele A11T (AF 0.006) sits at 0.0038.
- Distribution: recurrent GCK-MODY alleles stay high — G261R 0.90 [0.78, 0.97], G318R 0.93, E221K 0.91,
  G72R 0.93, R191W 0.80.
- **Discordant examples:** V455E is ClinVar P/LP with a feature prior of 0.91, but only 3 of 28
  genotype-first carriers (UK Biobank / Geisinger, PMID 36208030) are affected; the posterior is
  0.27 [0.15, 0.41] — the counts overturn both the classification and the prior. V199M, ClinVar conflicting
  at 6/6 affected, goes the other way at 0.90 [0.72, 0.99].
- Ranking: held-out feature prior AUC 0.784, ClinVar class 0.932 (59 variants with ≥5 carriers).

### Cost and runtime

Dated-rate API proxy per completed shard is in `cost_<GENE>_*.json` (GVF PROTOCOL_COST_EVAL rate
assumptions, no cached-input discount; not invoices). Summed over shards and both passes: HNF1A $60 (6,399
recorded calls, 1,000 papers), GCK $79 (7,357 calls, 1,050 papers), LDLR $182 (16,470 calls, 2,416 papers) —
$60–75 per 1,000 papers extracted including the source-recovery replay and the retry pass. Wall clock: discovery +
filtering ~1 h, download ~2.5 h per 340-paper shard, extraction ~4–6 h per shard under the shared Grok
quota, source recovery + replay ~3–5 h per shard, retry pass ~2 h, analysis ~4 min per arm.

### LDLR (familial hypercholesterolaemia) — `analysis/LDLR_protocol/REPORT.md`

- Collection → extraction: 20,062 candidates under the `hypercholesterol*` clause, 4,497 discovered, 3,105
  passed filtering, four shards × ~776 PMIDs with the top 625 extracted each (both passes): **2,416 papers
  extracted**, 2,431 in the union database, 433 with usable carrier counts, 1,554 observations, **871
  variants, 4,811 literature carriers** (95.5% affected; 66% single probands), 17,079 gnomAD carriers.
  3,777 count rows had no phenotype split and 813 were quarantined by the trust gate.
- Known-outcome check: all five common alleles at or below 0.045 — A391T (AF 0.080, 2/2 affected) → 0.0002,
  T726I (AF 0.0055, 7/8) → 0.0053, V827I (AF 0.0010, 12/12) → 0.045.
- Distribution: the well-observed pathogenic alleles keep their data — V408M (354/354) 0.98 [0.97, 0.99],
  G528D (182/182) 0.97, D200G (148/149, 110/111 of them in cascade screening) 0.92 [0.87, 0.95], W66G 0.89.
- Discordance: two variants. P476R, a ClinVar VUS at 16/16 affected, estimated 0.97 [0.87, 1.00]; and G335S,
  ClinVar P/LP with a feature prior of 0.13, 5/6 affected and 7 gnomAD carriers, estimated 0.27 [0.12, 0.47].
- Ranking: held-out feature prior AUC **0.881**, beating ClinVar class at 0.751 (231 variants with ≥5
  carriers) — the largest prior-over-classification margin of the four genes.
- **What the FH literature cannot do.** 95.5% of carriers are affected even in cascade-screening rows, so the
  counts carry almost no contrast. The earlier arms made this visible in the worst way: with a between-study
  random effect the all-affected rows were explained away and a handful of gnomAD carriers set the level
  (D200G at 0.14, V408M at 0.48). Under the protocol the estimates track the data, but the genotype-first
  hold-out still shows the count-informed posterior failing to beat a features-only prior on screening-
  ascertained carriers. Resolving FH penetrance needs an anchor with measured LDL-C, which is the biobank
  integration this proposal funds.

### BRCA2 (breast cancer) — `analysis/BRCA2_protocol/REPORT.md`

- Collection → extraction: 13,333 PubMed/PubMind candidates with the exact `"breast cancer"[Title/Abstract]`
  clause, discovery capped at 7,000, 5,875 passed relevance filtering, five shards × 1,175 PMIDs with the top
  400 extracted each (both passes, plus the replay extractions recovered from shards 2 and 3): **3,759 papers**
  in the union database, 575 with usable carrier counts, 8,182 variant × paper observations, **6,106 variants,
  15,689 literature carriers** (76% affected; 87% of observations are single probands), 2,006,331 gnomAD
  carriers counted as unaffected. 14,518 count rows had no usable phenotype split.
- Known-outcome check passes: all **31** common alleles sit at a posterior median of 0.002 with every 97.5%
  bound below 0.10 — N372H (gnomAD AF 0.28, 122/140 affected in the literature) at 0.000, K3326X (AF 0.008,
  140/182) at 0.012, I3412V (AF 0.036, 91/112) at 0.016. These are the cases a literature-only estimate calls
  penetrant and a population denominator correctly does not.
- Ranking: held-out feature prior AUC **0.885 [0.847, 0.916]** for in-silico features, 0.917 with ClinVar added,
  against ClinVar class alone at 0.896 and AlphaMissense alone at 0.708 (552 variants with ≥5 carriers).
  Held-out per-carrier negative log score 0.033 for the in-silico prior versus 0.040 for the intercept.
- Discordant examples, all four P/LP variants whose credible interval lies below 50%: Y1894X (22/22 affected in
  the literature but 45 gnomAD carriers) at 0.32 [0.22, 0.43]; W2586X (5/5, 16 gnomAD carriers) at 0.25;
  L2865X (2/12) at 0.21; R2336L (3/9) at 0.15. No non-P/LP variant is estimated above 50%.
- The prior's influence is visible where the literature is thin: S1982fs, a frameshift with no ClinVar record,
  has a truncating prior of 0.26 and 72 of 102 carriers affected, giving 0.67; K751Q, a missense with no record
  and 252 of 459 affected, stays near its data at 0.55.

## Stage log

| When (local) | Event |
| --- | --- |
| 2026-09-09 08:36 | HNF1A and GCK `gvf-run` launched (discovery → extraction, source recovery on) |
| 2026-09-09 08:53 | LDLR launched (`--extraction-top-n 2500`) |
| 2026-09-09 08:54 | BRCA2 launched (`--extraction-top-n 2000`; discovery hit the 7,000 PMID cap) |
| 2026-09-09 09:05 | Analysis stage validated end to end on the 12-paper HNF1A pilot database (5 hierarchical fits + 25 CV fits, 6 figures) |
| 2026-09-09 09:10 | Tier-2 relevance filtering runs at exactly the Azure request cap (50 calls/min per process); full-text download is sequential at ~10 s/paper. Projected completion: HNF1A and GCK this evening; LDLR and BRCA2 Thursday |
| 2026-09-09 09:20 | KCNQ1 control (existing 2,383-paper extraction) exposed strong between-paper heterogeneity in reported affected fractions; the model was upgraded to a variant × paper hierarchy (between-study SD τ) and the held-out target changed to genotype-first observations where the extractor records ascertainment. Synonymous alleles excluded from fits; in-sample posterior no longer used as a classifier (label leakage) |
| 2026-09-09 09:37 | Download stage measured at ~3 papers/min (sequential); HNF1A's 1,023 filtered PMIDs relaunched as 3 parallel `--pmid-file` shards of the same frozen code; GCK (3), LDLR (4, top-625 extraction each) and BRCA2 (5, top-400 each) queued to shard automatically when their Tier-2 filtering completes. See FREEZE.md "Sharded harvesting" |
| 2026-09-09 10:05 | KCNQ1 positive control complete (existing 2,383-paper extraction, 691 variants, 10,376 literature carriers). Literature-only arm: Spearman 0.43 with the curated Variant Browser penetrance; gnomAD-anchored arm: **Spearman 0.72** (134 variants with ≥5 carriers in both), benign common alleles G643S/P448R at ≈0.001 vs curated 0.003/0.006, G589D 0.13 [0.05, 0.24] vs 0.26, A341V 0.72 vs 0.83, V254M 0.73 vs 0.92; raw literature fractions for the same alleles are 0.75–0.91. Held-out feature prior (never saw the variant) identifies observed high-penetrance variants with AUC 0.80 [0.74, 0.86] in the anchored arm (ClinVar class 0.81). Caveat: in the anchored arm the models without AlphaMissense mix poorly (R̂ 1.3–1.5; mixture multimodality where a variant has only all-affected series and no anchor); warm-up lengthened to 1,500 iterations and non-converged models are flagged in every report |
| 2026-09-09 11:55 | All four genes sharded after Tier-2 filtering: HNF1A 3 shards × ~341 PMIDs (1023 filtered); GCK 3 shards × ~353 PMIDs (1060 filtered); LDLR 4 shards × ~776 PMIDs (3105 filtered); BRCA2 5 shards × ~1175 PMIDs (5875 filtered). Parent runs stopped after filtering; shard-aware analysis waiters armed (primary, sensitivity, gnomAD arms per gene) |
| 2026-09-09 15:45 | Azure Grok deployment quota (429, East US) fails ~8 extractions per shard per hour once ten shards extract concurrently (≈20% of papers so far), and swap reached 6.5 GB. Decision: pause the LDLR and BRCA2 shards (downloads kept on disk), let HNF1A/GCK finish, then give every gene a second `--resume-dir` pass at 2 workers per process, which re-extracts only the papers without an extraction JSON (the rate-limited ones) before the analysis runs (`orchestrate_gene.sh`); LDLR and then BRCA2 follow sequentially (`sequence_ldlr_brca2.sh`). Per-paper protocol unchanged |
| 2026-09-09 18:55 | Source-recovery stage measured at ~20 s per abstract-only paper (~45 min per shard). Second pass will skip it (`--no-source-recovery`) and the analysis unions each shard's pass-1 refreshed database with its pass-2 rebuild. LDLR and BRCA2 will run concurrently at 2 workers per process after HNF1A/GCK, so their non-LLM stages overlap |
| 2026-09-10 02:20–02:55 | HNF1A and GCK first-round analyses complete (3 shards each, both passes). HNF1A: 1,000 papers → 419 usable observations, 234 variants, 877 carriers (87% affected; 69% single-proband rows). GCK: 1,050 papers → 616 observations, 398 variants, 1,014 carriers (92% affected). Most excluded count rows are carrier totals without a phenotype split (GCK: 18,020 rows from the deep-mutational-scanning paper 37101203) |
| 2026-09-10 17:50 | Anchored-arm failure diagnosed: common polymorphisms brought ~180,000 gnomAD carriers into HNF1A and the global mixture weight collapsed (π ≈ 0.98), pushing known MODY variants to ≈0.04. Added the gnomAD AF > 0.1% exclusion (all arms) and a no-mixture anchored arm; KCNQ1 control re-run to choose the anchored model; HNF1A/GCK arms re-running. LDLR (4 shards, ~560 extractions each) and BRCA2 (5 shards, ~330 each) are finishing pass-1 extraction |
| 2026-09-10 18:05 | KCNQ1 control re-run with the frequency exclusion (8 common variants, 280 literature carriers removed): anchored model **without** the ascertainment mixture agrees best with the curated Variant Browser penetrance (Spearman 0.79, Pearson 0.83, median absolute difference 0.12 over 129 variants; e.g. Y111C 0.82 vs 0.91, G269S 0.62 vs 0.61, A370V 0.03 vs 0.08, G589D 0.43 vs 0.26); the anchored mixture model gives 0.68 (median difference 0.28) and under-estimates high-penetrance founder alleles. Decision: the anchored no-mixture arm (`<GENE>_gnomad_nomix`) is the headline per-variant estimate; the mixture arms are reported as ascertainment-aware sensitivity analyses |
| 2026-09-10 18:20 | Brett: do not exclude common variants — they are the known-outcome check (GWAS effect sizes or no effect). Frequency exclusion turned off in all arms; a common-variant check (literature fraction vs posterior, gnomAD anchors) added to every report; the mixture's ascertainment probability now also depends on the variant's allele frequency (not identified on proband-dominated HNF1A, so the anchored no-mixture arm stays the headline). KCNQ1 control and HNF1A/GCK arms re-running with common variants retained |
| 2026-09-11 04:02–06:43 | LDLR shards 3, 0 and 1 finished pass 1 (source-recovery replay re-extracted 247–261 papers per shard from recovered full text, ~7.3 h per shard; each ends with the expected non-fatal `fetch_paywalled.py exited 1`). Shard 2 (255/261 replays staged) and all five BRCA2 shards (46–56% of ~545 replays each) still in the replay stage; BRCA2 pass 1 projected 13:00–19:00 |
| 2026-09-11 06:45 | Zero-cost dry run of the analysis stage on the finished LDLR shards (3 of 4) and on the five pre-refresh BRCA2 pass-1 databases. LDLR: 1,825 papers → 1,128 usable observations, 701 variants, 3,748 carriers (97% affected; FH-ascertained case series), 387 with a ClinVar class (294 P/LP), 12 common (AF > 0.1%). BRCA2: 1,992 papers → 7,938 usable observations, 6,019 variants, 40,056 carriers (46% affected), 4,013 with a ClinVar class (488 P/LP, 872 B/LB, 2,047 VUS), 42 common. BRCA2 has 12,622 total-only rows, 6,959 of them from one population-frequency table (PMID 39779857), all excluded by the denominator rule |
| 2026-09-11 06:50 | Case-control check for the common-variant test: only 13 of the 166 BRCA2 common-variant rows (and none of the 19 LDLR rows) carry a `case_control` design, so the frozen contract is kept (no design exclusion). Common BRCA2 alleles are reported by affected-ascertained series (e.g. N372H 142/4,980, M784V 1,281/5,792, V2466A 146/148); the gnomAD anchor supplies 2,098–1,458,237 reference carriers per allele |
| 2026-09-11 06:56 | Scheduling change only: the LDLR/BRCA2 arms now run primary → gnomad_nomix (headline) → sens → gnomad, so the headline arm lands first (`auto_analyze_multi.sh`, `nomix_ldlr_brca2.sh`; arms and flags unchanged). Preliminary headline fits on the dry-run data started to measure BRCA2 fit time (6,019 variants vs 233–391 for HNF1A/GCK, whose arms took 3–4 min) |
| 2026-09-11 07:39 | LDLR shard 2 finished pass 1 (1,754 min); orchestrator started the LDLR retry pass at 07:40 (23–29 recovered full texts promoted per shard) |
| 2026-09-11 07:45 | Preliminary LDLR headline arm (3 of 4 shards, 692 variants, 286 s) exposes a failure mode of the population anchor for FH: the literature is 97% affected even in cascade-screening rows, the plain study effect grows to τ ≈ 5–6 logits (also 4.4 on HNF1A, 5.9 on BRCA2; 1.0 on GCK), and the gnomAD anchor (carriers assumed unaffected) then decides the variant-level estimate — e.g. D200G, 110/111 affected in cascade screening, is posted at 0.14. The genotype-first hold-out catches it (posterior anti-predictive, AUC 0.10 on 31 variants; ClinVar-class rate better), so LDLR's "discordant P/LP" list under the frozen arm must be read against the genotype-first evidence (H583Y 29/57 across 13 papers is a genuine reduced-penetrance example; D200G is an anchor artefact). Interpretation: hypercholesterolaemia is common and unmeasured in gnomAD, so "gnomAD carrier = unaffected" is invalid for FH; a biobank anchor with measured LDL-C is the remedy. No protocol change; the frozen arms run as specified |
| 2026-09-11 07:50 | BRCA2 preliminary headline fit (6,019 variants, 7,938 rows): ≈4–5 min per main NUTS fit, ≈65 s per CV fold fit, whole arm ≈1 h; every fit converged (R̂ ≤ 1.01, 0 divergences). Swap filled (17 GB machine, 9 gvf-run + 5 replay processes + fits), so the optional preliminary mixture arms were stopped; production processes unaffected |
| 2026-09-11 08:05 | Report fix (presentation only): anchored arms' evaluation set and labels include the gnomAD carriers assumed unaffected; the report now says so, adds a literature-only-label classification table, and lists discordant variants by literature evidence with gnomAD-added and genotype-first columns (FREEZE deviation 4). Regenerated: HNF1A, GCK, KCNQ1 arms. Under literature evidence the anchored no-mixture lists shrink (HNF1A_gnomad_nomix: none; GCK_gnomad_nomix: V455E and V199M; KCNQ1_gnomad_nomix: 4 + 3; LDLR prelim 31 → 11, all anchor-driven except where genotype-first data are absent) |
| 2026-09-11 08:10 | BRCA2 preliminary headline arm (pass-1 data, 6,019 variants, 1.98 M gnomAD carriers anchored, 55 min): held-out in-silico prior AUC 0.92 [0.89, 0.95] vs ClinVar 0.88 (anchored labels) and 0.87 vs 0.74–0.78 (literature labels); alternate-paper hold-out (122 variants, 14,350 carriers) NLL 0.71 for the count-informed posterior vs 1.56 for the ClinVar-class rate (Δ −0.85 [−1.55, −0.24]); common-variant check 31 alleles, 74% at ≈0 — the failures (I3412V AF 0.036, I2944F, D1420Y) have a gnomAD frequency but no allele count in the warehouse join, so no anchor was built and I3412V (B/LB, 339/1,230 in the literature) surfaces as a spurious 0.96 "discordant" example. Data-join gap, being checked before the BRCA2 arms run |
| 2026-09-11 08:25 | Data-join fix (FREEZE deviation 5): gnomAD max-AF record is now chosen among records that carry an allele count; AF-only `gnomad41_*` annotation records had been winning rounding ties and 476 BRCA2 variants (12 common) had a frequency but no anchor. BRCA2 dry run after the fix: all 42 common alleles anchored (I3412V 5,536 gnomAD carriers). HNF1A/GCK/KCNQ1/LDLR joins verified identical, no refit needed |
| 2026-09-11 09:50 | Deadline measure (FREEZE deviation 6): BRCA2's retry pass will list only the 397–400 PMIDs per shard that pass 1 extracted, sparing ~430 abstract-only re-download attempts per shard (~3 h) that could not change the top-400 extraction set (1–4 of them are inside it). LDLR's retry pass, already running unrestricted, re-attempted 77–~300 abstract-only papers per shard (35 min – 2 h) before extracting the 23–29 promoted full texts |
| 2026-09-11 10:15 | Corrected BRCA2 preliminary headline arm (pass-1 data, join fix, 48 min): the known-outcome check now passes for all 31 common alleles (posterior median 0.0005; every 97.5% bound < 0.10; I3412V 0.000 with 5,536 gnomAD carriers). Held-out in-silico prior AUC 0.93 [0.90, 0.95] vs ClinVar 0.89 (anchored labels), 0.84 vs 0.74–0.78 (literature labels). Alternate-paper hold-out (122 variants, 14,350 carriers): count-informed posterior NLL 0.58 vs ClinVar-class rate 1.56 (Δ −0.98 [−1.70, −0.26]). But τ = 7–10 logits and the six "discordant" P/LP variants by literature evidence (R2494X 80/81, R2318X 46/53, Y1894X 21/21, R2336H 18/20, E1308X 14/14, D2723H 6/6 affected) are all posted at 0.02–0.15 because gnomAD lists 8–45 carriers of each — the same anchor dominance as LDLR, now on BRCA2 truncating alleles whose literature is entirely case-ascertained. The genotype-first hold-out for BRCA2 is polluted by common polymorphisms typed in case series (N289H 138/139 "screening" carriers affected) and is not interpretable |
| 2026-09-11 12:16 | LDLR extraction complete (both passes): 2,416 papers extracted across four shards (606/602/597/611), 69–90 extraction failures logged per shard in pass 1 and retried. Analysis arms launched 12:14 (primary → gnomad_nomix → sens → gnomad) |
| 2026-09-11 13:40 | LDLR anchor provenance (Brett's question): the variantFeatures warehouse enumerates 7,748 LDLR coding SNVs on the MANE transcript; 984 of them have a gnomAD record, all from the v2.1 exome/genome datasets (~141k individuals; no v4 records for LDLR) — 640 missense, 312 synonymous, 32 stop-gained; indels and splice alleles are not in the warehouse, so they get neither anchor nor ClinVar class. Of the 871 modelled variants 221 received an anchor (146 of 362 P/LP); anchors are tiny except for common alleles (median 3 carriers, 75th percentile 7, A391T 12,192). A handful of gnomAD carriers assumed unaffected outweighs hundreds of affected literature carriers only because the between-study effect (τ ≈ 6) lets the literature rows be explained away |
| 2026-09-11 14:30 | **Headline arm changed to pooled counts** (FREEZE deviation 7, Brett's direction: "just regular counts", no anchor machinery). KCNQ1 control, 134 variants with ≥5 carriers vs curated Variant Browser penetrance — pooled: Spearman 0.746, Pearson 0.799, median abs. diff **0.103**; old anchored no-mixture: 0.816 / 0.832 / 0.118; anchored mixture 0.725 / 0.695 / 0.289; literature-only 0.374 / 0.438 / 0.206. Pooled is better calibrated in level, marginally worse in rank, and is the only arm that is not catastrophically wrong on LDLR. Common-variant check passes in all four genes under pooling (HNF1A median 0.0008, GCK 0.0051, KCNQ1 0.0163, LDLR 0.006; every 97.5% bound < 0.10). Held-out in-silico prior AUC is unchanged (HNF1A 0.914 vs 0.915; GCK 0.771 vs 0.776; KCNQ1 0.767 vs 0.770; LDLR 0.897 vs 0.898) because the prior model is the same. Pooled arms for HNF1A, GCK, LDLR and the KCNQ1 control are in `analysis/<GENE>_pooled`; BRCA2's is queued behind its primary arm |
| 2026-09-11 14:30 | Known limitation of pooling, recorded rather than engineered around: a variant reported only in affected probands with no gnomAD record is taken at face value, so the KCNQ1 discordance list grows from 7 to 13 (V417M, ClinVar B/LB, 7/7 affected across 4 papers, no gnomAD carriers → 0.84) while gaining genuine reduced-penetrance findings (G269S, P/LP, 52/136 → 0.40; R555C 20/63 → 0.36). The feature prior shrinks single-carrier variants (LDLR: 360 variants seen in one affected carrier, raw 1.00, posterior median 0.94, IQR 0.89–0.95) but cannot resolve a 7-of-7 proband series. The fix is a denominator (biobank carriers), not more model structure |
| 2026-09-11 16:20 | **Headline arm corrected to the lab's standing protocol** (FREEZE deviation 8): feature-conditioned prior + Beta posterior with prior strength 10, the Variant Browser `p_mean_w` construction from `iterations/kcnq1_mave_20260813`. My hierarchical arms had put the feature prior in as a diffuse mean (between-variant SD ≈ 2 logits) that a handful of observations overwhelms; the protocol puts it in as 10 pseudo-observations, which is what tempers a proband-only series. KCNQ1 V417M (B/LB, 7/7 probands, no gnomAD): prior 0.041 → posterior 0.44 [0.22, 0.67], vs 0.84 under plain pooling. LDLR D200G 0.92 [0.87, 0.95]; V408M 0.98; A391T 0.000. Common-variant check passes in all four genes (posterior medians 0.0008–0.016, every 97.5% bound < 0.10). Discordance lists collapse to what is defensible: LDLR 2 variants (P476R VUS 16/16 → 0.97; G335S P/LP 5/6 with 7 gnomAD carriers → 0.27), HNF1A 1, GCK 2, KCNQ1 13. Held-out prior AUC: LDLR 0.881, HNF1A 0.923, GCK 0.784, KCNQ1 0.803. Arms in `analysis/<GENE>_protocol`; BRCA2 queued |
| 2026-09-11 16:25 | KCNQ1 control under the protocol vs curated Variant Browser penetrance (134 variants): Spearman 0.734, Pearson 0.797, median abs. diff 0.104. The comparison is count-limited, not model-limited — our extraction supplies 14,182 matched carriers against the reference's 2,492, and the raw observed fractions themselves agree at Spearman 0.771, so the protocol posterior sits just under the ceiling the count difference imposes |
| 2026-09-11 18:00 | BRCA2 shard 2's source-recovery replay wedged and was cut loose. Diagnosis: the process sat at 99% CPU for 3h20m with no staged extraction after 14:38 and no LLM call after 14:43; a process sample showed the main thread entirely inside CPython's regex engine (`_sre_SRE_Scanner_search` / `sre_ucs2_match`) — catastrophic backtracking, not slow progress. The pending candidate list holds two enormous cleaned documents, PMID 37234870 at **207.6 MB** and PMID 40669753 at 107.5 MB. Killed the replay (`refresh_run_db.py`, pid 82580); `gvf-run` recorded the stage failure and continued into recovery layers. Cost: shard 2's 425 staged replay re-extractions are lost because the refresh lands them only at the end, so its 397 papers keep their pass-1 extractions; the retry pass still promotes and re-extracts shard 2's 47 recovered-source papers. The other four shards completed their replays normally (247–537 candidates each) |
| 2026-09-11 18:00 | Follow-up for the protocol backlog, not fixed inside the freeze: the variant regex pre-scan has no input-size ceiling and no timeout, so a single pathological document can consume a whole shard's recovery stage. A size cap with a recorded skip, or a scan timeout per document, belongs in `utils/variant_scanner.py`; the 207 MB and 107 MB BRCA2 files are reproducible test cases |
| 2026-09-12 06:25 | BRCA2 shard 3's replay wedged identically, on an **879 MB** cleaned document (PMID 40024972): 17 h of CPU, nothing staged for 14 h, main thread in the regex engine. Killed it; the shard resumed into recovery layers. Only 7 of its 544 candidates were still pending |
| 2026-09-12 06:26 | **Both wedged shards' replay work recovered.** `refresh_run_db.py` applies its acceptance gate per PMID *inside* the replay loop and writes the accepted extraction into `staged_extractions/`; only the database rebuild was still outstanding when the processes were killed. So the staged directories already hold gated output — originals for gated PMIDs, replays for accepted ones — plus papers the replay extracted that were never in the top-400 set. Promoted them into each run's `extractions/` (pass-1 originals kept in `pass1_pre_staged_backup/`): shard 2 397 → **661** extractions, shard 3 400 → **773**. The retry pass now rebuilds each shard's database from them and applies variantFeatures enrichment, the wrong-gene false-positive quarantine and the trust gate through the normal pipeline, rather than any hand-built database. Shards 0, 1 and 4 already carry their replay rows in their pass-1 refresh databases, which the analysis unions with the pass-2 rebuild, so all five shards end up on the same footing. Net BRCA2 corpus: 397+400+661+773+398 extractions across the shards, up from the 1,992 papers of the pass-1 dry run |
| 2026-09-12 10:18 | BRCA2 extraction complete (both passes). Final extractions per shard 397 / 400 / 655 / 743 / 398; the retry pass was restricted to each shard's extracted PMIDs (deviation 6) and logged 53–69 extraction failures per shard. Union corpus for the analysis: **3,759 papers**, 575 with usable carrier counts, 8,182 variant × paper observations, **6,106 variants, 15,689 carriers** (76% affected; 87% of observations are single probands), 14,518 count rows without a usable denominator. Analysis arms launched 10:18; the protocol (headline) arm runs after the primary |
