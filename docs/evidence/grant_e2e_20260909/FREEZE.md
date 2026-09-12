# Grant end-to-end freeze — 2026-09-09

Purpose: a frozen, reproducible run of the whole workflow (collect papers →
extract → estimate and refine penetrance with a full Bayesian model → report
distributions, classification performance, and discordant examples) for the
grant application. Genes in order: HNF1A, GCK, then LDLR, BRCA2.

## Frozen code

| Repository | State | Notes |
| --- | --- | --- |
| GeneVariantFetcher | tag `grant-freeze-20260909` = commit `76e42d33412787ce7ca8f9d6af077f56dc604a27` (clean tree) | extraction protocol; no code changes during the runs |
| BayesianPenetranceEstimator | commit `bd72a70e6b53010a9b7c3c9b68dceacceec301f3` + 14 uncommitted files from the 2026-09-08 pilot (hashes in `iterations/penetrance_pilot_20260908/reproducibility.json`) | analysis code for this freeze lives in `iterations/grant_e2e_20260909/` |
| variantFeatures | commit `3a2a37cd8475aa7ca071a9a2af0fe0145d9af0d5` + 34 uncommitted files (GCK/MAP4K2 identity repair) | `data/variants.db` read-only; sha recorded per gene at join time |

The uncommitted pilot changes in the two sibling repos should be committed by
Brett to make the freeze a pure commit triple; nothing in this iteration edits
those files.

## Model routing (from `.env`, names only)

Tier 2 relevance: `azure_ai/gpt-5.6-luna` (max reasoning). Table router:
`azure_ai/Kimi-K2.6-1`. Tier 3 extraction: `azure_ai/grok-4.3` (no cascade).
Claim verification / adjudication / arbiter: `azure_ai/gpt-5.6-sol` (xhigh).
Provider: Azure (recurring allocation, per TASKS.md budget policy).

## Launch contract

Launcher: `results/grant_e2e_20260909/launch_gene.sh` (records every launch in
`launches.jsonl`). Concurrency only: `MAX_WORKERS=3` per gene process (the
`.env` default of 1 was tuned for single-gene runs; concurrency does not change
extraction outputs). Flags: `--gold-free-run --no-publish-review
--allow-degraded-institutional` (VU proxy access is denied; publisher APIs and
PMC still run; the run is flagged degraded rather than halting).

Disease clause: the pipeline appends `"<term>"[Title/Abstract]` as an exact
phrase, so terms were chosen by PubMed hit count (2026-09-09):

| Gene | `--disease` | Gene query alone | With clause | Extra args |
| --- | --- | ---: | ---: | --- |
| HNF1A | diabetes | 3,056 | 1,094 | — |
| GCK | diabetes | 5,723 | 2,353 | `--max-pmids 3000` (avoid the sorted-PMID truncation at the 1,500 default) |
| LDLR | hypercholesterol* | 20,062 | 4,345 (2,293 with "familial hypercholesterolemia", which misses UK spelling) | `--max-pmids 5000` (4,497 PMIDs discovered incl. PubMind); sharded 4 ways after Tier-2 with `--extraction-top-n 625` per shard (≈2,500 extractions) |
| BRCA2 | breast cancer | 13,333 | 6,548 | `--max-pmids 7000` (discovery hit the 7,000 cap); sharded 5 ways after Tier-2 with `--extraction-top-n 400` per shard (≈2,000 extractions); kept at five shards to stay within the workstation's memory |

Cost calibration: the BMPR2 full run (974 Tier-2 calls, 756 extractions) was
$30.45 by the dated rate proxy; expect $25–60 per gene here.

## Sharded harvesting (deviation in orchestration only, 2026-09-09 09:37)

The full-text download stage of `gvf-run` is sequential (~3 papers/min), which
would have pushed LDLR and BRCA2 past the deadline. After each parent run
finished discovery and Tier-1/2 filtering, its `pmid_status/filtered_pmids.txt`
was split round-robin into N shards and each shard was run as a calibrated
`--pmid-file` run of the same frozen code (`results/grant_e2e_20260909/shard_launch.sh`;
`--no-corpus-sync` so parallel shards do not race on the corpus index; the
scope gate still applies, Tier 1/2 are not re-run because the parent already
applied them). Per-paper processing, models and prompts are identical; only the
orchestration is parallel. Shard databases are unioned at the observation level
by `build_observations.py --gvf-db <db1> <db2> ...` (database-local ids are
prefixed per shard). Parent run directories keep the discovery, abstracts and
filter decisions. Corpus folding of the new source was not run for shards. Because the shared
Grok deployment quota returned 429s that failed ~20% of extraction attempts
while ten shards extracted concurrently, every gene gets a second
`--resume-dir` pass (2 workers per process) that re-extracts only papers
without an extraction JSON before analysis; LDLR and BRCA2 run after HNF1A/GCK.
Papers that still fail after the retry pass are counted in each shard's
`extraction_failures.csv` and reported as acquisition/extraction loss.

## Deviations from the tagged code (documented, uncommitted)

1. `scripts/refresh_run_db.py` (one line): the override-CSV loader now accepts the
   `source_path` column that `scripts/fetch_linked_supplements.py` writes. Without
   it the whole source-recovery replay aborts (`ValueError: No source path column
   found`) whenever a run found linked supplements, so recovered paywalled full
   text was never re-extracted (observed on HNF1A shard 1 at 20:36, and the same
   failure mode is visible in the 2026-08-07 BMPR2 run). Applied 2026-09-09 20:40
   while the runs were live; shards whose refresh stage ran after that time used
   the fixed loader.
2. Pass-2 promotion (`promote_recovered.sh`, run by `resume_shards.sh` before a
   `--resume-dir --no-source-recovery` pass): for PMIDs listed in
   `source_qc/fetched_source_override.csv` with `action=refresh_replay`, the
   recovered full text replaces the abstract-only `pmc_fulltext` file and the
   pass-1 extraction JSON is cleared, so the normal extraction step re-reads the
   paper from full text. Pass-1 files are kept under `pass1_replaced/`. This is the
   same replacement the refresh stage performs; zero-variant full-text papers that
   the refresh would also re-read (`source_override.csv`) are not re-extracted here.
3. Arm scheduling (2026-09-11 06:56, `auto_analyze_multi.sh` and
   `nomix_ldlr_brca2.sh`): for LDLR and BRCA2 the analysis arms run primary →
   gnomad_nomix → sens → gnomad instead of primary → sens → gnomad → gnomad_nomix,
   so the headline (anchored, no-mixture) arm finishes first; each wait gives up
   after 8 h. Arms, flags and inputs are unchanged. HNF1A and GCK had already
   completed under the earlier order.
4. Report presentation (2026-09-11 08:05, `make_report.py` only; no fit, model or
   metric file changed): in the anchored arms the frozen evaluation set
   (`n ≥ 5`) and the "observed ≥ 50%" label count the gnomAD carriers assumed
   unaffected together with the literature carriers, but the report text called
   them "literature carriers". The report now states the definition, adds a
   second classification table scored against the literature fraction alone
   (≥ 5 literature carriers, common alleles excluded; written to
   `classification_metrics_literature_labels.csv`), and renders the discordance
   table for variants with ≥ 5 *literature* carriers, showing literature
   affected/carriers, the anchored gnomAD carriers and the genotype-first
   affected/carriers separately so anchor-driven entries are recognisable. The
   frozen-rule list stays in `fit/discordant_variants.csv`. All finished arms
   (HNF1A, GCK, KCNQ1) were re-rendered and re-copied; backup
   `make_report.py.bak_20260911`.
5. gnomAD record selection (2026-09-11 08:25, `join_features.py`): the warehouse
   stores AF-only annotation records (`gnomad41_exome`/`gnomad41_genome`, rounded
   AF, no AC/AN) beside the count-bearing gnomAD releases (`gnomad_exome`,
   `gnomad_genome`, `gnomad_r4_*_gene`). The join took the maximum-AF record
   regardless of counts, so whenever the rounded AF-only record won the tie the
   variant kept a frequency but lost its allele count and therefore its anchor
   (476 BRCA2 variants, 12 of them common, e.g. I3412V AF 0.036 with no anchor
   → posterior 0.96 for a B/LB polymorphism). The join now takes the maximum-AF
   record among those carrying an allele count and falls back to AF-only records
   only when none exists. HNF1A, GCK, KCNQ1 and LDLR joins are byte-identical on
   every gnomAD, ClinVar and AlphaMissense column (0 affected variants), so
   their fits stand; BRCA2's arms run with the corrected join. Backup
   `join_features.py.bak_20260911`.
6. BRCA2 retry-pass harvest scope (2026-09-11 09:50, `resume_shards.sh`, marker
   file `shards/BRCA2/PASS2_RESTRICT`): the retry pass lists only the PMIDs that
   pass 1 extracted (`extractions/` plus the promoted JSONs in
   `pass1_replaced/extractions/`, 397–400 per shard) instead of the full
   1,175-PMID shard list. Without it the harvester re-attempts the ~430
   abstract-only papers per shard (measured on the LDLR retry pass at ~30 s
   each, i.e. ~3 h per BRCA2 shard) although only 1–4 of them lie inside the
   top-400 extraction set. The extraction set is unchanged; the only difference
   from the HNF1A/GCK/LDLR procedure is that a retry success outside the top-400
   can no longer be promoted into it. LDLR's retry pass had already started under
   the unrestricted procedure. Backup `resume_shards.sh.bak_20260911`.
7. Headline arm changed to pooled counts (2026-09-11 14:30, Brett's direction;
   `fit_models.py` unchanged — the arm is `--gnomad-as-unaffected --no-mixture
   --tau-prior-sd 0.01`, i.e. the existing model with the between-study effect
   switched off). Each variant's literature affected/unaffected carriers and its
   gnomAD carriers enter one binomial; the only structure left is hierarchical
   shrinkage of the variant intercept toward the feature-informed prior
   (truncating class + AlphaMissense). Rationale: with a 5–7 logit study effect
   the all-affected literature rows were explained away and a handful of gnomAD
   carriers set the level (LDLR D200G: 148/149 affected, 110/111 of them in
   cascade screening, 13 gnomAD carriers → 0.14). Under pooling D200G is 0.91
   [0.87, 0.95] and the common-variant check still passes in every gene, so the
   ascertainment mixture and the study effect were not buying the known-outcome
   behaviour that motivated them — the gnomAD denominator was. Arms `<GENE>_pooled`
   are the new headline; `<GENE>_gnomad_nomix`, `<GENE>_gnomad`, `<GENE>_sens` and
   the literature-only `<GENE>` arm are retained as sensitivity analyses.
8. **Headline arm corrected to the lab's standing protocol** (2026-09-11 16:20,
   Brett; new script `scripts/fit_protocol.py`). Deviations 3–7 had replaced the
   lab's documented construction with a hierarchical model of my own. The
   protocol, implemented in
   `iterations/kcnq1_mave_20260813/scripts/run_analysis.py` and published as the
   Variant Browser's `p_mean_w`, is:

       logit(prior_v) = α + x_v·β      α ~ N(logit(pooled rate), 1.5), β ~ N(0, 1)
       p_v ~ Beta(prior_v·S + affected_v, (1 − prior_v)·S + unaffected_v),  S = 10

   The feature-conditioned prior enters as S = 10 pseudo-observations, so a
   variant reported as a handful of affected probands is tempered toward its
   feature-implied rate while a variant with hundreds of carriers overwhelms the
   prior. Carriers include gnomAD carriers as unaffected, as the public table
   does. This is the mechanism I had wrongly reported as absent: KCNQ1 V417M
   (ClinVar B/LB, 7/7 affected probands, no gnomAD record) has a feature prior of
   0.041 and a protocol posterior of 0.44 [0.22, 0.67], against 0.84 under plain
   pooling. Arms `<GENE>_protocol` are the headline; `_pooled`, `_gnomad_nomix`,
   `_gnomad`, `_sens` and the literature-only arm are retained as sensitivity
   analyses. `make_report.py` gained a protocol branch; no other script changed.
9. Recovery of two killed replays (2026-09-12 06:26, operational; no code
   changed). BRCA2 shards 2 and 3 had their source-recovery replays killed after
   they wedged in regex backtracking on a 207.6 MB and an 879 MB cleaned
   document. `refresh_run_db.py` applies its per-PMID acceptance gate inside the
   replay loop and writes the accepted extraction into `staged_extractions/`,
   leaving only `rebuild_db` for the end, so the staged directories already held
   gated output — originals where a replay was gated, replays where accepted —
   plus papers the replay extracted that had never been in the top-400 set.
   Those directories were copied into each run's `extractions/` (pass-1 originals
   retained in `pass1_pre_staged_backup/`), taking shard 2 from 397 to 661 and
   shard 3 from 400 to 773 extraction JSONs, and the retry pass then rebuilt each
   database and applied variantFeatures enrichment, the wrong-gene false-positive
   quarantine and the trust gate through the normal pipeline. No database was
   hand-assembled and nothing bypassed the acceptance gate.

## Runs

| Gene | Launched (UTC) | Run dir | Status |
| --- | --- | --- | --- |
| HNF1A | 2026-09-09 13:36:42 | `results/grant_e2e_20260909/HNF1A/20260909_083645` (discovery/filter) + `shards/HNF1A_{0,1,2}/` | complete (2 passes, 1,000 papers extracted); analyses in `docs/evidence/grant_e2e_20260909/analysis/HNF1A*` |
| GCK | 2026-09-09 13:36:44 | `results/grant_e2e_20260909/GCK/20260909_083647` + `shards/GCK_{0,1,2}/` | complete (2 passes, 1,050 papers extracted); analyses in `analysis/GCK*` |
| LDLR | 2026-09-09 13:53:21 | `results/grant_e2e_20260909/LDLR/20260909_085324` + `shards/LDLR_{0..3}/` | complete (2 passes, 2,416 papers extracted: 606/602/597/611 per shard; pass 1 ended 2026-09-11 07:39, pass 2 12:13); all four analysis arms complete 2026-09-11 12:22–12:49, in `analysis/LDLR*` |
| BRCA2 | 2026-09-09 13:54:27 | `results/grant_e2e_20260909/BRCA2/20260909_085430` + `shards/BRCA2_{0..4}/` | complete (2 passes; pass 1 ended 2026-09-12 07:07, pass 2 10:18). Final extractions per shard: 397 / 400 / 655 / 743 / 398. Shards 2 and 3 carry the promoted staged replays after their replays were killed (deviation 9); shards 0, 1 and 4 carry theirs in their pass-1 refresh databases. Analyses in `analysis/BRCA2*` |

## Analysis contract (BayesianPenetranceEstimator `iterations/grant_e2e_20260909/`)

1. `build_observations.py`: GVF DB → per (variant, PMID) carrier counts.
   Trusted `penetrance_data` rows only; cohort rows take precedence over
   `individual_records` within a PMID (GVF's aggregation rule); rows with both
   affected and unaffected are `explicit`, rows with affected + total are
   `total_derived`; variants keyed at protein level (1-letter), cDNA otherwise.
2. `join_features.py`: variantFeatures MANE-transcript consequences → per
   protein key: AlphaMissense, REVEL, CADD, BayesDel (mean over alleles), pLDDT,
   gnomAD AF, ClinVar classification (cDNA-matched allele preferred, else best
   review status; allele conflicts flagged).
3. `fit_models.py`: hierarchical logistic-binomial on variant × paper rows with
   an explicit ascertainment mixture (each row is an affected-ascertained series
   with probability π, informed by the paper's recorded ascertainment, or an
   informative observation), a between-variant effect σ and a between-study
   effect τ; single-carrier reports merged per variant; gnomAD anchor rows carry
   no study effect. Fitted by NUTS (PyMC). Fixed-effect sets: base; class
   (truncating); in-silico (class + AlphaMissense); ClinVar-only (class + P/LP,
   B/LB, no-record indicators); in-silico + ClinVar. Posterior per variant =
   the penetrance among carriers not selected for being affected. Five-fold
   cross-validation by variant scores the out-of-fold predictive on held-out
   rows (per-carrier negative log score, Brier, calibration slope). The
   in-sample posterior is reported only with a `_LEAKY` label.
4. `make_report.py`: figures and a Markdown report per gene; `SUMMARY.md`
   across genes.

4b. Split-paper validation (in `fit_models.py`): for every variant reported in
   two or more papers, papers are split in half (alternating by PMID); half A
   updates the feature prior into a variant posterior (importance-weighted
   draws), half B is scored. Comparators on the same held-out carriers: the
   features-only prior, the pooled rate, and the observed rate of the variant's
   ClinVar class learned from the other variants' half A. This is the
   leakage-free variant-level comparison of a count-informed penetrance
   estimate with classification; the in-sample posterior is reported only with
   a `_LEAKY` label. Synonymous and unkeyed rows are excluded from all fits.
4c. gnomAD-anchored arms: `<GENE>_gnomad` (`--gnomad-as-unaffected`, with the
   ascertainment mixture) and `<GENE>_gnomad_nomix` (`--no-mixture`, the lab's
   classic model: literature counts plus gnomAD carriers as unaffected, with a
   between-study effect only). The gnomAD allele count of the max-AF record enters
   as anchored unaffected carriers. The KCNQ1 control (agreement with the curated
   Variant Browser penetrance) decides which anchored model is headlined.
4d. Common variants are retained in every arm (Brett, 2026-09-10): alleles with
   gnomAD AF > 0.1% have a known outcome (at most a GWAS-scale effect), so they
   are the check the model must pass, not an exclusion. Every report carries a
   "known-outcome check": literature fraction (the ascertainment artifact) versus
   posterior penetrance for each common variant. In the mixture arms the
   ascertainment probability also depends on the variant's log10 allele
   frequency (case series of common alleles are ascertainment); with
   proband-dominated literature this term is not identified and the mixture
   still collapses when population anchors are present, which is why the
   anchored no-mixture arm is headlined. An optional `--max-gnomad-af` filter
   exists but is off.
5. Sensitivity arm (`<GENE>_sens`): identical pipeline with single-carrier
   observation rows removed (`--min-n-per-row 2`), i.e. without the
   proband-only case reports that inflate penetrance by ascertainment.
6. Positive control: the same analysis on the lab's existing KCNQ1 extraction
   database (`validation_runs/canonical_baseline/KCNQ1.db`), compared with the
   curated Variant Browser penetrance values
   (`compare_kcnq1_reference.py`).
7. Automation: `auto_analyze.sh` / `auto_analyze_sens.sh` / `auto_analyze_gnomad.sh` wait for each
   `RUN_STATUS.json`, read `active_db`, run `run_gene.sh`, and copy the report,
   run summaries and a dated-rate cost proxy into this directory.

Limitations carried into the report: literature counts pool across papers
(cohort overlap is not resolvable from the extraction), probands are included
(ascertainment inflates penetrance), and ClinVar classes are as of the
variantFeatures snapshot.
