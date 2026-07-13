# GVF Recall History

**Append-only.** This file is the durable, chronological record of recall
benchmarks and the changes that moved them. Add new dated entries at the **top of
the Timeline**; never delete or rewrite past entries (correct a past number only
by adding a new entry that supersedes it, with a note). Live/working status lives
in [`RECALL_STATUS.md`](RECALL_STATUS.md) (the single source of truth for the
*current* numbers) and the re-run procedure in
[`RECALL_REFRESH_RUNBOOK.md`](RECALL_REFRESH_RUNBOOK.md); this file is the *history*.

Grant target: **90% unique-variant recall**, submission **October 2026**.
Validation surface: four cardiac genes **KCNH2 / KCNQ1 / SCN5A / RYR2**
(gold standards in `gene_variant_fetcher_gold_standard/normalized/`). All
recall numbers below are **figures-skipped, DB-observed** scoring via
`scripts/run_recall_suite.py` unless noted, and depend on the Vanderbilt
Elsevier insttoken.

## Headline benchmark trajectory (four-gene aggregate)

| Date | Milestone | Unique-variant recall | Variant-row recall | carriers MAE |
|---|---|---:|---:|---:|
| 2025-11-18 | Project start (initial commit) | — | — | — |
| 2026-05-25 | First multi-gene scored artifact | ~1859/3013 (61.7%) | — | — |
| 2026-05-26 | Pre-SUA-sweep active DBs | 1884/3010 (62.6%) | 3714/6833 (54.4%) | — |
| 2026-05-26 | Post-SUA source-replay sweep | 2206/3010 (73.3%) | 4173/6833 (61.1%) | — |
| 2026-05-26 | Post-rollback + ClinVar/PubTator recovery | 2450/3010 (81.4%) | 5065/6833 (74.1%) | 1.055 |
| 2026-05-29 | cDNA↔protein matcher bridge + count-outlier guard | 2473/3010 (82.2%) | 5160/6833 (75.4%) | ~0.90 |
| **2026-06-05** | **Elsevier mmc supplement acquisition (landed)** | **2523/3010 (83.8%)** | **5350/6833 (78.3%)** | **0.882** |
| 2026-06-06 | refresh_recall re-extract (recall up, MAE regressed) | 2572/3010 (85.4%) | 5423/6833 (79.4%) | 1.274 |
| **2026-06-06** | **Duplicate-penetrance idempotency fix** | **2572/3010 (85.4%)** | **5423/6833 (79.4%)** | **0.614** |
| **2026-06-12** | **PDF-linearized table reconstruction + iter-2 quality gate + targeted KCNQ1 land** | **2590/3010 (86.0%)** | **5514/6833 (80.7%)** | **0.634** |
| **2026-06-12** | **+ targeted lands KCNH2/SCN5A/RYR2 (all four genes)** | **2591/3010 (86.1%)** | **5518/6833 (80.8%)** | **0.615** |
| **2026-07-12** | **Four-gene supplement reconciliation + gated SCN5A land** | **2596/3010 (86.2%)** | **5546/6833 (81.2%)** | **0.614** |
| — | **Target** | **2709/3010 (90.0%)** | — | → 0 |

Gap to the 90% unique-variant target: **113** variants (was 1126 at the 62.6%
starting point; ~89% of the original gap has been closed).

Per-gene unique-variant recall, latest canonical (2026-07-12):
KCNH2 **83.2%**, KCNQ1 **90.5%** (crossed the 90% target), SCN5A **86.8%**,
RYR2 **83.7%**.

---

## Timeline (newest first)

### 2026-07-12 — Four-gene idempotent supplement reconciliation and SCN5A gain
Closed the local supplement-fold gap across KCNH2/KCNQ1/RYR2/SCN5A from 289
papers to **0** (1,577/1,577 papers with convertible supplements folded).
Elsevier acquisition now reconciles individual `mmc` files rather than treating
"any supplement exists" as paper completion, reuses matching supplements across
sibling gene corpora, folds changed sources immediately, and keeps supplement
recovery separate from article-body recovery. The live pass recovered 64 missing
files across 49 papers and updated 427 contexts without re-downloading bodies.

Acceptance-gated re-extraction promoted only SCN5A (PMIDs 12193783, 20541041,
23973953, 29017927): SCN5A unique variants **1020→1027**, rows **2432→2461**,
and PMID coverage **619→622**, with carriers/affected MAE slightly improving and
unaffected MAE holding. KCNH2, KCNQ1, and RYR2 candidates regressed recall or
count coverage and were not promoted. Fresh aggregate: unique variants
**2596/3010 (86.2%)**, rows **5546/6833 (81.2%)**, carriers MAE **0.614**.
Grok and Claude CLI reviews informed the per-file completion model and the
supplement-only acquisition route.

### 2026-06-22 — Annotation-table guard: blank-cell wrong-gene + gnomAD-count-as-carrier (PMID 33013630)

**Garbage-source warning.** PMID **33013630** ("Are Variants Causing Cardiac
Arrhythmia Risk Factors in SUDEP?") is a **hypothesis/review**. Its Table 1 has
columns `Gene | Variant | gnomAD allele count | SIFT | PolyPhen-2 | References`
and **zero patient/carrier data** — it is variant *annotation*, never primary
carrier evidence (extraction rule 4). Treat any carrier counts sourced from this
PMID as invalid.

Two deterministic `pipeline/table_router.py` bugs were turning this annotation
table into carrier evidence (both hit *every* gene, including no-gold targets):
1. **Blank Gene cells (markdown rowspan) were not inherited.** Gene-grouped
   tables name the gene once and leave continuation rows blank; those rows were
   stamped with the *target* gene instead of the real one. So HCN4 `Val759Ile`,
   HCN2 `Pro802Ser`, SCN5A `Ala572Asp`, AKAP9 `Arg2607Gly` … all leaked into
   KCNH2. The flagged symptom was **KCNH2 `Val759Ile` = 870 affected / 0 / 870**
   — the 870 is HCN4's *gnomAD allele count*.
2. **`_looks_like_row_level_clinical_list` fired on the caption's clinical cue
   alone**, minting 1 inferred carrier per annotated row even with no subject
   column.

**Fix (landed in code, gold-free, generalized to the whole class):**
forward-fill blank Gene cells before the gene-filter; and an **annotation
table** — one whose only quantitative columns are variant annotation
(population-frequency allele counts like gnomAD/ExAC/TOPMed **or** in-silico
predictors like SIFT/PolyPhen/CADD/REVEL) with **no subject column** — no longer
infers a carrier per row, regardless of a clinical caption or "score" cue. The
current extractor returns **0** variants from this table (and from a
prediction-score-only table such as `REVEL score | CADD score`). 4 regressions
added in `tests/unit/test_table_router.py` (incl. a guard rail that a genuine
proband list with no annotation columns still infers); full offline suite green
(958 passed / 53 skipped). The literal 870 was already nulled in canonical
`02_strict` by the earlier count-outlier guard; the code fix stops the whole
class at the root.

**Did a good source do a good job? Yes — dropping 33013630 loses zero signal.**
The only genuine KCNH2 variants this review name-drops (`Arg744*`, `Gly924Ala`,
"identified in a SUDEP patient") are redundantly and far better captured by real
primary sources GVF already extracts well, e.g.:
`Arg176Trp` 25+ papers (incl. 19160088=128, 16754261=32/60, 21244686=88),
`Arg744*` (11802537=7, 22338672=1/3), `Gly924Ala` (**38370760**=65),
`Gly749Ala` (38370760=41), `Arg1047Leu` (38370760=67, +10 more). PMID
**38370760** is a standout high-quality cohort table GVF extracted cleanly with
affected/unaffected splits.

**Residual (not yet cleaned — see `TASKS.md`):** the May-29 canonical KCNH2
(21 vars from this PMID: 9 sole-source FP incl. `Val759Ile`, 12 with bogus
counts) and KCNQ1 (14 vars) DBs still carry this contamination; SCN5A/RYR2 are
clean. The live `Val759Ile=870/0/870` row also persists in the stale
`validation_runs/20260518` copies. Un-ingesting PMID 33013630 (the fixed
extractor yields 0 for it) is a precision/FP cleanup deferred to the next refresh.

### 2026-06-22 — Docs source-of-truth cleanup
Consolidated the recall docs so there is no competing "next-run plan" surface:
`TASKS.md` now owns the active Exact-Match Recovery Plan, `RECALL_STATUS.md` is a
short current metrics/failure-split snapshot, and this file owns the dated recall
trajectory. The embedded `RECALL_STATUS.md` session log and stale 2026-05-29
plan were removed from the live status file; their substantive benchmark history
is represented by the timeline entries below.

### 2026-06-12 — Targeted lands extended to the remaining three genes
After the KCNQ1 headline land, `scripts/targeted_land.py` was run on the other
three canonical genes. SCN5A PMID 19716085 promoted a cleaner re-extraction
(+1 unique variant, +4 variant rows; SCN5A **86.3%→86.4%**, and its carriers MAE
improved **0.489→0.454** as over-counted rows were replaced). KCNH2 PMID 32681117
and RYR2 PMID 24136861 each found a candidate whose re-extraction *held* recall —
promoted as cleaner, non-regressive (gold-gated, no recall/row/MAE regression).
Four-gene aggregate **86.0%→86.1%** unique (2591/3010), rows **80.7%→80.8%**
(5518/6833), carriers MAE **0.634→0.615**. Each gene backed up as
`<GENE>.db.before_targeted_land.db` before promotion.

### 2026-06-12 — PDF-linearized table reconstruction → iter-2 quality gate/selector → fast targeted land
PDF supplement tables that fold in linearized (one cell per line) are now
reconstructed into delimited rows before extraction (`ExpertExtractor.
_augment_pdf_linearized_tables`), with a post-extraction cDNA↔protein backfill.
The replay gate/selector were generalized to be **no-gold-safe**: the gene-scoped
deterministic table parse is the structural baseline — a fewer-row re-extraction
is accepted only if it does not lose paired variants and covers ≥85% of the
deterministic parse's positions (over-extracted prior → drop the excess; lossy
under-extraction → reject). Root insight: the blocker was never that extraction
dropped protein — the `refresh_run_db` count gate preferred a stale over-counted
cDNA-only extraction (KCNQ1 30758498: 182 rows, recovers 34/183 gold) over a
clean paired one (100 rows, recovers 85/183). Validated cross-gene (fires on
KCNQ1 + an untuned RYR2 paper; declines SCN5A 0/120). `scripts/targeted_land.py`
lands a single paper's win in **minutes not ~1h** (scan → bridge only candidates
→ re-extract+gate → surgical layer-preserving inject → gold-gated promote).
Landed KCNQ1 PMID 30758498 → **KCNQ1 87.9%→90.5% unique** (crosses target),
aggregate **85.5%→86.0% / rows 79.9%→80.7%**, gated (recall+row-recall+MAE),
backup `KCNQ1.db.before_targeted_land.db`. 787 unit tests incl. a non-cardiac
generalization test. Codex (gpt-5.5, read-only) reviewed the gate/selector design.

### 2026-06-12 — Source-layer backfill + cheap notation junk gate
Added a shared source-layer classifier, explicit `variant_papers.source_layer`
backfill for the four canonical DBs, and a scorer/migration reject for obvious
figure/regex-table junk: gene-symbol-as-variant, <=2-character protein strings,
and residue prose. Backups were written as
`.before_source_layer_20260612_093534` before any DB mutation. Re-scoring the
canonical DBs and the backups (fallback inference, no `source_layer` column)
produced identical aggregate recall/MAE/precision blocks, confirming the scorer
prefers the explicit column without changing semantics.

Recall and MAE were unchanged: uniqV **2572/3010 (85.4%)**, rows
**5423/6833 (79.4%)**, carriers MAE **0.614**, affected MAE **0.523**,
unaffected MAE **1.219**. The junk gate dropped 64 raw rows before scoring
(41 extra-on-gold-PMID rows, 8 counted extras). Headline counted precision is now
`5423/(5423+1660)` = **76.6%**; the loose raw gold-PMID upper bound is
`5423/(5423+13642)` = **28.4%**. The old `manual_or_legacy` bucket is now
`llm_text`.

### 2026-06-11 — Canonical rescore + precision decomposition
Fresh `scripts/run_recall_suite.py` run against the four canonical DBs confirmed
the current aggregate: uniqV **2572/3010 (85.4%)**, rows **5423/6833 (79.4%)**,
PMIDs **1274/1502 (84.8%)**, carriers MAE **0.614**, affected MAE **0.523**,
unaffected MAE **1.219**. Precision-vs-gold-PMIDs remains a false-positive
upper bound at **28.4%**, but decomposition shows only **1668/13683** extra
rows on gold PMIDs carry counts, giving `precision_vs_counted_gold_pmids`
**76.5%**. Added `variant_papers.source_layer` plus per-layer precision output
so ClinVar/PubTator/figure/linker rows are auditable without a manual sample.
Confirmed the MAE non-regression land gate already existed; fixed the stale
history note that claimed it was still missing.

### 2026-06-06 — Duplicate-penetrance idempotency fix (carriers MAE 1.274→0.614, recall preserved)
The 2026-06-05 `refresh_recall` re-extraction lifted recall (uniqV 83.8→85.4%,
rows 78.3→79.4%, patients 80.6→82.1%) but **spiked aggregate carriers MAE
0.882→1.274** — 94% of it KCNQ1, ~98% of that a single paper (PMID 32893267, the
Lahrouchi 2020 LQTS exome cohort). Root cause was **not** junk variants from
no-patient-data papers: the migration wrote one `penetrance_data` row per
supplement-table appearance of a variant (that paper lists each variant across
Table S4 + Table S14 × two ACMG schemes), and the scorer **sums** the rows linked
to each (pmid, variant) — so V254M scored 4× gold (100 vs 25). The count-outlier
guard can't catch it (uniform inflation lifts the per-paper median too, so nothing
trips the >10× rule) and the land gate checks only unique-variant recall, never
MAE, so it promoted silently.
- Fix is idempotency at the write path: exact-duplicate insert guards in
  `migrate_to_sqlite.insert_variant_data` for penetrance/individual/functional/
  phenotypes/variant_metadata, so re-migration or a variant repeated across table
  cells never writes a second identical row. Back-fill for DBs built before the
  guards: `dedup_existing_rows` (`scripts/dedup_db.py`).
- Canonical DBs deduped (backups `<GENE>.db.before_dedup.db`): aggregate carriers
  MAE **1.274→0.614** (below the 0.882 baseline), recall unchanged. Per-gene
  carriers MAE: KCNH2 0.860, SCN5A 0.489, KCNQ1 3.61→**0.897**, RYR2 0.323.
  Removed exact dups: penetrance 6883, phenotypes 7431 across the four DBs.
- Tests: `tests/unit/test_migrate_idempotent.py` (4) + the migration/scoring suite
  (125) green.
- Follow-up landed 2026-06-06: `refresh_recall` now has an MAE non-regression
  land gate (`_promotion_decision`) requiring recall hold plus no carriers,
  affected, or unaffected MAE regression before promotion.

### 2026-06-05 — Supplement acquisition: the real recall lever (+1.6pp uniqV)
Recon (`docs/SUPPLEMENT_ACQUISITION_PLAN.md`) re-framed the problem: the
Cloudflare-blocked publishers (Karger 0.3% / Sage 0.0% of the gap) are
near-irrelevant; the leak was that the **Elsevier full-text API fetches body
only**, dropping the `mmc` supplement mutation tables (~31% of the gap), and the
re-fold layer that makes on-disk supplements visible to Tier-3 was never wired
into the refresh path.
- Wired the supplement re-fold into `refresh_run_db` (Phase 0, no-network) + a
  nested-zip fix + a fold-gap QC counter.
- Added `ElsevierAPIClient.download_supplements` (mmc refs from the authenticated
  XML → open `ars.els-cdn.com` CDN) + `scripts/fetch_elsevier_supplements.py`.
- **Solved the #1 single blocker, SCN5A 29325976** — flagged 2026-06-01 as 87
  missing variants behind a "Cloudflare-blocked `mmc1.docx`". The supplement was
  in fact openly served on the ScienceDirect CDN; fetching + folding +
  re-extracting recovered **64 of 87** variants.
- Batch over the addressable Elsevier papers: 16 gained supplements; re-extract
  recovered 103/126 of their missing variants (per-paper proxy: +175 total).
- **Landed in the canonical DBs** via surgical injection (preserving
  clinvar/pubtator/figure layer rows): aggregate uniqV **82.2%→83.8%** (+50),
  rows **75.5%→78.3%** (+190), patients 77.8%→80.6%, carriers MAE 0.910→0.882.
  Nothing regressed. Backups: `{gene}.db.before_supplements_20260605.db`.
- Figures (the prior day's experiment) gave **0 recall lift** — the gaps are
  splice/intronic variants in supplement tables, not figures. Honest negative.
- T6 (Wiley/Springer supplements) investigated → access-blocked: Springer API key
  is off in `.env`; addressable Wiley DOIs return TDM 403. Reactivates when keys/
  access are restored.

### 2026-06-04 — Corpus consolidation, dashboard, figure fetch (0 recall)
- Consolidated all fetched source into `corpus/<GENE>/<PMID>/` + INDEX
  (6,382 papers), idempotent never-downgrade builder, corpus-as-cache for
  `gvf-run`.
- Static provenance/coverage/adjudication dashboard (`gvf dashboard`).
- Caption-triage à-la-carte figure fetch: 313 OA figures fetched; **0 gold
  recall lift** — confirmed the recall gap is supplements, not figures.

### 2026-06-01 — Acquisition replay & no-gold source QC (per-gene diagnostics)
Per-gene acquisition experiments on separate run dirs measured KCNH2 ~84.15% and
SCN5A ~84.19% uniqV (figures-skipped) after fetched-source replay + ClinVar/
PubTator. Diagnosis flagged **SCN5A 29325976** as the single largest gap (87
variants in a blocked supplement) — solved 2026-06-05. These were diagnostic runs
on `01_off`/turnkey dirs, not all consolidated into the canonical DBs.

### 2026-05-29 — Matcher bridge + count-outlier guard (+0.8pp uniqV, MAE 1.055→~0.90)
- `cli/compare_variants.py`: cDNA↔protein substitution bridge (matches `Y51X`
  gold to `c.153C>A` extracted by implied codon).
- `pipeline/extraction.py`: scope gene-column-less tables by caption (stops
  multi-gene leakage; PMID 26669661 161→24 SCN5A).
- Count-outlier guard cleared study-wide-N values (96 across 28 PMIDs): MAE win,
  no recall change. KCNQ1 was the dominant offender (carriers MAE 2.916).

### 2026-05-26 — SUA source-replay sweep + rollback + recovery (62.6%→81.4% uniqV)
The single biggest jump. Acceptance-gated source replay over the
`source_unbound_available` pool (147 PMIDs), then per-PMID rollback (restore the
higher-variant extraction), DB rebuild, and DB-observed ClinVar/PubTator
recovery. KCNQ1 moved 43.1%→84.2% (LQT supplement parser + native table shapes);
KCNH2 68.3%→82.6%.

### 2026-05-21 — Elsevier insttoken unblock (foundational acquisition win)
Vanderbilt institutional `X-ELS-Insttoken` installed; 242 of 246 previously-
paywalled Elsevier articles across the cardiac genes returned full text. This is
the credential the current cardiac-gene baselines depend on.

### 2025-11-18 → 2026-05 — Pipeline build-out (pre-systematic-benchmark)
Discovery → harvest → Tier1/2 filter → Tier3 extraction → migrate → recall
scoring pipeline; publisher routes (PMC, Elsevier/Wiley/Springer APIs, browser
recovery); gold-standard builders; the recall suite + MAE foundation. Systematic
four-gene recall benchmarking began ~2026-05-25.
