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
| — | **Target** | **2709/3010 (90.0%)** | — | → 0 |

Gap to the 90% unique-variant target: **186** variants (was 1126 at the 62.6%
starting point; ~83% of the original gap has been closed).

Per-gene unique-variant recall, latest canonical (2026-06-05):
KCNH2 **83.2%**, KCNQ1 **86.8%**, SCN5A **82.8%**, RYR2 **83.4%**.

---

## Timeline (newest first)

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
- Not yet done: an MAE non-regression check on the `refresh_recall` land gate
  (currently `if lm >= bm`, recall-only) so a count regression can't promote again.

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
