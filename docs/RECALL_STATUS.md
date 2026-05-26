# Recall Status

Last updated: 2026-05-25.

This note records the post-`19ae63f` state of the recall/generalization push so
the next run can resume from the same baseline without relying on chat history.

## Source Of Truth

Use this file as the current issue/status tracker for the recall push. Other
top-level handoff docs (`README.md`, `CLAUDE.md`, `CODEX.md`, and `TASKS.md`)
should link here instead of carrying independent live metric tables. If a metric
conflicts with this file, treat this file and the scored artifact below as
authoritative.

Historical recovery docs and scripts are still useful for debugging, but they
are not cold-start evidence unless they explicitly avoid gold-PMID-conditioned
inputs and KCNH2-only manual recovery.

## Current Git State

- Baseline-producing implementation commit: `19ae63f Improve deterministic variant recovery`
- Status-note commit before this cleanup pass: `4b69f50 docs: record current recall status`
- Branch: `main`
- Remote: `origin/main`
- Recent committed changes were intentionally scoped to deterministic extraction,
  acquisition recovery, comparison/normalization, migration safety, status docs,
  and tests.
- There are unrelated local dirty files in this checkout. Do not assume they are
  part of this cleanup unless they are explicitly staged/committed later.

## Current Scored Baseline

Use this artifact as the current scored baseline:

`recall_metrics/post_refresh_20260525/summary.json`

This is the honest DB-observed four-gene score after the 2026-05-25
source-driven refresh/rebuild pass. It includes KCNH2, KCNQ1, RYR2, and SCN5A;
ClinVar/PubTator recovery used DB PMIDs, not gold PMIDs.

| Metric | Matched / Gold | Recall |
| --- | ---: | ---: |
| PMIDs | 1175 / 1502 | 78.2% |
| Variant rows | 3620 / 6833 | 53.0% |
| Unique variants | 1859 / 3013 | 61.7% |
| Patients/carriers | 10911 / 18719 | 58.3% |
| Affected | 6617 / 12475 | 53.0% |
| Unaffected | 2855 / 3951 | 72.3% |

The target remains greater than 90% across the metrics that matter. The current
state is more reproducible and less gold-conditioned, but not close to solved.

Per-gene final recall from the same score:

| Gene | PMIDs | Variant rows | Unique variants | Patients | Affected | Unaffected |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| KCNH2 | 71.8% | 54.2% | 63.8% | 62.8% | 59.0% | 66.0% |
| KCNQ1 | 84.3% | 32.5% | 42.9% | 49.0% | 34.9% | 73.8% |
| RYR2 | 75.8% | 62.3% | 64.6% | 71.8% | 67.9% | 89.3% |
| SCN5A | 78.6% | 61.1% | 69.1% | 63.6% | 62.1% | 69.3% |

KCNE1 extraction completed in the four-gene run, but KCNE1 recall should not be
claimed until a normalized per-PMID gold input exists.

## 2026-05-25 Source-Driven Refresh Pass

Implemented and ran `scripts/refresh_run_db.py` as the safe alternative to
SQLite row patching. The refresh workflow selects stale/under-counted PMIDs from
source artifacts, rewrites canonical extraction JSON, rebuilds a fresh DB from
the full extraction directory, runs DB-observed recovery layers, then atomically
replaces the active DB with a backup.

| Gene | Replay candidates | Successful | Failed/restored | JSONs rebuilt | Final unique-variant recall |
| --- | ---: | ---: | ---: | ---: | ---: |
| KCNH2 | 241 | 232 | 9 | 1129/1129 | 63.8% |
| KCNQ1 | 127 | 106 | 21 | 2381/2381 | 42.9% |
| RYR2 | 101 | 96 | 5 | 638/638 | 64.6% |
| SCN5A | 301 | 242 | 59 | 1449/1449 | 69.1% |

Notable live findings:

- The refresh fixed many stale abstract-only and under-counted extractions, but
  most genes remain far below 90% variant-row and unique-variant recall.
- SCN5A high-yield deterministic table recoveries included PMIDs `26496715`
  (214 variants), `28341781` (53 variants), `30059973` (185 variants), and
  `32893267` (1080 parsed SCN5A rows). The `32893267` multi-gene consortium
  table is high-yield but should be audited for row-level precision and count
  semantics before using it as a headline improvement.
- Some known high-yield misses, including SCN5A `26921764`, still replayed to
  zero variants from available full text. Those are likely supplement/table-body
  or source-content failures, not stale SQLite rows.
- A large set of SCN5A 2022 PMIDs failed the circuit breaker with no target-gene
  variant signal; the script restored their previous JSONs and recorded them in
  `refresh_summary.json`.
- A paper-disagreement report bug with stale/missing context paths was fixed so
  missing local artifacts no longer abort aggregate scoring.

## 2026-05-20 Structure Probe

I rescored a four-gene probe with explicit DB paths for KCNH2, KCNQ1, RYR2, and
SCN5A under `validation_runs/structure_probe_20260520/`. KCNE1 remains excluded
from recall metrics because there is no normalized KCNE1 gold input.

The first KCNQ1 probe using `validation_runs/closeout_20260518_124343/dbs/KCNQ1_base.db`
was rejected because `multicohort_collapse_detector.py` surfaced a large SCN5A
leakage block in KCNQ1 PMID `32893267`. The usable four-gene probe instead uses
`validation_runs/20260517_203904/results/KCNQ1/20260517_204424/KCNQ1_after_32893267_gene_filter_patch.db`.

| Variant matching | PMIDs | Variant rows | Unique variants | Patients | Affected | Unaffected |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Exact | 1127 / 1502 (75.0%) | 4123 / 6833 (60.3%) | 2148 / 3013 (71.3%) | 12160 / 18719 (65.0%) | 7799 / 12475 (62.5%) | 2943 / 3951 (74.5%) |
| Fuzzy | 1136 / 1502 (75.6%) | 4178 / 6833 (61.1%) | 2186 / 3013 (72.6%) | 12299 / 18719 (65.7%) | 7890 / 12475 (63.2%) | 2959 / 3951 (74.9%) |

The exact-to-fuzzy gain is only about one percentage point for variant rows, so
the remaining deficit is not mainly notation matching. The next structure should
collect and score evidence in separate layers:

1. Per-PMID source acquisition status: not attempted, abstract-only, paywall
   stub, recovered full text, recovered browser text, supplement-only, and
   table-body-missing.
2. Paper-level disagreement summary: gold rows, matched rows, missing rows,
   count-mismatch rows, source status, table-reference count, and table-body
   count.
3. Row-level extraction claims: variant identifier, source location, raw count
   fields, derived total/affected/unaffected counts, and whether the count came
   from deterministic table parsing or LLM inference.
4. Small evidence-card claim verification for high-risk count fields, not broad
   second-pass full-paper extraction.
5. No-gold QC for new genes: variant-rich source with zero rows, missing
   supplements, suspicious study-wide count reuse, many null counts, and
   neighboring-gene leakage.

The refreshed paper-disagreement report shows the highest-yield blockers remain
source acquisition and missing table bodies, especially KCNQ1 `19716085`,
KCNQ1 `30758498`, KCNQ1 `19841300`, SCN5A `29325976`, SCN5A `26746457`,
RYR2 `19398665`, RYR2 `30170228`, KCNH2 `29650123`, KCNH2 `15840476`, and
KCNH2/SCN5A `24667783`. Existing claim-verification pilots show the useful
pattern is field-level validation: it fixed KCNH2 `19160088`
affected/unaffected direction when the card included the local table evidence,
while broader count-semantic cases often needed the verifier to withhold
autopopulation rather than invent a corrected value.

## 2026-05-21 Elsevier Insttoken Activation

Vanderbilt's institutional Elsevier token (`X-ELS-Insttoken`) was issued by
Elsevier Data Support (Jun Bautista) and installed into `.env` as
`ELSEVIER_INSTTOKEN`. File permissions on `.env` were tightened to user-only
(`chmod 600`). The token is sent only as a header in `harvesting/elsevier_api.py`
(never in URLs), and `.env` is gitignored.

### Unlock probe results (`scripts/test_insttoken_unlock.py`)

Probed every Elsevier-DOI row in each gene's `paywalled_missing.csv` from the
2026-05-18 turnkey run (plus the 2026-05-17 KCNQ1 run) by calling
`ElsevierAPIClient.fetch_fulltext()` and writing successful bodies to
`{PMID}_FULL_CONTEXT.md` directly inside the existing per-run
`pmc_fulltext/` directory. Pre-existing stub files were preserved as
`*.pre_insttoken_bak`.

| Gene | Probed | Full text | Metadata-only | Error | Saved |
| --- | ---: | ---: | ---: | ---: | ---: |
| KCNH2 | 70 | 69 | 0 | 1 | 69 |
| SCN5A | 78 | 75 | 1 | 2 | 75 |
| RYR2 | 51 | 51 | 0 | 0 | 51 |
| KCNE1 | 22 | 22 | 0 | 0 | 22 |
| KCNQ1 | 25 | 25 | 0 | 0 | 25 |
| **Total** | **246** | **242 (98.4%)** | **1** | **3** | **242** |

The four non-200 rows are not token failures:

- SCN5A PMIDs 15913580 and 17368591: `10.1016/j.cardiores.*` — *Cardiovascular
  Research* migrated from Elsevier to Oxford University Press in 2014. Old DOIs
  keep the `10.1016/` registrant prefix but the articles are now served by OUP
  and unreachable via the Elsevier API.
- SCN5A PMID 35577255: `10.1016/j.jogc.2022.04.013` — *Journal of Obstetrics
  and Gynaecology Canada*; likely outside Vanderbilt's ScienceDirect package.
- KCNH2: 1 unresolved candidate; see per-gene unlock CSV.

### Where the unlocked bodies live

All saved into the run's existing `pmc_fulltext/` (no separate subdirs), so the
standard extraction discovery path picks them up without per-PMID-dir plumbing:

- `validation_runs/turnkey_e2e_20260518_213934/results/KCNH2/20260518_213938/pmc_fulltext/*_FULL_CONTEXT.md`
- `validation_runs/turnkey_e2e_20260518_213934/results/SCN5A/20260518_213938/pmc_fulltext/*_FULL_CONTEXT.md`
- `validation_runs/turnkey_e2e_20260518_213934/results/RYR2/20260518_213938/pmc_fulltext/*_FULL_CONTEXT.md`
- `validation_runs/turnkey_e2e_20260518_213934/results/KCNE1/20260518_213938/pmc_fulltext/*_FULL_CONTEXT.md`
- `validation_runs/20260517_203904/results/KCNQ1/20260517_204424/pmc_fulltext/*_FULL_CONTEXT.md`

Each directory also contains an `insttoken_unlock_results.csv` with the
per-PMID outcome.

### Historical post-token PMID recall (before 2026-05-25 refresh)

The scoring path computes PMID recall as `matched_pmids / gold_pmids` where a
gold PMID counts as matched only if the SQLite DB has at least one extracted
variant row that matches a gold variant on that PMID (see
`cli/compare_variants.py:1949-1970`). Saving 242 full-text files into
`pmc_fulltext/` therefore does not move PMID recall on its own — a re-extraction
pass is required to convert the new text into DB rows.

The honest measurement against the current DBs is below
(`recall_metrics/post_insttoken_20260521/`):

| Gene | PMIDs (matched/gold) | PMID recall | vs 90% target |
| --- | ---: | ---: | ---: |
| KCNH2 | 230/262 | 87.8% | -2.2 |
| KCNQ1 | 90/112 | 80.3% | -9.7 |
| RYR2 | 345/508 | 68.0% | -22.0 |
| SCN5A | 752/1060 | 70.9% | -19.1 |
| **Aggregate (4-gene)** | **1133/1502** | **75.4%** | **-14.6** |

This table is retained for before/after context only. The current baseline is
the 2026-05-25 refresh score at the top of this file.

## Historical High-Yield Missing PMIDs

These are diagnostic targets from the 2026-05-21 gold-standard score. Do not
hard-code them into production logic. For current misses, use
`recall_metrics/post_refresh_20260525/*/missing_in_sqlite.csv` and the
paper-disagreement artifacts generated with that score.

### KCNH2

- `29650123`: 20 missing rows
- `15840476`: 19 missing rows
- `24667783`: 18 missing rows
- `11854117`: 9 missing rows
- `10973849`: 9 missing rows
- `27871843`: 8 missing rows

### RYR2

- `19398665`: 27 missing rows
- `30170228`: 24 missing rows
- `27452199`: 23 missing rows
- `22677073`: 15 missing rows
- `23595086`: 11 missing rows
- `29447731`: 9 missing rows

### SCN5A

- `29325976`: 87 missing rows, 109 carriers
- `26669661`: 27 missing rows, 299 carriers
- `20541041`: 26 missing rows
- `26921764`: 26 missing rows
- `26746457`: 25 missing rows
- `23631430`: 24 missing rows

## Current Main Barrier

As of 2026-05-25, the Elsevier insttoken source unlock has been consumed by the
source-driven refresh/rebuild pass. The dominant remaining barriers are no
longer stale SQLite rows:

1. **Missing or incomplete source bodies, especially supplements/tables.**
   Several refreshed papers now have real full text but still replay to zero or
   too few variants, which means the decisive variants are likely in missing
   supplements, table bodies, figures/pedigrees, or publisher content not
   represented in the saved markdown.
2. **Count semantics and table precision.** High-yield deterministic table
   recovery now finds many variants, but carrier, affected, and unaffected
   counts remain well below 90%. Multi-gene consortium tables such as
   `32893267` need row-level precision and count audits before being used as
   headline gains.
3. **Non-Elsevier paywalls still outstanding**:
   - Wiley: **resolved 2026-05-26.** Off-VPN probe with the current
     `WILEY_API_KEY` returned a 635 KB full-text PDF for `10.1002/humu.21126`
     (302 → CDN 200, `application/pdf`, 14 pages). The earlier "revoked"
     claim was stale; Human Mutation and other Wiley journals are reachable
     via `harvesting/wiley_api.py` with no network change.
   - Karger: no institutional agreement and Cloudflare blocks the Playwright
     stack. TDM request status unknown.
   - Sage/Liebert (PMID 23631430): Cloudflare fingerprint rejection.
4. **Notation and extraction-side generalization** still matter: affected/
   unaffected direction, cDNA-only normalization, multi-gene supplement
   filtering, and no-gold QC.

For ScienceDirect/Elsevier papers, paywall stubs / 403 / accepted-manuscripts
issues should now resolve at the API tier (Tier 2) instead of falling through
to Tier 3.5 browser scraping. The browser stack remains useful for
non-Elsevier paywalled sources.

## Downstream Extraction Issues Still Worth Fixing

Acquisition is the biggest blocker, but extraction is not finished. The important
remaining extraction problems are general, not PMID-specific:

1. Count semantics are still fragile.
   - Some papers report study-wide `N`, case/control counts, probands, families,
     and per-variant carriers in adjacent columns.
   - The extractor needs to preserve the raw count fields and classify them before
     converting to `patients.count`, `affected_count`, and `unaffected_count`.
   - This matters without gold standards because inflated carrier counts can look
     superficially complete.

2. Affected/unaffected direction can still be wrong.
   - The score shows many matched variants with count mismatches.
   - Add explicit parsing for column labels like affected, symptomatic, proband,
     case, control, unaffected, asymptomatic, and healthy controls.
   - Treat controls as unaffected only when the table context supports that.

3. cDNA-only and odd nucleotide notations need broader normalization.
   - RYR2 missing rows include raw cDNA-like forms without `c.` prefixes and splice
     style rows such as `40-2 A/G`.
   - This should be fixed as notation normalization, not as per-paper recovery.

4. Multi-gene supplement tables need robust row-level gene filters.
   - `19ae63f` removed the fixed cardiac-gene title filter and added non-cardiac
     regression coverage.
   - Still add more fixtures for multi-gene oncology/cardiomyopathy panels so new
     gene runs do not pull neighboring gene rows.

5. Figure/pedigree extraction remains a secondary route.
   - It can recover variants, but prior figure sweeps were slow and modest yield.
   - It should remain a fallback for papers where the source audit says variants
     are figure-only or pedigree-only.

6. New gene runs need no-gold QC.
   - For genes without gold standards, the pipeline should rank suspicious papers:
     gene mentioned plus no variants, supplement referenced but not downloaded,
     variant-rich full text but zero extracted rows, paywall marker present,
     many variants with null counts, and likely study-wide counts reused per row.

7. Oversized source contexts need bounded deterministic cleanup.
   - Data Scout now prefers an existing `*_CLEANED.md` sibling when raw
     `*_FULL_CONTEXT.md` exceeds the configured size guard.
   - Keep that preference in the no-gold path so large source artifacts do not
     waste model budget or bury useful tables.

## Generalization Audit

The `19ae63f` commit was checked for obvious PMID-specific production logic.
Committed production code does not branch on the high-yield recovery PMIDs such
as `29325976`, `28341781`, `33606749`, `27566755`, or `25904541`.

Remaining committed PMID mentions are examples or fixtures:

- `scripts/fetch_paywalled.py` usage examples
- `scripts/extract_figure_variants.py` usage example
- test fixtures and unit tests

The key generalization fixes already committed:

- No default LQTS PMID list in `scripts/fetch_paywalled.py`.
- No personal fallback NCBI email; `ENTREZ_EMAIL` or `NCBI_EMAIL` is required.
- Table-title gene filtering now uses generic HGNC-like symbols, not a fixed
  cardiac gene set.
- TP53/KRAS/BRAF/PIK3CA hotspot filtering now keeps those variants when that
  gene is the target.
- Institutional cookie domains are configurable through `GVF_COOKIE_DOMAINS` or
  `GVF_SSO_COOKIE_DOMAINS`.
- PMC proof-of-work supplement gates fail explicitly if unsolved rather than
  silently saving challenge HTML as source content.

## Bloat/Stale Audit

Current tracking surface to keep:

- This file: current metrics, current blockers, current next run plan.
- `docs/NEW_GENE_RUNBOOK.md`: no-gold workflow and QC expectations.
- `scripts/run_recall_suite.py` and `cli/compare_variants.py`: scored validation
  and gold-standard comparison entrypoints.
- `gvf gvf-run`, `scripts/fetch_paywalled.py`, `harvesting/pmc_pow.py`,
  `pipeline/extraction.py`, and `pipeline/table_router.py`: current recovery
  and extraction surfaces.

Historical or diagnostic surfaces:

- `scripts/recall_recovery/merge_v12_db.py`: KCNH2-only manual rescue; never
  cold-start evidence.
- `scripts/recall_recovery/ingest_clinvar.py` and `ingest_pubtator.py`: current
  only when run with DB-observed PMIDs. Gold-PMID enrichment is diagnostic only.
- `scripts/recall_recovery/recover_paywall_oa.py`: source-acquisition helper,
  now generalized to explicit input/output paths; not a default pipeline layer.
- `docs/TESTING.md`: reusable test prompt, not a current metrics source.
- `validation_runs/`, `results/`, `recall_metrics/`, and Python `__pycache__/`
  trees: generated/local artifacts and intentionally not the source of record
  unless a path is referenced from this file.

Already-addressed items that should not remain open:

- Default LQTS PMID list in `scripts/fetch_paywalled.py`.
- Personal fallback NCBI email in `scripts/fetch_paywalled.py`.
- Fixed cardiac-only table-title gene filter in `pipeline/extraction.py`.
- Target-gene hotspot filtering that dropped true TP53/KRAS/BRAF/PIK3CA variants.
- PMC proof-of-work challenge pages being treated as usable supplement content.

## Next Run Plan

1. **DONE 2026-05-21:** `ELSEVIER_API_KEY` + `ELSEVIER_INSTTOKEN` set in
   `.env`; 242/246 paywalled Elsevier candidates unlocked and saved as
   `_FULL_CONTEXT.md` into each gene's `pmc_fulltext/`. See
   "2026-05-21 Elsevier Insttoken Activation".
2. **DONE 2026-05-25:** Source-driven refresh/rebuild for KCNH2, KCNQ1, RYR2,
   and SCN5A using `scripts/refresh_run_db.py`; current score is
   `recall_metrics/post_refresh_20260525/summary.json`.
3. Use `recall_metrics/post_refresh_20260525/*/missing_in_sqlite.csv` plus the
   paper-disagreement report to separate remaining misses into: missing
   supplement/table body, figure/pedigree-only, source has variant text but
   extraction missed it, count-only mismatch, and notation mismatch.
4. Audit high-impact deterministic table recoveries, especially multi-gene
   consortium papers (`32893267`), for row-level gene filtering and count
   semantics before presenting them as reliable count improvements.
5. Address the residual non-Elsevier paywalls
   (Wiley key reissue, Karger TDM agreement, Sage CF) as smaller follow-up
   unblocks.
6. For any future source acquisition, rerun
   `python scripts/refresh_run_db.py --gene <GENE> --run-dir <RUN> --replace-db`
   instead of patching SQLite rows.
7. Promote only general parser/acquisition fixes to code. Keep per-PMID recovery
   artifacts in validation runs, not production branches.
8. Add no-gold QC output before using the pipeline for new gene-disease runs.
