# Recall Status

Last updated: 2026-05-21.

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

`validation_runs/turnkey_e2e_20260518_213934/targeted_additive_plus_scn5a_28341781_v10/recall_score/summary.json`

Aggregate recall from that score:

| Metric | Matched / Gold | Recall |
| --- | ---: | ---: |
| PMIDs | 891 / 1197 | 74.4% |
| Variant rows | 3451 / 5092 | 67.8% |
| Unique variants | 1804 / 2388 | 75.5% |
| Patients/carriers | 8055 / 10926 | 73.7% |
| Affected | 5867 / 8169 | 71.8% |
| Unaffected | 1944 / 2467 | 78.8% |

The target remains greater than 90% across the metrics that matter. The current
state is improved and more general, but not solved.

This scored baseline covers KCNH2, RYR2, and SCN5A. KCNE1 extraction completed
in the four-gene run, but KCNE1 recall should not be claimed until a normalized
per-PMID gold input exists.

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

### Post-token PMID recall (no re-extraction yet)

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

So the token unlock landed and 242 bodies are on disk, but the recall lift will
not be realized until a re-extraction pass consumes them. That is now the
single highest-leverage next step — see "Next Run Plan" below.

## Highest-Yield Remaining Missing PMIDs

These are diagnostic targets from the gold-standard score. Do not hard-code them
into production logic.

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

As of 2026-05-21, the Elsevier insttoken unblock has landed (see
"2026-05-21 Elsevier Insttoken Activation" above). The dominant remaining
barrier shifts from *acquisition* to *converting the now-available source into
DB rows*:

1. **Re-extraction pass** for KCNH2, KCNQ1, RYR2, and SCN5A so the 242 newly
   saved `_FULL_CONTEXT.md` files actually flow into the SQLite databases.
   Without this, PMID recall is stuck at the 2026-05-18 baseline.
2. **Non-Elsevier paywalls still outstanding**:
   - Wiley: `WILEY_API_KEY` in `.env` is revoked. Human Mutation
     (`10.1002/humu.*`) is unreachable until reissued.
   - Karger: no institutional agreement and Cloudflare blocks the Playwright
     stack. TDM request status unknown.
   - Sage/Liebert (PMID 23631430): Cloudflare fingerprint rejection.
3. **Extraction-side issues that recall still depends on** (unchanged from
   the 2026-05-20 audit below): count semantics, affected/unaffected direction,
   cDNA-only normalization, multi-gene supplement filtering, no-gold QC.

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
2. **Next:** Re-extract each gene against the consolidated `pmc_fulltext/`
   directories so the new bodies become SQLite rows. Recommended:
   `gvf gvf-run KCNH2 --email brett.m.kroncke.1@vanderbilt.edu --output <ts>`
   (and analogously for KCNQ1, RYR2, SCN5A). Honor existing `*_FULL_CONTEXT.md`
   files; do not redownload.
3. Rescore PMID/variant-row recall after each gene's re-extraction
   (`scripts/run_recall_suite.py --score --genes <GENE> --db <GENE>=...`).
   Expected: PMID recall rises substantially toward the 90% target, especially
   for SCN5A (currently 70.9%) and RYR2 (currently 68.0%) where the unlocked
   set is largest relative to gold PMID counts.
4. Confirm the recovered files are real article/supplement content, not paywall
   stubs. Use content-quality checks and inspect variant/table density. The
   `*.pre_insttoken_bak` companion files preserve the pre-token stubs for
   before/after comparison.
5. Address the residual non-Elsevier paywalls
   (Wiley key reissue, Karger TDM agreement, Sage CF) as smaller follow-up
   unblocks.
6. Promote only general parser/acquisition fixes to code. Keep per-PMID recovery
   artifacts in validation runs, not production branches.
7. Add no-gold QC output before using the pipeline for new gene-disease runs.
