# Recall Status

Last updated: 2026-05-29.

This note records the post-`19ae63f` state of the recall/generalization push so
the next run can resume from the same baseline without relying on chat history.

## 2026-05-29 Session — Applied Changes & Turnkey Reproduction

Committed on `codex/multigene-recall-improvement` (verify with `git log --oneline`):

- `457c70b feat: cDNA<->protein substitution bridge in recall matcher`
  (`cli/compare_variants.py`). Matches protein SNV/nonsense gold (e.g. `Y51X`) to
  cDNA-only extracted variants (`c.153C>A`) by implied codon, unique-candidate
  guarded. Recovers variants already in the DB that were unmatched on notation.
- `c5c4352 fix: scope gene-column-less tables by caption to stop multi-gene leakage`
  (`pipeline/extraction.py` + `tests/unit/test_extraction_table_parser.py`). A
  section table whose caption names a gene (e.g. "Supplemental Table 1. KCNQ1
  variants") no longer leaks its rows into a different gene's run. Fixes the
  PMID 26669661 case (161 -> 24 SCN5A). 3 regression tests added.

Verified four-gene state after this session (`scripts/run_recall_suite.py`):
PMIDs 83.7%, variant rows **75.4%** (was 74.1%), unique **82.2%** (81.4%),
patients **77.7%**, affected **76.6%**, carriers rows-mode MAE **0.90** (was ~1.06).

Live-DB operations applied this session (data artifacts are gitignored — a fresh
checkout reproduces them by RUNNING the steps, it does not inherit the DBs):

1. **Count-outlier guard** cleared 96 study-wide-N values across 28 PMIDs in all
   four extraction sets (backups: each run's `outlier_guard_backup_20260528/`) and
   in the live DBs (backups: `<GENE>.db.before_guard_20260528.db`). MAE win, no
   recall change. Reproduce on any run:
   `python scripts/apply_count_outlier_guard.py --extractions <run>/extractions --policy clear --backup-dir <run>/outlier_guard_backup --report-out /tmp/r.json`
   then rebuild + rescore.
2. **Supplement recovery** for single-gene supplements via re-extraction from an
   augmented FULL_CONTEXT + `scripts/refresh_run_db.py --replace-db` (RYR2 25844899
   5 -> 28). This also ran the figure-recovery layer (baseline had skipped it).

### Turnkey for a collaborator (fresh checkout)

The two code fixes above are in the branch, so they apply automatically. To run
independently and start testing:

1. Setup per the "Run On Another Computer Or VPN" section in `CLAUDE.md`
   (`pip install -e ".[dev]"`, `playwright install chromium`, `.env` with one LLM
   key + `NCBI_EMAIL`; `ELSEVIER_INSTTOKEN`/`WILEY_API_KEY` unlock paywalled text).
2. End-to-end on a gene: `gvf gvf-run <GENE> --email <you> --output results/`.
3. Score against gold: `python scripts/run_recall_suite.py --score --genes <GENE>
   --db <GENE>=<path-to.db> --outdir recall_metrics/<name>` (reports recall AND
   rows-mode MAE).
4. Local tests (no network): `.venv/bin/python -m pytest tests/ -q`. Note: the 5
   `tests/integration/test_triage.py` failures are pre-existing LLM-auth issues,
   not from these changes.

### Caveats / remaining levers (do not regress)

- Do NOT rebuild KCNQ1/KCNH2 via bare `migrate(extractions)+recovery` — it
  regressed KCNQ1 rows 77.7% -> 57.4% (their recall depends on build-process
  content the extractions alone don't regenerate). Use `refresh_run_db` with a
  real `--pmids-file` (an empty file + `--only-forced-pmids` is rejected).
- The gene-scoping fix covers the MARKDOWN table parser only. The FIXED-WIDTH
  parser (used for pdftotext'd supplements) still leaks multi-gene panels — the
  next general fix, and the reason 10 of 15 staged supplements were held back.
- Image-only supplements (`DownloadImage.aspx`, etc.) need vision/OCR.
- >90% on all metrics remains acquisition-bound (missing supplement/table bodies,
  Elsevier/NEJM, image tables), not extraction-bound. Audit precision before
  pushing recall harder (DBs hold many non-gold rows; figure recovery adds OCR
  noise; gold is incomplete so "not in gold" != wrong).

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
- Branch at verification: `codex/multigene-recall-improvement`
- Remote tracking branch: `origin/codex/multigene-recall-improvement`
- Recent committed changes were intentionally scoped to deterministic extraction,
  acquisition recovery, comparison/normalization, migration safety, status docs,
  and tests.
- There are unrelated local dirty files in this checkout. Do not assume they are
  part of this cleanup unless they are explicitly staged/committed later.

## Current Scored Baseline

Use this artifact as the current scored baseline:

`recall_metrics/post_rollback_recover_20260526_aggregate/summary.json`

This is the honest DB-observed four-gene score after the 2026-05-26 SUA source
replay, per-PMID rollback, DB rebuild, and recovery-layer pass. It includes
KCNH2, KCNQ1, RYR2, and SCN5A. ClinVar/PubTator recovery used DB-observed
PMIDs, not gold PMIDs; KCNH2 v12 manual recovery was not used; the figure layer
was skipped because no `--pmc-dir` was supplied.

The score was re-run on 2026-05-27 against the same four DB paths and
reproduced these values exactly.

| Metric | Matched / Gold | Recall | Gap to 90% |
| --- | ---: | ---: | ---: |
| PMIDs | 1253 / 1502 | 83.4% | 99 |
| Variant rows | 5065 / 6833 | 74.1% | 1085 |
| Unique variants | 2450 / 3010 | 81.4% | 259 |
| Patients/carriers | 14368 / 18719 | 76.8% | 2480 |
| Affected | 9418 / 12475 | 75.5% | 1810 |
| Unaffected | 3274 / 3951 | 82.9% | 282 |

The target remains greater than 90% across the metrics that matter. Unique
variant recall moved substantially, but variant rows, patient/carrier counts,
and affected counts are still well short of 90%.

Per-gene final recall from the same score:

| Gene | PMIDs | Variant rows | Unique variants | Patients | Affected | Unaffected |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| KCNH2 | 232/262 (88.5%) | 813/991 (82.0%) | 438/530 (82.6%) | 2258/2674 (84.4%) | 1429/1635 (87.4%) | 582/749 (77.7%) |
| KCNQ1 | 292/305 (95.7%) | 1317/1741 (75.6%) | 524/622 (84.2%) | 5943/7793 (76.3%) | 3256/4306 (75.6%) | 1258/1484 (84.8%) |
| SCN5A | 595/757 (78.6%) | 2213/3128 (70.7%) | 941/1183 (79.5%) | 4589/6219 (73.8%) | 3490/4876 (71.6%) | 1099/1343 (81.8%) |
| RYR2 | 134/178 (75.3%) | 722/973 (74.2%) | 547/675 (81.0%) | 1578/2033 (77.6%) | 1243/1658 (75.0%) | 335/375 (89.3%) |

Unique-variant gap to 90%:

| Gene | Current | At 90% | Need |
| --- | ---: | ---: | ---: |
| KCNH2 | 438/530 (82.6%) | 477/530 | +39 |
| KCNQ1 | 524/622 (84.2%) | 560/622 | +36 |
| SCN5A | 941/1183 (79.5%) | 1065/1183 | +124 |
| RYR2 | 547/675 (81.0%) | 608/675 | +61 |
| **Aggregate** | **2450/3010 (81.4%)** | **2709/3010** | **+259** |

The sum of per-gene 90% ceilings is +260 because each gene rounds up
independently. The aggregate unique-variant target is +259. SCN5A is the
largest single unique-variant blocker, while variant rows and count metrics have
larger aggregate gaps.

Trajectory:

| Stage | Aggregate unique variants |
| --- | ---: |
| Pre-SUA-sweep active DBs | 1884/3010 (62.6%) |
| Post-SUA sweep, before rollback/recovery | 2206/3010 (73.3%) |
| Post-rollback + DB-PMID recovery | 2450/3010 (81.4%) |
| Target | 2709/3010 (90.0%) |
| Gap | 259 unique variants, 8.6pp |

Do not confuse the `1884/3010` pre-SUA-sweep row with the older
`recall_metrics/post_refresh_20260525/summary.json` score, which was
`1859/3013`. The KCNQ1 gold denominator changed from 625 to 622 between those
artifacts.

KCNE1 extraction completed in the four-gene run, but KCNE1 recall should not be
claimed until a normalized per-PMID gold input exists.

### Rows-Mode MAE Baseline

Per-row mean absolute error (MAE) on matched-variant rows from the same
scored artifact, computed by taking `|gold_count - extracted_count|` for each
matched (PMID × variant) pair and averaging. Lower is better; target is
`MAE → 0`.

| Gene | carriers MAE | affected MAE | unaffected MAE | matched rows |
| --- | ---: | ---: | ---: | ---: |
| KCNH2 | 0.807 | 0.863 | 1.909 | 554 |
| KCNQ1 | 2.916 | 2.121 | 1.491 | 640 |
| SCN5A | 0.603 | 0.438 | 0.763 | 1302 |
| RYR2 | 0.292 | 0.204 | 1.815 | 610 |
| **Aggregate** | **1.055** | **0.744** | **1.216** | **3106** |

KCNQ1 is the dominant count-MAE offender (carriers MAE 2.916, ~10× RYR2's
0.292). Worst single outliers point at one bug class: study-wide N or
cohort/screened totals being assigned to per-variant rows.

| Gene | Top outlier (matched row) | Gold | Extracted | \|err\| |
| --- | --- | ---: | ---: | ---: |
| KCNQ1 | 29622001 G589D | 7 | 453 | 446 |
| KCNQ1 | 17192539 G589D | 50 | 243 | 193 |
| SCN5A | 16453024 S1103Y | 1 | 137 | 136 |
| KCNH2 | 19160088 L552S | 2 | 92 | 90 |
| SCN5A | 20470418 S1103Y | 85 | 26 | 59 |
| KCNQ1 | 23092362 T322A | 2 | 60 | 58 |

**Note on the gold-free path.** MAE is intrinsically gold-dependent. For
no-gold runs, substitute internal-consistency outlier detection: flag any
row whose carrier count is >10× the per-paper median, or where the same
large count is repeated across many distinct variants. The "refuse to reuse
study-wide N unless row-local evidence supports it" guard catches the
7-vs-453 class without needing a gold standard. All new validators are
required to run in both modes (gold-present and gold-absent); GVF must be
turnkey on new genes for which no curated answer exists.

**Count-outlier guard** is implemented as `pipeline/count_outlier_guard.py`
plus `scripts/apply_count_outlier_guard.py`. The detector flags rows whose
count is >10× the per-paper median AND >50 absolute (gold-free, internal
consistency). Policies:

- `--policy off` (or `--dry-run`): detect and report only.
- `--policy flag`: annotate variants with `count_outlier_flags` metadata,
  preserve raw value (default; safest).
- `--policy clear`: also zero the flagged count (raw preserved under flags).

Dry-run validation on production extractions confirms the detector catches
all known outliers: KCNQ1 PMID 29622001 G589D (value 453, 129× median),
SCN5A 16453024 S1103Y (137, 91× median), SCN5A 25904541 (max 3520, 3520×
median). To apply retroactively on a run, use the CLI and then rebuild the
DB with `harvesting.migrate_to_sqlite` and re-score.

## 2026-05-26 Post-Rollback And DB-PMID Recovery

After the SUA sweep, `/tmp/sua_sweep/rollback_and_recover.sh` compared each
forced PMID's current extraction JSON against the sweep backup. It restored the
backup when the backup had more raw variant rows and preserved the sweep output
when the sweep had more rows, then rebuilt each active DB and ran
`scripts/recall_recovery/run_all_layers.py`.

| Gene | Restored JSONs | Preserved sweep wins | Same count | Errors |
| --- | ---: | ---: | ---: | ---: |
| KCNH2 | 11 | 7 | 11 | 0 |
| KCNQ1 | 4 | 2 | 32 | 0 |
| SCN5A | 1 | 8 | 44 | 0 |
| RYR2 | 2 | 2 | 21 | 0 |

The KCNH2 restored JSONs had a raw backup-minus-current delta of 135 variant
rows. KCNH2 unique-variant recall moved from 68.3% pre-SUA-sweep to 64.3%
post-SUA-regression to 82.6% after rollback plus recovery. KCNQ1 moved from
43.1% pre-SUA-sweep to 84.2% after rollback plus recovery; the biggest lift
came from the LQT supplement parser and native table shapes.

The recovery pass was DB-PMID ClinVar/PubTator enrichment. Figures were skipped
because no `--pmc-dir` was supplied. Per-gene recovery additions recorded by
`run_all_layers.py`: KCNH2 ClinVar 1986 / PubTator 612, KCNQ1 3658 / 693,
SCN5A 3940 / 993, RYR2 1114 / 469.

## 2026-05-26 SUA Source Replay Sweep

Acceptance-gated source replay extended to the full `source_unbound_available`
(SUA) pool from the 2026-05-26 grafted baseline's
`paper_disagreement_report.csv`. 147 PMIDs (KCNH2 31 / KCNQ1 38 / SCN5A 53 /
RYR2 25), processed sequentially: `fetch_paywalled.py` (new PMC-supplement
fallback + Elsevier-insttoken path) → `refresh_run_db.py --only-forced-pmids
--replace-db --skip-recovery` → `run_recall_suite.py --score`. Pre-sweep DB
backups preserved per gene as `<GENE>.db.before_refresh_20260526_*.db`.
Fetch success rates: KCNH2 28/31, KCNQ1 35/38, SCN5A 42/53, RYR2 24/25.

Per-gene scored deltas (apples-to-apples, on the same production run dirs):

| Gene | unique_variants PRE | unique_variants POST | Δ unique | variant_rows PRE | variant_rows POST | Δ rows |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| KCNH2 | 362/530 (68.3%) | 341/530 (64.3%) | **-21** | 630/991 (63.6%) | 585/991 (59.0%) | -45 |
| KCNQ1 | 268/622 (43.1%) | 469/622 (75.4%) | **+201** | 566/1741 (32.5%) | 956/1741 (54.9%) | +390 |
| SCN5A | 818/1183 (69.1%) | 885/1183 (74.8%) | **+67** | 1912/3128 (61.1%) | 1985/3128 (63.5%) | +73 |
| RYR2  | 436/675 (64.6%) | 511/675 (75.7%) | **+75** | 606/973 (62.3%) | 647/973 (66.5%) | +41 |
| **Aggregate** | **1884/3010 (62.6%)** | **2206/3010 (73.3%)** | **+322 (+10.7pp)** | **3714/6833 (54.4%)** | **4173/6833 (61.1%)** | **+459 (+6.7pp)** |

KCNQ1 contributed the largest single-gene gain (+32.3pp on unique_variants),
driven by the LQT supplement parser hitting its native table shapes.

**KCNH2 regressed on net (-21 unique_variants)** because seven SUA PMIDs
(11844290, 14661677, 16969682, 19038855, 21185499, 21499742, 26937396)
returned 0 variants from the re-extracted thinner context — stub-as-newest.
`--only-forced-pmids` does not enforce per-PMID acceptance gating against the
prior extraction; it just rewrites. The fix is to per-PMID compare new vs
backup variant counts in `refresh_*/extraction_json_backup/` and restore the
backup for any regressor, then rebuild. See [[feedback-source-replay]].

**Recovery layers were skipped** (`--skip-recovery`) in the sweep score to keep
the source replay tight, so PMID-level recall regressed vs the 2026-05-25
baseline (1100/1502 73.2% vs 1175/1502 78.2%). The rollback/recovery pass above
then restored useful pre-sweep extraction JSONs and reran DB-PMID
ClinVar/PubTator recovery.

Sweep artifacts:

- Per-gene scored summaries: `recall_metrics/sua_sweep_20260526_<GENE>/<GENE>/`
- Per-gene refresh dirs: `<run_dir>/refresh_20260526_*/`
- Sweep run log: `/tmp/sua_sweep/sweep.log` (local, not tracked)

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

The honest measurement against the then-current DBs is below
(`recall_metrics/post_insttoken_20260521/`):

| Gene | PMIDs (matched/gold) | PMID recall | vs 90% target |
| --- | ---: | ---: | ---: |
| KCNH2 | 230/262 | 87.8% | -2.2 |
| KCNQ1 | 90/112 | 80.3% | -9.7 |
| RYR2 | 345/508 | 68.0% | -22.0 |
| SCN5A | 752/1060 | 70.9% | -19.1 |
| **Aggregate (4-gene)** | **1133/1502** | **75.4%** | **-14.6** |

This table is retained for before/after context only. The current baseline is
the post-rollback + DB-PMID recovery score at the top of this file.

## Historical High-Yield Missing PMIDs

These are diagnostic targets from the 2026-05-21 gold-standard score. Do not
hard-code them into production logic. For current misses, use
`recall_metrics/post_rollback_recover_20260526_aggregate/*/missing_in_sqlite.csv`
and the paper-disagreement artifacts generated with that score.

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

As of 2026-05-27, the Elsevier insttoken source unlock, SUA source replay,
per-PMID rollback, and DB-PMID ClinVar/PubTator recovery have been consumed by
the active four-gene DBs. The dominant remaining barriers are no longer stale
SQLite rows:

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

Sequenced after the 2026-05-26 SUA sweep + rollback + recovery pass. All new
work must be turnkey on genes-diseases that lack a gold standard
(`gvf gvf-run NEWGENE ...` should produce usable output plus QC flags
without comparison-based scoring).

**Done:**

1. **DONE 2026-05-21:** `ELSEVIER_API_KEY` + `ELSEVIER_INSTTOKEN` set in
   `.env`; 242/246 paywalled Elsevier candidates unlocked. See
   "2026-05-21 Elsevier Insttoken Activation".
2. **DONE 2026-05-25:** Source-driven refresh/rebuild for KCNH2, KCNQ1, RYR2,
   and SCN5A using `scripts/refresh_run_db.py`.
3. **DONE 2026-05-26:** SUA source replay, per-PMID rollback, DB rebuild, and
   DB-PMID ClinVar/PubTator recovery. Current score is
   `recall_metrics/post_rollback_recover_20260526_aggregate/summary.json`.

**Tier 1 — engineering foundation (no-gold-compatible by construction):**

4. **Acceptance gating in `refresh_run_db.py --only-forced-pmids`.** Before
   overwriting a per-PMID extraction JSON, compare new vs prior variant_count;
   if new < prior, restore from backup and skip. Codifies today's post-hoc
   rollback. Gold-free.
5. **Refuse-to-reuse-study-wide-N guard.** Post-extraction validator: any row
   whose carrier count is >10× the per-paper median is flagged and refused as
   a per-variant count unless row-local evidence confirms it. Highest single
   MAE payoff (catches the KCNQ1 29622001 G589D 7-vs-453 class). Internal
   consistency check, gold-free.
6. **Count provenance preservation + classification.** Preserve raw column
   names + raw values before writing final `patients.count` / `affected_count`
   / `unaffected_count`. Classify each count column as per-variant carrier,
   family count, proband count, cohort total, screened N, affected/case,
   unaffected/control. Refuse assignment when classification confidence is low.
7. **Evidence-card validator (dual-mode).** Two triggers: gold-error
   `|err|>=10` when gold is available; internal-consistency thresholds
   (>10× per-paper median; same large count repeated across rows; null-count
   clustering) when gold is absent. Pulls local table evidence; confirms,
   corrects, or withholds the count.
8. **Track recall AND rows-mode MAE per batch** in `scripts/run_recall_suite.py`
   output. Acceptance gate: accept changes only if recall improves without
   materially worsening MAE. No-gold variant: emit count-distribution
   statistics per PMID + outlier-row counts.

**Tier 2 — recall lever (missing-row attack, ~+200-500 unique variants):**

9. Investigate top missing-row PMIDs:
   `recall_metrics/post_rollback_recover_20260526_aggregate/*/missing_in_sqlite.csv`
   plus the paper-disagreement report. Highest-yield targets: SCN5A 29325976
   (87 rows), KCNQ1 30758498 (90), KCNQ1 17192539 (53). Separate misses into
   missing supplement/table body, figure/pedigree-only, source has variant
   text but extraction missed it, count-only mismatch, and notation mismatch.
10. Re-run recovery layers with `--pmc-dir` (figure layer was skipped on
    2026-05-26). Estimated +20-50 unique variants.
11. Diagnostic: walk the 7 zero-yield KCNH2 SUA papers (11844290, 14661677,
    16969682, 19038855, 21185499, 21499742, 26937396). Likely a
    clinical-case-list table shape the prompt does not match. Implement only
    general parser/prompt fixes uncovered.
12. cDNA-only and splice notation normalization (RYR2 raw `c.`-less forms,
    splice forms like `40-2 A/G`).

**Tier 3 — externally gated:**

13. Karger TDM agreement / Sage CF workaround (SCN5A 23631430 = 24 rows
    alone). Wiley TDM is verified working as of 2026-05-26.

**Operating rules:**

14. For any future source acquisition, rerun
    `python scripts/refresh_run_db.py --gene <GENE> --run-dir <RUN> --replace-db`
    (now acceptance-gated). Do not patch SQLite rows directly.
15. Promote only general parser/acquisition fixes to code. Keep per-PMID
    recovery artifacts in validation runs, not production branches.
16. New gene-disease runs (no gold) must surface no-gold QC output by default:
    gene mentioned with no variants, supplement referenced but not downloaded,
    variant-rich text with zero rows, paywall marker present, many null counts,
    suspicious study-wide-N reuse.
