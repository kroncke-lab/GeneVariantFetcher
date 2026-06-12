# Recall Status

Last updated: 2026-06-12.

This note records the current recall/generalization state so the next run can
resume from the same baseline without relying on chat history.

## How to read the numbers in this file (single source of truth)

Use the **2026-06-12 Current Canonical Baseline** immediately below for all live
planning. Older dated sections are retained as history and may contain numbers
that are superseded by the current baseline.

No other doc may restate a recall number; they link here.

## 2026-06-12 Current Canonical Baseline

Fresh run of `scripts/run_recall_suite.py` against the four canonical DBs:

- `results/KCNH2/e2e_working_20260529_full/02_strict/KCNH2.db`
- `validation_runs/20260517_203904/results/KCNQ1/20260517_204424/KCNQ1.db`
- `validation_runs/turnkey_e2e_20260518_213934/results/SCN5A/20260518_213938/SCN5A.db`
- `validation_runs/turnkey_e2e_20260518_213934/results/RYR2/20260518_213938/RYR2.db`

Four-gene aggregate:

| Metric | Matched / Gold | Recall | Gap to 90% |
| --- | ---: | ---: | ---: |
| PMIDs | 1274 / 1502 | 84.8% | 78 |
| Variant rows | 5423 / 6833 | 79.4% | 727 |
| Unique variants | **2572 / 3010** | **85.4%** | **137** |
| Patients/carriers | 15369 / 18719 | 82.1% | 1479 |
| Affected | 9929 / 12475 | 79.6% | 1299 |
| Unaffected | 3435 / 3951 | 86.9% | 121 |

Rows-mode MAE:

| Count field | Sum abs error / N | MAE |
| --- | ---: | ---: |
| Carriers | 2141 / 3485 | **0.614** |
| Affected | 1547 / 2956 | **0.523** |
| Unaffected | 323 / 265 | **1.219** |

Per-gene current recall:

| Gene | PMIDs | Variant rows | Unique variants | Patients | Affected | Unaffected | carriers MAE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| KCNH2 | 230/262 (87.8%) | 820/991 (82.7%) | 441/530 (83.2%) | 2256/2674 (84.4%) | 1404/1635 (85.9%) | 599/749 (80.0%) | 0.860 |
| KCNQ1 | 285/305 (93.4%) | 1408/1741 (80.9%) | 545/622 (87.6%) | 6472/7793 (83.0%) | 3407/4306 (79.1%) | 1313/1484 (88.5%) | 0.897 |
| SCN5A | 620/757 (81.9%) | 2429/3128 (77.7%) | 1021/1183 (86.3%) | 5016/6219 (80.7%) | 3832/4876 (78.6%) | 1184/1343 (88.2%) | 0.489 |
| RYR2 | 139/178 (78.1%) | 766/973 (78.7%) | 565/675 (83.7%) | 1625/2033 (79.9%) | 1286/1658 (77.6%) | 339/375 (90.4%) | 0.323 |

Headline precision is `precision_vs_counted_gold_pmids`, which restricts the
denominator to extra rows on gold PMIDs that carry at least one extracted count:
`5423 / (5423 + 1660) = 76.6%`. The looser raw proxy remains useful only as a
false-positive **upper bound**: `5423 / (5423 + 13642) = 28.4%`.

Why the raw proxy is pessimistic:

- 11,982 / 13,642 (88%) current extra-on-gold-PMID rows have zero patient counts
  and are ClinVar/PubTator-style linkage attributions rather than count-bearing
  paper extractions.
- Only 1,660 extra rows carry any carrier/affected/unaffected count.
- About 97% are well-formed variants absent from the count-curated gold packet,
  not malformed output.
- The scorer now rejects 64 obvious figure/regex-table junk rows before scoring
  (gene-symbol-as-variant, <=2-character protein notation, residue prose). This
  removed 41 extra-on-gold-PMID rows, including 8 counted extras, with recall and
  MAE unchanged.
- 53 structural/CNV rows are real biology but currently unmatchable by the
  variant matcher.

Interpretation: recall gains are mostly adding real signal; the 28.4% proxy
overstates true false positives by roughly 7x.

Current per-layer precision proxy from the scorer. The four canonical DBs now
have explicit `variant_papers.source_layer` values; the fallback-derived score
on the `.before_source_layer_20260612_093534` backups matches this block exactly.

| Source layer | Matched DB rows | Extra rows | Counted extra rows | precision_vs_gold_pmids | precision_vs_counted_gold_pmids |
| --- | ---: | ---: | ---: | ---: | ---: |
| clinvar | 464 | 2491 | 71 | 15.7% | 86.7% |
| figure | 236 | 465 | 39 | 33.7% | 85.8% |
| llm_table | 838 | 477 | 275 | 63.7% | 75.3% |
| llm_text | 455 | 823 | 168 | 35.6% | 73.0% |
| mixed | 1949 | 392 | 176 | 83.3% | 91.7% |
| pubtator | 12 | 159 | 0 | 7.0% | 100.0% |
| regex_table | 1256 | 3864 | 929 | 24.5% | 57.5% |
| regex_text | 213 | 4971 | 2 | 4.1% | 99.1% |

Current failure-mode split from `paper_disagreement_report.csv` and
`failure_taxonomy_report.py`:

| Failure mode | Missing rows | What it means |
| --- | ---: | --- |
| source_missing_or_stub | 568 | paper/source never landed or only a stub landed |
| source_abstract_only | 284 | abstract was available, but mutation tables/body were missing |
| available_source_underextraction | 248 | usable source exists but extraction missed rows |
| source_missing_table_bodies | 184 | supplement/full text landed without the relevant tables |
| partial_underextraction | 82 | some rows extracted, table not exhausted |
| count_semantics | 36 | variant present but carrier/affected/unaffected semantics wrong |

The next implementation lane is therefore source/table acquisition and binding
for the top PMIDs, then count semantics. For count semantics, point the
count-guard/evidence-card validation at `regex_table`: it is the dominant
count-bearing FP surface (`929` counted extras, 57.5% counted precision), while
figure rows are comparatively small and clean (`39` counted extras, 85.8%).
Broad model tweaking is lower yield.

## 2026-06-05 Session — Supplement Acquisition Landed In Canonical DBs

Recovered Elsevier `mmc` supplement mutation tables (the full-text API fetches
body only) + folded on-disk supplements, then surgically injected the 12
value papers into the canonical four-gene DBs (preserving clinvar/pubtator/figure
layer rows; verified no individual_records loss; backups at
`{gene}.db.before_supplements_20260605.db`). Method + full details:
`docs/SUPPLEMENT_ACQUISITION_PLAN.md`.

Four-gene aggregate (figures-skipped scoring, canonical DBs):

| Metric | Prior (2026-05-29) | Now (+supplements) | Δ |
|---|---|---|---|
| unique_variants | 2473/3010 (82.2%) | **2523/3010 (83.8%)** | +50 (+1.6pp) |
| variant_rows | 5160/6833 (75.5%) | **5350/6833 (78.3%)** | +190 (+2.8pp) |
| patients | (77.8%) | **80.6%** | +2.8pp |
| affected | (76.6%) | **78.7%** | +2.1pp |
| carriers MAE | 0.910 | **0.882** | better |

> **2026-06-06 — duplicate-penetrance idempotency fix (supersedes the carriers
> MAE above).** The 2026-06-05 `refresh_recall` re-extraction of multi-table
> cohort papers double-counted carriers (KCNQ1 PMID 32893267 V254M scored 4× gold
> = 100 vs 25), spiking aggregate carriers MAE to **1.274** (KCNQ1 3.61) while
> recall held. Root cause: migration wrote one `penetrance_data` row per
> supplement-table appearance of a variant and the scorer SUMS the rows linked to
> each (pmid, variant). Fixed by exact-duplicate insert guards in
> `migrate_to_sqlite.insert_variant_data` (migration is now idempotent) plus a
> `dedup_existing_rows` back-fill (`scripts/dedup_db.py`). Canonical DBs deduped
> (backups: `<GENE>.db.before_dedup.db`): **aggregate carriers MAE 1.274 → 0.614**
> (below the 0.882 baseline), **recall unchanged**. Per-gene carriers MAE now:
> KCNH2 0.860, SCN5A 0.489, KCNQ1 0.897, RYR2 0.323.

Per-gene unique_variants: KCNH2 82.8→83.2, KCNQ1 85.7→86.8, SCN5A 80.1→**82.8**
(largest gain — the SCN5A 29325976 mmc1.docx alone added 64 of its 87 missing
variants), RYR2 81.9→83.4. Aggregate gap-to-90% unique variants: 236 → **186**.
The two Cloudflare-blocked publishers (Karger 0.3% / Sage 0.0% of the gap) were
confirmed near-irrelevant to recall.

## 2026-06-01 Session — KCNH2/RYR2/SCN5A Acquisition Replay & No-Gold Source QC

KCNH2 missing-recall diagnosis showed the biggest immediate lever was source
acquisition, not another broad extraction run. The acquisition worklist now
reports PMID recall separately for PMIDs selected for full-text acquisition and
PMIDs that actually landed usable full text after `fetch_paywalled.py`.

KCNH2 acquisition experiment on
`results/KCNH2/e2e_working_20260529_full/01_off`:

- Gold-assisted acquisition worklist selected 39 PMIDs for fetch and 64 PMIDs
  for acquisition or source rebinding.
- `fetch_paywalled.py` landed 20 usable full-text PMIDs. This is
  `20/262 = 7.6336%` PMID recall for actual usable full text downloaded, versus
  `39/262 = 14.8855%` PMID recall for PMIDs merely selected for fetch.
- Staged refresh replayed those 20 sources without mutating canonical
  extractions; 17 passed the acceptance gate and 3 were rejected for variant
  regressions (`19038855`, `16399053`, `26022375`).
- Comparable KCNH2 score after ClinVar + PubTator, figures skipped:
  PMIDs **89.69%**, variant rows **84.96%**, unique variants **84.15%**.
  This improves the prior KCNH2 baseline of roughly 87.79% / 82.54% / 82.64%.

SCN5A and RYR2 were then audited with the same gold-assisted acquisition
worklist builder against the 2026-05-18 turnkey run:

- SCN5A selected 148 PMIDs for fetch (`148/757 = 19.5509%` selected-for-fetch
  PMID recall) and 191 PMIDs for acquisition or source rebinding
  (`191/757 = 25.2312%`). The fetch bucket represents 391 missing distinct
  variants, with another 242 in refresh-replay and 28 in manual/blocked source
  work.
- SCN5A existing-source replay first targeted the 41 already available
  refresh-replay PMIDs: 35 passed, 5 were regression-gated
  (`10973849`, `21964171`, `24606995`, `26173111`, `28491684`), and 1 was
  explosion-gated (`26669661`). After ClinVar + PubTator with figures skipped,
  this intermediate DB scored PMIDs **597/757 (78.86%)**, variant rows
  **2244/3128 (71.74%)**, unique variants **950/1183 (80.30%)**.
- SCN5A `fetch_paywalled.py` then ran on the 148 fetch-selected PMIDs. The run
  was interrupted before a native summary, but the output-dir summarizer
  recovered 143 attempted PMIDs and 44 usable full-text PMIDs
  (`44/757 = 5.8124%` actual usable-fulltext-downloaded PMID recall), versus
  `148/757 = 19.5509%` selected-for-fetch recall. The same 44 PMIDs were
  accepted by staged refresh; three 76-character placeholders (`10727647`,
  `15671436`, `17420262`) were rejected by the 500-character source gate.
- SCN5A after fetched-source replay plus DB-observed ClinVar + PubTator,
  figures skipped: PMIDs **629/757 (83.09%)**, variant rows
  **2331/3128 (74.52%)**, unique variants **996/1183 (84.19%)**, patients
  **4673/6219 (75.14%)**, affected **3546/4876 (72.72%)**, unaffected
  **1127/1343 (83.92%)**. Compared with the prior current scored SCN5A row
  below, this is +55 unique variants.
- Residual SCN5A diagnosis after that score: `29325976` is the largest single
  remaining gap (87 missing distinct variants). The main text explicitly says
  all SCN5A mutations are in Supplemental Table 2, but the supplement download
  is Cloudflare/redirect blocked (`mmc1.docx`, `supplements_downloaded=0`).
  Gold-assisted and no-gold worklists now route this pattern as
  `missing_variant_supplement` source acquisition instead of generic
  underextraction; the no-gold latest SCN5A source QC found 19 such PMIDs.
- Targeted `29325976` v3 experiment: `fetch_paywalled.py` now tries an
  authenticated Playwright/browser supplement fallback after requests-based
  supplement download fails. In this session Chrome cookie decryption loaded
  0 cookies; the supplement still failed (`Exceeded 30 redirects`, browser
  request 403, browser navigation timeout). Elsevier API body recovery succeeded
  and staged refresh accepted the source, but extraction found 0 variants. This
  confirms the 87 missing variants are in the blocked Supplemental Table 2, not
  in the API body. Residual post-fetch accounting for the current SCN5A DB:
  108 PMIDs selected for fetch (`108/757 = 14.2668%` selected-for-fetch PMID
  recall), 115 selected for source acquisition or rebinding (`15.1915%`), and
  the targeted run produced 1 usable-fulltext/source-refresh-successful PMID
  (`29325976`, `1/757 = 0.1321%`) with no variant recovery.
- Targeted SCN5A residual Wiley/Springer API-route batch on 2026-06-02:
  13 PMIDs were attempted from the 108-PMID residual fetch queue, covering 41
  missing distinct variants (`13/757 = 1.7173%` fetch-attempted PMID recall).
  `fetch_paywalled.py` landed 2 usable full-text PMIDs (`16643399`,
  `24667783`; `2/757 = 0.2642%` usable-fulltext-downloaded PMID recall),
  while 9 Wiley PMIDs stayed Cloudflare/TDM-403 blocked and 2 Springer PMIDs
  failed the content-quality gate. Corrected staged refresh on top of
  `refresh_20260601_154228/staged_extractions` attempted those 2 PMIDs,
  accepted 1 (`16643399`, 0 extracted variants; `1/757 = 0.1321%`
  source-refresh-successful PMID recall), and regression-gated `24667783`
  (`prior=25`, `new=1`). Corrected no-figure rescore remained unchanged:
  PMIDs **629/757 (83.09%)**, variant rows **2331/3128 (74.52%)**, unique
  variants **996/1183 (84.19%)**.
- Follow-up `24667783` diagnosis showed the fetched Springer/Nature source had
  the right `.doc` supplement, but antiword lost the SCN5A table rows. macOS
  `textutil` conversion recovered the Word table; scanner gene-context filtering
  removed the cross-gene `P73T` leak; the extraction artifact filter now rejects
  gene symbols such as `SCN5A` as malformed protein notations; and the scanner
  plus artifact guard now accept protein range deletions such as
  `p.Lys1505_Gln1507del` / `K1505_Q1507del`. The first forced one-PMID staged
  refresh recovered five previously missing gold rows (`N406S`, `E1784K`,
  `Y1795C`, `P1332L`, `M1498T`). The deletion-parser follow-up then recovered
  the remaining `24667783` gold row (`P.K1505_Q1507DEL`), leaving no missing
  rows for that PMID in the final PubTator-layer missing list. Final no-figure
  SCN5A score after DB-observed ClinVar + PubTator:
  PMIDs **629/757 (83.09%)**, variant rows **2337/3128 (74.71%)**, unique
  variants **996/1183 (84.19%)**, patients **4679/6219 (75.24%)**, affected
  **3552/4876 (72.85%)**, unaffected **1127/1343 (83.92%)**. Compared with the
  pre-`24667783` no-figure score, this is +6 matched variant rows and no
  PMID/unique movement.
- RYR2 selected 37 PMIDs for fetch (`37/178 = 20.7865%`) and 59 PMIDs for
  acquisition or source rebinding (`59/178 = 33.1461%`). `fetch_paywalled.py`
  was interrupted before writing a final summary, but the outcome summarizer
  recovered partial output-dir results: 17 usable full-text PMIDs landed
  (`17/178 = 9.5506%` actual usable-fulltext-downloaded PMID recall), covering
  80 missing distinct variants. Staged source refresh accepted 34/41 PMIDs
  (`34/178 = 19.1011%` source-refresh-successful PMID recall).
- RYR2 staged refresh replayed 41 source candidates: 34 passed, 4 failed source
  quality checks, and 3 were rejected by the regression gate
  (`31913406`, `33500567`, `33897349`). The rebuilt DB migrated 646/646
  extraction JSONs.
- RYR2 after DB-observed ClinVar + PubTator, figures skipped:
  PMIDs **143/178 (80.34%)**, variant rows **792/973 (81.40%)**, unique
  variants **573/675 (84.89%)**, patients **1669/2033 (82.10%)**, affected
  **1325/1658 (79.92%)**, unaffected **344/375 (91.73%)**. Compared with the
  prior current scored RYR2 row below, this is +26 unique variants.

New turnkey/no-gold tooling added in this session:

- `scripts/recall_audit/source_acquisition_audit.py` builds a source worklist
  without reading gold standards or recall discrepancies.
- `gvf gvf-run` now runs a `source-qc` step and writes `source_qc/` artifacts:
  source worklist, fetch queue, source override CSV, and summary JSON.
- `gvf gvf-run --source-recovery` is now the opt-in no-gold source recovery
  loop: fetch the source-QC queue, write selected-vs-successful acquisition
  outcome accounting, staged-refresh accepted source overrides, and score/report
  against the refreshed DB when recovery layers are not skipped.
- `scripts/recall_audit/summarize_acquisition_outcome.py` now supports both
  gold recall mode and no-gold worklist coverage mode. In gold mode it reports
  both selected-for-fetch PMID recall and actually usable-fulltext-downloaded
  PMID recall. When passed a `refresh_run_db.py` summary, it also reports
  `source_refresh_attempted` and `source_refresh_successful` PMID recall.
- Source/acquisition audits detect explicit prose pointers such as "All SCN5A
  mutations are listed in Supplemental Table 2" when the corresponding
  supplemental table body is absent. `fetch_input.csv` now includes these
  missing target-gene supplement cases, with DOI recovery from run-local
  `result.json` / artifact JSON where possible.
- `scripts/refresh_run_db.py` supports repeatable `--source-override-csv` and
  `--stage-extractions` so fetched sources can be replayed in a disposable
  extraction copy.
- `scripts/refresh_run_db.py` also supports repeatable/comma-separated
  `--replay-model` so a bounded replay can override a stuck or experimental
  default Tier 3 model without editing `.env`.
- `pipeline/source_quality.py` now shares the extractor's 500-character minimum
  full-text gate so tiny placeholder files are not counted as usable sources.
- `scripts/fetch_paywalled.py` preserves publisher supplement-link counts
  across Elsevier API body fallback and prints `supp_links` separately from
  `supp_downloaded`, so blocked supplement acquisition is visible in the run
  log and summary JSON.
- `scripts/fetch_paywalled.py` also threads an optional browser-backed
  supplement download fallback into the paywall enricher. If a supplement does
  download from a publisher page and the Elsevier API later replaces a stub main
  body, the recovered supplement markdown is appended to the API body instead of
  being overwritten.
- Publisher API fallback is now generalized for Elsevier, Wiley, and Springer:
  when a browser strategy is absent, empty, stubby, or supplement-only,
  `fetch_paywalled.py` tries the matching API client and writes successful
  bodies through the same quality gate and artifact schema.

Gold-free KCNH2 source QC smoke test on the same run produced 5,017 run PMIDs:
690 currently have usable full text, 4,327 are selected for fetch, and 89 are
selected for source refresh. This is broad operational coverage, not recall.
After tightening the denominator to exclude raw PubMed discovery lists by
default, gold-free source QC on the 2026-05-18 turnkey run produced:
SCN5A 1,496 actionable run PMIDs with 1,145 usable full text current
(`76.54%`), 345 selected for fetch (`23.06%`), and 6 routed to
manual/blocked; RYR2 638 actionable run PMIDs with 518 usable full text current
(`81.19%`) and 120 selected for fetch (`18.81%`). Use
`--include-discovery-pmids` only for intentionally broad diagnostics.

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
top-level handoff docs (`README.md`, `CLAUDE.md`, `AGENTS.md`, and `TASKS.md`;
`CODEX.md` is now a pointer stub) should link here instead of carrying
independent live metric tables. Companions: `docs/RECALL_HISTORY.md` is the
append-only benchmark/change history (the trajectory over time — never delete
its entries); `docs/RECALL_REFRESH_RUNBOOK.md` is the idempotent re-run procedure
(`scripts/refresh_recall.py`) for when new papers/permissions arrive. If a metric conflicts with this file, this file
is authoritative; within this file, the most recent dated per-gene session wins
over the older four-gene aggregate table (see "How to read the numbers" at the
top).

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

## Historical Scored Baseline (2026-05-26)

> **Dated 2026-05-26 — four-gene aggregate.** This is the latest *joint* score
> from the 2026-05-26 recovery pass. It is superseded by the 2026-06-11
> canonical baseline at the top of this file.

Aggregate scored-baseline artifact:

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

### Where the source bodies live

As of 2026-06-04, all fetched source (full text + supplements + figures) was
consolidated and deduplicated into a single discoverable home:

- **`corpus/<GENE>/<PMID>/`** — `{PMID}_FULL_CONTEXT.md`, `{PMID}_CLEANED.md`,
  `{PMID}_figures/`, `{PMID}_supplements/`, `{PMID}_artifacts.json`.
- **`corpus/INDEX.json`** / **`corpus/INDEX.csv`** — gene → PMID → paths,
  bytes, figure/supplement counts, and which run each best copy came from.

6,257 distinct (gene,PMID) papers (KCNQ1 2396, SCN5A 1496, KCNH2 1299, RYR2 638,
KCNE1 428), verified as a complete superset of every prior `pmc_fulltext/`
location. The old per-run `pmc_fulltext/` trees under `results/` and
`validation_runs/` were removed after consolidation; their run DBs and
`extractions/` remain in place. Build/refresh with the corpus builder
(`scripts/build_source_corpus.py` if promoted, else the one-off used on
2026-06-04).

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
the 2026-06-11 canonical score at the top of this file.

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
2. **Count semantics and regex-table precision.** High-yield deterministic
   table recovery now finds many variants, but carrier, affected, and unaffected
   counts remain well below 90%. With explicit source-layer provenance,
   `regex_table` is the count-bearing FP surface to validate first
   (`929` counted extras, 57.5% counted precision). Multi-gene consortium
   tables such as `32893267` need row-level precision and count audits before
   being used as headline gains.
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

## Next Run Plan (2026-05-29 — Codex-reviewed reconciliation)

This supersedes the prior tiering. The plan was investigated end-to-end (six
subsystem deep-dives) and independently reviewed by Codex (gpt-5.5, read-only,
file-verified). All new work must be turnkey on genes that lack a gold standard
(`gvf gvf-run NEWGENE ...` should produce usable output plus QC flags without
comparison-based scoring). Every DB-mutating step carries an explicit gate.

### Three corrections that shape the plan

1. **MAE machinery is not fully dormant.** The dormant pieces are the
   `count_classifier` / `count_outlier_guard` / `evidence_card` modules (built +
   unit-tested, only callable via standalone scripts; `count_provenance` is
   emitted by the prompt but dropped at `migrate_to_sqlite`). But the live
   extraction path ALREADY runs `_suppress_repeated_study_wide_counts`
   (`pipeline/extraction.py:2979`, called at :4111/:4132) and a default-on
   second-model adjudication `_maybe_adjudicate_extraction` (:4124,
   `enable_tier3_ensemble_qa` default True). Wiring the dormant modules must be
   reconciled with this existing suppression — flag-only first, never
   clear-by-default.
2. **The biggest clean recall win needs no network.** ~431 missing rows sit
   behind stale source binding (`source_desync`: a thin stub `_FULL_CONTEXT.md`
   is bound while a larger one is already on disk), and ~290 more sit in
   physical `{pmid}_supplements/` dirs that re-extraction never re-reads
   (discovery globs only `*_DATA_ZONES.md` / `*_CLEANED.md` / `*_FULL_CONTEXT.md`).
3. **Re-binding is gated on count, not yet on semantics.** The replay path
   (`refresh_run_db.py`) now has TWO internal-consistency gates: a regression
   gate (rejects new < prior) and a variant-explosion gate
   (`_is_variant_explosion`, default-on at :522 — rejects new > 10x prior AND
   >=400 absolute AND >=300 delta, the signature of garbage / wrong-paper /
   multi-gene leakage). `_largest_context_path` also filters abstract fallbacks
   and sub-500-character sources. What is still MISSING is a semantic / on-gene
   quality gate: none of these prove the recovered variants are real and
   on-target. C1 below should add that semantic gate — NOT lower the explosion
   thresholds, which would suppress legitimate supplement recoveries.

### Done

- **2026-05-21:** Elsevier insttoken; 242/246 paywalled Elsevier unlocked.
- **2026-05-25:** Source-driven refresh/rebuild (KCNH2/KCNQ1/RYR2/SCN5A).
- **2026-05-26:** SUA source replay, per-PMID rollback, DB rebuild, DB-PMID
  ClinVar/PubTator recovery. Baseline:
  `recall_metrics/post_rollback_recover_20260526_aggregate/summary.json`.
- **2026-05-27:** Acceptance gating in `refresh_run_db.py` (commit 8d51f79);
  rows-mode MAE tracked alongside recall in `run_recall_suite.py` (fbbeef2).
- **2026-05-28/29:** Count-outlier guard run manually (MAE 1.06→0.90);
  cDNA↔protein matcher bridge (457c70b); markdown caption gene-scoping
  (c5c4352). Classifier / guard / evidence-card modules committed but NOT yet
  wired into `gvf-run` (see Tier 1 C4 and Tier 2 B2).

### Tier 0 — Measure first (cheap, unblocks judgment on every other lever)

- **P1 — Precision metric + artifact in the scoring harness: DONE 2026-06-11.**
  `compute_precision_summary` is gold-PMID-restricted, frames the result as a
  false-positive upper bound, emits `unmatched_db_rows_on_gold_pmids.csv`, and
  now reports `counted_extra_on_gold_pmids`,
  `precision_vs_counted_gold_pmids`, and `by_source_layer` so recovery/linker
  rows are decomposable without manual sampling.

### Tier 1 — Cheap code wins, no-gold-safe

- **C3 — Fixed-width parser gene-scoping: DONE 2026-06-06.** The fixed-width
  (pdftotext) parser now scopes captions with explicit cardiac genes + LQT
  aliases, generic digit-bearing gene tokens via `_gene_symbol_tokens`, and
  contextual all-letter gene captions such as `LMNA mutations`. Regression tests
  cover BRCA1 and LMNA off-target leakage plus a no-gene prose caption that must
  remain claimable by the target.
- **R1a — Figure-reader precision gate: mostly done / lower priority.** Figures
  run on every `gvf-run` and now carry explicit `source_layer='figure'`; the
  figure bucket is comparatively clean after notation gating (`39` counted
  extras, 85.8% counted precision). Do not prioritize figure adjudication ahead
  of `regex_table` count validation.
- **C1 — Re-bind source discovery to the largest on-disk `_FULL_CONTEXT.md`**
  per PMID (~431 rows, no network), GATED: source-quality + title/PMID match +
  variant-explosion guard (reject suspicious count blow-ups), not just the
  existing under-count gate. Generalizes the manual per-gene re-bind.
  (`refresh_run_db.py` discovery + `_largest_context_path`.)
- **C2 — Supplement-fold reader.** Convert on-disk `{pmid}_supplements/` files
  (reuse `supplement_processing_service.process_supplement_files` with a no-op
  download callback) and append to `_FULL_CONTEXT.md` before extraction
  (~290 rows, no network). MUST land after C3 (multi-gene supplements leak
  otherwise).
- **C4 — Wire count guard + classifier into the live path, FLAG-ONLY, behind
  env knobs defaulting OFF** (`COUNT_GUARD_POLICY` / `COUNT_CLASSIFIER_POLICY`).
  Annotate only; reconcile with the existing
  `_suppress_repeated_study_wide_counts` (no double-clear). Then measure
  MAE+recall deltas; only after that consider selective `clear` with
  evidence/acceptance gates. (`pipeline/steps.py`, `config/settings.py`.)

### Tier 2 — Bounded builds

- **B1 — Reference-sequence cache + real translation + validation gate.** No
  reference CDS/protein exists in the repo today; validation is length-only and
  the cDNA↔protein bridge is integer-division codon math. Fetch + cache per-gene
  CDS/protein (biopython is already a dependency; reuse the Entrez pattern), add
  real translation in aggregation, and a position-in-range + reference-residue
  identity reject. Covers both the "reference CDS in aggregation" wish and the
  gold-free precision gate; generalizes to no-gold genes.
- **B2 — Evidence-card consumer with teeth** on high-uncertainty counts
  (confirm/correct/withhold). EXTENDS the existing default-on adjudication
  rather than building one; depends on C4 flags. (`pipeline/claim_verifier.py`.)
- **B3 — Affected/unaffected column-label parsing** (unaffected MAE ~1.2 is the
  worst; direction/label confusion: affected/symptomatic/proband/case vs
  unaffected/asymptomatic/control).

### Tier 3 — Multi-session / externally gated

- **R1b — Image-only / variant-rich-but-zero detector** to TARGET vision at the
  silent-zero-at-high-confidence class (the zero-variant QC today only re-runs a
  TEXT model). Feeds the R1a-gated vision reader.
- **R2 — Residual paywalls:** Karger (Cloudflare), Sage/Liebert (CF
  fingerprint). Elsevier + Wiley already covered.
- **R3 — cDNA-only / splice notation normalization** (RYR2 raw `c.`-less forms,
  splice forms like `40-2 A/G`).

### Recommended sequence

`P1 + C3 (+ R1a) → gated C1 → gated C2 → C4 flag-only → measure → selective C4
clear → B1 → B2 + B3 → R1b → R2 → R3`. C3 lands before C1/C2 (protects both from
multi-gene leakage); C4 stays flag-only until measured.

### Operating rules

- Every DB-mutating step carries a gate (under-count AND over-count/quality);
  no silent truncation — log what was dropped and why.
- For any source acquisition, rerun
  `python scripts/refresh_run_db.py --gene <GENE> --run-dir <RUN> --replace-db`
  (acceptance-gated). Do not patch SQLite rows directly.
- Promote only general parser/acquisition fixes to code. Keep per-PMID recovery
  artifacts in validation runs, not production branches.
- New gene runs (no gold) must surface no-gold QC by default: gene mentioned
  with no variants, supplement referenced but not downloaded, variant-rich text
  with zero rows, paywall marker present, many null counts, suspicious
  study-wide-N reuse.
