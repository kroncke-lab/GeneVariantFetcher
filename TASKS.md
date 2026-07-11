# GVF Tasks

## Current Focus
- Improve honest cold-start recall for KCNH2, RYR2, and SCN5A to 90% across PMIDs, variant rows, unique variants, patients, affected, and unaffected counts.
  - Current metrics and failure split are in `docs/RECALL_STATUS.md`.
  - The forward plan lives here. Do not duplicate live recall tables in this file.
- **Measurement loop is multi-gene now** — score KCNH2, RYR2, and SCN5A from `gene_variant_fetcher_gold_standard/normalized/*_recall_input.csv`; KCNQ1 scoring remains available when its gold input is in scope.
- **KCNE1 is extraction-only until a gold input exists.** There is no `KCNE1_recall_input.csv` in the current gold-standard package, so KCNE1 recall cannot be claimed yet.

## Autonomy at Scale (updated 2026-07-10) — trust-gate v1 landed; these are the open levers

Fleet-honesty (non-zero exit + `RUN_STATUS.json` on stage failure, fail-closed
regression gate) and the per-fact **trust gate v1** are MERGED (PR #140 + #142):
`pipeline/trust_gate.py` soft-quarantines gold-free-implausible counts into a
trusted/quarantine two-tier DB, default-on in `gvf-run` (Step 3.7). Full plan and
rationale in [`docs/AUTONOMY_ROADMAP.md`](docs/AUTONOMY_ROADMAP.md). The open
levers — what actually stands between here and "trust it unattended at scale":

- [ ] **Make the trusted tier the DEFAULT** the scorer, reports, and downstream
      tools read — they still read the WHOLE DB, so the tiering is informational
      today. Highest-leverage single step: add `trust_tier` filtering to
      `cli/compare_variants.py` (`extract_sqlite_data`) and report trusted-tier
      precision + the quarantine rate.
- [ ] **Count-role / evidence-type axis.** Promote count role (patient / cohort
      total / control / population) to a first-class field on `penetrance_data`
      and have the trust gate score role consistency (case-control ≠ carrier
      counts) — keeps BRCA case-control from looking like a cardiac penetrance
      error. Sharpens Lever 2 below.
- [ ] **Calibrate the trust gate per stratum** (cardiac + BRCA + one cold gene)
      with `scripts/precision_sample.py`; a pooled CI hides BRCA failure under
      cardiac volume. One-time human labeling per stratum, not per paper.
- [ ] **BRCA1/BRCA2 readiness** (see `project_brca_generalization_readiness` in
      memory): register BRCA2 in `PROTEIN_LENGTHS` / `gene_metadata`; handle
      BRCA-class notation BY CLASS (legacy BIC indels, IVS, exon CNVs), not a
      per-gene alias file; build a ClinVar silver standard to score against.
- [ ] **Fail-CLOSED validation on unknown genes.** `utils/variant_normalizer.py`
      `validate_position` returns True when a gene's length is unknown (MLH1,
      BRCA2) — quarantine those (reason `no_gene_length`) instead of trusting.
- [ ] **Fleet-scale nightly + acceptance metrics.** Wire `make regression-gate`
      into a real nightly/cron (needs the gitignored canonical DBs → local /
      self-hosted); define the gate's own acceptance metrics (trusted-tier
      precision on cardiac AND BRCA; quarantine rate not exploding on cold genes).
- [ ] **Fold the legacy guards into the trust record.** `carrier_guard` (NULLs
      counts) and `vf`-quarantine (DELETEs rows) should write `trust_tier` and
      preserve raw counts instead of destroying data.

The supplement-acquisition recall ceiling (Levers 1/3 below) is the other half —
finding more of the real evidence — and is orthogonal to the trust work.

### Stage 5 — study record + variant-class widening (landed 2026-07-10)

Teach extraction to model the study and stop dropping non-missense notations:

- [x] **C** — keep delins + IVS through `_filter_extraction_artifacts` (and
      migration regexes); unit tests in `test_extraction_table_parser.py`.
- [x] **A1–A2** — `study_design` / `ascertainment` / `cohort_source` /
      `population` / `study_summary` in prompts + `extraction_metadata` SQLite
      (also persists previously dropped `study_type`).
- [x] **A3–A5** — study context in `count_classifier` + `trust_gate`
      (`study_type_mismatch`, strengthened `population_count`; rule_version
      `tg2-*`) + distributions in `trust_report.py`.
- [x] **B** — `variant_class` + `structural_description` schema/persist;
      scanner exon/BIC/delins patterns; normalizer + `to_canonical_form`
      splice/delins/structural keys; synthetic gold seed
      `benchmarks/curated_extraction_eval/gold_overrides/STRUCTURAL_WIDENING_recall_input.csv`.

**Follow-ons:** re-extract or LLM backfill study fields on existing DBs; expand
structural/splice gold beyond the synthetic seed so B recall is measurable on
cardiac missense-heavy gold.

## Exact-Match Recovery Plan (2026-06-12) — START HERE

Root-caused decomposition of the exact-match gap to the manual gold curation.
This is the tracked checklist — update boxes as items land. **Both Claude and
Codex: resume the recall push from here.**

Gap by root cause after the 1B parser land: acquisition (PMID absent) ~426 ·
supplement/table not in source ~647 · in-source-but-not-extracted ~236 · matcher
~67 · count-mismatch is tracked in `docs/RECALL_STATUS.md`.

### Lever 1 — Supplement/table-body recovery on already-fetched papers (~883 rows; dominant)
Of the remaining approximate bucket: **~647 are NOT in the captured source text**
(supplement table never fetched/folded); **~236 ARE in the text but extraction
missed them** after the 1B parser land.

- [x] **1B — Table-reconstruction preprocessor (parser track, NO fetching).**
      Landed 2026-06-12. PDF-linearized supplement tables are reconstructed from
      one-cell-per-line text into delimited rows before deterministic table
      extraction. Strict no-MAE-regression land accepted the KCNQ1 replay that
      improved recall and carriers MAE; KCNQ1 `30758498` reconstructs but remains
      unpromoted until count semantics are fixed. Candidates that regressed
      unique/row recall or MAE were withheld. Current score is in
      `docs/RECALL_STATUS.md`.
- [ ] **1A — Supplement fetch+fold audit (acquisition/fold track, ~647 rows).**
      Many top PMIDs reference an eTable/Supplementary Table but have 0 supplement
      files on disk (e.g. KCNQ1 `17192539`: 21KB body, 48 missing variants absent
      from source). Audit "supplement referenced but absent/unfolded"; fetch
      (Elsevier mmc CDN works; Wiley/Springer gated); convert xlsx/docx/pdf; fold
      into FULL_CONTEXT before extraction. Ensure the CLEANED size-guard never drops
      a mutation-table region.
- [ ] **1C — Targeted re-extraction of the top extraction-gap PMIDs**,
      acceptance-gated via `scripts/refresh_recall.py`. Regenerate the ranked
      list from the post-1B score before running; the pre-1B list included KCNQ1
      `17192539` `30758498` `23631430` `19490272`; SCN5A `15840476` `20541041`
      `23631430` `21273195` `25163546` `24631775` `29325976`; RYR2 `19398665`
      `27452199`; KCNH2 `29650123` `16922724`.

### Lever 2 — Count-role attribution on matched rows
- [ ] Study-wide-N reuse is ~gone (~15). Residual is column-role confusion
      (affected/proband/case vs unaffected/asymptomatic/control vs study total).
      Point the count classifier / evidence-card validator at the `regex_table`
      layer; current counts and precision are in `docs/RECALL_STATUS.md`.
- [ ] **Un-ingest the PMID 33013630 annotation-table contamination** (precision
      FP cleanup). This review is a **garbage carrier source** — gnomAD allele
      counts + SIFT/PolyPhen, no patient data (see `docs/RECALL_HISTORY.md`
      2026-06-22). The root-cause `table_router.py` fix (blank-gene-cell inherit
      + no infer on population-annotation tables) is **landed**, so the current
      extractor yields 0 for it; the **May-29 canonical DBs are not yet cleaned**:
      KCNH2 21 vars (9 sole-source FP incl. `Val759Ile`, 12 bogus-count), KCNQ1
      14 vars (SCN5A/RYR2 clean). Drop this PMID's `variant_papers`/`penetrance_data`
      rows + orphaned nodes at the next refresh (reversible, backup first). All
      real KCNH2 variants it name-drops are already well-covered by primary
      sources (e.g. `Gly924Ala` 38370760=65), so this is pure precision gain, no
      recall loss.

### Lever 3 — Acquisition of absent-PMID rows (426; SCN5A has 96 absent PMIDs)
- [ ] Wiley/Springer supplements + remaining paywalls; restore Springer key; Wiley
      supplement via EZproxy route. Access-gated (see Blocked).

### Lever 4 — Matcher notation lanes (67 rows; cheap, already-extracted)
- [ ] indel `GxxxDel`/`LxxxIns`/`c.x_yDel` (27); same-codon-position substitutions
      that don't match — isoform offset / 1↔3-letter (27); compound `A + B` gold
      rows (8); frameshift (3); splice/IVS (2). Extend `cli/compare_variants.py`.

### Lever 5 — Structural/CNV matching lane (53 rows)
- [ ] Exon deletions, breakpoints, translocations: real biology, unmatchable today.
      Add a structural lane or exclude from the precision denominator.

## Completed
- [x] Springer API integration (full-text retrieval working)
- [x] Wire `gene_literature/supplements/` UnifiedSupplementFetcher into orchestrator
- [x] Fix all test failures and collection errors (18 → 0)
- [x] Canonical form comparison for variant matching
- [x] Remove 50-paper download cap (was silently dropping 350+ papers)
- [x] Fail-open Tier 1 filter (min_keywords 2→1, no-abstract papers pass through)
- [x] Add source-completeness report (Step 3.5: `source_completeness.json`)
- [x] Add zero-variant QA flagging and re-extraction (GVF_QA_MODEL env var)
- [x] Improve cohort/screening table carrier count extraction prompts
- [x] Fix pytest environment (unit tests passing locally)
- [x] **Paywall recovery pipeline (2026-05-11)**
  - [x] `harvesting/browser_html/cookie_loader.py` reads Chrome cookies → Playwright format
  - [x] `harvesting/browser_html/authenticated_pool.py` injects cookies + uses `channel="chrome"`
  - [x] `harvesting/browser_html/content_quality.py` quality gate (paywall phrases + body floor + padding)
  - [x] `harvesting/browser_html/dom_extract.py` DOM walker + CF detection + `pick_better_markdown`
  - [x] Modernized strategies — AHA (`#bodymatter` + CF wait + reload), Elsevier (`/abstract→/fulltext` + imprint domains), generic (DOM walker fallback + CF wait), Wiley (CF wait)
  - [x] `scripts/fetch_paywalled.py` CLI with DOI resolution, per-domain pacing, CF retry sweep, **PMC fallback**
  - [x] Verified against 10 audit cases: 10/10 correctly classified by the gate
- [x] **`cli/compare_variants.py` matcher refinements (2026-05-11)**
  - [x] `_positions_compatible(a, b)` with subset rule (handles `p.Pro926AlafsTer14` (926, 14) vs `P926fsX` (926,))
  - [x] `_cdna_indel_protein_positions(cdna)` codon math (c.842dupG → aa 281)
  - [x] `_protein_indel_position(variant)` protein indel parser
  - [x] `find_best_match(..., consumed: Set[str])` for 1-to-1 enforcement at match time
  - [x] cDNA/protein bridge pass wired into greedy 1-to-1 matcher
- [x] Cherry-picked from heuristic-chatelet worktree
  - [x] Scanner: parenthesized HGVS pattern `p.(Arg176Trp)` (commit `f6307da`)
  - [x] Harvest: log silent supplement-scrape failures with PMID context (commit `d3a639f`)
  - [x] Extraction: reject DATA_ZONES that override full text with garbage (commit `e6daa24`)
- [x] **Clinical mutation-list table extraction (2026-05-14)**
  - [x] Deterministically infer one carrier per row for clinical mutation-list tables with variant/proband/family/patient context but no explicit count column
  - [x] Preserve affected vs unaffected counts when the row text indicates asymptomatic/control/unaffected status
  - [x] Keep pure nomenclature/list tables out of extraction unless clinical row context is present
- [x] **Bulk KCNH2 stale-context reharvest (2026-05-14)**
  - [x] Reharvested all 213 stale/missed KCNH2 FULL_CONTEXT files in `results/KCNH2/20260506_102238/pmc_fulltext/`
  - [x] Recovered 17 real contexts; 196 still require paywall/manual/source access
  - [x] Wrote audit artifacts under `results/KCNH2/20260506_102238/reharvest_out_20260514/`
- [x] **KCNH2 v12 manual-recovery score + matcher patch (2026-05-15)**
  - [x] Scored `KCNH2_v12_manual_recovery_20260515.db` (PMID/row/unique-variant/patient recall recorded in `docs/RECALL_STATUS.md`)
  - [x] Added frameshift canonicalization for extraction spellings like `fsTer`, `fs/185`, `fs+*49`, and malformed `AlaX14`
  - [x] Extended cDNA/protein bridge matching to multi-base cDNA indel ranges; recovered `P926fsX` (PMID 26496715) and `P1034fsX` (PMID 29622001)
- [x] **Turnkey closeout and honest enrichment cleanup (2026-05-18)**
  - [x] Added `gvf gvf-run` as the one-command doctor -> extract -> recovery layers -> report path.
  - [x] Scored KCNH2/KCNQ1/RYR2/SCN5A under `validation_runs/closeout_20260518_124343/`.
  - [x] Changed ClinVar/PubTator recovery layers to default to DB-observed PMIDs instead of gold PMIDs; gold-PMID enrichment is now explicit diagnostic mode.
  - [x] Turnkey recovery layers now back up the SQLite DB before mutation.
- [x] **Four-gene turnkey long run and cleanup pass (2026-05-19)**
  - [x] Ran KCNH2, KCNE1, RYR2, and SCN5A end to end under `validation_runs/turnkey_e2e_20260518_213934/`.
  - [x] Confirmed KCNH2, RYR2, and SCN5A remain below 90% on all scored recall metrics except RYR2 unaffected.
  - [x] Confirmed KCNE1 extraction completes but cannot be scored without a gold recall input.
  - [x] Removed default KCNH2 v12 auto-merge, removed gold-PMID leakage from default recovery, added DB backups, broadened LLM provider checks, filtered non-article LinkOuts, added retry/challenge handling, and guarded Data Scout against oversized raw contexts.
- [x] **Wiley TDM API verified working (2026-05-26)** — off-VPN probe with the current `WILEY_API_KEY` in `.env` returned a real 635 KB PDF for `10.1002/humu.21126` (HTTP 302 → CDN 200, `application/pdf`, 14 pages). The earlier "revoked" claim was stale. Human Mutation and other Wiley journals are reachable now via `harvesting/wiley_api.py` without any network change.
- [x] **Source-driven DB refresh path (2026-05-25)**
  - [x] Added `scripts/refresh_run_db.py` as the safe alternative to SQLite row patching: select stale/under-counted source artifacts, rewrite canonical extraction JSON, rebuild the DB, then run recovery layers.
  - [x] Made recovery layers run in DB-PMID mode without requiring a gold standard, so no-gold genes still get ClinVar, PubTator, and figure recovery plus internal QC artifacts.
  - [x] Added source fingerprints and abstract-only fallback guards so reruns skip already-refreshed JSON unless the source changes.
- [x] **SUA replay rollback + DB-PMID recovery baseline (2026-05-26)**
  - [x] Rolled back per-PMID extraction JSON regressions from the source-unbound-available sweep while preserving sweep wins.
  - [x] Rebuilt active KCNH2/KCNQ1/SCN5A/RYR2 DBs and ran DB-observed ClinVar/PubTator recovery without gold-PMID enrichment.
  - [x] Re-scored the four-gene aggregate at `recall_metrics/post_rollback_recover_20260526_aggregate/summary.json`; current live metrics and gaps are in `docs/RECALL_STATUS.md`.
- [x] **KCNH2 acquisition replay + no-gold source QC (2026-06-01)**
  - [x] Added a gold-free `source_acquisition_audit.py` and wired it into `gvf-run` as `source-qc`.
  - [x] Added post-fetch outcome summarization that reports selected-for-fetch PMIDs separately from successfully usable-fulltext-downloaded PMIDs.
  - [x] Added `gvf-run` source recovery (source QC, paywall fetch, acquisition outcome summary, and staged refresh without requiring a gold standard), now default-on; pass `--no-source-recovery` for a fast PMC/free-text-only pass.
  - [x] Replayed 20 KCNH2 fetched sources through staged extraction; 17 passed acceptance gates after ClinVar+PubTator, figures skipped. Scored recall is in `docs/RECALL_STATUS.md`.
- [x] **RYR2 acquisition replay (2026-06-01)**
  - [x] Gold-assisted worklist selected 37 PMIDs for fetch and 59 PMIDs for source acquisition or rebinding.
  - [x] Post-fetch summarizer recovered partial interrupted-fetch output: 17 usable full-text PMIDs actually landed after the stricter 500-character source gate.
  - [x] Staged refresh replayed 41 candidates; 34 passed, 4 failed source-quality checks, and 3 were regression-gated.
  - [x] RYR2 scored after DB-observed ClinVar+PubTator, figures skipped; recall is in `docs/RECALL_STATUS.md`.
- [x] **SCN5A acquisition replay (2026-06-01)**
  - [x] Gold-assisted worklist selected 148 PMIDs for fetch and 191 PMIDs for source acquisition or rebinding.
  - [x] Existing-source replay accepted 35/41 PMIDs and improved SCN5A recall after DB-observed ClinVar+PubTator, figures skipped (see `docs/RECALL_STATUS.md`).
  - [x] Partial `fetch_paywalled.py` output produced 44 usable full-text PMIDs; the outcome summary now reports selected-for-fetch recall and actual refresh-successful recall separately (values in `docs/RECALL_STATUS.md`).
  - [x] Fetched-source replay accepted 44/47 candidates and scored after DB-observed ClinVar+PubTator, figures skipped (recall in `docs/RECALL_STATUS.md`).
  - [x] Residual audit now routes explicit missing target-gene supplement pointers to acquisition. SCN5A PMID `29325976` has 87 missing distinct variants behind a blocked Supplemental Table 2 download; latest no-gold source QC finds 19 `missing_variant_supplement` PMIDs.
  - [x] Added a Playwright/browser supplement download fallback and reran `29325976`: Elsevier API body recovery + staged refresh succeeded, but `mmc1.docx` still failed under 0-cookie browser context and extraction found 0 variants.
  - [x] Added generalized publisher API fallback for Elsevier/Wiley/Springer and ran a 13-PMID residual Wiley/Springer batch: 2 usable full texts landed, 1 accepted into refreshed extraction JSON, and corrected no-figure SCN5A recall held steady (see `docs/RECALL_STATUS.md`).
  - [x] Diagnosed `24667783` as a Word-supplement conversion miss, added `textutil` DOC fallback plus scanner/artifact guards, and accepted a forced one-PMID staged refresh that improved no-figure SCN5A row recall (see `docs/RECALL_STATUS.md`).
  - [x] Added SCN5A protein range-deletion scanning/artifact-filter support plus `refresh_run_db.py --replay-model`; recovered the remaining `24667783` `P.K1505_Q1507DEL` row (final no-figure SCN5A row recall in `docs/RECALL_STATUS.md`).

## Active Tasks
- [ ] **Adopt Azure-first / Anthropic-final model routing for the 101-paper
      staging loop.** Routine triage/table routing/extraction/debate should use
      Azure deployments (`gpt-5.4-nano`, `Kimi-K2.6-1`, `grok-4.3`, `gpt-5.4`,
      `DeepSeek-V4-Pro`). Reserve Anthropic for explicit final adjudication
      queues only: Sonnet 5 as the final screen and Opus 4.8 as the hard-case
      arbiter. Measure on the curated staging set before considering full-gene
      refreshes.
- [ ] **Close source/acquisition gaps to >90%** using the highest-yield PMIDs in
      the Exact-Match Recovery Plan above; SCN5A is now the largest remaining
      unique-variant blocker.
  - [ ] Resolve blocked supplement downloads for SCN5A `29325976` (`mmc1.docx`, Cloudflare/redirect blocked even after browser fallback in a 0-cookie session) and similar `missing_variant_supplement` PMIDs.
  - [ ] Resolve residual Wiley Cloudflare/TDM-403 and Springer content-gate failures from the 13-PMID SCN5A residual batch; current no-cookie run landed only `16643399` and `24667783`.
- [ ] **Investigate count semantics and cohort-table over-counting**. Preserve raw count columns and classify study-wide counts versus per-variant carriers before writing patient/affected/unaffected totals.
- [ ] **Create or import KCNE1 per-PMID gold input** before making KCNE1 recall claims.
- [x] **Obtain Elsevier INSTTOKEN** (done 2026-05-21). Token issued by Elsevier Data Support (Jun Bautista) against the `@vanderbilt.edu`-registered API key labeled `GeneVariantFetcher`. Installed into `.env` with user-only file perms; sent only as `X-ELS-Insttoken` header by `harvesting/elsevier_api.py`. Unlock probe across KCNH2/KCNQ1/RYR2/SCN5A/KCNE1 paywalled lists: 242/246 (98.4%) Elsevier candidates now return full text. Bodies saved into each run's `pmc_fulltext/` as `{PMID}_FULL_CONTEXT.md`. Details in `docs/RECALL_STATUS.md`.
- [x] **Re-extract KCNH2/KCNQ1/RYR2/SCN5A with consolidated insttoken full text** (done 2026-05-25). The 242 `_FULL_CONTEXT.md` insttoken bodies were consumed via `scripts/refresh_run_db.py` on the multigene suite, then graft-gated by the 2026-05-26 acceptance-gated source replay. Current scored state is in `docs/RECALL_STATUS.md`.
- [ ] **Integrate valid paywall recovery artifacts into canonical runs** before re-extraction. Keep paper-specific recovery lists in validation artifacts or `docs/RECALL_STATUS.md`, not in this checklist.
- [ ] **Manual PDF download for hard-blocked/stub PMIDs** only after the current-status source audit confirms the artifact is still missing. Feed real downloaded PDFs through `harvesting/format_converters.py`; avoid OneDrive on-demand placeholders.
- [ ] Set weekly recall cadence (Friday re-run + compare) — produces trajectory chart for grant.

## Blocked
- [ ] **Karger access** — Cloudflare interstitial doesn't clear even with 30+ cookies + `--no-headless`. TDM request status unknown.
- [ ] **Sage/Liebert (PMID 23631430)** — Sage's CF instance refuses Playwright fingerprint even with the right cookies.

## Backlog
- [ ] Create a KCNE1 gold test set and keep all non-KCNH2 gold inputs source-reconciled.
- [ ] Re-run extraction with regex disabled to measure regex vs LLM contribution.
- [ ] Expand quality gate test set beyond the 10 audit files.
- [x] Parameterize `scripts/recall_recovery/ingest_clinvar.py` and `ingest_pubtator.py` so cold-start genes can run the same recovery layers KCNH2 uses.

See `docs/RECALL_STATUS.md` for current measured status and `CLAUDE.md` for
recovery architecture and handoff details.
