# Current Recall Status - 2026-05-20

This note records the post-`19ae63f` state of the recall/generalization push so
the next run can resume from the same baseline without relying on chat history.

## Current Git State

- Current pushed head: `19ae63f Improve deterministic variant recovery`
- Branch: `main`
- Remote: `origin/main`
- Latest committed changes were intentionally scoped to deterministic extraction,
  acquisition recovery, comparison/normalization, migration safety, and tests.
- There are unrelated local dirty files in this checkout. Do not assume they are
  part of `19ae63f` unless they are explicitly staged/committed later.

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

The dominant blocker is still acquisition of the right source material, especially
full text and supplements from paywalled publisher pages. The Elsevier
institution token is likely high leverage:

- `ELSEVIER_API_KEY` is for API access.
- `ELSEVIER_INSTTOKEN` is the institutional entitlement token that should unlock
  subscription full text/supplements where Vanderbilt has access.
- After setting `ELSEVIER_INSTTOKEN`, rerun targeted paywall/supplement recovery
  first, then re-extract and rescore from the recovered `FULL_CONTEXT.md` files.

This is particularly important for ScienceDirect/Elsevier papers where we saw
paywall stubs, Cloudflare/403 behavior, or accepted manuscripts missing the
supplement.

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

## Next Run Plan

1. Set `ELSEVIER_API_KEY` and `ELSEVIER_INSTTOKEN`.
2. Run a targeted Elsevier recovery pass for the highest-yield missing
   ScienceDirect/Elsevier PMIDs, starting with SCN5A `29325976` and then the
   next highest-yield SCN5A/KCNH2 entries.
3. Confirm the recovered files are real article/supplement content, not paywall
   stubs. Use content-quality checks and inspect variant/table density.
4. Re-extract only the recovered PMIDs first.
5. Merge into a copied DB using preserve/targeted migration mode.
6. Rescore against the current baseline.
7. Promote only general parser/acquisition fixes to code. Keep per-PMID recovery
   artifacts in validation runs, not production branches.
8. Add no-gold QC output before using the pipeline for new gene-disease runs.
