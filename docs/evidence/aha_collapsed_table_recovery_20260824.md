# AHA collapsed-table recovery candidate — 2026-08-24

## Disposition

SHIP the scoped acquisition/parser fix and spend one final sealed gold-118
validation. This evidence authorizes the validation; it does not promote the
candidate or replace the accepted `20260824_postfix_gold118` headline.

## Sealed follow-through

The authorized validation subsequently completed and passed every promotion
gate. The promoted lock is
`benchmarks/codex_paper_eval/runs/20260824_aha_table_sourcebound_gold118`; its
554/283/78 result, count metrics, integrity hashes, and final budget are recorded
in `docs/evidence/gold118_aha_table_lock_20260824.md`. The remainder of this file
preserves the pre-lock authorization evidence.

## Root cause and general fix

AHA's current full-text DOM wraps cohort tables in
`<figure class="table">` and marks rows beyond the preview as `hidden`. The
browser snapshot contains the complete table, but the DOM renderer treated
every `figure` as caption-only. Its independent publisher-aware extraction kept
the prose and also omitted the table, so the length-only picker selected a
source with no patient rows.

The repair is publisher-general at the rendering layer and publisher-scoped at
the selection layer:

- render top-level tables nested in a figure, including hidden rows;
- bind rows/cells to their owning table so nested layout tables are not emitted
  twice;
- make table-completeness preference opt-in and enable it only for AHA;
- retain the original picker behavior for every other publisher.

Synthetic regressions pin collapsed rows, nested-table ownership, and the
opt-in/default split. The focused AHA/browser suite passes 30 tests; after the
canonical TASKS freshness guard was aligned with the already-active goal
heading/progress markers, the complete offline suite passes **2,478 tests**.

## Independent publisher evidence

The official AHA page for PMID 15466642 / DOI
`10.1161/01.CIR.0000144471.98080.CA` contains a patient table with 53 markdown
lines after conversion. The previously selected 19,574-character source had
zero table rows and none of the RYR2 identities. The repaired 22,211-character
DOM source contains all eight unique target-gene identities:

`P164S`, `R414L`, `I419F`, `A2403T`, `F4499C`, `A4510T`, `G4671R`, and
`I4848V`.

A second AHA page, PMID 27114410 / DOI
`10.1161/CIRCGENETICS.115.001370`, is a negative control: its accessible body
has zero table rows and none of the nine gold RYR2 identities. The fix does not
manufacture rows; that paper remains blocked on its access-controlled
supplement.

## Production-path micro replay

The repaired PMID 15466642 source was replayed through the normal production
extractor in an isolated temporary scaffold, with no gold used by selection or
extraction:

- table router: 8 RYR2 variants, no KCNQ1 rows from the same table;
- final extraction: exactly 8 variants, all source notations above;
- claim verification: 8/8 successful;
- off-gene or extra identities: 0;
- calls/tokens: 9 calls, 35,887 tokens (Grok 4.3: 1 call/15,542 tokens;
  GPT-5.6 Sol: 8 calls/20,345 tokens);
- public-list proxy: **$0.20972125**.

The eight identities exactly equal the frozen eight-row gold set for this
paper. If all other accepted predictions were held fixed, the deterministic
paper substitution would move 546/284/86 to **554/284/78**, recall 86.392% to
**87.658%**, raw precision 65.783% to **66.110%**, and F1 74.692% to
**75.374%**. This is a counterfactual diagnostic, not a locked headline. The
paper supplies no person-level counts under the current no-family-inference
contract, so conditional MAE and count coverage are unchanged by this replay.

## Requested external review

- AGY / Gemini 3.1 Pro at its installed maximum `high` first blocked the global
  picker heuristic. After it was made AHA-only and the second AHA control was
  checked, the final verdict was **SHIP**.
- Grok 4.6 at its installed maximum `xhigh` identified the nested-table edge
  during the first pass. After ownership filtering and a regression test were
  added, its final facts-only verdict was **SHIP** and explicitly endorsed one
  sealed validation rather than automatic promotion.

## Budget

Prior attributable spend was $33.06883430. Adding the isolated replay gives
**$33.27855555 used / $66.72144445 remaining** before the sealed run. CLI
Grok/AGY reviewer billing is unavailable and remains excluded, consistent with
the existing ledger.
