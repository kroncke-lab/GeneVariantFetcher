# Precision and cost levers (agent briefing)

Last reviewed: 2026-08-16.

Current-facing ranking of how to raise gold-120 **precision** and cut
**extraction cost**. This is not a second checklist and not a live metrics
table. `TASKS.md` owns what to do next; `docs/RECALL_STATUS.md` owns published
numbers. Dated evidence and paper-level tables live in the diagnostic cited
below.

If another note says “improve precision” without naming a denominator, ignore
it and use §1 here.

## Read with

| What | Where |
| --- | --- |
| Next actions | [`../TASKS.md`](../TASKS.md) gold-120 follow-through |
| Live / locked Gate 2 numbers | [`RECALL_STATUS.md`](RECALL_STATUS.md) |
| Paper-level extras, class split, Sol/Grok/Kimi shares | [`../benchmarks/codex_paper_eval/runs/20260813_gold120_verticalfix/diagnostics/current_gold_matcher_20260815/PRECISION_AND_COST.md`](../benchmarks/codex_paper_eval/runs/20260813_gold120_verticalfix/diagnostics/current_gold_matcher_20260815/PRECISION_AND_COST.md) |
| Affected/unaffected *value* exact-match | [`AFFECTED_UNAFFECTED_PRECISION.md`](AFFECTED_UNAFFECTED_PRECISION.md); diagnostic `PHENOTYPE_VALUE_PRECISION.md` in the same run dir |
| Remaining FNs | same dir, `NOTES.md` and `remaining_fn.tsv` |
| Gold label decisions | [`GOLD_CURATION_QUEUE_2026-08-14.md`](GOLD_CURATION_QUEUE_2026-08-14.md) |
| Measured $9.774 mix | [`PROTOCOL_COST_EVAL.md`](PROTOCOL_COST_EVAL.md) |

The 2026-08-15 diagnostic rescored the locked `20260813_gold120_verticalfix`
predictions against live gold and the current matcher. **No new extraction.**
The locked parent `report.json` was not rewritten and is not a new Gate 2 lock.

## 1. Which precision

There are two identity-precision numbers on that diagnostic:

- **Counted-extra (Gate 2):** 545/(545+8) = **98.55%** after the 2026-08-16
  rules (was 545/(545+10) = 98.20%). This is the
  repo-intended precision question: did we attach a patient number to a
  variant gold does not have? An extra counts against this number only if
  any of `carriers`, `affected`, `unaffected` is non-null. Matched identities
  with blank counts still help. Wrong counts on a *matched* row do not move
  this number (that is MAE).
- **Raw identity:** 545/1334 = **40.85%**. Every unmatched prediction is an
  FP, including count-null catalogue rows.

Do not collapse them. Dropping gold-incomplete extras raises raw precision
and hides true extractions; it leaves Gate 2 unchanged. The locked Gate 2
pass is still 95.70% (534/(534+24)) against a 77.3% floor.

Counted-extra is **not** precision of affected/unaffected values. Count
quality is the separate MAE / count-recall block (diagnostic carrier MAE
0.2875, supplied 160/633). Exact-match precision of emitted affected /
unaffected integers is a third denominator — see
[`AFFECTED_UNAFFECTED_PRECISION.md`](AFFECTED_UNAFFECTED_PRECISION.md).

## 2. What actually remains

- **Gate 2 FP budget is 8 rows**, down from 10. The two code-fixable ones
  are closed:
  - `c.693delCA` vs gold `c.692_693delCA` — a deletion-span bridge now treats
    the one-coordinate spelling as identical when the deleted bases match and
    the lone coordinate *is* a span endpoint. Both scorers, and twin identity
    too. Ambiguity refused; single-base deletions never bridge.
  - KCNH2 18808722 `L187P` (KCNQ1 gold on the same PMID) — the scanner's
    ±50-character conflict window could not see `…of KCNQ1, which results in
    … (L187P)`, so a document-level attribution pass now re-reads every
    occurrence with a ±240-character window. Measured across all 120 papers:
    **71 cross-gene tokens removed on 12 papers, 0 gold rows lost.** This
    cannot show up on a locked prediction set — it needs the paid re-extract.
  The remaining 8 are curator work, not code: two synonymous rows with counts
  (gold *does* carry synonymous variants with real counts — KCNQ1 `A344A` has
  23 carriers on 30758498 — so a blanket synonymous drop is wrong), four
  identity TPs against gold `0/0/0`, and the documented
  `GOLD_GAP_REAL_VARIANT` items. Full list in `PRECISION_AND_COST.md` §2a.
- **789 raw FPs.** Exclusive split: 70 other-gene gold (the scanner pass
  above is aimed at exactly this class, and lands on the next extraction), 34 on KCNH2 26746457
  eTable (24 already called correct, all count-null), 10 synonymous-shape, 1
  unmerged twin, 9 counted extras not in those rows, **667 remainder**.
- **584 of 789 extras sit on papers that already matched every gold row.**
  Largest: 19632626 (67), 32533946 (52), 29350269 (49), 33536282 (44),
  26746457 (34). Do not delete these to “fix” 40.82%.
- **88 FNs** on 22 papers (`remaining_fn.tsv`). Classes: 13 gold-side
  (editorial / wrong PMID / compound), 20 TIFF tables (29650123), 17
  caption-stub, 11 already recovered on KCNQ1 21956039 in the 4-paper check
  but not in the lock, 27 residual extraction.

## 3. Cost mix (do not shop models)

Locked gold-120 proxy **$9.774** / 527 calls
([`PROTOCOL_COST_EVAL.md`](PROTOCOL_COST_EVAL.md)):

| Share | Calls | Proxy | What to cut |
| --- | ---: | ---: | --- |
| GPT-5.6 Sol verification/vision | 393 (388) | $7.019 (72%) | Dominant. Locked traces: 244 figure-vision + 148 claim-verify + 1 adjudication. |
| Grok 4.3 extraction | 115 (114) | $2.587 (26%) | Already ~one extract per paper. |
| Kimi K2.6 routing | 19 (19) | $0.168 (2%) | Not a project. |

A successful deterministic table short-circuit *saves* a Grok call. The
2026-08-15 short-circuit *refuse* (census saw count columns, parse emitted
none) *adds* a Grok call. Do not advertise that refuse as a cost win.

Do not un-park the per-paper final check or enable `COUNT_RECOVERY_ENABLED`
as a cost/precision lever. Those are already parked / default-off for
measured reasons.

## 4. Ranked work (agents start at the top)

### $0 — do not spend the ~$10 re-extract for these

1. Curate the **10 counted extras** (only list that can still move 98.20%).
2. Do **not** delete the 780 identity-only extras to raise raw precision.
3. If raw precision is the target, **expand gold** on gold-exhausted
   catalogues, starting with 26746457 (24 source-confirmed).
4. Scorer only, generic: `c.693delCA` ↔ `c.692_693delCA`; RYR2
   `c.169-198_273+820del` ↔ EXON 3 (cds 169–273). No per-variant or
   per-gene aliases (`ΔKPQ`, hard-coded exon-3 maps). Structural matching
   already lives in `utils/structural_alleles.py` and both scorers use it;
   the 2026-08-15 rescore added **zero** structural TPs because those
   identities were not in the locked predictions (19216760 never extracted
   the exon event; 28087622 never extracted ΔKPQ).
5. Trace-only Sol audit: which of the 244 figure-vision calls added an
   identity the text extract already had.

### Needs extraction or acquisition — recall / MAE, not Gate 2 precision

6. Paid gold-120 re-extract (~$10 / ~40 min, live 119-attempt
   `tier2_gold_120.tsv`) realizes extraction-code already landed (KCNQ1
   21956039 11/11 on the 4-paper check). It is **not** required to rank or
   improve counted-extra precision. Next paid arm must use the 119-attempt
   manifest (KCNH2 10086972 removed).
7. Source: 29650123 TIFF tables in on-disk `mmc1.docx`; caption-stub bodies
   (19926015, 14678125, 31520628, 24667783).
8. `regex_table` **column binder** (one table IR, counts only from a bound
   column). The 26496715 class: 114 gold assertions in a table the pipeline
   read and flagged. That paper has **0 counted extras** today — binder is
   MAE / count-recall, not Gate 2. Large-table short-circuit refuse already
   landed; binder itself is not done. After the binder, not before: Luna-max
   post-layer count-ambiguity cards.

## 5. House constraints other agents keep hitting

- One checkout, `main`. No feature branches or worktrees unless Brett
  changes `CLAUDE.md`.
- Do not invent a review target on a clean tree. Do not rewrite locked
  runs or append-only history.
- Do not change gold without source verification. Do not collapse the
  19216760 exon-3 family pair (two families, 4/4/0 and 2/2/0) or edit gold
  counts without a curator.
- Do not add one-off gene-specific matcher aliases.
- Corpus writes go through `corpus/` (external volume) /
  `require_external_storage()`. Do not replace a broken symlink with a
  local directory.
- Offline tests: `.venv/bin/python -m pytest tests/unit -q`.

## 6. Already landed this iteration (do not redo)

Source-verified gold erratum (PMID remap 10086972→10086971, four
transcription renames, two true-duplicate drops). KCNH2 10086972 dropped
from live `gold_120`. Structural alleles module + NCBI exon maps. Large-table
short-circuit refuse. Eval matcher bridges (arrows, stop/fs, splice,
translation, twin-merge) — still eval-only in `run_eval.py`, not all ported
to `cli/compare_variants.py`. Always-on phenotype copy guard
(`pipeline/phenotype_count_guard.py`) — family/figure `affected == carriers`
with no counted split, plus (2026-08-16) two field-scoped refusals of a
non-closing partition and an unsourced `affected = 0`; see
[`AFFECTED_UNAFFECTED_PRECISION.md`](AFFECTED_UNAFFECTED_PRECISION.md).
cDNA deletion-span bridge in both scorers and in twin identity.
Document-level gene attribution in `utils/variant_scanner.py`.
`scripts/phenotype_value_precision.py` — the $0 rescore harness; run it before
proposing any new guard, and quote `kills / destroys`, not a bare percentage.
