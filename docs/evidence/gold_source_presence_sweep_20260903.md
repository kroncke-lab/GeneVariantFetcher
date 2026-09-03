# Gold-row source presence sweep (2026-09-03)

**Question.** How much of the named-variant gold is unreachable by *any*
reading protocol because the variant is not present in the source we hold on
disk — article text, converted supplements, or figure images?

**Method.** `scripts/recall_audit/gold_source_presence_sweep.py` classifies
every gold row of the mixed-gold inventory (`benchmarks/evaluation_tiers/
mixed_gold/inventory.tsv`, 1,533 non-quarantined gene-paper attempts, 7,316
gold rows) against everything under `corpus/<GENE>/<PMID>/`: the article
`FULL_CONTEXT`/`CLEANED` body, every supplement after conversion (`.docx`,
`.xlsx`, `.xls`, `.doc`, `.pptx`, `.zip` bundles through the production
`FormatConverter`; 626 files converted, 17 failed), and top-level PDFs via
`pdftotext`. Classification is blind to predictions and reads only the gold
variant string and the on-disk bytes. Nothing under `corpus/` was written.
Machine-readable outputs and the per-row table are in
[`gold_source_presence_sweep_20260903/`](gold_source_presence_sweep_20260903/).

The probe is deliberately wider than `tier_source_reachability.py`. Two blind
spots found and closed while building it:

- **Nonsense spellings.** `VariantNormalizer.get_all_forms("R222X")` yields
  `p.Arg222*` but papers write `Arg222Ter`/`Arg222Stop`; the scorer's
  `matches()` already bridges those, so the probe now does too (SCN5A 27232914
  was wrongly "absent").
- **Garbled PDF text.** KCNH2 15176425's 225 KB body is `pdftotext` output of
  a font-subset PDF: `(cid:NN)` glyph codes with no readable running text and
  zero variant tokens. A byte-count threshold calls that "full text". Rows in
  such bodies are now `text_absent_garbled_body`, an acquisition class.

## Classes

| class | meaning | owner | in recall denominator? |
|---|---|---|---|
| `present_in_body` | variant string in article body | reading protocol | yes |
| `present_in_supplement_only` | only in a converted supplement file | fold / source selection | yes |
| `present_in_pdf_only` | only in the article PDF text | source selection | yes |
| `source_absent` | nothing searchable on disk | acquisition | may exclude (hard ceiling) |
| `text_absent_stub_body` | < 6,000 searchable chars total (abstract/landing page/caption stubs) | acquisition | may exclude (hard ceiling) |
| `text_absent_garbled_body` | body dominated by unmapped glyph codes | acquisition (re-render/OCR) | may exclude (hard ceiling) |
| `text_absent_substitution` | plain substitution absent from every searchable byte, no figures | acquisition or gold error | may exclude (hard ceiling) |
| `text_absent_figures_present` | absent from text, figure images on disk | unknown (probe blind spot) | stays penalized |
| `text_absent_notation_inconclusive` | absent from text, indel/frameshift/structural notation | unknown (probe blind spot) | stays penalized |

"Hard ceiling" = the four acquisition classes. "Wide ceiling" adds the two
unknown classes and is an upper bound, not a claim.

## Results

All 1,533 attempts (including the 111 the inventory already marks
`source_unavailable`):

| stratum | gold rows | hard ceiling | wide ceiling | present in body |
|---|---:|---:|---:|---:|
| everything | 7,316 | 1,321 (18.1%) | 2,236 (30.6%) | 5,071 (69.3%) |
| the 1,422 runnable attempts | 7,107 | 1,122 (15.8%) | 2,037 (28.7%) | 5,061 (71.2%) |
| runnable, human-curated cardiac only | 5,772 | 1,037 (18.0%) | 1,821 (31.6%) | — |

So on the tranche suite the protocol is actually run on, **at most ~84% of
gold rows are reachable by any text-reading protocol** (hard ceiling), and the
honest range is 71–84% once the two undecidable classes are counted. Of the
6,442 gold rows with a non-zero carrier count, 1,134 sit behind the hard
ceiling and 2,005 behind the wide one: the count-coverage gap has the same
acquisition floor as identity recall.

By class on the runnable cohort: `present_in_body` 5,061;
`text_absent_figures_present` 582; `text_absent_substitution` 512;
`source_absent` 465; `text_absent_notation_inconclusive` 333;
`text_absent_stub_body` 122; `text_absent_garbled_body` 23;
`present_in_supplement_only` 9; `present_in_pdf_only` 0.

Three readings matter for the campaign:

1. **Binary supplements are not where the missing gold is.** Only 9 rows
   live in a converted supplement and not the body — the fold already lands
   supplement text. The 512 absent substitutions are absent from the
   supplements too.
2. **The absent-substitution rows are concentrated.** 512 rows over 112
   papers; 14 papers with ≥ 10 absent substitutions hold 275 of them. The
   largest are table-body acquisition failures where a `browser-html-generic`
   capture kept the prose and dropped every table: KCNQ1 14678125 (37 of 41
   rows absent, 13.8 KB body, zero table captions), KCNQ1 17192539 (51/57),
   SCN5A 21273195 (22/34), SCN5A 24631775 (25/27), SCN5A 11901046 (23/23).
   These are re-fetch targets (Wiley/Elsevier TDM), not extraction work.
3. **The wide-ceiling tail is figure-heavy.** 582 rows sit on papers with
   figure images and no text hit; SCN5A 20129283 alone contributes 76 (of 417
   gold rows) and KCNQ1 32893267 43. Prior à-la-carte figure-vision
   experiments recovered nothing; treat this class as a bound.

Per tranche (runnable rows / hard / wide / present): `mixed_gold_01` 242 /
75 / 93 / 149; `mixed_gold_02` 303 / 24 / 71 / 232; `mixed_gold_03` 399 / 39 /
94 / 305; `mixed_gold_04` 150 / 65 / 79 / 71. Tranche 1 is unusually
acquisition-bound (31% hard ceiling) because it holds KCNQ1 14678125 and all
three gene arms of the Finnish founder paper 15176425 (KCNH2 arm garbled,
SCN5A arm a 3 KB stub, KCNQ1 arm absent). Expect its baseline recall to look
worse than the accepted gold-118 lane for that reason alone; the paired delta
is unaffected.

## How to use these numbers

- Report tranche recall on **all gold** as the primary figure, and the
  reachable-gold figure (hard ceiling excluded) as a labeled diagnostic with
  the excluded row count. Never drop the two unknown classes from a
  denominator.
- The hard-ceiling paper list (`summary.md` § "Papers with the most
  text-absent gold rows") is the acquisition worklist. It is independent of
  any protocol change and should be worked through source recovery, which the
  frozen-source tranche runs deliberately do not do.
- Re-run the sweep after any corpus refresh; with a warm conversion cache it
  takes under a minute.
