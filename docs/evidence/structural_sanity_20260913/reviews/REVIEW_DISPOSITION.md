# Agy and Grok structural sanity review — 2026-09-13

Both requested CLIs completed independent headless consultations from the same
necessary aggregate scientific brief. They received no source code, local
paths, secrets or row-level research data. They were instructed to use no tools
or outside sources. The reviews are opinions based on supplied context;
local source inspection and numerical checks are identified separately below.

| Reviewer | Requested model | Returned identity | Elapsed | Reported usage | Reported cost |
|---|---|---|---:|---:|---:|
| Agy CLI | `gemini-3.1-pro-high`, high | Model not echoed; `SUCCESS`, one turn | 32.45 s | 18,661 total tokens | Unavailable |
| Grok CLI | `grok-4.6`, high | `grok-4.6-build`, one turn, `end_turn` | 78.07 s | 18,325 total tokens | $0.01817572 |

`agy_review.md` and `grok_review.md` preserve the full responses.
`agy_prompt.txt` and `grok_prompt.txt` preserve exact matching input text;
`shared_prompt.txt` is the common source. Structured raw responses, stderr,
execution status, timing, usage and CLI executable hashes remain beside them.
The Agy dollar charge is unknown, not zero. CLI costs are estimates, not
reconciled invoices. `run_cli_reviews.py` records the bounded launch method.

## Accepted conclusions and immediate checks

1. **BRCA2 polymer absence is an input omission.** The frozen manifest
   explicitly provides only the two experimental peptide frames, with 72
   COM-supported canonical positions and no IDR fallback. Its zero-polymer
   result does not mean that BRCA2 lacks disorder. Missing experimental
   coordinates alone do not justify marking an interval IDR.
2. **Use a provenance-backed canonical polymer layer.** Preserve real
   experimental frames and actual biological-copy identities. Add only
   explicitly supported contiguous IDR intervals, with canonical numbering
   and one sequence-polymer context per molecular sequence, not one repeated
   polymer context per PDB peptide copy. The formula remains
   `3.8 * sqrt(abs(i-j))` inside the same segment/chain. Distinct frames,
   ordered/IDR pairs and different IDR segments retain unavailable distances.
   Unsupported targets remain missing rather than density zero.
3. **Resolve state overlap explicitly.** If a bound experimental peptide
   overlaps a region supported as disordered in its unbound state, preserve
   named bound-3D and free-polymer alternatives or an explicit disjoint
   primary-state mask. Do not automatically turn equal averaging of assembly
   chain contexts into an equal mixture of alternative physical states.
   Test both canonical interval uniqueness and context duplication invariance.
4. **GCK's broad elevation follows the prior change.** The saved replay
   holds spatial weights exactly fixed, so its lift from median density
   18.24% to 37.41% is caused by the changed donor posterior values. The
   37.04% missense prior and 27.59% posterior for `A=0,U=1` are correct
   under the requested formula. This is not a result that every residue or
   every substitution causes MODY. Report the local deviation from the
   prior and the contribution of neighbor counts alongside the raw density.
5. **Endpoint interpretation and validation remain limited.** Pooled GCK
   affected counts contain hyperglycemia/MODY and activating/hypoglycemia
   evidence. A source-audited endpoint/mechanism ledger is the useful next
   correction; it must preserve uncertain or incompatible counts rather
   than assigning mechanism from a score. Tiny GCK MAE differences and the
   selected BRCA2 subset do not establish a robust incremental structural
   benefit or independent disease-outcome calibration. Keep sequence
   neighborhoods and fitted intercepts as comparators.

## Reviewer claims corrected or not adopted

- Both reviewers call the `A=0,U=1` mean a **floor**. It is the mean for that
  exact count pattern. It decreases further with additional unaffected
  carriers; there is no universal 27.6% floor.
- Agy calls strength 2.918755 **massive/highly rigid**. That overstates its
  absolute size. It contributes about 2.9 pseudo-observations and is material
  relative to a singleton, while larger observed counts can dominate it.
- Agy states that the small MAE gain proves a **genuine local signal**.
  The supplied internal comparison does not establish that conclusion;
  the sequence comparator is slightly better and the original clinical
  target set barely changes.
- Grok correctly distinguishes 246 population-only singletons from all
  singletons. An independent local recount confirms **332/634 total
  singletons (52.366%)**, comprising 246 population-only and 86 other-source
  units; see `local_gck_singleton_check.json`. There is no denominator error
  in the reported 52.4% figure.
- Grok's wording about residue collapsing inflating support is unsupported
  as a general directional claim. Variant-only exclusion and retention of
  distinct same-residue variants remain the user's fixed analysis question.
- Genuine experimentally observed biological-copy geometry is retained.
  The warning against inventing partner-filament or cross-fragment geometry
  must not be interpreted as deleting validated same-assembly distances.
- Neither reviewers' unrestricted count permutation nor new held-out
  outcome acquisition is an automatic acceptance gate. A future null
  diagnostic should preserve complete count-unit pairs and relevant source/
  count strata, with fixed predictions/evaluation semantics specified first.
  Current posterior-label metrics must continue to be described as internal
  diagnostics with fixed gene-by-class hyperparameters.

## Locally verified source and arithmetic

`local_source_receipt.json` pins the frozen extension README/runner, BRCA2
geometry documentation/summary, class-specific GCK reports and prior table.
`local_gck_singleton_check.json` recounts the saved 634-unit structural input.
The historical variance divisor is `M`; because `sum(w) < M`, this makes
variance smaller and strength larger than the normalized-weight sensitivity.
Alpha adds affected observations and beta adds all unaffected observations.
These rules, the missense-specific prior, distance-only positive-tail kernel,
same-segment polymer law and exact-variant exclusion remain unchanged.

No analysis, core pipeline or frozen evidence file was modified by the
consultations. BRCA2 geometry extension and any diagnostic additions require
their own separately inspectable outputs; agreement between reviewers is not
a substitute for the actual mapping and numerical checks.
