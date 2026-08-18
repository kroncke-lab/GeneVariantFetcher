# Cross-gene contamination in the collaborator-facing staged datasets

**Dated evidence, 2026-08-18.** Measurement only — this file states what is
wrong with already-published data and how much. The mechanism and the fixes are
in `docs/PROTOCOL_CHANGELOG.md` (2026-08-18, four rows); the remaining work is
`TASKS.md` §4c.

## What was wrong

A variant's gene was never grounded in the paper: it was stamped with the run's
target gene. A paper reporting two genes therefore had one gene's tables emitted
under the other. The sharpest demonstration is PMID 21232165, where the staged
BRCA1 and BRCA2 extractions are identical except the gene column — 67 variants
each, same tables, same order.

Fixed by `cc86a5e` and `eb6666b`.

## How much, per dataset

Measured by replaying the fixed source-grounded gene filter over the stored
extractions. Deterministic, no model calls. Per-paper counts and example
variants: `crossgene_contamination_20260818.json`.

| Dataset | Papers | Variants | Wrong-gene | Papers affected |
|---|---|---|---|---|
| BRCA1, staged 2026-08-14 | 50 | 7,346 | **259 (3.5%)** | 15 |
| BRCA2, staged 2026-08-14 | 46 | 1,426 | **113 (7.9%)** | 6 |
| BMPR2, staged 2026-08-07 | 758 | 5,838 | **0 (0.0%)** | 0 |

BMPR2 is single-gene PAH literature with no cross-gene tables, so it is clean
and needs no isolation. It is being re-extracted anyway as a **control**, to
show empirically that the fix is inert on uncontaminated data rather than
inferring it.

Worst affected: BRCA1 24549055 (49/217), BRCA1 27741520 (45/135), **BRCA2
26848529 (44/124)**, BRCA1 21232165 (37/67).

**PMID 26848529 is one of only two lead-approved BRCA2 collaborator papers**, so
this reaches the BRCA2 gold provenance, not only the staged evidence.

## Downstream exposure

Canonical scoring DBs carry rows for the affected PMIDs. These counts are total
rows on those papers, not wrong rows — they bound what a rescore would touch:

| DB | Rows on affected PMIDs | Papers | DB total |
|---|---|---|---|
| `validation_runs/canonical_baseline/BRCA1.db` | 1,783 | 12 | 89,248 |
| `validation_runs/canonical_baseline/BRCA2.db` | 479 | 6 | 6,840 |

`BRCA1.db` also stores BRCA2-only `c.5291C>G` for PMID 21232165. The
spaced-cDNA fix stops new occurrences but does not retroactively clean stored
rows.

## Isolation status

- The two contaminated run directories carry a local `QUARANTINE.md`
  (`benchmarks/codex_paper_eval/runs/20260814_experimental_146_rerun/production_runs/{BRCA1,BRCA2}/`).
  Those paths are gitignored, hence this tracked companion.
- **Nothing was deleted.** The source under `corpus/` is unaffected and remains
  valid — only the gene attribution on extracted variants is wrong — and the
  standing rule is that hard-won full-text extracts are never removed.
- **Variant Browser staging has NOT been touched.** The decision is to review
  the old-vs-new diff before isolating or replacing anything there, so the wrong
  calls are still live for collaborators as of this writing.
