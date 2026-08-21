# SCN5A Q1077 isoform identity diagnostic — 2026-08-21

## Decision

The constrained SCN5A lookup correction is suitable for trusted staging, but
the post-lock score is an in-sample diagnostic, not a new blind lock and not a
public precision claim. Public publication remains on hold.

## Root cause

The local VariantFeatures reference represents the common 2,015-aa SCN5A
Q1077del isoform. Much clinical literature reports against the 2,016-aa Q1077
isoform, so after Q1077 the same biological residue is commonly numbered one
higher in the paper. This is documented in source literature describing the
two isoforms and the mutation-study convention:

- [S1787N/Q1077del background study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4414567/)
- [SCN5A Q1077 and Q1077del isoforms](https://pmc.ncbi.nlm.nih.gov/articles/PMC6109095/)
- [SCN5A transcript/isoform comparison](https://pmc.ncbi.nlm.nih.gov/articles/PMC6636349/)

The pre-fix trusted projection held 53 true positives and 50 false positives as
`residue_offset_suspect`; 50 of those 53 true positives were SCN5A. A separate
parser defect treated bare uppercase proline in `P1730H` and `Pro871Leu` as an
HGVS `p` prefix and removed it.

## Fail-closed implementation

The paper's emitted notation and identity are never rewritten. The offset is a
VariantFeatures lookup key only, attempted after exact protein and exact cDNA
matching fail. It is allowed only when all of these conditions hold:

- gene is exactly SCN5A and the paper residue is 1078–2016;
- VariantFeatures proves the local reference is Q1077del: maximum real residue
  2015, optional stop row at 2016, and glutamate at 1077;
- the paper reference amino acid matches VariantFeatures at `N-1` and does not
  match at `N`;
- the call is a simple missense or stop call.

Structural alleles, generic offsets, unverified residues, cDNA/genomic calls,
other genes, and an unproven/long reference remain held. Exact same-coordinate
matches always win. Regression tests cover boundaries, exact precedence,
wrong-gene and long-reference cases, novel alternate residues, structural
alleles, and end-to-end preservation of the literature identity.

## Frozen evaluation and post-lock diagnostic

All rows use the same 119 attempts / 115 unique PMIDs and 633 gold identities.
Only the first row is the primary immutable lock.

| Projection | TP | FP | FN | Recall | Raw identity precision | Counted-extra precision proxy | Carrier supplied | Carrier MAE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Locked trusted primary | 489 | 438 | 144 | 77.25% | 52.75% | 97.80% | 178/633 | 0.247 |
| Post-lock trusted diagnostic | 546 | 452 | 87 | 86.26% | 54.71% | 98.03% | 188/633 | 0.250 |
| Post-lock raw diagnostic | 551 | 578 | 82 | 87.05% | 48.80% | 98.04% | 202/633 | 0.535 |

The lock-to-diagnostic change is a mixed correction bundle, not an isoform-only
A/B: +57 TP and +14 FP, or 80.3% precision on the net added identity matches.
That must not be confused with the scorer's 98.03% counted-extra proxy, whose
denominator answers a different question. The post-lock trusted projection
removes 126 raw false positives at the cost of 5 raw true positives.

The five remaining held true positives are KCNQ1 `A340E` and `S227del`, RYR2
`K4481R` and `R4037C`, and structural SCN5A `F1617del`. They are evidence for
targeted adjudication, not permission to generalize a ±1 rule or trust a
structural offset.

## Integrity

The original `LOCK.json` hash is unchanged. Diagnostic outputs live under
`benchmarks/codex_paper_eval/runs/20260821_current_changes_gold119/diagnostics/scn5a_q1077_fix/`.
A disjoint, pre-registered validation set is still required before promoting
the diagnostic result to a new acceptance lock or public quality claim.
