# Grok 4.6 xhigh adversarial review — SCN5A identity and collaborator staging

Grok 4.6 was run at `xhigh` effort on 2026-08-21 with tools and web search
disabled. It received the implemented invariants, regression coverage, frozen
119-paper arithmetic, and live BRCA2 staging invariants rather than permission
to edit the repositories.

## Verdict

- **GO:** merge and use the constrained changes for trusted staging.
- **HOLD:** public publication without more human adjudication.

The GO is conditional on preserving the literature identity at every export,
keeping the SCN5A offset lookup-only and missense/stop-only, leaving structural
and generic offsets held, retaining ambiguous BRCA1/2/BMPR2 identities outside
staging, and citing the original 489/438/144 lock as the only frozen result.

## Highest-risk challenges

1. The 98.03% counted-extra metric is not the precision of new matches. The
   mixed lock-to-diagnostic delta is 57 TP and 14 FP, or 80.3%.
2. The post-lock diagnostic combines the SCN5A lookup and proline-parser fixes;
   it cannot attribute the entire recall gain to one change or become a new
   blind lock.
3. Detaching/disputing BRCA2 gold is correct safety behavior but cannot be
   counted as extraction improvement. The 111 subjects require re-review.
4. No downstream Variant Browser/export path may replace the literature residue
   or cDNA with the VariantFeatures lookup coordinate.
5. The remaining KCNQ1, RYR2, and structural SCN5A true positives must not
   motivate a generic ±1 rule.
6. Overall trusted identity precision in the diagnostic is 54.7%; structural
   QC and the counted-extra proxy do not establish public paper-level quality.
7. A disjoint held-out SCN5A isoform set and source-grounded collaborator
   adjudication remain necessary.

This review influenced the release boundary: code and trusted staging may
advance, while public annotations remain unchanged.
