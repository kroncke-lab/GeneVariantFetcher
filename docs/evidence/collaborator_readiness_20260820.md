# Collaborator extraction readiness audit

Structural QC is not an empirical precision estimate. Publication remains on hold until a source-grounded sample is adjudicated and the Variant Browser imports trusted fields/quarantines.

| Gene | Gate | Papers | Full text | Variants | Nameless | Species-scope links | VF matched | VF held | VF quarantined |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BRCA1 | FAIL | 50/50 | 50/50 | 5348 | 0 | 0 | 2165/5348 (40.5%) | 346 | 324 |
| BRCA2 | FAIL | 45/45 | 44/45 | 1917 | 0 | 0 | 956/1917 (49.9%) | 129 | 137 |
| BMPR2 | FAIL | 50/50 | 45/50 | 554 | 0 | 0 | 300/554 (54.2%) | 56 | 32 |

## Gate failures

- **BRCA1:** mandatory VariantFeatures enrichment/quarantine gate leaves 346 ambiguous live variant identities
- **BRCA2:** mandatory VariantFeatures enrichment/quarantine gate leaves 129 ambiguous live variant identities
- **BMPR2:** mandatory VariantFeatures enrichment/quarantine gate leaves 56 ambiguous live variant identities

## Publication decision

**HOLD.** Passing these invariants is necessary but cannot establish paper-level precision. Keep public publication disabled pending a source-grounded sample and a trusted-field-only Variant Browser import.
