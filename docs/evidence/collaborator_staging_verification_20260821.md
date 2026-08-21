# BRCA1 / BRCA2 / BMPR2 trusted staging verification — 2026-08-21

## Outcome

All three exact collaborator manifests were imported into the private Variant
Browser staging/review database through the trusted identity projection. Public
annotations were not published. The raw GVF readiness audit remains FAIL by
design because it evaluates the unprojected databases; the trusted importer
held every ambiguous class listed below.

| Gene | Pinned papers | Staged evidence | Individuals | Table captions | Exact fact sources | Held identities |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| BRCA1 | 50 | 5,583 | 2,826 | 3 | 21,334 | 346 |
| BRCA2 | 45 | 1,918 | 513 | 24 | 3,824 | 129 |
| BMPR2 | 50 | 442 | 395 | 24 | 3,348 | 56 |

For all three snapshots, database verification found exact manifest order,
zero missing or extra papers, nonempty non-placeholder titles, and a nonempty
source location for every staged evidence row. The latest source-run metadata
records `trusted_projection=true` and the held VariantFeatures classes.

## BRCA2 canine/gold fail-closed repair

The corrected BRCA2 GVF database and 45-paper manifest contain no canine PMID
19944633. During staging verification, nine historical canine gold calls were
found still marked current even though their adjudications had been detached
and marked for re-review. The Variant Browser importer added changed subject
keys for fingerprint mismatches but omitted the vanished-evidence branch, so
its later gold invalidation query never saw those records.

The importer now adds every vanished reviewed subject to the invalidation set.
The current gold revision becomes non-current and `disputed`; history is kept
for audit. A regression test recreates a reviewed/gold subject, removes its
evidence upstream, reimports, and asserts the fail-closed result.

The live BRCA2 replay then verified:

- exact 45-paper manifest order;
- zero active evidence for PMID 19944633;
- zero current canine gold records;
- zero current BRCA2 gold records for any of the 111 detached/re-review
  subjects;
- zero current BRCA2 gold records overall; and
- disputed historical revisions retained for audit.

## Release boundary

The staging datasets are ready for authenticated collaborator adjudication, not
public publication. BRCA2's 111 detached subjects need human re-review, and all
three genes still need a source-grounded paper/identity precision sample. The
public `publish_annotations` command was deliberately not run.
