# Variant Browser collaborator-staging refresh
Date: 2026-08-14 (America/Chicago)
Dataset label: `collaborator_reextract_current_system_20260814`
Scope: collaborator-facing Azure `vb-curation` staging only. The public
`publish_annotations` command was not run.

## What changed
Re-extraction of the three fixed collaborator queues with the improved
matcher/extraction merged in GVF PR #182 (notation-twin merge, legacy-notation
matcher, circuit-breaker and gene-less mutation-table truncation fixes,
linkage codon-shadow projection gate). Same pinned manifests, disease
contexts, and calibrated `--pmid-file --no-source-recovery` protocol as the
2026-08-13 baseline refresh.

## Inputs and live result
| Gene | Pinned manifest | Papers | Evidence | Individuals | Exact fact sources |
| --- | --- | ---: | ---: | ---: | ---: |
| BMPR2 | `benchmarks/evaluation_tiers/reviewer_pmids_50_20260811/BMPR2.txt` | 50 | 518 (was 482) | 452 | 3,798 |
| BRCA1 | `benchmarks/curated_extraction_eval/review_pmids_50/BRCA1.txt` | 50 | 7,430 (was 7,260) | 3,668 | 26,915 |
| BRCA2 | `benchmarks/curated_extraction_eval/review_pmids_20260811_brca2_provenance/BRCA2.txt` | 46 | 2,353 (was 2,346) | 604 | 4,956 |

Preflight for every gene: papers N → N, add 0, remove 0 — live membership and
review order match the permanent manifests exactly.

## Adjudication preservation
All 111 historical BRCA2 adjudications were preserved and are marked
`needs_re_review` (adjudications-relinked=0, adjudications-needing-rereview=111).
BMPR2 and BRCA1 carried no attached adjudications at import time.

## Identity reconciliation
The first BRCA2 attempt aborted before mutation: canonical
`gvf:BRCA2:label:p.e1107*` would merge two legacy identities. A live audit
showed `gvf:p.E1107*` carries state (1 gold record, 1 evidence row) and
`gvf:E1107*` carries none (0 feedback / 0 evidence / 0 adjudications / 0 gold).
The state-less key was moved transactionally to
`archive:merged-20260814:gvf:E1107*`; `gvf:p.E1107*` was retained unchanged.
No identity was deleted, so the reconciliation is reversible.

## Operational note
The very first publish attempt failed with an ODBC `HYT00` login timeout —
the Azure SQL staging instance was resuming from auto-pause; an immediate
retry succeeded. Source runs: `production_runs/<GENE>/20260814_1238*` (logs
`bmpr2.log`, `brca1.log`, `brca2.log`, publish logs `vb_publish*.log`).
