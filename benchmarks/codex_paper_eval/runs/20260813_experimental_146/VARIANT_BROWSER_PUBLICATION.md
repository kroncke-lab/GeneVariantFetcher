# Variant Browser collaborator-staging refresh

Date: 2026-08-13 (America/Chicago)

Dataset label: `collaborator_reextract_current_system_20260813`

Scope: collaborator-facing Azure `vb-curation` staging only. The public
`publish_annotations` command was not run.

## Inputs and live result

| Gene | Manifest | Manifest SHA-256 | DB SHA-256 | Papers | Evidence before -> after | Individuals | Exact facts | SourceRun |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| BMPR2 | `benchmarks/evaluation_tiers/reviewer_pmids_50_20260811/BMPR2.txt` | `b459f5b97d6caadf39ea25562d5521af74f2c23cbde2002f07dcd6d2cfd0abda` | `1c4da69704bdb49440bfcdcb6163519402d9fe2424650b85e8d746707ce71ac9` | 50 | 528 -> 482 | 470 | 3,871 | 110 |
| BRCA1 | `benchmarks/curated_extraction_eval/review_pmids_50/BRCA1.txt` | `241e80798c787e09fc0382559dffe13f21be007903e1e1be826871940c84327e` | `bf29dc2fddabee146a0701184b62a4e6183d4f18947c2e3cd51804fda7de5c2e` | 50 | 6,299 -> 7,260 | 3,663 | 27,172 | 111 |
| BRCA2 | `benchmarks/curated_extraction_eval/review_pmids_20260811_brca2_provenance/BRCA2.txt` | `32d9dc822106da9a78aeef469f7bbc38428754e35bbf8da3362e11e5e41f5cd0` | `06b4a45015469404188168a3d23df6cd6457f863e3833b55c1ca79fded5c225b` | 46 | 2,735 -> 2,346 | 591 | 4,920 | 112 |

For every gene, live PMID membership and `review_order` match the permanent
manifest exactly. The BRCA2 manifest retains the six documented
non-collaborator provenance exclusions.

## BRCA2 adjudication preservation

Before import, BRCA2 held 111 adjudications (3 confirm, 8 correct-counts,
23 wrong-paper, 77 wrong-variant) and 111 current accepted gold records. Audit
exports were taken before mutation. After import:

- 111/111 adjudication stable keys and 111/111 gold record keys are preserved.
- All 111 adjudications have `needs_re_review=true` and no current evidence FK.
- 98 subjects matched a current subject key but their evidence fingerprint
  changed; their gold records are disputed/non-current.
- 13 subjects no longer exist in the imported evidence. Their historical gold
  revisions remain detached and are excluded by the normal current-snapshot
  export contract.
- Default exportable BRCA2 adjudications: 0. Default exportable BRCA2 gold: 0.

Audit CSV hashes:

| Export | Rows | SHA-256 |
| --- | ---: | --- |
| Adjudications before | 111 | `54ea0e67a445fad7576887dce5e8c68fb35d3d7d6e7da744e3605b1306e43b8b` |
| Gold before | 111 | `b7dc21f183ee54d78cfe36659f429c8a7122c837daf45cdeee0d1b142a494741` |
| Adjudications after | 111 | `8ddd96dd03c28986b2bd58daee1c89f398d7a4051d75c5ff03610d3f85829451` |
| Gold after | 111 | `9d0e019ea74df4712a2df4fa2c4f9421cae5a0d0127ef0827db73e1b97ef1de6` |

The CSVs are retained in `variant_browser_pre_refresh/` beside this report.

## Identity reconciliation

The first BRCA2 attempt aborted before mutation because six legacy identities
would collapse into canonical `gvf:BRCA2:cdna:c.2808_2811delacaa`. Live audit
showed that only `gvf:c.2808_2811delACAA` carried state: one wrong-variant call
and its gold record for PMID 26833046. It was retained unchanged. The other five
identities had zero evidence, staging variants, reviewed validation, feedback,
adjudications, or gold and were moved transactionally to:

- `archive:merged-20260813:gvf:A938fs`
- `archive:merged-20260813:gvf:c=c.2808_2811delACAA|p=p.A938Pfs*21`
- `archive:merged-20260813:gvf:c=c.2808_2811delACAA|p=p.Ala938Profs*21`
- `archive:merged-20260813:gvf:c=c.2808_2811delACAA|p=p.Ala938ProfsX21`
- `archive:merged-20260813:gvf:p.Ala938ProfsX21`

No identity was deleted, so this reconciliation is reversible and the reviewed
subject key remained stable.

## Count-recovery boundary

A carrier-only dry-run used `azure_ai/gpt-5.6-sol` at `high` reasoning on BRCA2
PMIDs 12942367 and 22382806. It found 26 gaps, accepted 0, rejected four proposed
counts whose quotes did not uniquely bind the count to the normalized variant,
and wrote 0 rows. Burn was two calls, 33,136 tokens, 10.79 provider-seconds, and
a $0.188 public-price proxy. Legacy BIC notation was not manually inferred.

## Verification

- BMPR2: 50 papers, exact manifest/order, 482 evidence, 0 adjudications.
- BRCA1: 50 papers, exact manifest/order, 7,260 evidence, 0 adjudications.
- BRCA2: 46 papers, exact manifest/order, 2,346 evidence, 111 preserved calls,
  all 111 requiring re-review, and zero default-exportable stale calls.
- Source run label for all three: `20260813_174842`.
- Dataset label for all three: `collaborator_reextract_current_system_20260813`.
- Public annotation publication: not run.
