# Collaborator curated-150 staging readiness — 2026-08-24

## Disposition

The BRCA1, BRCA2, and BMPR2 50-paper queues are ready for private collaborator
review in Variant Browser. Each live queue exactly matches its frozen manifest,
contains 50 full-text papers, has resolved bibliographic titles, and uses the
trusted import projection. Public annotation publication remains on hold and
was not run.

This is a reviewer-readiness decision, not a claim that automated extraction is
fully curated. Ambiguous variant identities and untrusted count fields were
withheld so that the remaining uncertainty is visible as manual review work
rather than disguised as accepted evidence.

## Paper curation

The original fixed cohorts were screened against the paper text and PubMed
metadata. Reviews and compilations, wrong-gene interaction papers, somatic-only
studies, explicitly negative target-gene cohorts, off-disease studies, cell-only
studies, and abstract-only sources were excluded. Mixed studies were retained
only when current-study human observations for the target gene were separable.

Twenty-nine papers were replaced: 12 in BRCA1, 8 in BRCA2, and 9 in BMPR2. The
frozen manifests and the complete add/remove ledger are in
`benchmarks/curated_extraction_eval/review_pmids_50_20260824_curated/`.

## Final source and staging checks

| Gene | Exact cohort | Full text | LLM calls / decisions | Live evidence | Individuals | Exact fact sources | Held identities |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| BRCA1 | 50/50 | 50/50 | 598 / 622 | 4,470 | 1,803 | 14,387 | 151 |
| BRCA2 | 50/50 | 50/50 | 439 / 464 | 1,241 | 632 | 3,770 | 21 |
| BMPR2 | 50/50 | 50/50 | 306 / 329 | 456 | 411 | 2,956 | 29 |

All three upstream runs completed with exit code zero, no stage failures or
warnings, no impossible count partitions, no wrong gene symbols, and
write-time-verified trace integrity. They were produced from clean commit
`086e90bbc3136d6a3fddfa76c6f1a9b91181f228`, with gold access disabled.

The raw upstream readiness audit reports placeholder titles and ambiguous live
identities. The Variant Browser trusted importer closes those two gaps before
review: it resolves every title through PubMed and withholds the ambiguous
identity classes. The final live queues have zero unresolved titles. Held
classes are:

| Gene | No notation | Residue offset suspect | Residue unverified |
| --- | ---: | ---: | ---: |
| BRCA1 | 15 | 41 | 95 |
| BRCA2 | 8 | 4 | 9 |
| BMPR2 | 20 | 4 | 5 |

## Source-line sample

Thirty trusted-projection records were checked directly against their local
full-text source, ten per gene. All 30 were supported.

- BRCA1: PMIDs 12872263 (`c.5545T>C`, two cases), 28664449
  (`c.1069A>T`), 29088781 (`c.181T>G`), 30040829 (`c.5266dupC`),
  30606148 (`c.67_75delGAGTGTCCC`), 30675319 (`c.131G>T`, five
  patients/seven families), 31341521 (`c.2368A>G`), 31771539
  (`c.4343delG`), 32879886 (`c.3257del`, three unrelated validation
  cases), and 33403015 (`c.5266dupC`).
- BRCA2: PMIDs 12552570 (`K3326X`), 12655567 (`5146_5149delTATG`),
  14517958 (`p.G2044A`), 22977638 (`c.2403insA`, four unrelated
  patients), 26271414 (`c.8633_8754del`), 26848529 (`c.865A>C`, 13
  carriers), 15876480 (`c.353A>G`), 18431501 (`c.865A>C`), 21232165
  (`c.1151C>T`, one case and zero controls), and 25948282
  (`c.9118-2A>G`, three patients).
- BMPR2: PMIDs 11484688 (`c.354T>G`, 41 carriers with 18 affected and
  23 unaffected), 15146475 (`c.15_19del`), 15358693 (`p.Gln42Arg`),
  18503968 (`p.Arg491Trp`, two cases), 21801371 (`p.Arg491Trp`, three
  cases), 22632830 (`p.Q403*`), 30894412 (`p.Arg491Gln`, 22 carriers
  with 8 affected and 14 unaffected), 32502478 (`p.Arg491Trp`),
  32966279 (`p.Arg491Trp`), and 38473983 (`p.Arg332*`).

## Live Variant Browser verification

The current snapshots share dataset label `collaborator_curated50_20260824`
and status `in_review`. Server-side rendering of each live paper queue returned
HTTP 200, contained all 50 manifest PMIDs and the dataset label, and showed no
metadata warning.

- BRCA1: pair 11, `https://variantbrowser.org/review/11/papers/`
- BRCA2: pair 8, `https://variantbrowser.org/review/8/papers/`
- BMPR2: pair 16, `https://variantbrowser.org/review/16/papers/`

BRCA2 has 111 historical adjudications. They remain preserved, but none is
linked to the fresh snapshot and all 111 are marked for re-review because the
evidence fingerprints changed. BRCA1 and BMPR2 have no prior adjudications.

## Reproducibility

The three run directories are under
`results/vb_curated150_20260824/<GENE>/20260824_181400/`. Their trusted-source
preflights bind the active database, database digest, clean producer commit,
extractor code digest, gold-free status, and verified trace. The manifest
membership test is `tests/unit/test_review_pmid_cohort.py`.

No `publish_annotations` command was run. The public site remains unchanged.
