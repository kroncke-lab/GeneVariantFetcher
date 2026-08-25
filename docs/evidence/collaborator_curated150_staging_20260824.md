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

## Metric status and preregistered gold

Only three of the 150 papers overlap an approved human-curated fixture, and
BMPR2 has no approved gold rows. Precision, recall, F-score, and MAE are
therefore undefined for this cohort; reviewer readiness is not a quality score.

The exact 50-paper rosters are now frozen into separate 30-paper calibration and
20-paper holdout packets per gene at
`benchmarks/curated_extraction_eval/gold150_preregistered_20260824/`. Source
hashes, exhaustive `NONE` semantics, family-versus-individual count roles, and
the one-score holdout contract are pinned. Human source curation is the next
step; the extraction output may not construct or approve the answer key.

## Reproducibility

The three run directories are under
`results/vb_curated150_20260824/<GENE>/20260824_181400/`. Their trusted-source
preflights bind the active database, database digest, clean producer commit,
extractor code digest, gold-free status, and verified trace. The manifest
membership test is `tests/unit/test_review_pmid_cohort.py`.

No `publish_annotations` command was run. The public site remains unchanged.

---

## Correction — 2026-08-25: BRCA2 alleles were staged in the BRCA1 queue

The readiness disposition above was independently re-verified against live
`vb-curation` on 2026-08-25. Every membership, title, count, adjudication and
render claim reproduced exactly: 50/50 exact manifest match per gene, zero
unresolved titles, `review_order` 1-50 with no nulls, evidence/individual counts
matching the table, BRCA2's 111 adjudications preserved and all flagged for
re-review, all three pages HTTP 200, and no public publication.

One claim did not survive. **The BRCA1 queue contained 152 BRCA2 variants**,
spread across **13 of the 50 BRCA1 papers** — every one a combined
*BRCA1 and BRCA2* study, where a joint variant table was attributed wholesale to
the target gene. All 152 sat beyond `c.5711`, which is BRCA1's ceiling under
both RefSeq and legacy BIC numbering, so no numbering convention explains them.
They included the BRCA2 founder alleles `c.5946delT` and `c.6275_6276delTT`.

This is a floor, not a total: leakage at positions at or below `c.5592` is
indistinguishable from a real BRCA1 coordinate by this test, so the true count is
higher. The 30/30 source-line sample above was not wrong, only underpowered
against a defect of this rate and not targeted at cross-gene attribution.

Two further corrections to the reading of the table above:

- **BRCA2's low hold count (21 vs BRCA1's 151) did not mean BRCA2 was cleaner.**
  BRCA2 was absent from `utils.variant_normalizer.PROTEIN_LENGTHS`, so
  `variant_scanner` fell back to its 9999 default and BRCA2 had no protein bound
  at all. It was unmeasured, not clean. BRCA2 is now registered at 3418.
- **BRCA2's seven out-of-CDS cDNA rows are legacy BIC numbering**, not
  mis-attribution (`c.10462` is BIC for `I3412V`). They are a notation-vintage
  display issue and were left in place.

### Cause and fix

`classify_unmatched` in `scripts/enrich_from_variantfeatures.py` has always
bounded *protein* residues against the variantFeatures-observed `max_pos`
(`misparse_out_of_range`), but its cDNA-only branch returned
`cdna_only_unmatched` with no bound at all — and that class is on the Variant
Browser trusted importer's admit-list (`TRUSTED_UNMATCHED_VF_CLASSES`). A BRCA2
allele written only in `c.` notation therefore entered a BRCA1 run unchallenged.

Fixed at GVF commit `683ad33`: a new `cdna_out_of_range` class whose ceiling is
derived from the same `max_pos` the protein branch uses (three bases per residue
plus the stop codon), so a new gene needs no registration. UTR coordinates are
excluded rather than measured, intronic offsets are judged on their base
position, and `legacy_source_notation` still resolves first. Because the
importer's admit-list is an allowlist, the new class is held fail-closed with no
Variant_Browser change.

### Applied to the live BRCA1 queue

Enrichment was re-derived on the frozen 20260824 BRCA1 run and pair 11 re-imported
under dataset label `collaborator_curated50_20260825_cdna_range_guard`. The
extraction itself was not re-run and did not change; the re-derivation is recorded
in that run's `run_manifest.json` under `provenance.enrichment_rederived`, and the
prior database is kept at `BRCA1.before_cdna_range_guard_20260825.db`.

| | before | after |
| --- | ---: | ---: |
| Papers | 50 | 50 (add 0, remove 0) |
| Evidence | 4,470 | 4,315 |
| Individuals | 1,803 | 1,713 |
| Exact fact sources | 14,387 | 13,933 |
| Held identities | 151 | 312 (`cdna_out_of_range`=161) |
| Out-of-CDS cDNA rows | 164 | **2** |

Exactly 161 identities were reclassified, with every other class unchanged
(matched 1525, `novel_in_range` 1260, `residue_unverified` 95,
`residue_offset_suspect` 41, `legacy_source_notation` 40,
`no_notation_suspect` 15). The two surviving out-of-CDS rows are genuine BRCA1
variants in legacy BIC numbering that carry a protein notation and matched
variantFeatures (`c.5623`=R1835Q, `c.5622`=Arg1835*); they are correctly kept.
About nine of the 161 sit in the same legacy-BIC band and are held anyway — the
intended trade for a queue whose contract is to withhold ambiguous identities
rather than present them as reliable.

BRCA2 (pair 8) and BMPR2 (pair 16) were not re-imported and are unchanged.
`publish_annotations` was still not run; the public site remains unchanged.

### Not a metric statement

These 150 papers are not human-adjudicated, so they are not an evaluation
denominator and no precision, recall or MAE figure moves as a result of this
change. The human-curated cardiac gold-120 remains the only scoring standard.
