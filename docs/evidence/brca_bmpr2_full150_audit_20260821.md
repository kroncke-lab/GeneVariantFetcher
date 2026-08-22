# BRCA1/BRCA2/BMPR2 exact-150 blind extraction audit — 2026-08-21

## Decision

The complete local candidate is ready for source-grounded review, but it is
**not approved for public publication**. The pipeline processed exactly 50
BRCA1, 50 BRCA2, and 50 BMPR2 papers without using gold values for selection,
extraction, replay matching, or acceptance. Pinned manifests, staged JSON, and
final SQLite paper membership agree exactly for each gene. PMID 19944633, the
canine BRCA2 paper that exposed the earlier harness defect, is absent from every
manifest, staged payload, and final database.

No Variant Browser public annotations were changed. The fail-closed importer
and publisher changes are code-only safeguards for a later reviewed promotion.

## Final local candidates

| Gene | Final database | Papers | Live variants | Paper links | Fact provenance | Penetrance rows | Individuals | VF matched | Quarantined variants | Trusted / held penetrance | Family links / raw values |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| BRCA1 | `results/vb_full150_candidatefixed_20260821/BRCA1/20260821_153608/BRCA1.final6_20260821.db` | 50 | 3,582 | 4,612 | 18,727 | 2,072 | 4 | 1,566 | 187 | 1,885 / 187 | 168 / 147 |
| BRCA2 | `results/vb_full150_candidatefixed_20260821/BRCA2/20260821_153608/BRCA2.final5_20260821.db` | 50 | 722 | 854 | 4,304 | 534 | 58 | 344 | 52 | 383 / 151 | 154 / 142 |
| BMPR2 | `results/vb_full150_candidatefixed_20260821/BMPR2/20260821_153608/BMPR2.final5_20260821.db` | 50 | 260 | 301 | 2,348 | 212 | 97 | 157 | 6 | 205 / 7 | 3 / 3 |

“Family links” counts links whose stored provenance contains a
`family_count` role. “Raw values” is the subset with an associated raw count.
These observations are preserved for review but masked from individual-carrier
totals by the trust gate.

The BRCA1 quarantine consists of 54 out-of-range misparses and 133 wrong-gene
rows. BRCA2 quarantines 52 wrong-gene/residue-mismatch rows; BMPR2 quarantines 6
wrong-gene rows. Quarantined identities are not included in the live-variant
column.

## Blind execution contract

1. Each gene used its exact pinned 50-paper manifest. Acquisition verified
   source identity before accepting a body or cache entry.
2. Extraction and staged replay ran with gold discovery disabled. The replay
   gate uses only source notation, source location, typed count role, and
   compatible variant identity.
3. A replayed count can attach only to one unique compatible identity. Missing,
   ambiguous, or conflicting untyped evidence restores the pre-replay backup.
   BIC digits are never treated as cDNA and codon position alone never aliases
   variants.
4. Table/LLM and SQLite reconciliation share the same ref+position+alt
   frameshift rule. Thus `V340G` may fold into `p.Val340Glyfs*6`, but an
   alternate-residue mismatch, position mismatch, or alt-less codon-only
   frameshift cannot.
5. Mandatory VariantFeatures enrichment and high-confidence quarantine ran
   before a database became a final candidate. Publication additionally
   requires the trusted projection and remains a separate explicit action.

This is a general extraction contract: none of the above steps receives the
expected answer for a paper.

## Validation and adversarial review

The locked, gold-free cardiac validation is
`benchmarks/codex_paper_eval/runs/20260821_candidate_local_gold119`. Gold was
disabled throughout extraction and trace generation; scoring occurred only
after the prediction hash lock. It produced 546 TP, 290 FP, and 87 FN: 86.26%
recall, 65.31% raw precision, 74.34% F1, 98.03% counted-extra precision, and
94.44% count-bearing-only precision. Carrier values were supplied for 184/633
gold rows, with 38.27% nonzero recall and 0.255 conditional MAE. The run records
578 calls and 2.552 million tokens over 119 attempts; the live tier is now 118
attempts / 114 unique PMIDs after unrelated papers were quarantined.

Grok 4.6 at `xhigh` challenged the implementation in several rounds. It found
an ambiguity waiver in replay, insertion-order-dependent migration, unsafe
partial penetrance merges, loss of predicted-protein parentheses, and a live
table path that retained the rich frameshift only in `source_notation`. The
code now fails closed at those boundaries, shares one alias helper, and has
direct tests for both merge directions and the exact live regex-to-merge path.

The requested Claude desktop adversarial pass could not be completed because
macOS was locked when control was retried. No stale Claude answer was reused or
represented as a fresh review; the completed Claude exact-150 staging was
preserved and audited rather than rerun.

## Structural invariants

Across all three final candidates:

- manifest set = staged extraction set = database paper set;
- there are zero placeholder titles, wrong-gene live rows, nameless variants,
  negative counts, impossible unquarantined partitions, or duplicate
  penetrance strata;
- PMID 19944633 is absent;
- BRCA1 PMID 32068069 legitimately has no extracted variants because its cohort
  is BRCA-negative/off-target;
- two incomplete `c.5266_5267ins` mentions in BRCA1 PMID 27376475 remain
  excluded because the source payload omits the inserted sequence.

The remaining duplicate-cDNA groups are explicit conflicting protein claims or
insufficiently proven aliases, not rows safe to collapse: 22 in BRCA1, 4 in
BRCA2, and 3 in BMPR2. None has multiple penetrance strata. The BRCA2 groups are
`c.2808_2811delACAA`, `c.5851_5854delAGTT`, `c.5972C>T`, and
`c.9098_9099insA`; the BMPR2 groups are `c.1021G>A`, `c.633A>G`, and
`c.981T>C`. The BRCA1 groups are `c.1016delA`, `c.1360_1361delAG`,
`c.1934delC`, `c.2302delA`, `c.2389_2390delGA`, `c.2679_2682delGAAA`,
`c.3294delT`, `c.3418_3419insTGACTACT`, `c.3432G>C`, `c.3450delT`,
`c.3756_3759delGTCT`, `c.4065_4068delTCAA`, `c.4120_4121delAG`,
`c.4165_4166delAG`, `c.4228delG`, `c.442_444delCAG`, `c.4484G>T`,
`c.4886_4887delinsC`, `c.5035delC`, `c.5177_5180delGAAA`, `c.5434C>G`, and
`c.981_982delAT`. They are not blindly normalized away.

## Promotion blockers

- Draw and source-adjudicate a precision sample from these 150 papers. The
  structural audit is not an empirical precision estimate.
- Re-review the 111 BRCA2 subjects detached by the corrected exact manifest.
- Run a disjoint, pre-registered validation rather than tuning again on the
  locked cardiac cohort.
- Complete the fresh Claude adversarial pass when the desktop session is
  available.

Until those gates are satisfied, the correct disposition is **local candidate
complete; public publication HOLD**.
