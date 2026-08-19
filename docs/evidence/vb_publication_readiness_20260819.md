# Variant Browser publication readiness — BRCA1, BRCA2, BMPR2

**Dated evidence, 2026-08-19. Verdict: HOLD all three genes.**

This is an accuracy assessment of the re-extracted output, not a diff against
the previous run. It exists because a change-based comparison cannot answer
"is what we are about to publish correct?".

## Why BMPR2 is not a control

An earlier pass used BMPR2 as a control on the reasoning that the cross-gene fix
changed nothing there. That is circular: nobody had ever validated BMPR2, so "no
change" is equally consistent with "uniformly wrong in a way this fix does not
touch". BMPR2 is assessed here on its own merits.

The audit vindicated the objection twice over:

- BMPR2's `wrong_gene` bucket was **empty**, which reads as clean. The audit then
  found 6 of 14 cross-gene in the adjacent bucket. The empty bucket was
  *blindness*: BMPR2's cross-gene errors hide in blank rowspan cells rather than
  in populated conflicting ones.
- variantFeatures flagged **zero** BMPR2 rows as `misparse_out_of_range`, so the
  pre-existing classifier reported nothing gene-related for BMPR2 at all. The
  residue check added 2026-08-19 flags **64**.

## Independent evidence: variantFeatures

Step 3.6 (`vf-enrich`) validates a variant's gene against a 36GB reference
warehouse. It was gated behind `full_coverage`, which defaults off, so **it never
ran on any `--pmid-file` run — which is how every collaborator review set was
built**. The "SKIPPED" log line sat inside the same branch, so nothing reported
its absence. `VARIANTFEATURES_DB` was also unset.

Ungated and run against the 150-paper re-extraction:

| Gene | GVF variants | matched | `wrong_gene_residue_mismatch` | `misparse_out_of_range` | total flagged |
|---|---:|---:|---:|---:|---:|
| BRCA1 | 4,462 | 38.8% | **504** | 59 | 594 |
| BRCA2 | 1,427 | 32.9% | **230** | 18 | 304 |
| BMPR2 | 458 | 48.0% | **64** | 0 | 85 |

Per-row lists: `<run>/<GENE>_vf_false_positives.csv`.

`wrong_gene_residue_mismatch` means the reference residue at that position does
not match this gene. That is an upper bound on cross-gene: a residue-numbering
or isoform difference produces the same signal. The confident cross-gene subset
is the part that *also* matches another gene at that position — 171 BRCA1,
181 BRCA2, 3 BMPR2 on a direct check.

## Sampled adjudication against the source papers

150 rows read against corpus full text, every defect claim adversarially
re-checked: **63 upheld, 0 overturned**.

| Gene | Est. precision | 95% interval |
|---|---:|---|
| BRCA1 | ~52% | 29–74% |
| BRCA2 | ~56% | 32–77% |
| BMPR2 | ~72% | 50–85% |

The intervals are wide (n≈14 per bucket) and must not be dropped when quoting
these. They cover sampling error only; three factors push the true rate up, not
down — the audit checked three axes only, 53 rows fell outside every bucket, and
auditors recorded "uncertain" rather than guessing. Treat these as floors.

**Even the optimistic end of the cleanest gene is one bad row in seven.**

## Verified defects (checked by hand, not relayed on trust)

1. **PMID 23469205 — BRCA1 founder alleles published under BRCA2.** The source
   row reads `| ID_2017 | 29 | (+) | BRCA1 | c.5382insC | ...`. `c.5382insC` is
   the Ashkenazi BRCA1 founder allele; also `C61G` / `c.300T>G` (canonical BRCA1
   RING domain), `R1751X`, `A1669V`. Someone searching BRCA2 for their family's
   5382insC would find it there.
2. **PMID 19944633 — canine BRCA2 staged as human.** Title: "Single nucleotide
   variation in exon 11 of **canine** BRCA2 in healthy and cancerous mammary
   tissue" (27 mentions of "canine"). 16 variants on dog residue numbering that
   does not map to human BRCA2. No gene-attribution fix addresses this; it is a
   paper-relevance failure upstream.
3. **PMID 30283497 is NOT the same problem** and must not be removed with it: a
   yeast assay *of human BRCA1 variants* (C61G, A1708E, M1775R). The
   discriminator is whether the species adjective modifies the gene ("canine
   BRCA2") or names the assay system ("BRCA1 ... in Yeast").
4. **Fabricated penetrance.** Both non-null BMPR2 penetrance values sampled were
   invented, both in the frightening direction — 100% and 43.9% against papers
   reporting 0/1 affected and 10–20%. Penetrance is the entire clinical question
   for a BMPR2 carrier.
5. **53 rows carry no variant identifier at all** (BMPR2 30, BRCA1 14, BRCA2 9)
   — nameless rows on a public page.

## Corrections to earlier claims in this iteration

- "Contamination is 3.5% BRCA1 / 7.9% BRCA2" — **undercount**. It replayed only
  the final filter and missed router-level gene-column filtering.
- "429 variants not present in the source" — **inflated by the audit tool**,
  which searched `185delAG` while papers write `185 delAG`. Whitespace-insensitive
  matching cut it to 374; the remainder is still unresolved.
- "The prose scanner is the dominant remaining cross-gene source" — **wrong**.
  The audit found most cross-gene errors arriving through structured table
  paths, which are the paths hardened first.

## Not covered

Count/penetrance accuracy beyond the sampled rows; the ~1,100 BRCA1 rows the
audit estimates carry a wrong carrier count; the residue mismatches that match
no gene (5.9% BRCA1, 6.3% BRCA2, 10.7% BMPR2), which are not cross-gene and have
no diagnosis yet; and every gene other than these three.
