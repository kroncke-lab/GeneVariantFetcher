# Population source inventory

`population_variants.csv.gz` contains **30,101 actual gnomAD gene-query records**
for the five genes, including indels and alleles absent from the literature.
These records came from complete public gnomAD 4.1.1 gene responses cached on
2026-09-11. Fresh gene and canonical-transcript annotation queries on 2026-09-12
return exactly the same genomic allele inventory and identical exome, genome,
and joint AC, AN, homozygote counts, and filter lists for every allele.

| Gene | All source rows | Observed, QC pass | Canonical coding/splice, QC pass | Canonical missense, QC pass |
| --- | ---: | ---: | ---: | ---: |
| HNF1A | 3,326 | 2,389 | 1,302 | 723 |
| GCK | 3,542 | 2,292 | 956 | 468 |
| LDLR | 5,497 | 4,218 | 2,216 | 1,230 |
| BRCA2 | 10,769 | 9,037 | 6,811 | 4,231 |
| KCNQ1 | 6,967 | 3,912 | 1,716 | 833 |

AC0 records and QC failures remain inspectable but have `population_eligible=False`.
Every canonical protein reference parsed from these annotations matches the
corresponding gene-verified UniProt sequence. The collector does not decide
literature overlap or assign affected counts; the union stage owns those steps.

## Scope and annotation identity

The gene API retrieves variants in the union of CDS intervals across the gene's
transcripts, padded by **75 bases** on each side. Overlapping intervals are merged
and results are paginated internally. The full genomic gene span, deep intronic
regions, remote regulatory elements, CNVs, and structural variants are not
enumerated by this endpoint. The corresponding transcript endpoint restricts the
same procedure to the selected transcript. This boundary is implemented in the
[official browser source, pinned to the inspected commit](https://github.com/broadinstitute/gnomad-browser/blob/5c14bfea3298f46854df14781adc1ccead58da3b/graphql-api/src/queries/variant-datasets/gnomad-v4-variant-queries.ts).

The gene summary can choose another transcript. Consequently, primary `hgvsc`,
`hgvsp`, `consequence`, and `transcript_id` fields come exclusively from an exact
canonical-transcript query. Other-transcript gene summaries are retained in
separate `gene_annotation_*` columns. Unavailable canonical annotation stays
blank. These population alleles can remain in the broad source inventory, but
they cannot silently become canonical structural donors.

All five selected transcript versions match the archived literature features:
HNF1A ENST00000257555.11; GCK ENST00000403799.8; LDLR ENST00000558518.6;
BRCA2 ENST00000380152.8; KCNQ1 ENST00000155840.12. Protein sequences and gene
assignment are independently checked against the corresponding official
UniProt records listed in `population_provenance.json`.

## Counts and QC

For these autosomal genes, `gnomad_carriers = joint_ac - joint_hom`: each
homozygote contributes two alternate alleles but one carrier. Use the official
joint cohort once. The table preserves all assay-specific counts; adding exome,
genome, and joint counts would repeat observations. Population carriers are
treated as unaffected under the user's specified assumption. AC remains
available for sensitivity analysis.

`qc_pass` requires every present sequencing assay to pass its filters. A missing
assay is permitted. This reconstructs the documented joint VCF categories:
PASS, EXOMES_FILTERED, GENOMES_FILTERED, or BOTH_FILTERED.
The gnomAD team explains these semantics in its
[joint-filter clarification](https://discuss.gnomad.broadinstitute.org/t/is-joint-combined-genome-exome-faf-unreliable-if-either-genome-or-exome-fails-filters/88/3).

The API's `joint.filters` is retained as source data, with a separate
`joint_api_pass` field. It is unsuitable as the joint VCF filter: the inspected
summary code reads singular `joint.filter` even though its data selection asks
for plural `joint.filters`, and can append AC0 based on the exome count. Thus
some API joint filter lists are empty despite assay QC failures, or show AC0
despite a positive joint count. The derived policy uses the explicit assay
lists and positive joint AC.

## Reproduce and extend

`fetch_annotations.py` copies the earlier complete source snapshots from the
local VariantFeatures cache and downloads the bounded canonical annotation and
UniProt queries. Requests are saved alongside the parser; raw responses remain
under ignored `results/population_inclusive_penetrance_20260912/raw/`.
`build_population.py` performs identity, count, QC, and WT checks and writes the
compact CSV and provenance. The downstream union can run directly from the
committed compact snapshot without an API call or warehouse database.

`full_locus_availability.csv` records the five GRCh38 gene spans for a possible
separate full-locus analysis. The official region endpoint supports interval
queries, with limits of less than 2.5 Mb and no more than 30,000 records per
request; larger responses need subdivision. No full-locus expansion is mixed
into this CDS-and-flank snapshot.
