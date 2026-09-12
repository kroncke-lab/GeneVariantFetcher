# Full genomic gene-span population extension

This additive source expands the earlier CDS-and-flank inventory to the exact
GRCh38 genomic start/stop interval of each of the five genes. Every returned
small-variant allele is retained, including deep intronic and other noncoding
positions. Only records with positive adjusted joint AC and passing QC are
eligible for the population-inclusive prior. No potential-but-unobserved allele
is inserted as an unaffected observation.

The five intervals and boundaries come from the official gene metadata frozen
in the parent population provenance. `full_locus_summary.csv` reports the
completed inventories. Use each `GENE_provenance.json` file's `output_files`
list to load its CSV parts; every compressed part is at most 1.1 MB.

The collector queries the official gnomAD region endpoint in nonoverlapping
20 kb intervals. It subdivides a request if the service reports a size limit.
For each gene it checks complete interval coverage, genomic allele identities,
autosomal count consistency, duplicate concordance, and exact agreement with
every earlier CDS-footprint allele whose POS lies inside the gene span. This
includes matching all exome, genome, and joint AC/AN/homozygote/filter fields.
The region query implementation and bounds are inspectable in the
[official gnomAD browser source](https://github.com/broadinstitute/gnomad-browser/blob/5c14bfea3298f46854df14781adc1ccead58da3b/graphql-api/src/graphql/resolvers/variants.ts).

`gnomad_carriers = joint_ac - joint_hom` counts each autosomal carrier once per
allele. The user assumes these carriers are unaffected. The three sequencing
summaries are preserved separately; the calculation uses the official joint
cohort once. QC requires every present assay to pass, with absent assays
permitted. The raw API joint-filter field remains visible but is not used as
the joint VCF PASS flag; see the parent inventory's documented API filtering
issue and the
[gnomAD team's joint-filter definition](https://discuss.gnomad.broadinstitute.org/t/is-joint-combined-genome-exome-faf-unreliable-if-either-genome-or-exome-fails-filters/88/3).

Region-only records receive no invented protein or transcript annotation.
The union joins the previous canonical annotations by exact genomic allele ID.
Remaining variants stay explicitly unannotated for structural donor purposes.
Membership in a gene's genomic span does not establish a functional effect on
that gene, and overlapping genes can share a genomic region. These are
small-variant records selected by POS; structural variants/CNVs and deletions
anchored outside the interval require separate inventories.

Raw requests and responses remain in the ignored
`results/population_inclusive_penetrance_20260912/raw/full_locus/` directory.
Every CSV row identifies its source interval and request/response digest.
Per-gene provenance records raw file hashes, retrieval times, overlap checks,
and output part hashes. `gnomad_r4` is a rolling selector, so these outputs are
frozen snapshots of the current 4.1.1 dataset, checked against the earlier
4.1.1 source counts. See the official
[4.1.1 release description](https://gnomad.broadinstitute.org/news/2026-03-gnomad-v4-1-1).

From the GVF root, rerun the collector using its existing raw cache:

```sh
.venv/bin/python docs/evidence/population_inclusive_penetrance_20260912/population/full_locus/fetch_full_locus.py
```

The parent CDS-and-flank snapshots remain unchanged and provide the narrower
universe sensitivity analysis.
