# Population-allele predictor snapshot

This is a read-only snapshot of existing VariantFeatures annotations for every
row in the five-gene population inventory. It does not run annotation jobs or
download additional predictors. Join `population_predictors.csv.gz` on
`gene,variant_id`, where `variant_id` is the population genomic allele ID.

The identity check requires exact GRCh38 chromosome, position, reference allele,
alternate allele, and warehouse gene membership. Only a chromosome `chr` prefix
is removed. Indels are not normalized or proxied, and genomic alleles are not
collapsed to protein substitutions. Every input row remains present, including
alleles with no warehouse match and rows the population analysis later excludes
for consequence, counts, transcript, or QC. Missing scores remain missing.

The three scalar features are AlphaMissense, GPN-Star M447 signed calibrated
LLR, and AlphaGenome Atlas AVI. All three are stored in the warehouse's
`annotations_pathogenicity` table. GPN entropy is in conservation annotations;
it is not the requested signed M447 LLR and is not substituted here.

- AlphaMissense increases with predicted pathogenicity; it is not penetrance.
- More negative signed GPN-Star LLR indicates greater predicted impact or
  constraint; it is not a probability.
- Higher AVI indicates greater predicted variant impact. AVI includes
  AlphaMissense information, so these are not independent predictors.

An available score requires exactly one distinct finite `(version, score)` pair.
Multiple versions are a conflict even if their values agree. No arbitrary
latest-version rule is used. AVI additionally must agree exactly with the
authoritative scalar `AVI_SCORE` context's `(dataset_version, raw_score)` pair.
The conflict status is saved alongside each missing score.

`predictor_annotation_rows.csv.gz` preserves queried score rows and their source
versions/timestamps. `avi_scalar_contexts.csv.gz` preserves the scalar context
records used for validation. `coverage.csv` reports availability by gene and
predictor; `predictor_provenance.json` records the logical read snapshot, input
and output hashes, collector hash, version distribution, and query intervals.
The SQLite source was opened with `mode=ro`, `query_only=ON`, and one read
transaction. No source database write occurred.

The historical literature-only analysis's archived AlphaMissense values remain
in `../../gck_structural_pilot_20260912/eligibility/archived_alphamissense.csv`.
This new genomic snapshot does not silently replace or relabel those archived
values. Any comparison combining historical protein keys and genomic alleles
must state its identity and score-source policy.

Reproduce from the GeneVariantFetcher repository with its Python environment:

```bash
.venv/bin/python docs/evidence/population_inclusive_penetrance_20260912/predictors/collect_population_predictors.py
.venv/bin/python docs/evidence/population_inclusive_penetrance_20260912/predictors/validate_population_predictors.py
```

The external source volume and both VariantFeatures storage symlinks must be
mounted. The script fails if the existing source database is absent. Because
the warehouse can receive new annotations, a later rerun is a new dated
snapshot; use the saved rows and hashes to reproduce this snapshot's analysis.
The validator checks complete inventory preservation and reconstructs score
resolution from the saved annotation/context rows, without reopening SQLite.
