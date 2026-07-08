# Variant Browser Paper Test Cases

This directory contains a 50-paper, source-verified exploration slice drawn from
the curated extraction benchmark. It is intended for browser/review QA and human
toy testing, not as a replacement for the full scored benchmark.

`cases.tsv` includes:

- 24 public-browser cardiac papers: `KCNH2`, `KCNQ1`, `SCN5A`, `RYR2`
- 6 review-staging papers: `BRCA2`
- 20 candidate/non-public papers: `APOE`, `BRCA1`, `MYBPC3`

Use `browser_scope` to filter to the current public browser genes when needed.
Every row has cached source and source-verified gold in the 101-paper curated
fixture, so this set can be used without running extraction mode.

Coverage summary:

| dimension | coverage |
| --- | --- |
| genes | APOE 6, BRCA1 7, BRCA2 6, MYBPC3 7, KCNH2/KCNQ1/SCN5A/RYR2 6 each |
| eras | pre-2000 11, 2000s 12, 2010s 20, post-2020 7 |
| strategies | table 20, text 13, mixed 10, figure 7 |
| size buckets | tiny 6, small 7, medium 21, large 6, very_large 10 |

Columns:

- `gene`, `pmid`, `pubmed_url`: paper identity.
- `browser_scope`: `public_browser`, `review_staging`, or `candidate_not_public`.
- `disease_area`: broad domain label for cross-domain sampling.
- `year`, `era`: PubMed publication year and era bucket.
- `strategy`: primary evidence layout from the curated benchmark registry.
- `size_bucket`: rough human-review size bucket based on gold rows/carriers.
- `gold_*`: source-verified benchmark gold counts.
- `source_status`: cached-source status from the fixture manifest.
- `why_selected`: the intended manual-QA failure mode or coverage reason.
- `title`: PubMed title.
