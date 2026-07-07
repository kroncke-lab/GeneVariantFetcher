# benchmarks/

Runnable evaluation sets for the GeneVariantFetcher (GVF) extraction pipeline.

A "benchmark" here is a **small, fixed, curated set of gold-standard papers** plus
a one-command runner that reports extraction quality (recall + MAE) on just that
set — so pipeline changes can be measured **cheaply and repeatably**, without
re-running the entire ~1,500-paper gold standard. These are the fast inner loop;
`scripts/run_recall_suite.py` + `docs/RECALL_STATUS.md` remain the authoritative
full-gold scorer and metrics.

## Sets

- **[`curated_extraction_eval/`](curated_extraction_eval/README.md)** — 101
  hand-picked, strategy-diverse gold papers across KCNH2, KCNQ1, SCN5A, RYR2,
  BRCA1, BRCA2, MYBPC3, and APOE, spanning tables, figures, in-text evidence,
  negative cases, and false-positive guard cases. It scores existing DBs or runs
  the fixed PMID set through the regular default `gvf-run` post-selection
  pipeline. **Start here** — its README is self-contained.

  ```bash
  python benchmarks/curated_extraction_eval/run_benchmark.py        # score canonical DBs (free)
  python benchmarks/curated_extraction_eval/run_benchmark.py --db KCNH2=/path/to.db
  python benchmarks/curated_extraction_eval/run_benchmark.py --mode extract --email you@example.com
  ```

## Adding a new set

Copy the structure of `curated_extraction_eval/`: a `registry.tsv` listing the
papers, a `build_fixture.py` that regenerates all derived files from it, a frozen
`gold/` subset (+ `gold_overrides/` for papers the repo gold doesn't cover), a
`run_benchmark.py`, an `add_paper.py`, and a self-contained `README.md`. Keep
paper full text out of git (gitignore `sources/`), matching repo policy.
