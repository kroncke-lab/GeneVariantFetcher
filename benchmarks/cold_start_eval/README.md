# Cold-Start Benchmark (clean room)

**Question this answers:** can GVF cold-start an *unseen* gene end to end —
discover → acquire → extract — not just extract well on a pre-fetched, PMID-fixed
corpus?

This is deliberately the complement of
[`../curated_extraction_eval/`](../curated_extraction_eval/README.md). The curated
benchmark freezes the corpus and the PMID list to isolate *extraction* quality.
This one measures the stages the curated benchmark skips: **discovery, filtering,
and acquisition**, starting from nothing.

## Why a separate benchmark

The four cardiac gold genes (KCNH2, KCNQ1, SCN5A, RYR2) and the eight curated
registry genes have been used repeatedly for gold-gated tuning. Strong numbers
there are partly *in-distribution* — they don't demonstrate that an arbitrary
gene works. A clean-room run on a gene that appears in **neither** the gold
standard nor the registry, with an **empty corpus**, is the honest test of the
"give me any gene" claim.

## Isolation guarantees

The harness enforces four isolations so the result reflects a true cold start
(see `run_cold_start.py` for the code and citations):

| Isolation | Mechanism | Prevents |
|---|---|---|
| Empty corpus | `GVF_CORPUS_DIR` → a fresh empty dir | Warm-starting from cached full text |
| No corpus write-back | `--no-corpus-sync` | Polluting the repo `corpus/` with the cold gene |
| Fresh output dir | unique `--output` | Scavenging a prior run's downloads |
| No gold PMIDs | **omit** `--pmid-file` | Leaking gold into discovery (live discovery only) |

Because an unseen gene has no gold CSV, recovery runs DB-observed and per-layer
recall scoring is skipped by design. Discovery/acquisition metrics come from the
gold-free source-QC bundle (`source_qc/source_acquisition_summary.json`).

## Genes to use

Pick a gene in neither the gold standard nor the curated registry, non-cardiac,
not already cached in `corpus/`. Suggested (disease in parentheses):

- **LDLR** (familial hypercholesterolemia)
- **MLH1 / MSH2** (Lynch syndrome)
- **TP53** (Li-Fraumeni syndrome)
- **CFTR** (cystic fibrosis)
- **PAH** (phenylketonuria) · **PTEN** · **VHL**

The harness refuses genes that already have coverage unless you pass
`--allow-non-cold-gene`.

## Run it

Dry run (default — prints the exact hermetic command, creates the isolated
corpus dir, launches nothing):

```bash
source .venv/bin/activate
python benchmarks/cold_start_eval/run_cold_start.py LDLR
```

Live run (long — full discovery + acquisition + extraction):

```bash
python benchmarks/cold_start_eval/run_cold_start.py LDLR --run
```

After a live run, acquisition metrics print automatically and are also at:

```
benchmarks/cold_start_eval/runs/LDLR/<timestamp>/source_qc/source_acquisition_summary.json
```

## Scoring extraction (optional, after the run)

Cold start has no gold, so it is not auto-scored. To measure extraction recall
without leaking gold into acquisition:

1. Curate a gold answer for the cold gene **after** the run — e.g. a
   `gene_variant_fetcher_gold_standard/normalized/<GENE>_recall_input.csv` or a
   `curated_extraction_eval/gold_overrides/<GENE>_recall_input.csv`.
2. Score the produced DB against it:

   ```bash
   python scripts/run_recall_suite.py --score \
     --gold-dir gene_variant_fetcher_gold_standard/normalized \
     --results-dir benchmarks/cold_start_eval/runs
   ```

3. For a true-precision read on the cold gene's extras, use
   [`scripts/precision_sample.py`](../../scripts/precision_sample.py).

Report the result explicitly as **cold-start, arbitrary-gene** validation —
distinct from the in-distribution cardiac gold metrics.
