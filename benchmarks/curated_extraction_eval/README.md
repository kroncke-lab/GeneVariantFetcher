# Curated Extraction Eval — a small, fixed benchmark for prompt/harness/guardrail changes

**What this is, in one sentence:** a hand-picked set of **24 gold-standard papers**
(across KCNH2, KCNQ1, SCN5A, RYR2) that the GeneVariantFetcher (GVF) pipeline
already extracts *very well*, so you can measure whether a change to the
**prompts, harness, guardrails, or variant matcher** helped or hurt — on
**recall** and **mean absolute error (MAE)** — **without re-running the entire
gold standard and burning a huge number of LLM tokens.**

If you are GPT‑5.5, Gemini, Codex, or a human who just opened this folder: this
README is self-contained. You do not need any other context to use it.

---

## Why it exists

The full GVF gold standard is ~1,500 papers / ~6,800 variant rows. Scoring a
pipeline change against all of it is slow and, if you re-extract, expensive. This
benchmark is the **fast inner loop**: a *fixed*, *curated*, *strategy-diverse*
subset you can run repeatedly to catch regressions and confirm improvements
cheaply. It is deliberately small (24 papers, ~740 gold variant rows) and
deliberately **easy for the current pipeline** — so any drop in the numbers is a
real signal that a change regressed something, not just noise from hard papers.

It complements (does not replace) the full scorer
(`scripts/run_recall_suite.py`) and the authoritative metrics in
`docs/RECALL_STATUS.md`. Use this for quick iteration; confirm big changes on the
full gold standard before claiming a headline number.

---

## Quick start

```bash
# From the repo root. Uses the project venv if present.

# 1) SCORE MODE (default, NO LLM tokens, takes seconds):
#    Score the four canonical DBs on the 24-paper subset.
python benchmarks/curated_extraction_eval/run_benchmark.py

# 2) Score a DB *you* produced (e.g. after a matcher/scorer/normalizer change):
python benchmarks/curated_extraction_eval/run_benchmark.py \
    --db KCNH2=/path/to/your/KCNH2.db --db SCN5A=/path/to/your/SCN5A.db

# 3) EXTRACT MODE (opt-in, COSTS LLM TOKENS): run the fixed PMID set through
#    the regular default gvf-run pipeline after selection, then score it.
#    Needs a working .env.
python benchmarks/curated_extraction_eval/run_benchmark.py \
    --mode extract --email you@example.com
```

The runner prints an aggregate scorecard (recall + MAE), a per-paper table
grouped by extraction strategy, and — once a baseline exists — the **delta vs the
committed baseline** so a regression is obvious at a glance.

---

## The two modes (read this before trusting a number)

| | `score` (default) | `extract` |
|---|---|---|
| **Cost** | Free. No LLM calls. Seconds. | **Costs LLM tokens.** Minutes. |
| **What it does** | Scores an *existing* SQLite DB against the frozen gold subset. | Runs every gene in the registry through default `gvf-run` from the fixed PMID lists -> builds DBs -> scores them. |
| **Use it to test** | Matcher, scorer, normalizer, or any DB you built elsewhere. | Prompt, extraction harness, guardrail, count-classifier changes. |
| **How** | Wraps `scripts/run_recall_suite.py` pointed at `./gold`. | Wraps `python -m cli gvf-run <GENE> --pmid-file pmids/<GENE>.txt` with default source recovery, recovery layers, corpus sync, and report behavior. |

`extract` mode starts after paper selection: the PMID files are fixed and
discovery/Tier-1/Tier-2 filtering are intentionally skipped. Everything after
that uses the regular driver: source cache/download, preprocessing, Data Scout,
Tier-3 extraction, migration, recovery layers, source QC/source recovery, corpus
sync, and report generation. If `sources/` has been materialized by
`build_fixture.py`, the runner uses it as a local mini corpus; otherwise it falls
back to the repo's normal `corpus/` cache. `--fast` is available only for a
cheap diagnostic pass and adds `--no-source-recovery`; do not use it for the
benchmark number.

---

## What the metrics mean

- **Recall** = matched / gold. Reported for: `pmids`, `variant_rows`,
  `unique_variants` (the headline — distinct variants found), `patients`,
  `affected`, `unaffected`. Higher is better; the program target is 90%.
- **Rows-mode MAE** = mean of `|gold_count − extracted_count|` over matched
  (paper × variant) rows, for carriers / affected / unaffected counts. Lower is
  better; target → 0. This benchmark intentionally includes papers with
  non-trivial count disagreement (e.g. PMIDs 10973849, 22949429) so MAE has
  signal to move — a count-semantics regression will show up here.
- **Precision proxy** (in the full scorer artifacts): an extra-rows-vs-gold rate,
  a *loose false-positive upper bound*, not clean precision. See
  `docs/RECALL_STATUS.md` for why it overstates false positives.

---

## Current baseline (the numbers to hold or beat)

Measured by `score` mode against the four **canonical DBs** in
`docs/RECALL_STATUS.md`. The machine-readable copy is **`expected_baseline.json`**
(the runner diffs against it automatically); regenerate with `--write-baseline`.

| Metric | Subset recall | (Full gold standard, for contrast) |
|---|---:|---:|
| PMIDs | 100% (24/24) | ~85% |
| Variant rows | ~99.6% | ~81% |
| **Unique variants** | **~99.5%** | ~86% |
| Patients / Affected / Unaffected | ~99.9% / ~99.8% / 100% | ~85% / ~84% / ~87% |
| Carriers MAE | ~0.16 | ~0.61 |

The subset scores far higher than the full gold standard **on purpose** — these
are papers the pipeline already nails, so the benchmark is a sensitive regression
detector. A change that drops subset unique-variant recall below ~99% or pushes
carriers MAE above ~0.2 is doing real damage and should be investigated before it
ships.

---

## What's in this folder

```
curated_extraction_eval/
├── README.md            ← you are here
├── registry.tsv         ← THE editable list of papers (gene, pmid, strategy, note). Append a row to grow the set.
├── add_paper.py         ← one-command add: appends to registry, checks gold+source, rebuilds
├── manifest.md          ← human-readable table of the papers (auto-generated)
├── manifest.csv         ← machine-readable: gene, pmid, strategy, gold counts, gold source, title, corpus path, PubMed URL
├── expected_baseline.json ← committed baseline the runner diffs against
├── run_benchmark.py     ← the runner (score + extract modes)
├── build_fixture.py     ← regenerator: reads registry.tsv + gold and rebuilds every derived file
├── pmids/<GENE>.txt      ← PMID lists for `gvf gvf-run --pmid-file` (auto-generated)
├── gold/normalized/<GENE>_recall_input.csv ← FROZEN gold SUBSET (auto-generated; same schema as the full gold, so the scorer runs unmodified)
├── gold_overrides/<GENE>_recall_input.csv  ← curator-supplied gold for papers the repo gold standard doesn't cover (see its README)
└── sources/<GENE>/<PMID>/... ← local mini corpus snapshot (gitignored; same layout as corpus/<GENE>/<PMID>/)
```

`.last_run/` and `.extract_runs/` (gitignored) hold the most recent scorer and
extraction artifacts: `summary.json`, `report.md`, and per-gene
`discrepancies.csv` (one row per gold/DB variant with match status and count
diffs — the place to look when a paper regresses).

---

## Adding papers to the set (it is built to grow)

The set is meant to **accrete problem cases**: when you find a gene/variant paper
that extracts badly — a missed supplement table, a figure the vision layer drops,
a count it gets wrong — add it here so future changes are measured against it and
can't silently re-break it. Adding is one command:

```bash
# 1) Add the paper. --strategy is how its variants are encoded (table|text|figure|mixed).
python benchmarks/curated_extraction_eval/add_paper.py SCN5A 12345678 \
    --strategy table --note "supplement mutation table dropped by the API body fetch"
```

`add_paper.py` appends a row to `registry.tsv`, tells you whether the paper has a
**gold answer** and **cached source**, and rebuilds the fixture. (Or just append a
tab-separated `gene<TAB>pmid<TAB>strategy<TAB>note` line to `registry.tsv` by hand
and run `python build_fixture.py`.)

Every paper needs a **gold answer** to be scoreable. Two cases:

- **Gene already has a gold standard (KCNH2/KCNQ1/SCN5A/RYR2) and the PMID is in
  it** → nothing else to do; the gold is pulled automatically.
- **A new gene (BRCA1, APOE, MYBPC3, …) or a paper the curator never scored** →
  supply the expected variants in
  `gold_overrides/<GENE>_recall_input.csv` (same schema as the repo gold; only
  `variant,pmid,carriers,affected,unaffected` are required). See
  [`gold_overrides/README.md`](gold_overrides/README.md). This is how the set
  grows beyond the four cardiac genes. `add_paper.py` and `build_fixture.py` both
  tell you exactly when an override is needed.

The paper's **source** must be cached at `corpus/<GENE>/<PMID>/` (true after any
`gvf-run` that fetched it) so `build_fixture.py` can materialize it under
`sources/<GENE>/<PMID>/`. `score` mode can run without source; `extract` mode
needs either a complete local `sources/` mini corpus or the normal repo `corpus/`
cache.

> New problem papers will (correctly) score **below** the well-extracted originals
> until the pipeline is fixed. That is the point of a regression set — the
> per-paper table and the baseline delta show exactly which ones still fail.
> After a real improvement lands, re-run `run_benchmark.py --write-baseline` to
> move the bar up.

## How the original 24 papers were chosen

Goal: papers the pipeline **already extracts well**, spanning the **strategies**
a curator uses to pull variants out of a paper, balanced across the four genes.

1. **Candidate pool:** every gold PMID, scored against the canonical DBs.
2. **"Well-extracted" filter:** ≥6 gold variants, ≥85% row recall, and accurate
   counts (summed |count error| ≤ 1.0 per matched row). Every pick here was then
   **re-validated live** against the canonical DBs (≥95% row recall; most 100%).
3. **Strategy tag** (from the dominant `variant_papers.source_layer`):
   - **table** — variants in mutation-list tables (`regex_table`/`llm_table`)
   - **text** — variants in running prose (`llm_text`/`regex_text`)
   - **figure** — variants read from pedigrees/figures by the vision layer
   - **mixed** — table + figure both contribute (exercises routing)
4. **Curation:** 24 picked to balance gene × strategy, span small (7 variants) to
   large (90), and include count-disagreement cases so MAE has signal.

Final mix: **9 table · 8 text · 5 figure · 2 mixed**, across **KCNH2 6 · KCNQ1 5
· SCN5A 7 · RYR2 6**. See `manifest.md` for the full list with titles and links.

> KCNH2 has no pure-*figure* paper in the set: in this gold standard its
> well-extracted papers are table/text, and the clean figure-sourced papers are
> in SCN5A and RYR2. That is a property of the data, recorded here so it is not
> mistaken for an omission.

### Reproduce the selection / regenerate the fixture

`registry.tsv` is the source of truth for *which* papers are in the set and
*why*; `build_fixture.py` reads it (plus the repo gold + `gold_overrides/` +
`corpus/`) and regenerates the gold subset, PMID lists, manifest, and local
mini-corpus source snapshot deterministically:

```bash
python benchmarks/curated_extraction_eval/build_fixture.py            # rebuild everything
python benchmarks/curated_extraction_eval/build_fixture.py --no-sources  # skip the full-text copy
```

To change the set, edit `registry.tsv` (or use `add_paper.py`) and re-run. The
original candidate pool can be regenerated by scoring the full gold against the
canonical DBs and filtering as described above.

---

## Gotchas

- **Canonical DB paths** are hard-coded in `run_benchmark.py` and `build_fixture.py`
  to match `docs/RECALL_STATUS.md`. If those DBs move, update both (one constant
  each). On a fresh checkout without them, pass `--db GENE=/path` explicitly.
- **The gold subset is frozen.** It is a row-filtered copy of
  `gene_variant_fetcher_gold_standard/normalized/<GENE>_recall_input.csv`. If the
  full gold is re-curated, re-run `build_fixture.py` to refresh the subset.
- **`sources/` is gitignored** (paper full text is never committed, matching repo
  policy). It exists locally after `build_fixture.py` and mirrors the
  `corpus/<GENE>/<PMID>/` layout. `score` mode does not need it; `extract` mode
  uses it as `GVF_CORPUS_DIR` when complete, otherwise it falls back to the
  normal repo `corpus/`.
