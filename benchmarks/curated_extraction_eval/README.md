# Curated Extraction Eval — a small, fixed benchmark for prompt/harness/guardrail changes

**What this is, in one sentence:** a hand-picked set of **101 gold-standard papers**
across **8 gene-disease pairs** — the 4 cardiac channelopathy genes (KCNH2,
KCNQ1, SCN5A, RYR2) plus hereditary cancer (BRCA1, BRCA2), hypertrophic
cardiomyopathy (MYBPC3), and Alzheimer's/lipid disorders (APOE) — spanning
**publication years 1989–2025**, all four extraction strategies, cohort sizes
from single case reports to 100+ variant screens, and **deliberately-included
failure-mode and false-positive-guard cases**, so you can measure whether a change
to the **prompts, harness, guardrails, or variant matcher** helped or hurt — on
**recall** and **mean absolute error (MAE)** — **without re-running the entire
gold standard and burning a huge number of LLM tokens.**

If you are GPT‑5.5, Gemini, Codex, or a human who just opened this folder: this
README is self-contained. You do not need any other context to use it.

---

## Why it exists

The full GVF gold standard is ~1,500 papers / ~6,800 variant rows. Scoring a
pipeline change against all of it is slow and, if you re-extract, expensive. This
benchmark is the **fast inner loop**: a *fixed*, *curated*, *diversity-spanning*
subset you can run repeatedly to catch regressions and confirm improvements
cheaply. It is small (101 papers, ~3,000 gold variant rows) but **deliberately
diverse** — it pairs papers the pipeline already nails (the high-recall floor
that makes regressions obvious) with genes, eras, and failure modes the pipeline
does *not* yet handle well, so the aggregate carries real headroom (~90%) for a
protocol change to move. Per-paper scores tell you which is which: a
high-recall paper dropping is a regression; a known-hard paper rising is a win.

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

Measured by `score` mode against the eight **canonical DBs** wired into
`run_benchmark.py` (the 4 cardiac DBs from `docs/RECALL_STATUS.md` plus the
BRCA1/BRCA2/MYBPC3/APOE run DBs). The machine-readable copy is
**`expected_baseline.json`** (the runner diffs against it automatically);
regenerate with `--write-baseline`.

| Metric | Subset recall (101 papers) | (Full cardiac gold standard, for contrast) |
|---|---:|---:|
| PMIDs | 99.0% (100/101) | ~85% |
| Variant rows | 90.1% (2728/3029) | ~81% |
| **Unique variants** | **92.3% (1863/2019)** | ~86% |
| Patients / Affected / Unaffected | 94.6% / 93.9% / 97.7% | ~85% / ~84% / ~87% |
| Carriers MAE | 0.300 | ~0.61 |

This set sits near ~90% **on purpose**: it blends papers the pipeline nails
(the regression floor) with deliberately-hard new-gene and failure-mode papers
(the headroom). So read the **per-paper table**, not just the aggregate — a
high-recall paper dropping is a real regression; the known-hard papers (e.g.
KCNH2 29650123, RYR2 19398665/27452199, SCN5A 24144883, KCNQ1 24667783) are
extraction-gap targets that *should* rise when a protocol change works. The
4 `negative_cases/` guards (run `pytest benchmarks/curated_extraction_eval/negative_cases`)
must stay green — they assert the protocol does NOT mint carriers from
annotation/predictor/wrong-gene tables.

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
├── fixtures/            ← deterministic synthetic cases used by unit tests; never folded into curated gold
├── negative_cases/      ← PRECISION/guard cases: assert the protocol does NOT mint carriers from annotation/predictor/wrong-gene data (gold-free, self-contained pytest). See its README.
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

## How the 101 papers were chosen

The set was built in two arms to span **four axes**: gene-disease pair, era,
extraction strategy, cohort size — and to weight in **failure modes**, not just
papers the pipeline nails.

**Cardiac arm (73 papers · KCNH2/KCNQ1/SCN5A/RYR2).** Drawn from the repo gold
standard (free to score). Per gold PMID, the canonical DB gave strategy (dominant
`variant_papers.source_layer`), era (`papers.publication_date`), and approximate
recall; picks were chosen to add era anchors (e.g. KCNH2 7889573 / 1995, RYR2
34661651 / 2022), size anchors (KCNQ1 32893267 / 249 variants down to 3-variant
case reports), strategy coverage, and in-source failure cases the pipeline misses
(KCNH2 29650123, RYR2 19398665/27452199, SCN5A 24144883, KCNQ1 24667783).

**Non-cardiac arm (28 papers · BRCA1/BRCA2/MYBPC3/APOE).** These genes have no
repo gold standard, so each paper's gold answer was **curated by hand into
`gold_overrides/` from the cached full text/tables** (the extraction DB was used
only as a row locator, never trusted as gold). Counts the source did not state
were left out rather than inferred; APOE allele-frequency-only papers were dropped
in favor of ones with concrete coding variants + patient counts. The arm spans
1989 (APOE 2539388) to 2025 (MYBPC3 40453736) and includes count-role-confusion
cases (BRCA2 21356067, BRCA1 33468216) on purpose.

**Strategy tags** (`table | text | figure | mixed`, from the dominant source
route) — current mix: **54 table · 24 text · 12 figure · 11 mixed**, across
**KCNH2 15 · KCNQ1 19 · SCN5A 19 · RYR2 20 · BRCA1 8 · BRCA2 6 · MYBPC3 8 ·
APOE 6**. Plus 4 `negative_cases/` guards (gold-free FP assertions). See
`manifest.md` for the full list with titles and links.

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
