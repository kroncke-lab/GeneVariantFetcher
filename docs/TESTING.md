# GVF Thorough Test Prompt

Paste this into any agent (Claude Code, codex, etc.) standing in front of a
fresh checkout of GeneVariantFetcher. It exercises every layer of the
pipeline and the post-extraction recovery stack, then scores against the
gold standard.

---

## Prompt

> Run a thorough test of GeneVariantFetcher end-to-end. Use the project
> `.venv`. Do not modify CLAUDE.md, .env, or anything in `results/` that
> isn't a generated artifact. Treat any failure as a stop-and-investigate
> signal — don't paper over errors.
>
> ### 1. Environment + static checks
> 1. Verify `.venv/` exists and is Python 3.11+; if missing, recreate with
>    `python3.11 -m venv .venv && .venv/bin/pip install -e ".[dev]"`.
> 2. Confirm required env vars: `NCBI_EMAIL`, `NCBI_API_KEY`,
>    `OPENAI_API_KEY` *or* `AZURE_AI_API_KEY`+`AZURE_AI_API_BASE`,
>    `ELSEVIER_API_KEY`, `WILEY_API_KEY`. Report missing ones.
> 3. Run linters: `.venv/bin/python -m ruff check . --no-fix` and
>    `.venv/bin/python -m ruff format --check .`. Report any issues.
>
> ### 2. Unit tests
> Run `.venv/bin/python -m pytest tests/unit -q`. Expect 390+ passing,
> at most 1 skipped. Capture the count and any failures.
>
> ### 3. Module-import smoke
> Verify every importable module loads cleanly:
> ```bash
> .venv/bin/python -c "
> import importlib, pkgutil
> for pkg in ['cli', 'pipeline', 'harvesting', 'utils', 'gene_literature', 'gui']:
>     for _, name, _ in pkgutil.walk_packages([pkg], prefix=pkg + '.'):
>         try: importlib.import_module(name)
>         except Exception as e: print(f'IMPORT FAIL {name}: {e}')
> print('All modules imported')"
> ```
> No `IMPORT FAIL` lines should appear.
>
> ### 4. End-to-end pipeline for KCNH2 (small)
> Run a bounded e2e to confirm the pipeline produces a valid DB. Use
> `--max-pmids 25` so it finishes in minutes, not hours:
> ```bash
> .venv/bin/python main.py --gene KCNH2 --max-pmids 25 --output results/test_kcnh2
> ```
> Expect `results/test_kcnh2/KCNH2/<timestamp>/KCNH2.db` with non-zero
> papers and variants. Print the row counts.
>
> ### 5. Recovery scripts smoke test
> Each of these should import + parse `--help` without error:
> ```bash
> for s in merge_v12_db ingest_clinvar ingest_pubtator recover_paywall_oa; do
>   .venv/bin/python scripts/recall_recovery/$s.py --help || echo "FAIL $s"
> done
> .venv/bin/python scripts/extract_figure_variants.py --help
> .venv/bin/python scripts/run_recall_suite.py --help
> ```
>
> ### 6. Figure-reader on existing data
> If `results/KCNH2/20260517_074737/pmc_fulltext/*_figures/` exists, pick
> one PMID with figures and run:
> ```bash
> .venv/bin/python scripts/extract_figure_variants.py \
>   --gene KCNH2 --pmid 19038855 \
>   --pmc-dir results/KCNH2/20260517_074737/pmc_fulltext \
>   --out /tmp/figure_test
> ```
> Expect a JSON report under `/tmp/figure_test/19038855.json` with a
> non-empty `distinct_variants` list (this PMID's Table 2 yielded 34
> variants in our last run).
>
> ### 7. Score against gold
> ```bash
> .venv/bin/python scripts/run_recall_suite.py --score --genes KCNH2 \
>   --db KCNH2=results/KCNH2/20260517_074737/KCNH2.db \
>   --outdir recall_metrics/test_$(date +%Y%m%d_%H%M%S)
> ```
> The current baseline (post-cleanup, post-figure-reader) is:
>
> | Metric | Recall | Target |
> |---|---|---|
> | pmids | 90.8% | 90% ✓ |
> | variant_rows | 82.2% | 90% |
> | unique_variants | 83.0% | 90% |
> | patients | 85.8% | 90% |
> | affected | 87.5% | 90% |
> | unaffected | 81.8% | 90% |
>
> Anything materially below this is a regression; investigate.
>
> ### 8. Report
> Produce a single-page summary:
> * Test step → pass/fail + 1-line note
> * Final recall numbers vs baseline above
> * Any modules that failed to import
> * Any failing unit tests
> * One paragraph: "if this branch shipped today, would you be comfortable?"
