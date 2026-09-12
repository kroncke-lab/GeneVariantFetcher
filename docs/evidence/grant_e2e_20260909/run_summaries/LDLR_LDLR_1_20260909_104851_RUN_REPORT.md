# GVF Run Report — LDLR

- Started: 2026-09-11T07:40:22.483199
- Duration: 208.8 min
- Run dir: `/Users/kronckbm/GitRepos/GeneVariantFetcher/results/grant_e2e_20260909/shards/LDLR_1/LDLR/20260909_104851`
- DB: `/Users/kronckbm/GitRepos/GeneVariantFetcher/results/grant_e2e_20260909/shards/LDLR_1/LDLR/20260909_104851/LDLR.db`
- Stage status: ✓ core stages ok; best-effort warnings recorded

## Best-effort Warnings

These quality/metadata stages did not remove core extracted evidence, so they do not change the process exit code:

  - ⚠️ source accounting disagreement: on-disk scan says 700 full text / 77 abstract-only, finalized extraction records say 602 / 0

## Paper Final Check

_Skipped: disabled (paper_final_check_enabled=false)_

Final-check findings remain recorded in `paper_final_check`; no actionable trust composition was reported for this run.

## Doctor

Required:
  - NCBI_EMAIL: ✓
  - LLM provider key: ✓
Recommended:
  - NCBI_API_KEY: ✓
  - ELSEVIER_API_KEY: ✓
  - WILEY_API_KEY: ✓
  - SPRINGER_API_KEY: ✓
LLM provider keys:
  - OPENAI_API_KEY: –
  - AZURE_AI_API_KEY: ✓
  - ANTHROPIC_API_KEY: ✓
Subscription unlocks (not required, but lift recall):
  - ELSEVIER_INSTTOKEN: ✓
  - GVF_EZPROXY_PREFIX: ✓
External storage:
  - corpus/: ✓ (linked)
- NCBI reachable: True

## Recovery Progression

| Layer | ClinVar added | PubTator added | Figures added |
|---|---:|---:|---:|
| 0_baseline | — | — | — |
| 1_clinvar | None | — | — |
| 2_pubtator | — | None | — |
| 3_figures | — | — | 21 |

_No gold standard CSV was available, so recovery ran without recall scoring._

## Source Acquisition QC (Gold-Free)

| Signal | PMIDs | Coverage |
|---|---:|---:|
| Usable full text now | 463/777 | 59.6% |
| Selected for fetch | 313/777 | 40.3% |
| Selected for supplement-only fetch | 6/777 | 0.8% |
| Selected for source refresh | 227/777 | 29.2% |
| Manual or blocked | 4/777 | 0.5% |
| Zero-variant usable full text | 229/777 | 29.5% |

- Worklist: `/Users/kronckbm/GitRepos/GeneVariantFetcher/results/grant_e2e_20260909/shards/LDLR_1/LDLR/20260909_104851/source_qc/source_acquisition_worklist.csv`
- Fetch queue: `/Users/kronckbm/GitRepos/GeneVariantFetcher/results/grant_e2e_20260909/shards/LDLR_1/LDLR/20260909_104851/source_qc/fetch_input.csv`
- Supplement-only queue: `/Users/kronckbm/GitRepos/GeneVariantFetcher/results/grant_e2e_20260909/shards/LDLR_1/LDLR/20260909_104851/source_qc/supplement_input.csv`
- Refresh source overrides: `/Users/kronckbm/GitRepos/GeneVariantFetcher/results/grant_e2e_20260909/shards/LDLR_1/LDLR/20260909_104851/source_qc/source_override.csv`

## Next steps
- Inspect per-layer outputs under `/Users/kronckbm/GitRepos/GeneVariantFetcher/results/grant_e2e_20260909/shards/LDLR_1/LDLR/20260909_104851/layers`.
- For no-gold source recovery, run fetch_paywalled.py on the source QC fetch queue, then refresh_run_db.py with `--source-override-csv` and `--stage-extractions`, or rerun `gvf-run` with `--source-recovery`.
- The scored DB is at `/Users/kronckbm/GitRepos/GeneVariantFetcher/results/grant_e2e_20260909/shards/LDLR_1/LDLR/20260909_104851/LDLR.db` for ad-hoc queries.
