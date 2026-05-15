# GeneVariantFetcher Handoff

GVF extracts genetic variants, carrier counts, and phenotype data from
biomedical literature for the Kroncke Lab variant interpretation pipeline. The
grant target is 90% unique-variant recall by June 2026.

## Current Measured State

Latest scored KCNH2 run: 2026-05-15 manual-recovery v12.

| Metric | Matched / gold | Recall |
|---|---:|---:|
| PMIDs | 184 / 262 | 70.2% |
| Variant rows | 542 / 991 | 54.7% |
| Unique variants | 323 / 530 | 60.9% |
| Patients/carriers | 1758 / 2674 | 65.7% |
| Affected | 1095 / 1635 | 67.0% |
| Unaffected | 461 / 749 | 61.5% |

Local artifacts:
- DB: `results/KCNH2/20260506_102238/end_to_end_20260515_manual_recovery/KCNH2_v12_manual_recovery_20260515.db`
- Recall output: `recall_metrics/kcnh2_manual_recovery_after_matchfix_20260515/`
- Gold input: `gene_variant_fetcher_gold_standard/normalized/KCNH2_recall_input.csv`

The old 59.1% KCNH2 number from 2025-12-11 is now historical only. Use the
2026-05-15 score above as the current local baseline.

## What Changed In This Branch

- Source recovery stack: authenticated Playwright recovery in
  `scripts/fetch_paywalled.py`, publisher strategies under
  `harvesting/browser_html/`, quality gating, per-domain pacing, and PMC
  fallback.
- Clinical table extraction: `pipeline/table_router.py` and
  `pipeline/extraction.py` preserve clinical mutation-list tables and infer one
  carrier per clinical row when there is patient/proband/family context but no
  explicit count column.
- Recall runner: `scripts/run_recall_suite.py` scores normalized per-gene
  recall inputs from `gene_variant_fetcher_gold_standard/normalized/`.
- Matcher: `cli/compare_variants.py` has positional-digit guards, greedy
  one-to-one assignment, cDNA/protein bridging, and 2026-05-15 fixes for
  frameshift spellings such as `fsTer`, `fs/185`, `fs+*49`, and malformed
  `AlaX14`.
- Gold-standard builders: `scripts/build_gold_standard_from_varbrowser.py` and
  `scripts/build_ryr2_gold_standard_from_xlsx.py`.

## Manual Acquisition Queue

The manual acquisition queue is the set of passed-filter PMIDs that do not have
usable assembled full text after the automatic harvesters and fallback routes.
They are not "assigned to a human" by code; they are papers GVF could not fetch
with the currently available credentials and network position.

Current local KCNH2 status:
- 170 PMIDs were in the manual-queue subset.
- 18 were recovered and integrated into the v12 DB.
- 152 still lack usable full context.
- Current missing-in-SQLite gap after the matcher patch: 449 variant rows across
  107 PMIDs. The top 10 PMIDs account for 227 rows.

Top current KCNH2 losses:

| PMID | Missing rows | Main blocker |
|---|---:|---|
| 15840476 | 86 | Elsevier/Heart Rhythm subscription full text |
| 14661677 | 24 | Source/extraction coverage |
| 29650123 | 21 | Full-ish context, extraction miss |
| 24667783 | 20 | Supplement/table detail not fully present |
| 16922724 | 20 | Wiley subscription/Cloudflare |
| 23098067 | 16 | Open access body, extraction/table miss |
| 23631430 | 12 | Sage/Liebert Cloudflare |
| 17905336 | 10 | Extraction miss |
| 12402336 | 9 | Wiley Human Mutation |
| 27871843 | 9 | Source/extraction coverage |

## Highest ROI Blocker

PMID 15840476 is still the largest single blocker. Elsevier full-text retrieval
needs both headers:

```text
X-ELS-APIKey:    <ELSEVIER_API_KEY>
X-ELS-Insttoken: <ELSEVIER_INSTTOKEN>
```

`.env` already supports `ELSEVIER_INSTTOKEN` through the Elsevier API client.
The API key alone returns metadata/abstract only. On a Vanderbilt VPN or
institutional machine, first try `scripts/fetch_paywalled.py`; if the API still
does not unlock the article, manually download a real PDF in the browser and
feed it through `harvesting/format_converters.py`. Do not use OneDrive
on-demand placeholders; they were checked on 2026-05-14 and were invalid PDFs.

## Run On Another Computer Or VPN

From a clean checkout:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
pip install -r gui/requirements.txt
python -m playwright install chromium
```

Create `.env` with at least:

```bash
OPENAI_API_KEY=...
NCBI_EMAIL=...
NCBI_API_KEY=...
ELSEVIER_API_KEY=...
SPRINGER_API_KEY=...
WILEY_API_KEY=...
# Optional but highest value:
ELSEVIER_INSTTOKEN=...
```

Default local tests avoid live network calls:

```bash
.venv/bin/python -m pytest tests/ -q
```

Live network/institutional checks are marked separately:

```bash
GVF_TEST_OUTPUT_DIR=/tmp/gvf_tests .venv/bin/python -m pytest -m requires_network tests/integration -q
```

Recover queued or paywalled papers after VPN/login access is available:

```bash
.venv/bin/python scripts/fetch_paywalled.py \
  --input results/KCNH2/20260506_102238/pmc_fulltext/paywalled_missing.csv \
  --output results/KCNH2/20260506_102238/pmc_fulltext \
  --no-headless
```

Then rerun extraction/migration and score:

```bash
.venv/bin/python scripts/run_recall_suite.py --score --genes KCNH2 \
  --db KCNH2=results/KCNH2/20260506_102238/end_to_end_20260515_manual_recovery/KCNH2_v12_manual_recovery_20260515.db \
  --outdir recall_metrics/kcnh2_after_vpn_recovery
```

## Active Work

1. Close KCNH2 source/extraction coverage to 90% recall. The current
   variant-row gap to 90% is 350 rows.
2. Unblock PMID 15840476 with `ELSEVIER_INSTTOKEN`, VPN/IP access, or a valid
   manual PDF.
3. Re-run focused extraction for top missing-but-present sources, especially
   29650123, 24667783, 23098067, and 17905336.
4. Investigate `count_mismatches=117`; many top mismatches look like SQLite
   over-counting carrier totals from cohort tables, separate from recall.
5. Run fresh KCNQ1 and SCN5A extraction DBs so the multi-gene recall runner can
   score them locally. RYR2 has normalized gold input but still needs source
   reconciliation before treating per-PMID recall as final.

## Files To Know

- `cli/compare_variants.py` - gold-standard matcher and recall summary.
- `scripts/run_recall_suite.py` - multi-gene recall scoring.
- `scripts/fetch_paywalled.py` - canonical paywall recovery entry point.
- `harvesting/browser_html/` - authenticated browser recovery strategies.
- `harvesting/browser_html/content_quality.py` - paywall/stub quality gate.
- `pipeline/table_router.py` and `pipeline/extraction.py` - clinical table
  extraction.
- `gene_variant_fetcher_gold_standard/` - normalized recall inputs and source
  exports.

## House Rules

- Use the project `.venv` when available.
- Do not commit agent scratch files: `HEARTBEAT.md`, `IDENTITY.md`, `SOUL.md`,
  `TOOLS.md`, `USER.md`, or `.openclaw/`.
- Do not commit `.env`, local `results/`, SQLite DBs, or generated
  `recall_metrics/`.
- Pre-commit hooks may reformat staged files. Never use `--no-verify`; re-stage
  after hooks run.
