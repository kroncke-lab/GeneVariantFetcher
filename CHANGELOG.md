# Changelog

All notable changes to GeneVariantFetcher will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - 2026-05-21

### Added
- `scripts/test_insttoken_unlock.py` probes paywalled Elsevier DOIs through
  `ElsevierAPIClient` with the institutional token and writes unlocked bodies
  as `{PMID}_FULL_CONTEXT.md` into each run's existing `pmc_fulltext/`
  directory. Pre-existing stub files preserved as `*.pre_insttoken_bak`.

### Changed
- `ELSEVIER_INSTTOKEN` is now active in `.env` (issued 2026-05-21 by Elsevier
  Data Support against the Vanderbilt-registered API key). Unlock probe across
  KCNH2/KCNQ1/RYR2/SCN5A/KCNE1 paywalled lists: **242/246 (98.4%) full text**.
- `docs/RECALL_STATUS.md` adds a "2026-05-21 Elsevier Insttoken Activation"
  section with the probe table, file-landing inventory, and a measured
  pre-re-extraction PMID-recall baseline (KCNH2 87.8%, KCNQ1 80.3%, RYR2 68.0%,
  SCN5A 70.9%, aggregate 75.4%). The "Next Run Plan" now leads with
  re-extraction rather than acquisition.
- `docs/vanderbilt_api_access.md`, `CLAUDE.md`, and `TASKS.md` updated to
  reflect the insttoken activation and the new highest-ROI step (re-extraction).

## [Pre-insttoken] - 2026-05-20

### Added
- `gvf gvf-run` turnkey driver for doctor checks, extraction, deterministic
  recovery layers, and recall/report handoff.
- `gvf audit-paywalls` for acquisition blocker review.
- Deterministic recall-recovery layers under `scripts/recall_recovery/` for
  ClinVar, PubTator, paywall/OA recovery, and historical KCNH2 v12 merging.
- Current recall source of truth:
  `docs/RECALL_STATUS.md`.
- Recall-audit scripts under `scripts/recall_audit/` for PMID side-by-side
  review, failure taxonomy, gene-leakage checks, study-wide N detection,
  multi-cohort collapse review, PMID replay, and acquisition status summaries.

### Changed
- Browser/Playwright paywall recovery is now part of the documented acquisition
  stack.
- Elsevier documentation now distinguishes `ELSEVIER_API_KEY` from the
  institutional `ELSEVIER_INSTTOKEN`, the current highest-yield missing
  credential.
- Recall docs now defer live metrics to `docs/RECALL_STATUS.md` instead of
  repeating potentially stale extraction-statistics artifacts.

### Fixed
- Table-router parsing now drops rows that explicitly name an off-target gene
  when the router forgot to map the gene column.
- Extraction artifact filtering has test coverage for bare numeric variant
  names such as `378`.

## [2.1.0] - 2026-02-10

### Added
- Comprehensive variant scanner (regex across all text, not just tables)
- KCNH2 variant alias dictionary (4,766 aliases → 856 canonical variants)
- PDF table extraction from supplements (pdfplumber + openpyxl)
- Springer HTML downloader (bypasses OA-only API via institutional proxy)
- Fuzzy variant matching (frameshift, nonsense, deletion, duplication, position ±1)
- `--scout-first` CLI flag for improved extraction context
- FULL_CONTEXT rebuild script (intelligent source selection with quality scoring)
- Artifact filtering (rejects p.XXX, p.unknown, position validation)
- Broader table parser headers (cDNA, protein, HGVS, AAChange, etc.)

### Changed
- Consolidated variant normalizers into single canonical module
- Excel processing limit increased from 100 to 10,000 rows
- Table regex re-enabled with normalizer-based validation
- Aggregation now normalizes variants before grouping

### Fixed
- Excel truncation bug (silently dropping data after 100 rows)
- Springer API returning metadata-only JATS files
- FULL_CONTEXT.md garbled/truncated from poor source selection
- Extraction JSON parsing (protein_notation field mismatch)
- SyntaxError in vandy_access.py
- Competing normalizer imports across pipeline

### Removed
- `build/` artifacts from git tracking
- Dead test scripts (test_playwright_v*, etc.)
- Stray files (=, 1, egg-info)
- Duplicate normalizer (improved_variant_normalizer.py)

### Performance
- Unique variant recall: 19.4% → 59.1% (+205%)
- Carrier recall: 71.3%
- Affected individual recall: 67.6%

## [2.0.0] - 2026-01-28

### Added
- SQLite-based paper storage (replacing flat file outputs)
- Multi-format harvesting (PDF, DOCX, XLSX, HTML, XML)
- Publisher API integrations (Elsevier, Springer, Wiley, PMC)
- LLM-powered extraction pipeline with structured outputs
- Pydantic models for variant data validation
- CLI with Typer framework
- FastAPI web interface for browsing results

### Changed
- Complete architecture rewrite from single-script to modular pipeline
- Moved from regex-only extraction to hybrid LLM + regex approach

## [1.0.0] - 2025-12-01

### Added
- Initial release
- Basic PubMed search and paper discovery
- Simple regex-based variant extraction
- CSV output format
