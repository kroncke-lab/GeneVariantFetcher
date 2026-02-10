# Changelog

All notable changes to GeneVariantFetcher will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
