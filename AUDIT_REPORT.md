# GeneVariantFetcher - Pipeline Architecture Audit Report

**Date:** 2026-01-14
**Branch:** claude/audit-pipeline-architecture-kSf9k

---

## Executive Summary

This codebase is a sophisticated 6-step pipeline for extracting genetic variant data from biomedical literature. The primary interface is a **FastAPI-based web GUI** (not a desktop GUI). The code is generally well-architected but has notable areas of bloat, duplication, and some broken/unused modules.

**Key Findings:**
- ~12,000 lines of core Python code across `pipeline/`, `harvesting/`, `gene_literature/`
- **1 broken file** that should be fixed or removed
- **2 redundant entry points** that should be consolidated
- **3 modules** with limited/no real usage in production workflow
- **1 minimal CLI stub** that adds little value

---

## 1. Actual Pipeline Architecture

### How It Runs

The application runs as a **web-based GUI** using FastAPI + WebSockets:

```bash
python main.py          # Launches FastAPI server at localhost:8000
```

This opens a browser-based interface that provides:
- Visual pipeline configuration
- Real-time progress monitoring via WebSocket
- Job management (pause/resume/stop)
- Settings configuration

### Pipeline Steps (6 Total)

| Step | Module | Description |
|------|--------|-------------|
| **0** | `gene_literature/synonym_finder.py` | Optional: discover gene synonyms |
| **1** | `gene_literature/discovery.py` | Fetch PMIDs from PubMind/PubMed/EuropePMC |
| **1.5** | `harvesting/abstracts.py` | Fetch paper metadata & abstracts |
| **1.6** | `pipeline/filters.py` | Tier 1+2 filtering before download |
| **2** | `harvesting/orchestrator.py` | Download full-text + supplements |
| **2.5** | (internal) | Identify abstract-only extraction candidates |
| **3** | `pipeline/extraction.py` | LLM extraction (with Data Scout) |
| **4** | `pipeline/aggregation.py` | Penetrance validation |
| **5** | `harvesting/migrate_to_sqlite.py` | JSON → SQLite migration |

### Filter Tiers

| Tier | Class | Purpose | Model |
|------|-------|---------|-------|
| 1 | `KeywordFilter` | Fast regex-based filtering | None |
| 2 | `InternFilter` | LLM classification | gpt-4o-mini |
| 2 (alt) | `ClinicalDataTriageFilter` | Clinical data triage | gpt-4o-mini |
| 3 | `ExpertExtractor` | Full extraction | gpt-4o-mini → gpt-4o cascade |

---

## 2. Identified Issues

### 2.1 BROKEN FILE: `gene_literature/relevance_checker.py`

**Severity: HIGH** - This file has multiple syntax/logic errors and is non-functional.

```python
# Line 7 - Wrong import
from dataclasses import dataclasses  # Should be: from dataclasses import dataclass

# Line 13 - Wrong decorator
@dataclasses  # Should be: @dataclass

# Lines 27-30, 50-54, 167-170 - Broken type hints
api_key: Optional[check_relevance] = None  # Should be: Optional[str]
model: check_relevance = "..."  # Should be: str

# Lines 136-142 - Completely broken code
import json_match  # Not a real module
import result      # Not a real module
json_match = result.search(...)  # Variable shadowing import
```

**Recommendation:** Delete this file or fix it completely. It is only imported by `collect_literature.py` which itself is barely used.

---

### 2.2 REDUNDANT ENTRY POINTS

**File 1:** `main.py` (213 lines)
**File 2:** `run_gui.py` (173 lines)

Both files do essentially the same thing:
1. Check GUI dependencies (fastapi, uvicorn, websockets)
2. Load environment variables
3. Check for incomplete jobs
4. Launch uvicorn server

**Duplication:** ~80% code overlap

**Recommendation:** Delete `run_gui.py` entirely. `main.py` is the canonical entry point and is already documented as such.

---

### 2.3 MINIMAL/STUB CLI: `pipeline/cli.py`

This 46-line file only provides a single command:

```python
@app.command()
def version() -> None:
    """Print the installed package version."""
    # ...
```

**Also has a bug** on line 16:
```python
def _get_version() -> Optional[version]:  # Should be: Optional[str]
```

**Recommendation:** Either expand this to be a real CLI alternative to the GUI, or remove it and document `python main.py --cli` as the CLI interface.

---

### 2.4 UNUSED/LOW-VALUE MODULES

#### A. `gene_literature/collect_literature.py` (159 lines)

Standalone CLI tool that duplicates discovery functionality already in `discovery.py`. Only imports the broken `relevance_checker.py`.

**Usage:** Zero imports from any other module
**Recommendation:** Delete or mark as deprecated

#### B. `gene_literature/relevance_checker.py` (187 lines)

As noted above, this file is broken. Uses Anthropic API (not OpenAI like rest of codebase) and has syntax errors.

**Usage:** Only imported by `collect_literature.py` (which itself is unused)
**Recommendation:** Delete

#### C. `gene_literature/synonym_relevance_checker.py`

Only imported by `synonym_finder.py` for optional LLM-based synonym filtering.

**Status:** Optional feature, works correctly
**Recommendation:** Keep (it's properly integrated)

---

### 2.5 DUPLICATE FUNCTION DEFINITIONS

**File:** `gene_literature/pubmind_fetcher.py`

Contains both:
- `PubMindFetcher.fetch_pmids_for_gene()` - class method (line 134)
- `fetch_pmids_for_gene()` - module-level function (line 329)

The module-level function just wraps the class method. This pattern is intentional for convenience but creates redundancy.

**Recommendation:** Keep as-is (this is a valid convenience pattern)

---

## 3. Code Bloat Analysis

### 3.1 Verbose Files (Top 5 by Lines)

| File | Lines | Notes |
|------|-------|-------|
| `tests/compare_variants.py` | 1,450 | Testing utility - OK |
| `gui/server.py` | 1,222 | Main server - complex but necessary |
| `harvesting/migrate_to_sqlite.py` | 1,115 | SQLite migration - complex but necessary |
| `harvesting/wiley_api.py` | 1,003 | Publisher API - conditional feature |
| `harvesting/orchestrator.py` | 949 | Download coordinator - necessary |

**Verdict:** Most large files are justified by their complexity. The test utility (`compare_variants.py`) could be refactored but isn't blocking.

### 3.2 Overly Verbose Patterns

**Automated Workflow Logging:** `automated_workflow.py` uses heavy emoji-decorated logging:
```python
logger.info("📚 STEP 1: Discovering relevant papers...")
logger.info("✓ Found %d PubMind PMIDs, %d PubMed PMIDs...")
```

This is fine for user visibility but adds clutter to the code.

**Excessive Step Numbering:** The workflow uses sub-steps (1.5, 1.6, 2.5) which makes the "6-step pipeline" claim in documentation misleading. It's really ~9 discrete operations.

---

## 4. Proposed Fixes

### Priority 1: Delete Broken/Unused Files

```bash
# Files to delete
rm gene_literature/relevance_checker.py     # Broken, unused
rm gene_literature/collect_literature.py    # Unused, imports broken file
rm run_gui.py                               # Redundant with main.py
```

### Priority 2: Fix Type Hint Bug

**File:** `pipeline/cli.py`, line 16

```python
# Before
def _get_version() -> Optional[version]:

# After
def _get_version() -> Optional[str]:
```

### Priority 3: Update Documentation References

Remove mentions of `run_gui.py` from:
- `gui/README.md`
- Any other documentation

### Priority 4: Consider Consolidating CLI

Either:
1. **Expand `pipeline/cli.py`** to include full workflow commands (replacing `python main.py --cli`)
2. **Or delete `pipeline/cli.py`** and keep `main.py --cli` as the canonical CLI interface

---

## 5. README Revision Recommendations

The current README is generally good but could be improved:

### Current Issues

1. **Misleading step count:** Says "6 steps" but actual workflow has ~9 operations
2. **Missing `run_gui.py` deprecation:** If we delete it, remove any references
3. **Inconsistent entry points:** Documents both `main.py` and `automated_workflow.py`

### Proposed README Changes

```markdown
## Quick Start

```bash
# Install
pip install -e .
pip install -r gui/requirements.txt

# Set required environment variables (or use .env file)
export OPENAI_API_KEY="your-key"
export NCBI_EMAIL="your@email.com"

# Launch the GUI (recommended)
python main.py
```

## Pipeline Overview

The pipeline executes these stages:

1. **Discovery** - Fetch PMIDs from PubMind/PubMed/EuropePMC
2. **Metadata** - Fetch abstracts and paper metadata
3. **Filtering** - Tier 1 (keywords) + Tier 2 (LLM triage)
4. **Download** - Full-text articles + supplemental materials
5. **Extraction** - LLM-based variant/patient data extraction
6. **Aggregation** - Penetrance validation
7. **SQLite** - Normalize to database

## Running Modes

| Mode | Command | Use Case |
|------|---------|----------|
| **GUI** | `python main.py` | Interactive use (recommended) |
| **CLI** | `python main.py --cli GENE --email EMAIL` | Scripting/automation |
```

### Remove from README

- References to `run_gui.py` (redundant)
- The separate `automated_workflow.py` examples (consolidate to `main.py --cli`)

---

## 6. Summary Table

| Category | Item | Action |
|----------|------|--------|
| **Broken** | `relevance_checker.py` | DELETE |
| **Unused** | `collect_literature.py` | DELETE |
| **Redundant** | `run_gui.py` | DELETE |
| **Bug** | `pipeline/cli.py` line 16 | FIX type hint |
| **Docs** | README.md | UPDATE per recommendations |
| **Minimal** | `pipeline/cli.py` | EXPAND or DELETE |

---

## 7. Files Safe to Keep (Not Bloat)

These files may look unused but are legitimate optional features:

- `harvesting/wiley_api.py` - Conditional publisher integration
- `harvesting/elsevier_api.py` - Conditional publisher integration
- `gene_literature/synonym_finder.py` - Used when `--auto-synonyms` flag set
- `gene_literature/synonym_relevance_checker.py` - Used by synonym_finder
- `pipeline/abstract_extraction.py` - Integrated into main extraction pipeline
- All `tests/` files - Testing infrastructure

---

## Appendix: File Inventory

### Core Pipeline (~12,000 lines)
- `pipeline/` - 7 files, ~4,000 lines
- `harvesting/` - 10 files, ~6,000 lines
- `gene_literature/` - 9 files, ~3,500 lines (includes 2 to delete)

### Supporting Code (~4,500 lines)
- `utils/` - 9 files, ~2,000 lines
- `gui/` - 3 files, ~2,500 lines

### Entry Points (~1,200 lines)
- `main.py` - 213 lines (primary)
- `automated_workflow.py` - 830 lines (CLI backend)
- `run_gui.py` - 173 lines (DELETE)

### Tests (~7,000 lines)
- `tests/` - 15 files

**Total: ~55 Python files, ~25,000 lines**
