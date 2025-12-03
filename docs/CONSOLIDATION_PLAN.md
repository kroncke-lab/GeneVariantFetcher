# Code Consolidation Plan

## Overview
This document identifies redundancies in the GeneVariantFetcher repository and provides a plan to consolidate duplicate code, reducing maintenance burden and improving code clarity.

## Identified Redundancies

### 1. Paper Sourcing Classes (HIGH PRIORITY)
**Location:**
- `sourcer.py` (root) - `PaperSourcer` class
- `pipeline/sourcing.py` - `PaperSourcer` class

**Issue:** Nearly identical implementations with minor differences:
- `sourcer.py` has simpler config handling
- `pipeline/sourcing.py` has more sophisticated config integration with `get_settings()`
- Both query PubMind, PubMed, and EuropePMC

**Consolidation Strategy:**
- Keep: `pipeline/sourcing.py` (more feature-complete, better config integration)
- Remove: `sourcer.py` (standalone file, less integrated)
- Update: All imports of `sourcer.PaperSourcer` → `pipeline.sourcing.PaperSourcer`

**Files to Update:**
- Check for imports of `from sourcer import` or `import sourcer`
- Update `automated_workflow.py` if it uses `sourcer.py`

---

### 2. Filter Modules (HIGH PRIORITY)
**Location:**
- `pipeline/filters.py` - Contains `KeywordFilter`, `InternFilter`, `ClinicalDataTriageFilter`
- `harvesting/filters.py` - Appears to be duplicate (needs verification)
- `gene_literature/clinical_data_triage.py` - CLI wrapper + uses `ClinicalDataTriageFilter` from `pipeline.filters`

**Issue:**
- `ClinicalDataTriageFilter` exists in both `pipeline/filters.py` and potentially `harvesting/filters.py`
- `gene_literature/clinical_data_triage.py` is a CLI wrapper that imports from `pipeline.filters`

**Consolidation Strategy:**
- Keep: `pipeline/filters.py` as the canonical filter implementation
- Keep: `gene_literature/clinical_data_triage.py` as CLI wrapper (it's a different purpose)
- Remove: `harvesting/filters.py` if it's truly duplicate
- Verify: Check if `harvesting/filters.py` has unique functionality

---

### 3. Retry Utilities (MEDIUM PRIORITY)
**Location:**
- `utils/retry_utils.py` - Uses `tenacity` library, comprehensive retry decorators
- `pipeline/utils/retry_utils.py` - Custom implementation with `functools.wraps`, simpler

**Issue:**
- Two different implementations of retry logic
- `utils/retry_utils.py` is more feature-rich (uses tenacity)
- `pipeline/utils/retry_utils.py` is simpler but less flexible

**Consolidation Strategy:**
- Keep: `utils/retry_utils.py` (more robust, uses industry-standard library)
- Remove: `pipeline/utils/retry_utils.py`
- Update: All imports from `pipeline.utils.retry_utils` → `utils.retry_utils`
- Note: May need to update function names if they differ (`llm_retry` vs custom decorator)

---

### 4. Example Script Wrapper (LOW PRIORITY)
**Location:**
- `tests/example_automated_workflow.py` - Just calls `automated_workflow.main()`

**Issue:**
- Trivial wrapper with no added value
- According to CLAUDE.md, example scripts should not be used in production

**Consolidation Strategy:**
- Keep: For now (documentation/example purposes)
- Consider: Adding deprecation warning or removing if truly unused
- Note: Low priority - minimal impact

---

### 5. Rerun Scripts Overlap (MEDIUM PRIORITY)
**Location:**
- `rerun_extraction.py` - Re-runs extraction on existing markdown files
- `rerun_aggregate_helper.py` - Aggregates penetrance from rerun extractions
- `automated_workflow.py` - Contains extraction and aggregation logic

**Issue:**
- `rerun_extraction.py` duplicates extraction logic from `automated_workflow.py`
- `rerun_aggregate_helper.py` duplicates aggregation logic
- Both scripts have similar parallel processing code

**Consolidation Strategy:**
- Extract: Common extraction logic into `pipeline/extraction.py` (already exists as `ExpertExtractor`)
- Extract: Common aggregation logic into `pipeline/aggregation.py` (already exists)
- Refactor: `rerun_extraction.py` to use shared functions from pipeline modules
- Refactor: `rerun_aggregate_helper.py` to use shared aggregation function
- Keep: Scripts as CLI entry points, but make them thin wrappers

---

### 6. Literature Collection (LOW PRIORITY)
**Location:**
- `gene_literature/collector.py` - `LiteratureCollector` class
- `gene_literature/collect_literature.py` - CLI wrapper for `LiteratureCollector`

**Issue:**
- `collect_literature.py` is just a CLI wrapper (appropriate separation)
- However, `collector.py` may overlap with `pipeline/sourcing.py` functionality

**Consolidation Strategy:**
- Keep: Both files (appropriate separation of concerns)
- Verify: If `LiteratureCollector` and `PaperSourcer` have overlapping functionality
- Consider: Whether `LiteratureCollector` should use `PaperSourcer` internally

---

### 7. PubMind Fetcher Functions (LOW PRIORITY)
**Location:**
- `gene_literature/pubmind_fetcher.py` - `PubMindFetcher` class + `fetch_pmids_for_gene()` function
- `pipeline/sourcing.py` - Uses `PubMindFetcher` internally

**Issue:**
- `fetch_pmids_for_gene()` is a standalone function that wraps `PubMindFetcher`
- May be redundant if `PaperSourcer` already handles this

**Consolidation Strategy:**
- Keep: Both (function is convenience wrapper, class is full implementation)
- Verify: If standalone function is used elsewhere

---

## Consolidation Priority Order

1. **HIGH PRIORITY** (Do First):
   - Consolidate `sourcer.py` → `pipeline/sourcing.py`
   - Verify and consolidate filter modules

2. **MEDIUM PRIORITY** (Do Second):
   - Consolidate retry utilities
   - Refactor rerun scripts to use shared pipeline modules

3. **LOW PRIORITY** (Do Third):
   - Review example scripts
   - Verify literature collection overlap

---

## Implementation Steps

### Step 1: Paper Sourcing Consolidation
1. Search for all imports of `sourcer`:
   ```bash
   grep -rE 'from sourcer|import sourcer|sourcer\.' .
   ```
2. Update all imports to use `pipeline.sourcing`
3. Delete `sourcer.py`
4. Test that `automated_workflow.py` still works

### Step 2: Filter Module Consolidation
1. Read `harvesting/filters.py` to verify if it's duplicate
2. If duplicate, delete `harvesting/filters.py`
3. Verify `gene_literature/clinical_data_triage.py` correctly imports from `pipeline.filters`
4. Test filter functionality

### Step 3: Retry Utilities Consolidation
1. Search for imports of `pipeline.utils.retry_utils`:
   ```bash
   grep -r "from pipeline.utils.retry_utils\|import pipeline.utils.retry_utils" .
   ```
2. Update imports to `utils.retry_utils`
3. Verify function names match (may need adapter)
4. Delete `pipeline/utils/retry_utils.py`
5. Test retry functionality

### Step 4: Rerun Scripts Refactoring
1. Extract common extraction logic from `rerun_extraction.py` into shared function
2. Extract common aggregation logic from `rerun_aggregate_helper.py` into shared function
3. Refactor scripts to call shared functions
4. Test rerun functionality

---

## Testing Checklist

After each consolidation step:
- [ ] Run `automated_workflow.py` with a test gene
- [ ] Run `rerun_extraction.py` on existing data
- [ ] Run `rerun_aggregate_helper.py` on existing extractions
- [ ] Verify filter functionality works
- [ ] Check that all imports resolve correctly
- [ ] Run any existing tests

---

## Risk Assessment

**Low Risk:**
- Removing `sourcer.py` (if `pipeline/sourcing.py` is feature-complete)
- Consolidating retry utilities (if function signatures match)

**Medium Risk:**
- Filter consolidation (need to verify no unique functionality)
- Rerun script refactoring (need to ensure backward compatibility)

**Mitigation:**
- Keep backup of files before deletion
- Test thoroughly after each step
- Use git branches for each consolidation step

---

## Expected Benefits

1. **Reduced Maintenance:** Fewer files to maintain
2. **Clearer Architecture:** Single source of truth for each component
3. **Easier Testing:** Centralized logic easier to test
4. **Better Code Reuse:** Shared functions reduce duplication
5. **Improved Documentation:** Less confusion about which module to use

---

## Notes

- Some "redundancies" may be intentional (e.g., CLI wrappers vs. classes)
- Always verify that "duplicate" code doesn't have unique functionality
- Consider backward compatibility for any public APIs
- Update documentation (CLAUDE.md, README.md) after consolidation
