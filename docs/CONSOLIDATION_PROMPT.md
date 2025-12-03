# Prompt: Consolidate Redundant Code in GeneVariantFetcher

## Context
You are working on the GeneVariantFetcher repository, a biomedical literature extraction pipeline. A consolidation plan has been created at `docs/CONSOLIDATION_PLAN.md` identifying redundant code that should be consolidated.

## Your Task
Execute the consolidation plan step-by-step, starting with HIGH PRIORITY items, then MEDIUM, then LOW priority.

## Instructions

### Phase 1: HIGH PRIORITY Consolidations

#### Task 1.1: Consolidate Paper Sourcing Classes
**Goal:** Remove `sourcer.py` and migrate all usage to `pipeline/sourcing.py`

**Steps:**
1. Search for all imports/references to `sourcer.py`:
   ```bash
   grep -r "from sourcer\|import sourcer\|sourcer\." --include="*.py" .
   ```

2. Compare `sourcer.py` and `pipeline/sourcing.py`:
   - Identify any unique functionality in `sourcer.py` that's missing in `pipeline/sourcing.py`
   - If unique functionality exists, merge it into `pipeline/sourcing.py`
   - Ensure `pipeline/sourcing.py` has all features from `sourcer.py`

3. Update all imports:
   - Change `from sourcer import PaperSourcer` → `from pipeline.sourcing import PaperSourcer`
   - Change `from sourcer import query_papers_for_gene` → `from pipeline.sourcing import query_papers_for_gene`
   - Update any other imports as needed

4. Delete `sourcer.py` after confirming all imports are updated

5. Test:
   - Run `python automated_workflow.py BRCA1 --email test@example.com --max-pmids 5 --max-downloads 2` (or similar test)
   - Verify no import errors occur

#### Task 1.2: Consolidate Filter Modules
**Goal:** Verify filter modules and remove any true duplicates

**Steps:**
1. Read `harvesting/filters.py` and compare with `pipeline/filters.py`:
   - Are they identical or do they have unique functionality?
   - If identical, delete `harvesting/filters.py`
   - If different, document the differences

2. Verify `gene_literature/clinical_data_triage.py`:
   - Ensure it correctly imports `ClinicalDataTriageFilter` from `pipeline.filters`
   - If it imports from elsewhere, update the import

3. Search for any other imports of filter classes:
   ```bash
   grep -rE 'from.*filters import' --include="*.py" .
   ```

4. Test filter functionality:
   - Run a test that uses filters (if test exists)
   - Verify no import errors

### Phase 2: MEDIUM PRIORITY Consolidations

#### Task 2.1: Consolidate Retry Utilities
**Goal:** Remove `pipeline/utils/retry_utils.py` and migrate to `utils/retry_utils.py`

**Steps:**
1. Compare the two retry utility files:
   - `utils/retry_utils.py` uses `tenacity` (more robust)
   - `pipeline/utils/retry_utils.py` uses custom implementation
   - Identify function name differences

2. Search for imports:
   ```bash
   grep -r "from pipeline.utils.retry_utils\|import pipeline.utils.retry_utils" --include="*.py" .
   ```

3. Update imports:
   - Change to `from utils.retry_utils import ...`
   - If function names differ, create adapter functions or update call sites

4. Delete `pipeline/utils/retry_utils.py`

5. Test:
   - Verify retry functionality still works
   - Check for any runtime errors

#### Task 2.2: Refactor Rerun Scripts
**Goal:** Make rerun scripts use shared pipeline modules instead of duplicating logic

**Steps:**
1. Analyze `rerun_extraction.py`:
   - Identify code that duplicates `automated_workflow.py` extraction logic
   - Extract common extraction logic into a shared function (or use existing `ExpertExtractor`)

2. Analyze `rerun_aggregate_helper.py`:
   - Verify it uses `pipeline.aggregation.aggregate_penetrance` (it should already)
   - If it duplicates aggregation logic, remove duplication

3. Refactor `rerun_extraction.py`:
   - Use `pipeline.extraction.ExpertExtractor` directly
   - Use `pipeline.aggregation.aggregate_penetrance` directly
   - Keep only CLI-specific logic in the script

4. Test:
   - Run `rerun_extraction.py` on existing data
   - Run `rerun_aggregate_helper.py` on existing extractions
   - Verify outputs are correct

### Phase 3: LOW PRIORITY Consolidations

#### Task 3.1: Review Example Scripts
**Goal:** Document or remove truly redundant example scripts

**Steps:**
1. Review `tests/example_automated_workflow.py`:
   - It's a thin wrapper - consider adding deprecation notice or removing
   - Check if it's referenced in documentation

2. Review other example scripts in `tests/`:
   - Identify which are truly redundant vs. useful examples
   - Document findings

#### Task 3.2: Verify Literature Collection Overlap
**Goal:** Ensure `LiteratureCollector` and `PaperSourcer` don't have unnecessary overlap

**Steps:**
1. Compare `gene_literature/collector.py` and `pipeline/sourcing.py`:
   - Do they serve different purposes?
   - Could `LiteratureCollector` use `PaperSourcer` internally?

2. Document findings (no changes needed if separation is appropriate)

## General Guidelines

1. **Before Each Change:**
   - Read the relevant files completely
   - Understand the differences between "duplicate" code
   - Verify no unique functionality will be lost

2. **During Changes:**
   - Make one consolidation at a time
   - Test after each change
   - Commit changes incrementally (if using git)

3. **After Changes:**
   - Run basic smoke tests
   - Verify imports work
   - Check for any broken functionality

4. **If Uncertain:**
   - Document the uncertainty
   - Keep the code if it's not clearly redundant
   - Add comments explaining why both exist

## Testing Commands

After consolidations, run these to verify:

```bash
# Test main workflow
python automated_workflow.py BRCA1 --email test@example.com --max-pmids 5 --max-downloads 2

# Test rerun scripts (if you have existing data)
python rerun_extraction.py automated_output/BRCA1/YYYYMMDD_HHMMSS BRCA1

# Check for import errors
python -c "from pipeline.sourcing import PaperSourcer; from pipeline.filters import ClinicalDataTriageFilter; from utils.retry_utils import standard_retry; print('Imports OK')"
```

## Success Criteria

- [ ] All HIGH PRIORITY consolidations complete
- [ ] All MEDIUM PRIORITY consolidations complete
- [ ] No import errors
- [ ] Main workflow still works
- [ ] Rerun scripts still work
- [ ] Documentation updated if needed

## Notes

- Some files may appear redundant but serve different purposes (e.g., CLI wrappers)
- Always verify functionality before deleting code
- Keep backups or use version control
- Update this prompt file with any findings or deviations from the plan
