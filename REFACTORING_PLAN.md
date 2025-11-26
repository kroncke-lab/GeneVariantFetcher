# GeneVariantFetcher Refactoring Plan: Simplified Approach

## Executive Summary

This plan implements three key improvements to the automated variant extraction workflow:
1. **Graceful Degradation**: Make workflows resilient to missing API keys
2. **Parallelism Controls**: Expose and configure concurrent processing
3. **SQLite Integration**: Unified workflow including database migration

## Design Principles

- **Simplicity over Abstraction**: Leverage existing LiteLLM abstraction instead of building new layers
- **Backwards Compatibility**: All new features are opt-in via flags
- **Progressive Enhancement**: Each phase adds value independently
- **YAGNI**: Only build what's immediately needed

---

## Phase 1: Graceful Degradation (Simplified Abstraction)

### Problem
Current state: Hard exit on missing API key prevents harvest-only workflows

```python
# automated_workflow.py:332-337
if not (os.getenv("AI_INTEGRATIONS_OPENAI_API_KEY") or os.getenv("OPENAI_API_KEY")):
    logger.error("⚠️  ERROR: OpenAI API key not found!")
    sys.exit(1)  # ❌ Hard exit
```

### Solution
Soft warning + skip flag instead of full abstraction layer

#### Changes to `automated_workflow.py`:

1. **Replace hard API key check** (lines 332-337):
```python
# Check for API keys (soft warning with escape hatch)
has_api_key = bool(os.getenv("AI_INTEGRATIONS_OPENAI_API_KEY") or os.getenv("OPENAI_API_KEY"))

if not has_api_key:
    if args.skip_extraction:
        logger.warning("⚠️  No OpenAI API key found - extraction will be skipped")
    else:
        logger.error("⚠️  ERROR: OpenAI API key not found!")
        logger.error("Please set AI_INTEGRATIONS_OPENAI_API_KEY or OPENAI_API_KEY")
        logger.error("Or use --skip-extraction to run harvest-only workflow")
        sys.exit(1)
```

2. **Add CLI flags**:
```python
parser.add_argument("--skip-extraction", action="store_true",
                   help="Skip LLM extraction (harvest papers only)")
parser.add_argument("--model", type=str, default=None,
                   help="LLM model to use (e.g., 'gpt-4o', 'claude-3-opus', 'anthropic/claude-3-5-sonnet-20241022')")
```

3. **Conditional extraction** in workflow function:
```python
# STEP 3: Extract Variant and Patient Data (Optional)
if not skip_extraction:
    logger.info("\n🧬 STEP 3: Extracting variant and patient data using AI...")
    # ... existing extraction code ...
else:
    logger.info("\n⏭️  STEP 3: Skipping extraction (--skip-extraction flag set)")
    extractions = []
```

4. **Pass model parameter to extractor**:
```python
from config.settings import get_settings
settings = get_settings()

# Override model if specified
if args.model:
    extractor = ExpertExtractor(models=[args.model], tier_threshold=tier_threshold)
else:
    extractor = ExpertExtractor(tier_threshold=tier_threshold)
```

#### Changes to `rerun_extraction.py`:

Same pattern - add soft warning and `--skip-extraction` flag for consistency.

#### Why This Approach?

- ✅ **Leverages LiteLLM**: Already supports OpenAI, Anthropic, local models via model strings
- ✅ **No new abstractions**: Reuses existing `BaseLLMCaller` and `ExpertExtractor`
- ✅ **Immediate value**: Solves the "no API key" problem today
- ✅ **Future-proof**: Can add full abstraction layer later if needed

---

## Phase 2: Parallelism Controls

### Problem
Hard-coded parallelism with no visibility:
```python
max_workers = min(8, len(markdown_files))  # Line 180
```

### Solution
Expose parallelism as configurable parameters with enhanced logging

#### Changes to `automated_workflow.py`:

1. **Add CLI flags**:
```python
parser.add_argument("--max-workers", type=int, default=None,
                   help="Number of parallel workers (default: auto-detect based on CPU)")
parser.add_argument("--max-workers-cap", type=int, default=8,
                   help="Maximum cap for auto-detected workers (default: 8)")
parser.add_argument("--rate-limit-delay", type=float, default=0.0,
                   help="Delay in seconds between API requests to avoid rate limits (default: 0)")
```

2. **CPU-based heuristic** (replace line 180):
```python
import os
import time
import threading

# Determine worker count
if args.max_workers:
    max_workers = args.max_workers
else:
    # Auto-detect based on CPUs
    cpu_count = os.cpu_count() or 4
    max_workers = min(args.max_workers_cap, cpu_count)

max_workers = min(max_workers, len(markdown_files)) if markdown_files else 1

logger.info(f"🔧 Parallelism configured: {max_workers} workers (cap: {args.max_workers_cap})")
```

3. **Enhanced logging** in `process_paper_file`:
```python
def process_paper_file(md_file):
    """Process a single paper file (for parallel execution)"""
    thread_id = threading.current_thread().name
    start_time = time.time()

    # Extract PMID from filename
    pmid_match = md_file.stem.split('_')[1] if '_' in md_file.stem else None

    if not pmid_match:
        logger.warning(f"[{thread_id}] Could not extract PMID from filename: {md_file.name}")
        return None

    logger.info(f"[{thread_id}] Processing PMID {pmid_match}...")

    # Rate limiting
    if args.rate_limit_delay > 0:
        time.sleep(args.rate_limit_delay)

    # ... existing processing code ...

    elapsed = time.time() - start_time
    logger.info(f"[{thread_id}] ✓ Completed PMID {pmid_match} in {elapsed:.1f}s ({num_variants} variants)")
```

4. **Progress tracking** in results collection:
```python
# Collect results as they complete
completed = 0
total = len(markdown_files)

for future in as_completed(future_to_file):
    md_file = future_to_file[future]
    completed += 1

    try:
        result = future.result()
        if result:
            extractions.append(result)

        # Progress indicator
        logger.info(f"📊 Progress: {completed}/{total} papers ({completed/total*100:.1f}%)")
    except Exception as e:
        logger.error(f"⚠ Failed to process {md_file.name}: {e}")
```

#### Changes to `rerun_extraction.py`:

Apply same parallelism controls and logging enhancements.

---

## Phase 3: SQLite Integration

### Problem
Migration is separate manual step - users want unified workflow

### Solution
Add optional `--migrate-to-sqlite` flag to workflow

#### Changes to `automated_workflow.py`:

1. **Add CLI flags**:
```python
parser.add_argument("--migrate-to-sqlite", action="store_true",
                   help="Migrate extraction results to SQLite database after workflow")
parser.add_argument("--db-path", type=str, default=None,
                   help="SQLite database path (default: auto-detect as {GENE}.db)")
parser.add_argument("--cleanup-after-migrate", action="store_true",
                   help="Archive PMC files after successful migration")
```

2. **Add Step 5: Migration** (after aggregation):
```python
# ============================================================================
# STEP 5: Migrate to SQLite (Optional)
# ============================================================================
migration_stats = None

if migrate_to_sqlite:
    logger.info("\n🗄️  STEP 5: Migrating to SQLite database...")

    try:
        from migrate_to_sqlite import (
            create_database_schema,
            migrate_extraction_directory as migrate_dir,
            cleanup_data_directory
        )

        # Determine database path
        if db_path:
            db_file = db_path
        else:
            db_file = f"{gene_symbol}.db"

        logger.info(f"Database: {db_file}")

        # Create schema and migrate
        conn = create_database_schema(db_file)
        migration_stats = migrate_dir(conn, extraction_dir)

        logger.info(f"✓ Migrated {migration_stats['successful']}/{migration_stats['total_files']} papers to {db_file}")

        # Optional cleanup
        if cleanup_after_migrate:
            logger.info("🧹 Running cleanup and archival...")
            cleanup_results = cleanup_data_directory(
                output_path,
                delete_empty_dirs=True,
                archive_pmc=True,
                delete_pmc_after_archive=True,
                dry_run=False
            )
            logger.info(f"✓ Archived {len(cleanup_results['archives_created'])} directories")

        conn.close()

    except Exception as e:
        logger.error(f"⚠️  Migration failed: {e}", exc_info=True)
        migration_stats = {"error": str(e)}

else:
    logger.info("\n⏭️  STEP 5: Skipping SQLite migration (use --migrate-to-sqlite to enable)")
```

3. **Update summary** (lines 242-289):
```python
summary = {
    "gene_symbol": gene_symbol,
    "workflow_timestamp": datetime.now().isoformat(),
    "statistics": {
        "pmids_discovered": len(pmids),
        "papers_downloaded": num_downloaded,
        "papers_extracted": len(extractions),
        "total_variants_found": total_variants,
        # ... existing stats ...
    },
    "extraction": {
        "skip_extraction": skip_extraction,
        "model_used": args.model or "default (from config)",
        "parallel_workers": max_workers,
    },
    "output_locations": {
        # ... existing paths ...
    },
    "penetrance_validation": {
        # ... existing validation ...
    }
}

# Add migration stats if performed
if migration_stats:
    summary["migration"] = {
        "database_path": db_file if migrate_to_sqlite else None,
        "papers_migrated": migration_stats.get("successful", 0),
        "migration_errors": migration_stats.get("failed", 0),
    }
```

#### Changes to `migrate_to_sqlite.py`:

**Already well-designed for import** - no changes needed! The script already exposes:
- `create_database_schema(db_path)` → returns connection
- `migrate_extraction_directory(conn, extraction_dir)` → returns stats dict
- `cleanup_data_directory(...)` → returns cleanup results

Just need to ensure it's importable:
```python
# Add to top of file if not present
if __name__ == "__main__":
    sys.exit(main())
```

---

## Implementation Order

### Sprint 1: Graceful Degradation (2-3 hours)
1. ✅ Modify API key check in `automated_workflow.py`
2. ✅ Add `--skip-extraction` and `--model` flags
3. ✅ Add conditional extraction logic
4. ✅ Apply same changes to `rerun_extraction.py`
5. ✅ Test harvest-only workflow

### Sprint 2: Parallelism Controls (2-3 hours)
1. ✅ Add parallelism flags to both scripts
2. ✅ Implement CPU-based heuristic
3. ✅ Enhance logging with thread IDs and timing
4. ✅ Add progress tracking
5. ✅ Test with different worker counts

### Sprint 3: SQLite Integration (1-2 hours)
1. ✅ Add migration flags to `automated_workflow.py`
2. ✅ Implement Step 5 migration logic
3. ✅ Update summary format
4. ✅ Test end-to-end unified workflow

### Sprint 4: Documentation & Polish (1 hour)
1. ✅ Update CLI help text
2. ✅ Add usage examples
3. ✅ Document new flags in README
4. ✅ Create migration guide

---

## Testing Plan

### Test 1: Harvest-Only Workflow (No API Key)
```bash
# Unset API keys
unset OPENAI_API_KEY
unset AI_INTEGRATIONS_OPENAI_API_KEY

# Should complete successfully (harvest only)
python automated_workflow.py SCN5A --email test@example.com \
    --max-pmids 10 --max-downloads 5 --skip-extraction
```

**Expected:** Downloads papers, skips extraction, no crash

### Test 2: Custom Parallelism
```bash
# Limit to 2 workers with rate limiting
python automated_workflow.py TTR --email test@example.com \
    --max-pmids 10 --max-downloads 5 \
    --max-workers 2 --rate-limit-delay 1.0
```

**Expected:** Processes 2 papers at a time, logs thread IDs, 1s delay between requests

### Test 3: Unified Workflow with Migration
```bash
# Complete workflow: harvest → extract → aggregate → migrate
python automated_workflow.py BRCA1 --email test@example.com \
    --max-pmids 20 --max-downloads 10 \
    --migrate-to-sqlite --db-path BRCA1_test.db
```

**Expected:** Full workflow completes, BRCA1_test.db created with data

### Test 4: Alternative LLM Provider
```bash
# Use Anthropic Claude instead of OpenAI
export ANTHROPIC_API_KEY="your-key-here"

python automated_workflow.py TP53 --email test@example.com \
    --max-pmids 10 --max-downloads 5 \
    --model "anthropic/claude-3-5-sonnet-20241022"
```

**Expected:** Uses Claude via LiteLLM, extracts variants successfully

---

## Benefits Summary

### For Users
- ✅ **No more crashes** on missing API keys
- ✅ **Harvest-only mode** for paper collection
- ✅ **One-command workflow** including database migration
- ✅ **Rate limiting control** to avoid API throttling
- ✅ **Progress visibility** with detailed logging
- ✅ **Provider flexibility** via LiteLLM

### For Maintainers
- ✅ **Simple codebase** - no unnecessary abstractions
- ✅ **Backwards compatible** - all changes are opt-in
- ✅ **Well-tested** - clear test scenarios
- ✅ **Extensible** - can add full abstraction layer later if needed

---

## Future Enhancements (Out of Scope)

These can be added later if needed:

1. **Resume Capability**: Save `.workflow_state.json` for crash recovery
2. **Progress Bars**: Use `tqdm` for visual progress indicators
3. **Config File**: Support `.workflow-config.yaml` for complex setups
4. **Full Extractor Abstraction**: If multiple extractor backends are needed
5. **Caching Layer**: Cache LLM responses to reduce API costs
6. **Batch Processing**: Process multiple genes in one command

---

## Migration Guide for Users

### Before (Old Workflow)
```bash
# Step 1: Run workflow
python automated_workflow.py TTR --email me@example.com

# Step 2: Manually migrate (separate command)
python migrate_to_sqlite.py --data-dir automated_output/TTR/20251125_114028
```

### After (New Workflow)
```bash
# Single unified command
python automated_workflow.py TTR --email me@example.com \
    --migrate-to-sqlite --cleanup-after-migrate

# Or harvest-only mode (no API key needed)
python automated_workflow.py TTR --email me@example.com \
    --skip-extraction --max-downloads 100
```

---

## Success Criteria

- ✅ Workflow completes without API key when `--skip-extraction` is used
- ✅ Parallelism is configurable and logged clearly
- ✅ SQLite migration is optional and integrated
- ✅ All existing functionality still works
- ✅ No breaking changes to existing scripts
- ✅ Test suite passes
- ✅ Documentation is updated

---

## Conclusion

This simplified approach achieves the goals of the original plan while:
- **Reducing complexity** by leveraging existing abstractions
- **Delivering value faster** through incremental improvements
- **Maintaining flexibility** for future enhancements

Total implementation time: **6-9 hours** vs. **20+ hours** for full abstraction approach.
