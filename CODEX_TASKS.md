# CODEX_TASKS.md - GVF Recall Improvement

## Goal
Improve variant extraction recall from 5% to 90% by fixing supplementary table downloads.

## Root Cause
Papers reference variants in "online suppl. table 1" but GVF doesn't:
1. Parse these references from main text
2. Know where Karger (and other publishers) host supplement files
3. Follow supplement URLs that aren't on the same page as the article

## Test Case
**PMID 26496715** (DOI: 10.1159/000440608)
- Paper says: "54 KCNH2 mutations in online suppl. table 1"
- GVF extracted: 4 variants (from main text only)
- Missing: 50+ variants in supplementary table

---

## Task 1: Audit supplement_scraper.py
**Time estimate:** 30 min

**Goal:** Document current supplement detection logic and identify gaps.

**Steps:**
1. Read `harvesting/supplement_scraper.py`
2. Document which publishers have specific handlers (currently: Nature, Elsevier)
3. Document what patterns the generic scraper looks for
4. Identify why Karger supplements aren't being found

**Output:**
Add a docstring block at the top of the file documenting:
- Supported publishers and their URL patterns
- Generic fallback patterns
- Known gaps (Karger, Oxford, Springer, etc.)

**Acceptance criteria:**
- The file has a comprehensive docstring explaining coverage
- Gap list includes at least: Karger, Oxford Academic, Springer, AHA Journals

---

## Task 2: Create Karger publisher handler
**Time estimate:** 45 min

**Goal:** Add `scrape_karger_supplements()` method to find Karger supplement URLs.

**Context:**
- Karger DOIs: `10.1159/XXXXXX`
- Karger article URLs: `karger.com/crd/article/...` or `karger.com/doi/...`
- Karger supplement pattern: `karger.com/doi/suppl/10.1159/XXXXXX/`
- Alternative: Look for "Supplementary Material" section on article page

**Steps:**
1. Add `scrape_karger_supplements(self, html: str, base_url: str) -> List[Dict]` method
2. Look for:
   - Links in "Supplementary Material" section
   - Direct download links matching patterns like `/suppl/`, `/supplement/`
   - Data-attribute links (some publishers use JS)
3. Extract filename and download URL for each supplement
4. Update the domain routing in `doi_resolver.py` to use this handler for `karger.com`

**Test:**
```bash
# After implementation, this DOI should find supplements:
curl -sL "https://doi.org/10.1159/000440608" | grep -i suppl
```

**Acceptance criteria:**
- New `scrape_karger_supplements()` method exists
- `doi_resolver.py` routes karger.com to this handler
- Method finds at least the supplementary table link for test DOI

---

## Task 3: Add supplement reference parser
**Time estimate:** 30 min

**Goal:** Parse main text for "online suppl. table N" references and extract expected table count.

**Location:** Create new file `harvesting/supplement_reference_parser.py`

**Steps:**
1. Create regex patterns to find supplement references in text:
   - "online suppl. table 1"
   - "supplementary table S1"
   - "see supplement"
   - "additional file N"
   - "(online suppl. material)"
2. Extract:
   - Number of expected supplement tables
   - Any variant counts mentioned (e.g., "54 mutations in online suppl. table 1")
3. Return structured data: `{"expected_tables": 3, "mentioned_variants": 54}`

**Interface:**
```python
def parse_supplement_references(text: str) -> Dict[str, Any]:
    """
    Parse main text for references to supplementary materials.
    
    Returns:
        Dict with keys:
        - expected_tables: int (number of distinct tables referenced)
        - table_refs: List[str] (specific references found)
        - mentioned_variants: int (variant count if explicitly mentioned)
    """
```

**Test with:**
```python
text = "We identified 54 KCNH2 mutations (online suppl. table 1, www.karger.com/doi/...)..."
result = parse_supplement_references(text)
assert result["expected_tables"] >= 1
assert result["mentioned_variants"] == 54
```

**Acceptance criteria:**
- New module exists at `harvesting/supplement_reference_parser.py`
- Finds "online suppl. table N" pattern
- Extracts variant counts when mentioned
- Has unit tests in `tests/test_supplement_reference_parser.py`

---

## Task 4: Integrate supplement URL extraction from text
**Time estimate:** 30 min

**Goal:** Extract explicit supplement URLs from main text and follow them.

**Context:**
Papers often include the actual URL: "online suppl. material, www.karger.com/doi/10.1159/000440608/suppl"

**Location:** `harvesting/orchestrator.py`

**Steps:**
1. In `_process_supplements()`, before calling `get_supplemental_files()`:
   - Parse the main markdown for URLs containing "suppl" or "supplement"
   - Use regex: `https?://[^\s<>"]+(?:suppl|supplement)[^\s<>"]*`
2. Add any found URLs to the supplement list
3. Deduplicate with URLs found by the scraper

**Test:**
Add test that parses:
```
"available at www.karger.com/doi/10.1159/000440608/suppl"
```
and returns that URL.

**Acceptance criteria:**
- Orchestrator extracts supplement URLs from markdown text
- Found URLs are added to the supplement download list
- No duplicate downloads

---

## Task 5: Add post-download validation
**Time estimate:** 30 min

**Goal:** Warn when expected supplements weren't downloaded.

**Location:** `harvesting/orchestrator.py`

**Steps:**
1. After downloading supplements, use `supplement_reference_parser` to check:
   - Did we download as many tables as were referenced?
   - If text says "54 mutations in suppl table" but we found 4 variants, flag it
2. Add warning log when expected != actual
3. Add field to manifest: `supplement_gap_detected: bool`

**Output format in logs:**
```
⚠️ PMID 26496715: Main text references 1 supplementary table(s) but 0 were downloaded
⚠️ PMID 26496715: Paper mentions 54 variants but extraction found 4 - likely missing supplements
```

**Acceptance criteria:**
- Post-download check compares expected vs actual supplements
- Warning logged when mismatch detected
- Manifest updated with gap flag

---

## Task 6: End-to-end test on PMID 26496715
**Time estimate:** 30 min

**Goal:** Verify the full pipeline now downloads the Karger supplement.

**Steps:**
1. Create integration test in `tests/test_karger_integration.py`
2. Test should:
   - Process PMID 26496715
   - Verify supplement file is downloaded
   - Verify supplement content includes variant table data
   - Verify extraction finds >50 variants (vs. 4 before)

**Test skeleton:**
```python
def test_karger_supplement_download():
    """PMID 26496715 should download Karger supplementary table with 54 KCNH2 mutations."""
    harvester = PMCHarvester(output_dir="test_output", gene_symbol="KCNH2")
    success, result, content = harvester.download_pmid("26496715")
    
    assert success, f"Download failed: {result}"
    
    # Check supplements directory
    supp_dir = Path("test_output/26496715_supplements")
    assert supp_dir.exists(), "Supplements directory not created"
    assert len(list(supp_dir.iterdir())) > 0, "No supplements downloaded"
    
    # Check unified content mentions KCNH2 variants
    assert "KCNH2" in content
    # After fix, should find reference to supplementary table content
```

**Acceptance criteria:**
- Integration test passes
- At least one supplement file is downloaded for PMID 26496715
- Test can be run with `pytest tests/test_karger_integration.py -v`

---

## Execution Order

Run tasks in this order (each builds on previous):

1. **Task 1** (Audit) - Understand current state
2. **Task 3** (Reference Parser) - Create the utility first
3. **Task 2** (Karger Handler) - Add publisher-specific logic
4. **Task 4** (URL Extraction) - Extract URLs from text
5. **Task 5** (Validation) - Add sanity checks
6. **Task 6** (E2E Test) - Verify everything works

## Running with Codex CLI

```bash
cd /mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher

# Task 1: Audit (review/documentation)
codex exec "Complete Task 1 from CODEX_TASKS.md - audit supplement_scraper.py and document coverage gaps"

# Task 2: Karger handler
codex exec "Complete Task 2 from CODEX_TASKS.md - create scrape_karger_supplements method"

# Task 3: Reference parser
codex exec "Complete Task 3 from CODEX_TASKS.md - create supplement_reference_parser.py"

# Task 4: URL extraction integration
codex exec "Complete Task 4 from CODEX_TASKS.md - extract supplement URLs from main text"

# Task 5: Post-download validation  
codex exec "Complete Task 5 from CODEX_TASKS.md - add validation warnings for missing supplements"

# Task 6: End-to-end test
codex exec "Complete Task 6 from CODEX_TASKS.md - create integration test for PMID 26496715"
```

## Success Metrics

After all tasks complete:
- PMID 26496715 downloads supplementary table
- Extraction recall improves from ~5% to target 50%+ on reprocessing
- New infrastructure supports adding more publishers easily
