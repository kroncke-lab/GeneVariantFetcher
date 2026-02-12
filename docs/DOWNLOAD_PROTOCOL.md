# GVF Paper Download Protocol

**Last updated:** 2026-02-12
**Summary of all fixes and learnings from Feb 2026 optimization sprint.**

## Download Pipeline Flow

```
PMID → PMCID lookup → PMC XML → Markdown
                    ↓ (fail)
                    Publisher API fallback (Elsevier/Springer/Wiley)
                    ↓ (fail)
                    Unpaywall OA lookup → PDF → Markdown
                    ↓ (fail)
                    DOI resolver → Web scraping
                    ↓ (fail)
                    Free URL from PubMed → Web scraping
                    ↓ (fail)
                    Mark as paywalled
```

## Key Components

### 1. Pre-Download Consolidation (`pipeline/steps.py`)
Before ANY downloads, `_consolidate_prior_downloads()` scans all prior `gvf_output/` directories for existing `FULL_CONTEXT.md` files. Copies the largest version found into the current harvest directory. This avoids re-downloading papers from previous runs.

**Why it matters:** The Jan 28 run had 249 PMC papers. The Feb 11 fresh run only got 225 because of API intermittency. Without consolidation, we lose prior work.

### 2. Resume Support (`harvesting/orchestrator.py`)
`harvest()` checks for existing `FULL_CONTEXT.md` files before starting and skips already-downloaded PMIDs. Each `download_pmid()` call is wrapped in try/except so one failure doesn't kill the batch.

### 3. PMC XML Fallback
When PMC has a PMCID but XML is unavailable (restricted, rate-limited, or embargoed), the pipeline now falls through to publisher APIs via DOI instead of immediately marking as paywalled.

### 4. Publisher APIs

| API | DOI Prefixes | Config Key | Notes |
|-----|-------------|------------|-------|
| **Elsevier** | `10.1016/` + others | `ELSEVIER_API_KEY` | XML namespace stripping fixed (case-insensitive regex) |
| **Springer** | `10.1007/`, `10.1038/`, `10.1186/` + others | `SPRINGER_API_KEY` | Key was missing from settings.py until Feb 12 fix |
| **Wiley** | `10.1111/`, `10.1002/` + others | `WILEY_API_KEY` | TDM API + web scraping fallback |

**All keys loaded via Pydantic BaseSettings from `.env`** — do NOT rely on `os.environ`.

**Last-resort:** If no DOI prefix matches a known publisher, all available APIs are tried sequentially.

### 5. Supplement Handling
- PMC papers: Supplements fetched via EuropePMC API + PMC supplementary files
- Publisher papers: DOI resolver scrapes publisher pages for supplement links
- `resolve_and_scrape_supplements()` returns `List[Dict]` with `url`/`name` keys (NOT file paths)
- Supported formats: Excel (.xlsx/.xls), Word (.docx/.doc), PDF, text, CSV
- Excel supplements are 17x more effective than PDFs for variant data

### 6. Figure Extraction
- PMC papers: Figures extracted from PMC article HTML page
- Publisher papers: Figures extracted from PDF supplements (when `EXTRACT_FIGURES=true`)
- Only PMC papers get standalone figure image files — publisher papers rely on text descriptions
- Pedigree analysis requires GPT-4o vision (post-processing step, not part of download)

## Configuration (`.env`)

```
ELSEVIER_API_KEY=your_key
WILEY_API_KEY=your_key
SPRINGER_API_KEY=your_key
NCBI_EMAIL=your_email
EXTRACT_FIGURES=true
EXTRACT_PEDIGREES=true
```

## Commits (Feb 12, 2026)

| Commit | Description |
|--------|-------------|
| `f30ea28` | Pre-download consolidation from prior runs |
| `603eec0` | Springer API key in settings.py + Elsevier XML namespace fix |
| `f339f58` | PathLike crash fix + resume support |
| `60f68a3` | PMC XML fallback to publisher APIs + broader DOI matching |

## Performance Results

| Run | Full Text | Coverage | Notes |
|-----|----------|----------|-------|
| Jan 28 | 249 | - | Original run, PMC only |
| Feb 11 | 225/600 | 37.5% | Fresh run, Springer broken, no consolidation |
| **Feb 12** | **829/1059** | **78.3%** | All fixes applied |

## Known Limitations

1. **Unpaywall** returns invalid JSON for ~95% of queries (likely rate-limited or blocked)
2. **Wiley TDM** fails for many older articles ("article not found")
3. **~230 papers** remain genuinely paywalled (no open access path via any API)
4. **Figure extraction** only works for PMC-sourced papers (publisher figures require different approach)
5. **OOM risk** on large batches with figure extraction enabled — disable `EXTRACT_FIGURES` for bulk runs
