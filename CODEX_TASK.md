# Task: Save Raw Assets in browser_fetch.py

Modify `cli/browser_fetch.py` to save ALL raw assets when scraping papers.

## Requirements

### 1. SAVE RAW HTML
In `_scrape_article_html()`, before converting to markdown, save the original HTML to `{PMID}_raw.html` in the downloads directory.

### 2. EXTRACT FIGURES
- Find all `<img>` tags that look like figures:
  - class contains 'figure'
  - inside `<figure>` or `<picture>` elements
  - src contains 'figure' or 'fig'
- Download each image to `{PMID}_figures/` subdirectory
- Name files as `figure_001.jpg`, `figure_002.png` etc (preserve original extension)
- Also save any figure captions found near the images to a `{PMID}_figures/captions.json`

### 3. KEEP EXISTING MARKDOWN
The existing markdown generation should still work unchanged. The `{PMID}_FULL_CONTEXT.md` output continues as before.

### 4. OUTPUT STRUCTURE
After scraping PMID 12345678, we should have:
```
12345678_raw.html           # Original HTML (NEW)
12345678_figures/           # Folder with images (NEW)
  figure_001.jpg
  figure_002.png
  captions.json
12345678_FULL_CONTEXT.md    # Existing markdown (unchanged)
```

## Implementation Notes

- Focus on `_scrape_article_html()` method in the `BrowserFetcher` class
- Use the existing `page` object (Playwright) to download images
- Handle errors gracefully - if a figure download fails, log warning but continue
- Use `requests` or playwright's `page.goto()` for image downloads
- Check if image URL is absolute or relative and handle accordingly

## Testing
After changes, test with:
```bash
python cli/browser_fetch.py ~/gvf_output/abstract_only_papers.csv --pmids 19490382 --headless
```

Commit changes when done.
