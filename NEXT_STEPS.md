# Supplement Handling Improvements

**Problem:** The Gene Variant Fetcher pipeline fails to extract variants from ~94.6% of papers because most variant data is in supplementary files, not the main article text.

**Proposed Improvements:**

1.  **Integrate Browser Fallback More Widely:**
    *   Modify `doi_resolver.py` to use `browser_supplement_fetcher.py` when direct HTTP requests fail or when specific domains (e.g., Karger, Oxford, Wiley) are encountered.
    *   Create a more generic browser-based supplement fetcher that can be used for other publishers as needed.
2.  **Improve Generic Scraper:**
    *   Add more keywords and patterns to the generic scraper to identify supplements, including "ESM," "MOESM," "electronic supplementary," "supporting information," `/downloadSupplement`, `/MediaObjects/`, and `/figures/` paths.
    *   Implement logic to follow links that say "Supplementary" even without a file extension.
3.  **Enhance Publisher-Specific Scrapers:**
    *   Create handlers for publishers like Wiley and AHA Journals.
    *   Improve existing handlers to handle JavaScript-rendered content and API-based downloads.
4.  **Playwright for All Publishers:**
    *   Consider using Playwright as a primary scraping mechanism for all publishers to address the limitations of traditional HTML parsing. However, be mindful of the performance and resource implications of this approach.
5.  **Address Wiley Challenges:**
    *   Fully integrate the Wiley supplement handler.

**Next Steps:**

1.  Integrate the Karger browser fallback into the main pipeline within `doi_resolver.py`.
2.  Enhance the generic scraper with additional keywords and patterns.