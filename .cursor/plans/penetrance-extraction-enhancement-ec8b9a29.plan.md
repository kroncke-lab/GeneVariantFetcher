# Plan: Enhance Supplemental File Extraction with DOI Resolution and Web Scraping

**Goal:** Reliably extract supplemental files for articles when the PMC API does not provide them. This involves a fallback mechanism that uses the article's DOI to find and scrape the publisher's website.

**Strategy:** Implement a "waterfall" fallback strategy. The system will first attempt to retrieve files via the PMC API. If that fails, it will resolve the article's DOI to the publisher's page and use a domain-specific scraper to find the supplemental file links.

---

### Revised Plan Steps

#### 1. Refactor `get_supplemental_files` to be a Controller

The existing `get_supplemental_files` function will be modified to orchestrate the data retrieval process.

- **Current Logic:** Relies solely on the PMC API.
- **New Logic:**

  1.  Call the existing PMC API method.
  2.  If the PMC method returns no files, call a new internal method, `_get_supplemental_files_from_doi`, passing the article's DOI and PMID.
  3.  Return the combined results.

#### 2. Implement the DOI Resolver and Scraper Router (`_get_supplemental_files_from_doi`)

This new private method will be the core of the fallback logic. It will resolve the DOI and route to the correct scraper.

- **Action:** Create a new method `_get_supplemental_files_from_doi(self, doi: str, pmid: str) -> List[Dict]`.
- **Implementation Details:**

  1.  Initialize a `requests.Session` object to persist connections and cookies.
  2.  Define a standard set of browser headers (e.g., `User-Agent`, `Accept`, `Referer`) to spoof a real browser and avoid being blocked by anti-bot measures.
  3.  Perform a GET request to `https://doi.org/{doi}` with `allow_redirects=True`.
  4.  Capture the final URL from the response (`response.url`).
  5.  Parse the network location (domain) from the final URL using `urllib.parse.urlparse`.
  6.  Use an `if/elif/else` block to route to the appropriate scraper function based on the domain (e.g., `"nature.com"`, `"gimjournal.org"`, `"sciencedirect.com"`).
```python
# Pseudocode for the router
def _get_supplemental_files_from_doi(self, doi: str, pmid: str) -> List[Dict]:
    # 1. Resolve DOI to final URL with session and headers
    try:
        response = self.session.get(f"https://doi.org/{doi}", allow_redirects=True)
        final_url = response.url
        domain = urlparse(final_url).netloc
    except Exception as e:
        # Log error and return empty
        return []

    # 2. Route to specific scraper based on resolved domain
    if "nature.com" in domain:
        return self._scrape_nature_supplements(response.text, final_url)
    elif "gimjournal.org" in domain or "sciencedirect" in domain:
        return self._scrape_elsevier_supplements(response.text, final_url)
    else:
        return self._scrape_generic_supplements(response.text, final_url)
```


#### 3. Implement Publisher-Specific Scrapers

We will create separate methods for each publisher to handle their unique HTML structures.

- **A. Nature Scraper (`_scrape_nature_supplements`)**
  - **Input:** HTML content and the final URL.
  - **Action:** Use `BeautifulSoup` to parse the HTML. Look for sections related to "Supplementary Information" and extract `href` attributes from `<a>` tags.

- **B. Elsevier/GIM Scraper (`_scrape_elsevier_supplements`)**
  - **Input:** HTML content and the final URL.
  - **Action:** This requires a multi-step approach due to dynamic content loading.

    1.  First, use `BeautifulSoup` to look for simple `<a>` tags containing supplemental file links.
    2.  If that fails, search the HTML for `<script type="application/json">` blocks that might contain article metadata, including file URLs.
    3.  As a final fallback, use regular expressions to search for patterns like `mmc1.pdf`, `mmc2.xlsx`, which are common in Elsevier's "multimedia components" (MMCs).

- **C. Generic Scraper (`_scrape_generic_supplements`)**
  - **Input:** HTML content and the final URL.
  - **Action:** This will be a best-effort scraper. It will use `BeautifulSoup` to find all `<a>` tags and filter them based on keywords in the link text or `href`, such as "supplement", "supporting", "appendix", ".pdf", ".docx", ".xlsx".

#### 4. Integrate and Test

- **Action:**

  1.  Wire up all the new methods within the `GeneVariantFetcher` class.
  2.  Create a set of test cases with DOIs from various publishers (Nature, GIM/Elsevier, and others) to validate that the correct scrapers are called and that they successfully extract file URLs.
  3.  Ensure that failures (e.g., DOI resolution error, scraping failure) are handled gracefully and do not crash the program.

This revised plan directly addresses the technical challenges of web scraping and provides a more resilient and effective solution for extracting supplemental files.