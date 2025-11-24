# Plan: Enhance Supplemental File Extraction with DOI Resolution and Web Scraping

**Goal:** Reliably extract supplemental files for articles when the PMC API does not provide them.

**High-Level Strategy:** Implement a "waterfall" fallback strategy. The system will first attempt to retrieve files via the PMC API. If that fails, it will resolve the article's DOI to find the definitive article page and then use a domain-specific scraper to find the supplemental file links.

**Key Technical Challenges:** This plan is designed to overcome two primary obstacles in web scraping academic publishers:
1.  **Anti-Bot Protections:** Publishers like Nature and Elsevier employ strict measures that block standard `requests` calls, resulting in "403 Forbidden" errors.
2.  **Dynamic Content & Non-Standard URLs:** Many sites load supplemental file links with JavaScript *after* the initial page load. Furthermore, URLs often use internal IDs (like PIIs) instead of the article's DOI, making URL "guessing" unreliable.

---

### Plan Steps

#### 1. Refactor `get_supplemental_files` to be a Controller

The existing `get_supplemental_files` function will be modified to orchestrate the data retrieval process.

- **Current Logic:** Relies solely on the PMC API.
- **New Logic:**
  1.  Call the existing PMC API method.
  2.  If the PMC method returns no files, call a new internal method, `_get_supplemental_files_from_doi`, passing the article's DOI and PMID.
  3.  Return the combined results.

#### 2. Implement the DOI Resolver and Scraper Router

This is a critical step to get the correct, final URL of the article page. **Do not construct or guess URLs.** We must resolve them.

- **Action:** Create a new method `_get_supplemental_files_from_doi(self, doi: str, pmid: str) -> List[Dict]`.
- **Implementation Details:**

  1.  **Spoof Browser Headers:** To avoid being blocked (the "403 Forbidden" wall), we must mimic a real browser. Use a `requests.Session` to persist headers and cookies across redirects.
     ```python
     headers = {
         'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
         'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
         'Referer': 'https://pubmed.ncbi.nlm.nih.gov/'
     }
     session = requests.Session()
     session.headers.update(headers)
     ```
  2.  **Resolve the DOI:** Perform a GET request to `https://doi.org/{doi}` with `allow_redirects=True`. This will follow the redirect chain to the final publisher page.
  3.  **Capture Final URL and Route:** Get the final URL from `response.url`, parse its domain, and route to the appropriate scraper.

- **Pseudocode for the Router:**
  ```python
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

Create separate, targeted methods for each publisher to handle their unique HTML structures and content loading mechanisms.

- **A. Nature Scraper (`_scrape_nature_supplements`)**
  - **Input:** HTML content and the final URL.
  - **Action:** Use `BeautifulSoup` to parse the HTML. Look for sections related to "Supplementary Information" and extract `href` attributes from `<a>` tags.

- **B. Elsevier/GIM Scraper (`_scrape_elsevier_supplements`)**
  - **Input:** HTML content and the final URL.
  - **Challenge:** These sites often load supplemental links dynamically via JavaScript. `BeautifulSoup` alone may not be sufficient as it only sees the initial HTML.
  - **Action:** Implement a multi-step approach.
    1.  **Check Static HTML:** First, use `BeautifulSoup` to look for simple `<a>` tags containing supplemental file links.
    2.  **Inspect JSON Scripts:** If that fails, search the HTML for `<script type="application/json">` blocks. These often contain article metadata, including file URLs, in a structured format.
    3.  **Regex for MMCs:** As a final fallback, use regular expressions. Elsevier often calls supplements "Multimedia Components" (MMCs) and links to them with patterns like `mmc1.pdf`, `mmc2.xlsx`. The links themselves may have a structure like `https://www.gimjournal.org/cms/attachment/{internal-id}/{filename}`.

- **C. Generic Scraper (`_scrape_generic_supplements`)**
  - **Input:** HTML content and the final URL.
  - **Action:** This will be a best-effort scraper. It will use `BeautifulSoup` to find all `<a>` tags and filter them based on keywords in the link text or `href`, such as "supplement", "supporting", "appendix", ".pdf", ".docx", ".xlsx".

#### 4. Integrate and Test

- **Action:**
  1.  Wire up all the new methods within the `GeneVariantFetcher` class.
  2.  Create a set of test cases with DOIs from various publishers (Nature, GIM/Elsevier, and others) to validate that the correct scrapers are called and that they successfully extract file URLs.
  3.  Ensure that failures (e.g., DOI resolution error, scraping failure) are handled gracefully and do not crash the program.
