#!/usr/bin/env python3
"""Fetch paywalled biomedical articles using the local Chrome session.

Drives the Tier 3.5 browser stack with an ``AuthenticatedBrowserPool`` that
inherits cookies from the user's Chrome profile, so publishers that gate
content behind institutional SSO (e.g. VUMC's Microsoft sign-in) serve full
text instead of a paywall stub.

Usage::

    # Test on the four hardest paywalled PMIDs from our missing-variants gap.
    python -m scripts.fetch_paywalled --output ./paywall_test

    # Custom PMID list.
    python -m scripts.fetch_paywalled --pmid 15840476 --pmid 10973849 \\
        --output ./out --no-headless

    # Pre-resolved DOI (skips the NCBI roundtrip).
    python -m scripts.fetch_paywalled --pmid-doi 15840476=10.1161/01.CIR.0000164255.06478.96 \\
        --output ./out

The script writes one ``{PMID}/FULL_CONTEXT.md`` per success, mirroring what
the main extraction pipeline expects to find. A summary line per PMID is
printed at the end (publisher / chars / paywall-marker count / outcome).
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Path bootstrap so `python scripts/fetch_paywalled.py` works without -m.
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import requests  # noqa: E402

from harvesting.browser_html import (  # noqa: E402
    AuthenticatedBrowserPool,
    BrowserHTMLFetcher,
    cookie_domain_summary,
    load_chrome_cookies,
    validate_article_content,
)
from harvesting.browser_html.dom_extract import extract_body_markdown  # noqa: E402
from harvesting.browser_html.strategies import (  # noqa: E402
    find_strategy,
    registered_names,
)
from harvesting.format_converters import FormatConverter  # noqa: E402
from harvesting.supplement_scraper import SupplementScraper  # noqa: E402


# Body selectors used by the PMC HTML fallback. PMC articles render the
# body inside <article>, <main>, or a #mc/.tsec container depending on the
# theme.
_PMC_BODY_SELECTORS = (
    "article",
    "main",
    "div.tsec",
    "#mc",
    "div#maincontent",
    "div.article-page",
    "body",
)


def europepmc_lookup_pmcid(pmid: str, session: requests.Session) -> Optional[str]:
    """Return PMCID for a PMID if Europe PMC has one, else None.

    Europe PMC's search API returns the PMCID directly when an article has
    been deposited in PMC. NIH-funded papers usually do; subscription papers
    usually don't.
    """
    url = (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
        f"?query=EXT_ID:{pmid}%20AND%20SRC:MED&format=json&resultType=lite"
    )
    try:
        r = session.get(url, timeout=20)
        if r.status_code != 200:
            return None
        results = r.json().get("resultList", {}).get("result", [])
        if not results:
            return None
        pmcid = results[0].get("pmcid")
        if pmcid and not pmcid.startswith("PMC"):
            pmcid = f"PMC{pmcid}"
        return pmcid
    except Exception:
        return None


def try_pmc_fallback(
    pmid: str,
    output_dir: Path,
    session: requests.Session,
) -> Optional[Dict]:
    """Attempt to recover a stub via the PMC HTML deposit.

    Many subscription papers (especially NIH-funded ones) have a PMC version
    that's free to read after the embargo. We look up the PMCID via Europe
    PMC, fetch the PMC HTML, and run our DOM extractor on it.

    Returns a result row dict on success, or None if no PMC version exists
    or the extraction failed quality gates. On success the row's ``path``
    points at the rewritten FULL_CONTEXT.md.
    """
    pmcid = europepmc_lookup_pmcid(pmid, session)
    if not pmcid:
        return None

    pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
    try:
        r = session.get(
            pmc_url,
            timeout=30,
            headers={
                "User-Agent": (
                    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                    "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36"
                ),
                "Accept": "text/html,application/xhtml+xml",
            },
        )
    except Exception as e:
        LOG.warning("PMC fallback fetch failed for PMID %s: %s", pmid, e)
        return None
    if not r.ok or len(r.content) < 5000:
        return None

    html = r.text
    markdown = extract_body_markdown(html, _PMC_BODY_SELECTORS)
    if not markdown:
        return None

    ok, reason = validate_article_content(markdown)
    if not ok:
        # PMC content didn't pass the gate either — could be an abstract-
        # only PMC record. Don't overwrite the stub from Tier 3.5.
        LOG.info("PMC fallback for %s extracted but failed gate: %s", pmid, reason)
        return None

    # Write a PMC FULL_CONTEXT.md and source marker.
    pmid_dir = output_dir / pmid
    pmid_dir.mkdir(parents=True, exist_ok=True)
    (pmid_dir / "FULL_CONTEXT.md").write_text(markdown, encoding="utf-8")
    (pmid_dir / "page.html").write_text(html, encoding="utf-8")
    (pmid_dir / "source.txt").write_text(
        f"Recovered via PMC fallback: {pmc_url}\n", encoding="utf-8"
    )
    return {
        "pmid": pmid,
        "strategy": "pmc_fallback",
        "outcome": "success",
        "reason": reason,
        "chars": len(markdown),
        "supp_files": 0,
        "final_url": pmc_url,
        "path": str(pmid_dir / "FULL_CONTEXT.md"),
        "pmcid": pmcid,
    }


LOG = logging.getLogger("fetch_paywalled")

# Highest-value paywalled PMIDs from the May 2026 gap analysis.
DEFAULT_PMIDS: Tuple[str, ...] = (
    "15840476",  # Tester 2005, Circulation — 86 gold variants
    "10973849",  # Splawski 2000, Circulation — 59 variants
    "26496715",  # Heart Rhythm 2015 — 53 variants
    "11854117",  # Splawski 2002, Circulation — 44 variants
)

# Override map for PMIDs whose DOI we already know (avoids an NCBI lookup).
# Filled in lazily by parse_pmid_doi_overrides() from CLI args.


def _ncbi_email() -> str:
    return (
        os.environ.get("ENTREZ_EMAIL")
        or os.environ.get("NCBI_EMAIL")
        or "brett.kroncke@gmail.com"
    )


def pubmed_resolve_doi(
    pmid: str, session: requests.Session, max_attempts: int = 4
) -> Optional[str]:
    """Resolve a PMID to its DOI via NCBI esummary.

    NCBI E-utilities rate-limits unauthenticated callers to 3 req/s; we retry
    with exponential backoff on 429. Returns None when no DOI is registered
    or the API repeatedly fails.
    """
    import time

    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        f"?db=pubmed&id={pmid}&retmode=json"
    )
    delay = 0.5
    for attempt in range(1, max_attempts + 1):
        try:
            r = session.get(url, timeout=20)
            if r.status_code == 429:
                time.sleep(delay)
                delay *= 2
                continue
            r.raise_for_status()
            item = r.json().get("result", {}).get(str(pmid), {})
            for artid in item.get("articleids", []) or []:
                if (artid.get("idtype") or "").lower() == "doi":
                    return artid.get("value")
            return None
        except Exception as e:
            if attempt == max_attempts:
                LOG.warning("DOI resolution failed for PMID %s: %s", pmid, e)
                return None
            time.sleep(delay)
            delay *= 2
    return None


def make_session() -> requests.Session:
    s = requests.Session()
    s.headers.update(
        {
            "User-Agent": (f"GVF-PaywalledFetch/1.0 (mailto:{_ncbi_email()})"),
        }
    )
    return s


def parse_pmid_doi_overrides(raw: List[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for entry in raw or []:
        if "=" not in entry:
            continue
        k, v = entry.split("=", 1)
        k, v = k.strip(), v.strip()
        if k and v:
            out[k] = v
    return out


def write_outputs(pmid: str, result, output_dir: Path) -> Tuple[Path, Optional[str]]:
    """Write FULL_CONTEXT.md + raw HTML; return (path, body_text_for_validation)."""
    pmid_dir = output_dir / pmid
    pmid_dir.mkdir(parents=True, exist_ok=True)

    body = result.main_markdown or ""
    full_ctx = pmid_dir / "FULL_CONTEXT.md"
    full_ctx.write_text(body, encoding="utf-8")

    if result.main_html:
        (pmid_dir / "page.html").write_text(result.main_html, encoding="utf-8")

    meta = {
        "pmid": pmid,
        "publisher": result.publisher,
        "final_url": result.final_url,
        "supp_files": result.supp_files,
        "figure_paths": [str(p) for p in (result.figure_paths or [])],
        "notes": result.notes,
        "error": result.error,
        "markdown_chars": len(body),
    }
    (pmid_dir / "result.json").write_text(
        json.dumps(meta, indent=2, default=str), encoding="utf-8"
    )

    return full_ctx, (body or None)


def fetch_one(
    fetcher: BrowserHTMLFetcher,
    pmid: str,
    doi: Optional[str],
    output_dir: Path,
    pmc_session: Optional[requests.Session] = None,
) -> Dict:
    """Run Tier 3.5 for a single PMID; return a result row.

    On stub/empty outcomes from Tier 3.5, attempt the PMC fallback (Europe
    PMC PMCID lookup + PMC HTML fetch) so subscription papers with NIH PMC
    deposits are recovered automatically.
    """
    row: Dict = {
        "pmid": pmid,
        "doi": doi or "",
        "strategy": "",
        "outcome": "skipped",
        "reason": "",
        "chars": 0,
        "supp_files": 0,
        "final_url": "",
    }
    if not doi:
        row["reason"] = "no DOI"
        return row

    # Tell the user which strategy will run, before we burn time on the fetch.
    strategy = find_strategy(doi=doi, allowlist=None)
    row["strategy"] = strategy.NAME if strategy else "(none)"
    if strategy is None:
        row["reason"] = "no matching strategy"
        return row

    result = fetcher.fetch(pmid=pmid, doi=doi, pub_date=None)
    if result is None:
        row["outcome"] = "skipped"
        row["reason"] = "fetcher returned None"
        return row

    row["final_url"] = result.final_url or ""
    row["supp_files"] = len(result.supp_files or [])

    full_ctx_path, body = write_outputs(pmid, result, output_dir)
    row["chars"] = len(body or "")
    row["path"] = str(full_ctx_path)

    # Apply quality gate (the fetcher itself also runs the validator, but we
    # surface the reason here for the per-PMID summary table).
    ok, reason = validate_article_content(body or "")
    if result.error:
        row["outcome"] = "error"
        row["reason"] = result.error
    elif not body:
        row["outcome"] = "empty"
        row["reason"] = ";".join(result.notes) or "no markdown"
    elif not ok:
        row["outcome"] = "paywall_or_stub"
        row["reason"] = reason
    else:
        row["outcome"] = "success"
        row["reason"] = reason

    # PMC fallback. If Tier 3.5 left us with a stub or empty body, try the
    # Europe PMC route — many NIH-funded subscription papers have a free
    # PMC deposit. The fallback overwrites FULL_CONTEXT.md only when its
    # extraction passes the quality gate, so a failed PMC attempt doesn't
    # mask the Tier 3.5 stub for diagnostics.
    if row["outcome"] in ("paywall_or_stub", "empty") and pmc_session is not None:
        pmc_row = try_pmc_fallback(pmid, output_dir, pmc_session)
        if pmc_row is not None:
            print(
                f"  PMC fallback succeeded: {pmc_row['pmcid']} -> {pmc_row['chars']} chars"
            )
            row["outcome"] = "success_via_pmc"
            row["strategy"] = f"{row['strategy']}+pmc_fallback"
            row["reason"] = f"PMC fallback ({pmc_row['pmcid']}): {pmc_row['reason']}"
            row["chars"] = pmc_row["chars"]
            row["final_url"] = pmc_row["final_url"]
            row["path"] = pmc_row["path"]
            row["pmcid"] = pmc_row["pmcid"]

    return row


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--pmid",
        action="append",
        default=None,
        help="PMID to fetch (repeat for multiple). Defaults to the 4 top-value test PMIDs.",
    )
    parser.add_argument(
        "--pmid-doi",
        action="append",
        default=[],
        help="PMID=DOI override (repeat). Skips the NCBI esummary lookup.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path("./paywall_test"),
        help="Output directory (default: ./paywall_test)",
    )
    parser.add_argument(
        "--profile",
        default=None,
        help="Chrome profile name (e.g. 'Default'). Default: merge all profiles.",
    )
    parser.add_argument(
        "--headless",
        dest="headless",
        action="store_true",
        default=True,
    )
    parser.add_argument(
        "--no-headless",
        dest="headless",
        action="store_false",
        help="Run with a visible browser window (useful for debugging SSO challenges).",
    )
    parser.add_argument(
        "--allow",
        action="append",
        default=None,
        help="Strategy NAME to allow (repeat). Default: all registered strategies.",
    )
    parser.add_argument(
        "--no-chrome-channel",
        dest="use_chrome_channel",
        action="store_false",
        default=True,
        help="Use bundled Chromium instead of the locally-installed Chrome binary.",
    )
    parser.add_argument(
        "--timeout-s",
        type=int,
        default=120,
        help="Per-paper timeout in seconds (default: 120).",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    pmids: List[str] = list(args.pmid) if args.pmid else list(DEFAULT_PMIDS)
    overrides = parse_pmid_doi_overrides(args.pmid_doi)
    output_dir: Path = args.output.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"PMIDs:       {', '.join(pmids)}")
    print(f"Output dir:  {output_dir}")
    print(f"Strategies:  {registered_names()}")
    print(f"Headless:    {args.headless}")
    print()

    # ---- 1. Load cookies from local Chrome ----
    print("Loading Chrome cookies…")
    cookies = load_chrome_cookies(profile_name=args.profile)
    summary = cookie_domain_summary(cookies)
    print(f"  total cookies: {len(cookies)}")
    for d, n in list(summary.items())[:20]:
        print(f"    {d:40s} {n}")
    print()

    # ---- 2. Build the authenticated pool ----
    pool = AuthenticatedBrowserPool(
        cookies=cookies,
        headless=args.headless,
        use_chrome_channel=args.use_chrome_channel,
    )

    # ---- 3. Build minimal settings for the fetcher ----
    class _S:
        enable_browser_html_fallback = True
        browser_html_publisher_allowlist = args.allow or registered_names()
        browser_html_headless = args.headless
        browser_html_max_per_run = max(50, len(pmids) * 2)
        browser_html_per_paper_timeout_s = args.timeout_s
        browser_html_min_embargo_months = None

    scraper = SupplementScraper()
    converter = FormatConverter()
    session = make_session()

    fetcher = BrowserHTMLFetcher(
        scraper=scraper,
        converter=converter,
        session=session,
        output_dir=output_dir,
        settings=_S(),
        validate_content_quality=validate_article_content,
        pool=pool,
        bypass_embargo=True,
        enabled_override=True,
    )

    # ---- 4. Resolve DOIs ----
    import time as _time

    print("Resolving DOIs…")
    doi_map: Dict[str, Optional[str]] = {}
    for pmid in pmids:
        if pmid in overrides:
            doi_map[pmid] = overrides[pmid]
            print(f"  {pmid}  (override) -> {doi_map[pmid]}")
            continue
        doi = pubmed_resolve_doi(pmid, session)
        doi_map[pmid] = doi
        print(f"  {pmid}  -> {doi or '(no DOI)'}")
        _time.sleep(0.4)  # NCBI E-utilities: 3 req/s cap without API key.
    print()

    # ---- 5. Fetch ----
    import time as _t

    def _domain_for(doi: Optional[str]) -> str:
        # Map DOI prefixes to the publisher domain so we can pace per-host.
        if not doi:
            return ""
        if doi.startswith("10.1161"):
            return "ahajournals.org"
        if doi.startswith("10.1016"):
            return "sciencedirect.com"
        if doi.startswith("10.1093"):
            return "academic.oup.com"
        if doi.startswith("10.1002") or doi.startswith("10.1111"):
            return "onlinelibrary.wiley.com"
        if doi.startswith("10.1159"):
            return "karger.com"
        return "other"

    last_domain_at: Dict[str, float] = {}
    rows: List[Dict] = []

    def _run_with_pacing(pmid: str, gap_seconds: float) -> Dict:
        dom = _domain_for(doi_map[pmid])
        now = _t.time()
        wait_for = gap_seconds - (now - last_domain_at.get(dom, 0.0))
        if wait_for > 0 and last_domain_at:
            print(f"  pacing: sleeping {wait_for:.1f}s before {dom}")
            _t.sleep(wait_for)
        r = fetch_one(fetcher, pmid, doi_map[pmid], output_dir, pmc_session=session)
        last_domain_at[dom] = _t.time()
        return r

    try:
        for i, pmid in enumerate(pmids):
            print(f"--- PMID {pmid} ---")
            row = _run_with_pacing(pmid, gap_seconds=15.0)
            rows.append(row)
            print(
                f"  strategy={row['strategy']:14s} outcome={row['outcome']:18s} "
                f"chars={row['chars']:>7d} supp={row['supp_files']:>2d} "
                f"url={row['final_url'][:80]}"
            )
            if row.get("reason"):
                print(f"  reason: {row['reason']}")
            print()

        # Retry sweep — anything that failed with a Cloudflare interstitial
        # often resolves with a longer cool-down. CF escalates difficulty on
        # rapid same-domain hits, so wait substantially longer before retrying.
        retry_idxs = [
            i
            for i, r in enumerate(rows)
            if r["outcome"] != "success"
            and "cloudflare" in (r.get("reason") or "").lower()
        ]
        if retry_idxs:
            print(f"\nRetrying {len(retry_idxs)} CF-blocked PMID(s) after cool-down…\n")
            for i in retry_idxs:
                pmid = rows[i]["pmid"]
                print(f"--- retry PMID {pmid} ---")
                new_row = _run_with_pacing(pmid, gap_seconds=45.0)
                if new_row["outcome"] == "success":
                    rows[i] = new_row
                    print(
                        f"  retry succeeded: chars={new_row['chars']} "
                        f"reason={new_row['reason']}"
                    )
                else:
                    rows[i]["notes"] = (
                        rows[i].get("notes") or ""
                    ) + " retry_also_failed"
                    print(f"  retry failed: {new_row.get('reason')}")
                print()
    finally:
        try:
            pool.close()
        except Exception:
            pass

    # ---- 6. Final summary ----
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    by_publisher: Dict[str, Dict[str, int]] = {}
    for r in rows:
        p = r["strategy"] or "(none)"
        by_publisher.setdefault(p, {"success": 0, "fail": 0})
        if r["outcome"] in ("success", "success_via_pmc"):
            by_publisher[p]["success"] += 1
        else:
            by_publisher[p]["fail"] += 1
    for p, c in sorted(by_publisher.items()):
        print(f"  {p:14s}  success={c['success']}  fail={c['fail']}")

    success = sum(1 for r in rows if r["outcome"] in ("success", "success_via_pmc"))
    print(f"\nOverall: {success}/{len(rows)} PMIDs succeeded.")

    # Persist a CSV-ish JSON for downstream inspection.
    (output_dir / "summary.json").write_text(
        json.dumps(rows, indent=2, default=str), encoding="utf-8"
    )
    print(f"Wrote {output_dir / 'summary.json'}")

    return 0 if success == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
