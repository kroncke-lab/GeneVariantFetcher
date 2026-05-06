"""Per-strategy embargo eligibility check.

Each strategy declares its own ``EMBARGO_MONTHS``. This module just provides
the comparison logic and a helper to fetch a publication date from PubMed.
"""

from __future__ import annotations

import datetime
import logging
from dataclasses import dataclass
from typing import Optional, Tuple

logger = logging.getLogger(__name__)


@dataclass
class EmbargoChecker:
    """Compares a paper's pub date against a strategy's embargo policy."""

    today: Optional[datetime.date] = None

    def __post_init__(self):
        if self.today is None:
            self.today = datetime.date.today()

    def is_eligible(
        self,
        pub_date: Optional[datetime.date],
        embargo_months: Optional[int],
    ) -> Tuple[bool, str]:
        """Has the embargo lapsed?

        - ``embargo_months is None`` → always eligible (we just don't know,
          so let the strategy try and fail fast on access denial)
        - ``embargo_months == 0`` → free immediately
        - ``pub_date is None`` and ``embargo_months > 0`` → ineligible (we
          can't prove the embargo has lapsed)
        - otherwise → eligible iff pub_date + embargo_months <= today
        """
        if embargo_months is None:
            return True, "no embargo policy declared"
        if embargo_months <= 0:
            return True, "no embargo"
        if pub_date is None:
            return False, "embargo present but pub_date unknown"

        # Approximate month math: 30.44 days/month average is plenty for
        # a 12-month embargo that's actually enforced loosely by publishers.
        delta = self.today - pub_date
        elapsed_months = delta.days / 30.44

        if elapsed_months >= embargo_months:
            return (
                True,
                f"{elapsed_months:.1f} mo elapsed >= {embargo_months} mo embargo",
            )
        return (
            False,
            f"{elapsed_months:.1f} mo elapsed < {embargo_months} mo embargo",
        )


def get_pub_date_from_pmid(pmid: str, pmc_api) -> Optional[datetime.date]:
    """Fetch publication date for a PMID via PubMed Entrez.

    Returns None on any failure — callers must handle that gracefully.
    The PMC API client already rate-limits and retries.
    """
    try:
        # Reuse the Entrez session managed by pmc_api. We piggyback on
        # the existing rate limiter rather than introducing a second one.
        from Bio import Entrez

        if hasattr(pmc_api, "_rate_limit"):
            pmc_api._rate_limit()
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        article = records["PubmedArticle"][0]["MedlineCitation"]["Article"]
        # Try Journal/JournalIssue/PubDate first (most reliable)
        journal_pubdate = (
            article.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
        )
        date = _parse_pubdate_dict(journal_pubdate)
        if date is not None:
            return date

        # Fallback: ArticleDate
        for ad in article.get("ArticleDate", []):
            try:
                return datetime.date(
                    int(ad["Year"]), int(ad["Month"]), int(ad.get("Day", 1))
                )
            except Exception:
                continue
        return None
    except Exception as e:
        logger.debug(f"Could not fetch pub_date for PMID {pmid}: {e}")
        return None


_MONTHS = {
    "jan": 1,
    "feb": 2,
    "mar": 3,
    "apr": 4,
    "may": 5,
    "jun": 6,
    "jul": 7,
    "aug": 8,
    "sep": 9,
    "oct": 10,
    "nov": 11,
    "dec": 12,
}


def _parse_pubdate_dict(pd: dict) -> Optional[datetime.date]:
    """PubMed PubDate fields can be year-only, year+month, or full date."""
    year = pd.get("Year")
    if not year:
        # MedlineDate fallback (e.g. "2014 Jan-Feb")
        ml = pd.get("MedlineDate", "")
        if ml:
            try:
                year = int(str(ml).split()[0])
            except Exception:
                return None
        else:
            return None
    try:
        year_i = int(str(year))
    except Exception:
        return None

    month_raw = str(pd.get("Month", "1")).strip().lower()
    month = _MONTHS.get(month_raw[:3], None)
    if month is None:
        try:
            month = int(month_raw)
        except Exception:
            month = 1

    day_raw = pd.get("Day", "1")
    try:
        day = int(str(day_raw))
    except Exception:
        day = 1

    try:
        return datetime.date(year_i, max(1, min(12, month)), max(1, min(28, day)))
    except Exception:
        return None
