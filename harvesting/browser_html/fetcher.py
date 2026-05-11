"""Tier 3.5 entrypoint — BrowserHTMLFetcher.

Orchestrates strategy lookup, embargo gating, browser pool lifecycle,
content validation, and per-attempt logging to ``browser_html_log.csv``.
"""

from __future__ import annotations

import csv
import datetime
import logging
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List, Optional

from .base import FetchContext, FetchResult, PublisherStrategy
from .browser_pool import BrowserPool
from .embargo import EmbargoChecker
from .strategies import find_strategy, registered_names

logger = logging.getLogger(__name__)


class BrowserHTMLFetcher:
    """Public Tier 3.5 entrypoint.

    Lifecycle:
        - constructed once per ``PMCHarvester`` instance, no browser yet
        - first ``fetch()`` call lazily starts the browser pool
        - ``close()`` shuts the browser down at end of harvest

    Disabled by default. Driven by config:
        - ``enable_browser_html_fallback`` (bool)
        - ``browser_html_publisher_allowlist`` (list of strategy NAMEs)
        - ``browser_html_min_embargo_months`` (cap to override per-strategy)
        - ``browser_html_headless`` (bool)
        - ``browser_html_max_per_run`` (int)
        - ``browser_html_per_paper_timeout_s`` (int)
    """

    def __init__(
        self,
        scraper: Any,
        converter: Any,
        session: Any,
        output_dir: Path,
        settings: Optional[Any] = None,
        validate_content_quality: Optional[Any] = None,
    ):
        self.scraper = scraper
        self.converter = converter
        self.session = session
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.settings = settings
        self._validate = validate_content_quality

        self._enabled = bool(getattr(settings, "enable_browser_html_fallback", False))
        self._allowlist = list(
            getattr(
                settings,
                "browser_html_publisher_allowlist",
                ["aha", "oxford", "wiley", "elsevier_open", "generic"],
            )
        )
        self._headless = bool(getattr(settings, "browser_html_headless", True))
        self._max_per_run = int(getattr(settings, "browser_html_max_per_run", 50))
        self._timeout_s = int(getattr(settings, "browser_html_per_paper_timeout_s", 90))
        self._use_profile = bool(getattr(settings, "browser_html_use_profile", False))
        self._profile_path = getattr(settings, "browser_html_profile_path", None)
        self._channel = getattr(settings, "browser_html_channel", "chrome")
        self._slow_mo_ms = int(getattr(settings, "browser_html_slow_mo_ms", 0))
        # Optional global override; if set and >= a strategy's own embargo,
        # the global wins. None means "respect each strategy's own value."
        self._min_embargo_override: Optional[int] = getattr(
            settings, "browser_html_min_embargo_months", None
        )

        self._pool: Optional[BrowserPool] = None
        self._embargo = EmbargoChecker()
        self._attempts_made = 0
        self.log_path = self.output_dir / "browser_html_log.csv"
        self._init_log()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def is_enabled(self) -> bool:
        return self._enabled

    def fetch(
        self,
        pmid: str,
        doi: str,
        pub_date: Optional[datetime.date] = None,
    ) -> Optional[FetchResult]:
        """Run Tier 3.5 for a single PMID. Returns None when ineligible."""
        if not self._enabled:
            return None
        if self._max_per_run > 0 and self._attempts_made >= self._max_per_run:
            self._log_row(pmid, doi, "", "skipped", "max_per_run cap reached", 0, 0, 0)
            return None

        strategy = find_strategy(doi=doi, allowlist=self._allowlist)
        if strategy is None:
            self._log_row(pmid, doi, "", "skipped", "no strategy matched", 0, 0, 0)
            return None

        # Embargo gate
        embargo_months = self._effective_embargo_months(strategy)
        eligible, embargo_reason = self._embargo.is_eligible(pub_date, embargo_months)
        if not eligible:
            self._log_row(
                pmid,
                doi,
                strategy.NAME,
                "skipped",
                f"embargo: {embargo_reason}",
                0,
                0,
                0,
            )
            return None

        if self._pool is None:
            self._pool = BrowserPool(
                headless=self._headless,
                slow_mo=self._slow_mo_ms,
                use_profile=self._use_profile,
                profile_path=self._profile_path,
                channel=self._channel,
            )

        ctx = FetchContext(
            pmid=pmid,
            doi=doi,
            output_dir=self.output_dir,
            scraper=self.scraper,
            converter=self.converter,
            session=self.session,
            timeout_s=self._timeout_s,
        )

        start = time.time()
        result: Optional[FetchResult] = None
        try:
            with self._pool.page() as page:
                page.set_default_timeout(self._timeout_s * 1000)
                result = strategy.fetch(page, ctx)
        except Exception as e:
            elapsed = time.time() - start
            self._attempts_made += 1
            logger.warning("Tier 3.5 strategy %s crashed: %s", strategy.NAME, e)
            self._log_row(pmid, doi, strategy.NAME, "error", str(e), 0, 0, elapsed)
            return FetchResult(publisher=strategy.NAME, error=str(e))

        elapsed = time.time() - start
        self._attempts_made += 1

        if result is None:
            self._log_row(
                pmid,
                doi,
                strategy.NAME,
                "error",
                "strategy returned None",
                0,
                0,
                elapsed,
            )
            return None

        # Surface the silent failure mode where the strategy "succeeded" but
        # produced no markdown — observed on cohort-paper URLs whose page
        # loaded but the publisher-aware scraper couldn't find body content.
        if result.main_markdown is None and result.error is None:
            result.notes.append(
                "strategy returned empty main_markdown after extraction"
            )

        # Optional content validation — refuse to count "404 page" content.
        if result.is_usable() and self._validate is not None:
            try:
                ok, why = self._validate(result.main_markdown)
                if not ok:
                    result.main_markdown = None
                    result.notes.append(f"content validation failed: {why}")
            except Exception as e:
                result.notes.append(f"validator crashed: {e}")

        outcome = "success" if result.is_usable() else "failed"
        self._log_row(
            pmid,
            doi,
            strategy.NAME,
            outcome,
            result.error or ";".join(result.notes),
            len(result.main_markdown or ""),
            len(result.supp_files or []),
            elapsed,
        )
        return result

    def close(self) -> None:
        """Shut the browser down. Safe to call multiple times."""
        if self._pool is not None:
            try:
                self._pool.close()
            except Exception:
                pass
            self._pool = None

    # ------------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------------

    @property
    def attempts_made(self) -> int:
        return self._attempts_made

    @property
    def registered_strategies(self) -> List[str]:
        return registered_names()

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _effective_embargo_months(self, strategy: PublisherStrategy) -> Optional[int]:
        """Resolve per-strategy embargo, with optional global override."""
        own = strategy.EMBARGO_MONTHS
        if self._min_embargo_override is None:
            return own
        if own is None:
            return self._min_embargo_override
        return max(own, self._min_embargo_override)

    def _init_log(self) -> None:
        if self.log_path.exists():
            return
        with open(self.log_path, "w", newline="") as f:
            csv.writer(f).writerow(
                [
                    "timestamp",
                    "pmid",
                    "doi",
                    "strategy",
                    "outcome",
                    "reason",
                    "markdown_chars",
                    "supplement_count",
                    "elapsed_s",
                ]
            )

    def _log_row(
        self,
        pmid: str,
        doi: str,
        strategy_name: str,
        outcome: str,
        reason: str,
        markdown_chars: int,
        supp_count: int,
        elapsed_s: float,
    ) -> None:
        try:
            with open(self.log_path, "a", newline="") as f:
                csv.writer(f).writerow(
                    [
                        datetime.datetime.now().isoformat(timespec="seconds"),
                        pmid,
                        doi,
                        strategy_name,
                        outcome,
                        (reason or "")[:500],
                        markdown_chars,
                        supp_count,
                        f"{elapsed_s:.2f}",
                    ]
                )
        except Exception as e:
            logger.warning("Failed to write browser_html_log row: %s", e)
