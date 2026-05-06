"""Lazy single-browser Playwright pool for Tier 3.5.

One headless browser per harvester instance, opened on first use, closed
when the pool's context manager exits. Stealth patches mask common
Playwright fingerprints so publishers don't immediately serve us a CAPTCHA.
"""

from __future__ import annotations

import logging
from contextlib import contextmanager
from typing import Any, Optional

logger = logging.getLogger(__name__)


_STEALTH_SCRIPT = """
Object.defineProperty(navigator, 'webdriver', { get: () => undefined });
Object.defineProperty(navigator, 'plugins', {
    get: () => [
        { name: 'Chrome PDF Plugin', filename: 'internal-pdf-viewer' },
        { name: 'Chrome PDF Viewer', filename: 'mhjfbmdgcfjbbpaeojofohoefgiehjai' },
        { name: 'Native Client', filename: 'internal-nacl-plugin' },
    ],
});
Object.defineProperty(navigator, 'languages', { get: () => ['en-US', 'en'] });
window.chrome = window.chrome || { runtime: {} };
"""


class BrowserPool:
    """Manages a single Playwright browser context, lazily started.

    Use as a context manager OR call ``close()`` explicitly. Calling
    ``new_page()`` before the pool is started raises RuntimeError —
    callers should prefer ``page()`` which auto-starts.
    """

    def __init__(self, headless: bool = True, slow_mo: int = 0):
        self.headless = headless
        self.slow_mo = slow_mo
        self._playwright = None
        self._browser = None
        self._context = None
        self._started = False

    @property
    def started(self) -> bool:
        return self._started

    def start(self) -> None:
        """Start Playwright and open a browser context."""
        if self._started:
            return

        try:
            from playwright.sync_api import sync_playwright
        except ImportError as e:
            raise RuntimeError(
                "Playwright is not installed. Install with "
                "`pip install playwright && playwright install chromium` "
                "or `pip install -e '.[browser]'`."
            ) from e

        self._playwright = sync_playwright().start()
        launch_args = [
            "--disable-blink-features=AutomationControlled",
            "--disable-dev-shm-usage",
            "--no-first-run",
            "--no-default-browser-check",
            "--window-size=1920,1080",
        ]
        self._browser = self._playwright.chromium.launch(
            headless=self.headless,
            slow_mo=self.slow_mo,
            args=launch_args,
        )
        self._context = self._browser.new_context(
            user_agent=(
                "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                "AppleWebKit/537.36 (KHTML, like Gecko) "
                "Chrome/131.0.0.0 Safari/537.36"
            ),
            viewport={"width": 1920, "height": 1080},
            locale="en-US",
            timezone_id="America/New_York",
            color_scheme="light",
        )
        # Apply stealth patches to every new page automatically.
        self._context.add_init_script(_STEALTH_SCRIPT)
        self._started = True
        logger.info("BrowserPool started (headless=%s)", self.headless)

    def new_page(self) -> Any:
        """Return a fresh Playwright Page. Auto-starts the pool if needed."""
        if not self._started:
            self.start()
        return self._context.new_page()

    @contextmanager
    def page(self):
        """Yield a page and close it after use."""
        if not self._started:
            self.start()
        page = self._context.new_page()
        try:
            yield page
        finally:
            try:
                page.close()
            except Exception:
                pass

    def close(self) -> None:
        """Close the browser and Playwright runtime."""
        if not self._started:
            return
        for closer in (
            lambda: self._context.close() if self._context else None,
            lambda: self._browser.close() if self._browser else None,
            lambda: self._playwright.stop() if self._playwright else None,
        ):
            try:
                closer()
            except Exception:
                pass
        self._context = None
        self._browser = None
        self._playwright = None
        self._started = False
        logger.info("BrowserPool closed")

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False
