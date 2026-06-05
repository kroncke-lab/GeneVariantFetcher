"""Playwright pool that inherits the local Chrome session.

Same lifecycle as ``BrowserPool``, with two differences:

1. Cookies loaded from the user's Chrome profile are injected into the
   Playwright ``BrowserContext`` before the first navigation. This is how we
   reuse publisher SSO sessions (e.g. VUMC institutional access) without
   asking the user to copy-paste cookies or rerun any login.
2. User agent is forced to the major-version Chrome currently installed (we
   pin against a recent stable, no need to dynamic-probe), so the publisher
   sees a fingerprint consistent with the cookies it issued.

The Chrome profile itself is NOT opened (Chrome locks its user-data-dir while
running), so this pool is safe to use while Chrome is open.

When ``persistent_profile_path`` is provided, the pool launches a dedicated
Playwright/Chrome user-data directory instead. That path is intended for
authorized interactive publisher access: run visible once, complete SSO or
publisher access checks in that profile, and reuse the same browser state for
later article and supplement fetches.
"""

from __future__ import annotations

import logging
from contextlib import contextmanager
from pathlib import Path
from typing import Any, List, Optional

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


class AuthenticatedBrowserPool:
    """Playwright pool that injects Chrome cookies into the context.

    Mirrors ``BrowserPool`` but with cookie inheritance. Keep API surface
    compatible (``start``, ``new_page``, ``page``, ``close``) so it can be
    swapped in wherever a BrowserPool is expected.
    """

    def __init__(
        self,
        cookies: Optional[List[dict]] = None,
        headless: bool = True,
        slow_mo: int = 0,
        use_chrome_channel: bool = True,
        persistent_profile_path: Optional[str] = None,
        use_stealth: bool = False,
    ):
        self._cookies = cookies or []
        self.headless = headless
        self.slow_mo = slow_mo
        self.use_chrome_channel = use_chrome_channel
        self.persistent_profile_path = persistent_profile_path
        # Stealth backend: use patchright (a patched Playwright drop-in) instead
        # of vanilla Playwright. Empirically (2026-06-05) patchright + the real
        # Chrome channel + HEADFUL clears Wiley's Cloudflare managed challenge;
        # the headless-shell and new-headless are still detected. So enabling
        # stealth forces channel=chrome and headful, and skips the custom UA /
        # headers / init-script (those re-introduce a detectable fingerprint;
        # patchright manages fingerprinting itself).
        self.use_stealth = use_stealth
        self._playwright = None
        self._browser = None
        self._context = None
        self._started = False

    @property
    def started(self) -> bool:
        return self._started

    @property
    def cookie_count(self) -> int:
        return len(self._cookies)

    def start(self) -> None:
        if self._started:
            return
        if self.use_stealth:
            try:
                from patchright.sync_api import sync_playwright
            except ImportError as e:
                raise RuntimeError(
                    "Stealth backend requested but patchright is not installed. "
                    "Install with `pip install patchright && patchright install chromium`."
                ) from e
            # CF managed challenge only clears with real Chrome, headful.
            self.use_chrome_channel = True
            if self.headless:
                logger.warning(
                    "use_stealth forces headful Chrome (headless is CF-detected); "
                    "overriding headless=True -> False."
                )
                self.headless = False
        else:
            try:
                from playwright.sync_api import sync_playwright
            except ImportError as e:
                raise RuntimeError(
                    "Playwright is not installed. Install with "
                    "`pip install playwright && playwright install chromium`."
                ) from e

        self._playwright = sync_playwright().start()
        launch_args = [
            "--disable-blink-features=AutomationControlled",
            "--disable-dev-shm-usage",
            "--no-first-run",
            "--no-default-browser-check",
            "--window-size=1920,1080",
        ]

        # Prefer the installed Chrome binary over the bundled Chromium when
        # available — the fingerprint matches the cookies we just imported
        # (issued to the user's real Chrome). Falls back automatically if
        # Chrome isn't installed.
        launch_kwargs = dict(
            headless=self.headless,
            slow_mo=self.slow_mo,
            args=launch_args,
        )
        if self.use_chrome_channel:
            launch_kwargs["channel"] = "chrome"

        if self.use_stealth:
            # Minimal context: let real Chrome supply its native UA / client
            # hints. A custom UA or extra headers here are a detectable
            # inconsistency that defeats the whole point.
            context_kwargs = dict(
                viewport={"width": 1920, "height": 1080},
                locale="en-US",
                timezone_id="America/New_York",
            )
        else:
            context_kwargs = dict(
                user_agent=(
                    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                    "AppleWebKit/537.36 (KHTML, like Gecko) "
                    "Chrome/131.0.0.0 Safari/537.36"
                ),
                viewport={"width": 1920, "height": 1080},
                locale="en-US",
                timezone_id="America/New_York",
                color_scheme="light",
                # Some publishers gate on Accept-Language and Sec-CH-UA hints — let
                # the channel default to whatever Chrome would send, but make sure
                # Accept-Language is set explicitly.
                extra_http_headers={
                    "Accept-Language": "en-US,en;q=0.9",
                },
            )

        if self.persistent_profile_path:
            profile_dir = Path(self.persistent_profile_path).expanduser()
            profile_dir.mkdir(parents=True, exist_ok=True)
            persistent_kwargs = {**launch_kwargs, **context_kwargs}
            try:
                self._context = self._playwright.chromium.launch_persistent_context(
                    str(profile_dir),
                    **persistent_kwargs,
                )
            except Exception as e:
                if self.use_chrome_channel:
                    logger.warning(
                        "Chrome channel persistent launch failed (%s); falling back to bundled Chromium.",
                        e,
                    )
                    persistent_kwargs.pop("channel", None)
                    self._context = self._playwright.chromium.launch_persistent_context(
                        str(profile_dir),
                        **persistent_kwargs,
                    )
                else:
                    raise
            self._browser = None
            logger.info(
                "AuthenticatedBrowserPool started with persistent profile %s",
                profile_dir,
            )
        else:
            try:
                self._browser = self._playwright.chromium.launch(**launch_kwargs)
            except Exception as e:
                if self.use_chrome_channel:
                    logger.warning(
                        "Chrome channel launch failed (%s); falling back to bundled Chromium.",
                        e,
                    )
                    launch_kwargs.pop("channel", None)
                    self._browser = self._playwright.chromium.launch(**launch_kwargs)
                else:
                    raise

            self._context = self._browser.new_context(**context_kwargs)
        # patchright handles fingerprint hardening itself; our hand-rolled init
        # script can re-introduce detectable patches, so skip it under stealth.
        if not self.use_stealth:
            self._context.add_init_script(_STEALTH_SCRIPT)

        # Inject cookies before any navigation.
        if self._cookies:
            accepted = self._inject_cookies(self._cookies)
            logger.info(
                "AuthenticatedBrowserPool: injected %d/%d cookies",
                accepted,
                len(self._cookies),
            )

        self._started = True
        logger.info(
            "AuthenticatedBrowserPool started (headless=%s, channel=%s, persistent_profile=%s)",
            self.headless,
            "chrome" if self.use_chrome_channel else "chromium",
            bool(self.persistent_profile_path),
        )

    def _inject_cookies(self, cookies: List[dict]) -> int:
        """Inject cookies, retrying individuals if a batch is rejected.

        Playwright rejects an entire add_cookies() call if any cookie is
        malformed. We try the bulk path first and fall back to per-cookie
        injection so a single bad cookie doesn't lose us the rest of the jar.
        """
        try:
            self._context.add_cookies(cookies)
            return len(cookies)
        except Exception as e:
            logger.warning("bulk add_cookies failed (%s); retrying per cookie", e)

        accepted = 0
        for c in cookies:
            try:
                self._context.add_cookies([c])
                accepted += 1
            except Exception as e:
                logger.debug(
                    "rejected cookie %s on %s: %s",
                    c.get("name"),
                    c.get("domain"),
                    e,
                )
        return accepted

    def new_page(self) -> Any:
        if not self._started:
            self.start()
        return self._context.new_page()

    @contextmanager
    def page(self):
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
        logger.info("AuthenticatedBrowserPool closed")

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False
