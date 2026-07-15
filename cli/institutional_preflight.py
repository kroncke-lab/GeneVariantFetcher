"""Institutional-access preflight guard + post-run source-integrity audit for gvf-run.

Two independent defenses against the silent failure where a full gene harvest
runs at abstract-only coverage but still reports success. (When a paywalled
paper is unreachable, ``harvesting/orchestrator.py`` writes an ABSTRACT-ONLY
stub and returns ``success=True``; the doctor check is advisory and never blocks;
and the source-recovery subprocess treats a nonzero exit as non-fatal. So a run
with a dead EZproxy session churns through and looks like it worked.)

1. :func:`probe_institutional_access` — a LIVE precondition check run *before* a
   full-dataset harvest. Unlike ``gvf_run.institutional_auth_status`` (which is
   presence-only: "is an env var set / is browser_cookie3 importable?"), this
   actually routes a Wiley DOI through the configured EZproxy with the loaded
   Chrome session cookie and verifies the response is real licensed full text —
   NOT an SSO/login redirect (the tell of an expired cookie) and NOT a Cloudflare
   wall. ``gvf-run`` BLOCKS a full run when this fails, unless the operator opts
   out with ``--allow-degraded-institutional`` (or ``GVF_PREFLIGHT_SKIP=1``).

   Why a live probe: the standalone ``scripts/check_ezproxy.py`` probe (and the
   presence checks) FALSE-PASS on an expired session, because an expired cookie
   302-redirects to a login page that is a large HTTP 200 with no Cloudflare
   marker — which their "200 + >4KB + not-CF" test accepts. This probe adds the
   final-URL check (a login/IdP host ⇒ FAIL) that closes that hole.

2. :func:`audit_source_integrity` — a post-run GROUND-TRUTH check. Counts the
   fraction of a run's resolved papers that landed as ABSTRACT-ONLY stubs vs.
   full text and flags the run degraded when that fraction is extreme. This
   catches what a point-in-time preflight structurally cannot: a cookie that
   expires mid-run, or a per-publisher entitlement gap the single Wiley probe
   never exercised. The threshold is deliberately HIGH (near-total failure),
   because a healthy run's full-text success is only ~30-50% (see
   docs/QUICKSTART.md) — so a normal run is already ~50-70% abstract-only and a
   low threshold would false-alarm on every healthy run.

Calibration note: the preflight block is tuned for the VUMC/Vanderbilt setup
(EZproxy + Elsevier insttoken + publisher keys), where a missing/dead EZproxy IS
a real coverage gap worth stopping for. A deliberate non-EZproxy run (pure
open-access, or API-token-only) should use ``--allow-degraded-institutional`` or
``--no-source-recovery``.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

logger = logging.getLogger("gvf_run.preflight")


def _int_env(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name) or default)
    except (TypeError, ValueError):
        return default


def _float_env(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name) or default)
    except (TypeError, ValueError):
        return default


# The default canary DOI is a Wiley (Cloudflare-fronted, subscription) article —
# the publisher class EZproxy exists to unlock. Overridable for when it is
# retracted/relocated/flipped to open access. Env parsing is DEFENSIVE: this
# module is imported inside a guarded run path, so a malformed override value
# must degrade to the default, never raise at import time.
_DEFAULT_PROBE_DOI = os.environ.get("GVF_PREFLIGHT_DOI") or "10.1111/jce.14865"
_DEFAULT_TIMEOUT_S = _int_env("GVF_PREFLIGHT_TIMEOUT_S", 45)

# A response whose FINAL (post-redirect) URL host/path matches one of these is a
# login/SSO landing page, not article full text — the signature of an expired or
# missing institutional session. This is the check scripts/check_ezproxy.py lacks.
_LOGIN_URL_MARKERS = (
    "login",
    "signin",
    "sign-in",
    "idp",
    "shibboleth",
    "/sso",
    "wayf",
    "openathens",
    "authenticate",
    "adfs",
    "okta",
    "microsoftonline",
    "auth.",
)

_CF_MARKERS = ("just a moment", "cf-chl", "challenge-platform", "/cdn-cgi/challenge")

# The header the orchestrator's abstract-only fallback writes as the first line
# of the stub FULL_CONTEXT.md (harvesting/orchestrator.py:_build_abstract_only_fallback).
_ABSTRACT_ONLY_HEADER = "# ABSTRACT-ONLY FALLBACK"

# Default: only flag a run degraded when full-text acquisition has essentially
# collapsed. See the module docstring for why this is high, not 0.5.
_DEFAULT_MAX_ABSTRACT_ONLY_RATIO = _float_env("GVF_MAX_ABSTRACT_ONLY_RATIO", 0.85)


# ---------------------------------------------------------------------------
# Preflight: live institutional-access probe
# ---------------------------------------------------------------------------
@dataclass
class AccessReport:
    """Result of the live institutional-access probe."""

    viable: bool  # a live full-text path was confirmed
    should_block: bool  # a full run should stop unless the operator overrides
    ezproxy_configured: bool = False
    ezproxy_env_var: Optional[str] = None
    cookies_loaded: int = 0
    ez_cookies: int = 0
    cookie_error: Optional[str] = None
    live_probe: str = "skipped"  # "PASS" | "blocked" | "skipped" | "error"
    probe_detail: str = ""
    insttoken: bool = False
    elsevier_key: bool = False
    wiley_key: bool = False
    springer_key: bool = False
    reason: str = ""
    lines: list[str] = field(default_factory=list)


def _looks_cf(text: str, status: int, headers) -> bool:
    low = (text or "")[:4000].lower()
    server = str(headers.get("server", "")).lower() if headers else ""
    if status == 403:
        return True
    if any(m in low for m in _CF_MARKERS):
        return True
    # A Cloudflare-served challenge can also surface as a 429/503 whose
    # interstitial isn't in the first 4KB; treat CF-served throttle statuses as a
    # wall so a routed request isn't mistaken for real licensed content.
    return "cloudflare" in server and status in (429, 503)


def _looks_login(url: str) -> bool:
    return any(m in (url or "").lower() for m in _LOGIN_URL_MARKERS)


def _detect_ezproxy_env() -> Optional[str]:
    for k in (
        "GVF_EZPROXY_PREFIX",
        "GVF_EZPROXY_HOST",
        "PROXY_LOGIN_PREFIX",
        "PROXY_HOST",
    ):
        if os.environ.get(k):
            return k
    return None


def _load_cookies_with_timeout(timeout_s: int):
    """Load Chrome cookies without ever hanging the run.

    Chrome cookie decryption needs the macOS Keychain "Chrome Safe Storage" key; a
    non-interactive / unauthorized context can stall INDEFINITELY there. Delegate
    to ``load_chrome_cookies(timeout_seconds=...)``, which loads each domain in a
    short-lived child process and *terminates* a stuck one — a plain worker thread
    can't be killed, so a ThreadPoolExecutor would still join-and-hang on exit.
    Returns ``(cookies, error)``.
    """
    try:
        from harvesting.browser_html.cookie_loader import load_chrome_cookies

        # Per-domain cap so a stalled Keychain lookup is killed, not awaited.
        per_domain = min(timeout_s, 20)
        return list(load_chrome_cookies(timeout_seconds=per_domain) or []), None
    except Exception as exc:  # noqa: BLE001
        return [], f"cookie load failed: {str(exc)[:160]}"


def _is_ez_cookie(c: dict) -> bool:
    name = (c.get("name") or "").lower()
    dom = (c.get("domain") or "").lower()
    return name.startswith("ezproxy") or "proxy.library." in dom


def probe_institutional_access(
    *, timeout_s: Optional[int] = None, doi: Optional[str] = None
) -> AccessReport:
    """Live check of whether authenticated full-text recovery actually works.

    Returns an :class:`AccessReport`. ``should_block`` is True when a full-dataset
    run would silently degrade to abstract-only. All imports are lazy so importing
    this module stays cheap and offline-safe.
    """
    timeout_s = timeout_s or _DEFAULT_TIMEOUT_S
    doi = doi or _DEFAULT_PROBE_DOI

    insttoken = bool(os.environ.get("ELSEVIER_INSTTOKEN"))
    elsevier_key = bool(os.environ.get("ELSEVIER_API_KEY"))
    wiley_key = bool(os.environ.get("WILEY_API_KEY"))
    springer_key = bool(os.environ.get("SPRINGER_API_KEY"))

    rpt = AccessReport(
        viable=False,
        should_block=True,
        insttoken=insttoken,
        elsevier_key=elsevier_key,
        wiley_key=wiley_key,
        springer_key=springer_key,
        ezproxy_env_var=_detect_ezproxy_env(),
    )

    try:
        from harvesting.browser_html import ezproxy

        rpt.ezproxy_configured = ezproxy.is_configured()
    except Exception as exc:  # noqa: BLE001
        rpt.ezproxy_configured = False
        rpt.probe_detail = f"ezproxy import failed: {str(exc)[:120]}"

    # 1. EZproxy not wired at all → no automated Wiley/Karger/Sage bypass.
    if not rpt.ezproxy_configured:
        rpt.live_probe = "skipped"
        rpt.reason = (
            "EZproxy is not configured (no GVF_EZPROXY_PREFIX/HOST) — the "
            "automated Cloudflare bypass for Wiley/Karger/Sage is unavailable."
        )
        rpt.should_block = True
        rpt.lines = _summary_lines(rpt)
        return rpt

    # 2. Load the session cookie (with a hard timeout so a Keychain stall can't
    #    hang the run). No cookie ⇒ no live session ⇒ block.
    cookies, cookie_error = _load_cookies_with_timeout(timeout_s)
    rpt.cookies_loaded = len(cookies)
    rpt.ez_cookies = sum(1 for c in cookies if _is_ez_cookie(c))
    rpt.cookie_error = cookie_error

    if cookie_error:
        rpt.live_probe = "error"
        rpt.reason = (
            f"could not read Chrome cookies ({cookie_error}). On macOS this is "
            "usually a Keychain access problem — run gvf-run in your own logged-in "
            "Terminal and grant access to 'Chrome Safe Storage'."
        )
        rpt.should_block = True
        rpt.lines = _summary_lines(rpt)
        return rpt

    if rpt.ez_cookies == 0:
        rpt.live_probe = "blocked"
        rpt.reason = (
            "no EZproxy session cookie found in Chrome — log into the library "
            "proxy once in Chrome (VUMC users: sign in with your @vumc.org email)."
        )
        rpt.should_block = True
        rpt.lines = _summary_lines(rpt)
        return rpt

    # 3. LIVE probe: route a Wiley DOI through EZproxy and confirm real full text.
    try:
        from scripts.fetch_paywalled import (
            make_session,
            hydrate_session_with_browser_cookies,
        )

        session = make_session()  # EZproxy rewriter already installed
        hydrate_session_with_browser_cookies(session, cookies)
        url = f"https://onlinelibrary.wiley.com/doi/full/{doi}"
        r = session.get(url, timeout=timeout_s, allow_redirects=True)
        final_url = str(getattr(r, "url", "") or "")
        cf = _looks_cf(r.text, r.status_code, r.headers)
        login = _looks_login(final_url)
        ok = (r.status_code == 200) and not cf and not login and len(r.text) > 4000
        rpt.live_probe = "PASS" if ok else "blocked"
        rpt.viable = ok
        rpt.should_block = not ok
        rpt.probe_detail = (
            f"status={r.status_code} cf={cf} login_redirect={login} "
            f"bytes={len(r.text)} final={final_url[:120]}"
        )
        if not ok:
            if login:
                rpt.reason = (
                    "the EZproxy request redirected to a login/SSO page — the "
                    "session cookie is expired. Re-log into the proxy in Chrome."
                )
            elif cf:
                rpt.reason = "the request hit a Cloudflare challenge (routing not clearing the wall)."
            else:
                rpt.reason = "the routed request did not return licensed full text."
    except Exception as exc:  # noqa: BLE001
        rpt.live_probe = "error"
        rpt.viable = False
        rpt.should_block = True
        rpt.reason = f"live EZproxy probe errored: {str(exc)[:160]}"
        rpt.probe_detail = rpt.reason

    rpt.lines = _summary_lines(rpt)
    return rpt


def _summary_lines(rpt: AccessReport) -> list[str]:
    yn = lambda b: "yes" if b else "NO"  # noqa: E731
    return [
        f"EZproxy configured ({rpt.ezproxy_env_var or 'GVF_EZPROXY_*'}) : {yn(rpt.ezproxy_configured)}",
        f"Institutional session cookies loaded          : {rpt.ez_cookies}"
        + (f"  (of {rpt.cookies_loaded} total)" if rpt.cookies_loaded else "")
        + (f"  [{rpt.cookie_error}]" if rpt.cookie_error else ""),
        f"Live EZproxy full-text probe                  : {rpt.live_probe}"
        + (
            f"  ({rpt.probe_detail})"
            if rpt.probe_detail and rpt.live_probe != "PASS"
            else ""
        ),
        f"ELSEVIER_INSTTOKEN present                    : {yn(rpt.insttoken)}",
        f"WILEY_API_KEY / SPRINGER_API_KEY present      : {yn(rpt.wiley_key)} / {yn(rpt.springer_key)}",
    ]


def format_block_message(
    rpt: AccessReport, *, allow_flag: str = "--allow-degraded-institutional"
) -> str:
    """The operator-facing block message (also used, softened, for the override warning)."""
    body = "\n".join(f"  - {ln}" for ln in rpt.lines)
    return (
        "\n"
        "==================================================================\n"
        "  GVF PREFLIGHT: INSTITUTIONAL ACCESS DEGRADED — RUN HALTED\n"
        "==================================================================\n"
        "A full-dataset harvest was requested, but institutional/paywall access\n"
        "is NOT live. Continuing would SILENTLY produce an abstract-only run\n"
        "(paywalled papers degrade to abstract stubs and the run still reports\n"
        "success), losing most full-text coverage — Wiley, Karger, Sage,\n"
        "ScienceDirect.\n"
        f"\nWhy: {rpt.reason}\n"
        "\nDetected:\n"
        f"{body}\n"
        "\nTo fix (institutional path):\n"
        "  1. In .env set GVF_EZPROXY_HOST=login.proxy.library.vanderbilt.edu\n"
        "     (use the login. subdomain — the apex cert is not valid).\n"
        "  2. Log into the library EZproxy ONCE in Chrome via Vanderbilt SSO\n"
        "     using your @vumc.org identity. Non-@vumc.org logins are denied.\n"
        "  3. Ensure Chrome cookies are readable by GVF (grant Keychain access\n"
        "     to 'Chrome Safe Storage'; update browser-cookie3 if needed).\n"
        "  4. Verify:  .venv/bin/python scripts/check_ezproxy.py   (must exit 0)\n"
        "\nThen re-run. To proceed ANYWAY with reduced coverage (not recommended),\n"
        f"re-run with {allow_flag}  (or set GVF_PREFLIGHT_SKIP=1).\n"
        "For a deliberate fast/calibration run, use --no-source-recovery or\n"
        "--pmid-file, which skip this check.\n"
        "==================================================================\n"
    )


# ---------------------------------------------------------------------------
# Post-run: source-integrity audit (ground truth)
# ---------------------------------------------------------------------------
@dataclass
class IntegrityReport:
    """How much of the run actually landed as full text vs. abstract-only stubs."""

    full_text: int = 0
    abstract_only: int = 0
    total: int = 0
    ratio: float = 0.0  # abstract_only / total
    threshold: float = _DEFAULT_MAX_ABSTRACT_ONLY_RATIO
    degraded: bool = False
    message: str = ""


def _is_abstract_only_stub(path: Path) -> bool:
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:  # first non-empty line
                s = line.strip()
                if s:
                    return s.startswith(_ABSTRACT_ONLY_HEADER)
        return False
    except OSError:
        return False


def audit_source_integrity(
    run_dir: Path, *, threshold: Optional[float] = None
) -> IntegrityReport:
    """Scan a run dir's ``*_FULL_CONTEXT.md`` files and measure the stub ratio.

    Counts each paper's resolved full-text file and classifies it as full text or
    an ABSTRACT-ONLY stub (by its header). ``degraded`` is True only when the stub
    ratio meets/exceeds ``threshold`` (default 0.85 — near-total full-text
    failure), so a normal ~30-50%-full-text run is not flagged.
    """
    threshold = _DEFAULT_MAX_ABSTRACT_ONLY_RATIO if threshold is None else threshold
    rpt = IntegrityReport(threshold=threshold)
    try:
        files = list(Path(run_dir).rglob("*_FULL_CONTEXT.md"))
    except OSError as exc:
        rpt.message = f"integrity audit could not scan {run_dir}: {exc}"
        return rpt

    # Dedupe by PMID (the filename stem before ``_FULL_CONTEXT.md``): the same
    # paper can appear under more than one subdir (e.g. pmc_fulltext plus a
    # recovery dir). Count each PMID once and let full text win over a stub copy —
    # we are measuring "did this paper ultimately land full text?", per paper.
    suffix = "_FULL_CONTEXT.md"
    stub_by_pmid: dict = {}
    for f in files:
        pmid = f.name[: -len(suffix)] or f.name
        stub = _is_abstract_only_stub(f)
        stub_by_pmid[pmid] = (
            stub if pmid not in stub_by_pmid else (stub_by_pmid[pmid] and stub)
        )

    for is_stub in stub_by_pmid.values():
        if is_stub:
            rpt.abstract_only += 1
        else:
            rpt.full_text += 1
    rpt.total = rpt.full_text + rpt.abstract_only
    rpt.ratio = (rpt.abstract_only / rpt.total) if rpt.total else 0.0
    rpt.degraded = rpt.total > 0 and rpt.ratio >= threshold

    pct = rpt.ratio * 100
    if rpt.total == 0:
        rpt.message = "source integrity: no resolved full-text files found to audit."
    elif rpt.degraded:
        rpt.message = (
            f"source integrity DEGRADED: {rpt.abstract_only}/{rpt.total} "
            f"({pct:.0f}%) of resolved papers are ABSTRACT-ONLY stubs "
            f"(>= {threshold * 100:.0f}% threshold) — full-text acquisition "
            "essentially failed. Institutional access was likely down during "
            "recovery; fix access (scripts/check_ezproxy.py) and re-run."
        )
    else:
        rpt.message = (
            f"source integrity OK: {rpt.full_text}/{rpt.total} "
            f"({100 - pct:.0f}%) resolved papers have full text "
            f"({rpt.abstract_only} abstract-only)."
        )
    return rpt
