#!/usr/bin/env python3
"""Verify institutional EZproxy routing clears the Cloudflare wall (Wiley etc.).

Mission-critical sanity check: confirms that, with ``GVF_EZPROXY_PREFIX``/``HOST``
set and a valid EZproxy session cookie loaded, a Wiley Online Library request
returns licensed full text instead of the Cloudflare 403 "Just a moment…" page
OR an expired-session redirect to the SSO login page.

    GVF_EZPROXY_HOST=login.proxy.library.vanderbilt.edu \\
        python scripts/check_ezproxy.py --doi 10.1111/jce.14865

This is a thin CLI over the SAME live probe ``gvf-run`` uses as its institutional
preflight (``cli/institutional_preflight.probe_institutional_access``), so its
verdict cannot disagree with what actually gates a run. In particular it does NOT
false-pass a login/SSO redirect (a large HTTP 200 with no Cloudflare marker) the
way a naive "200 + big body + not-CF" check would.

Exit codes: 2 = EZproxy not configured; 0 = live full-text access confirmed;
1 = configured but the routed request did not return licensed full text.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--doi",
        default=None,
        help="Wiley DOI to probe (default: the preflight canary DOI).",
    )
    args = ap.parse_args()

    # Load .env the same way gvf-run does, so `GVF_EZPROXY_HOST` set in .env is
    # honored when this is run standalone (otherwise it would report "not
    # configured" despite a configured .env).
    try:
        from utils.bootstrap import initialize_runtime

        initialize_runtime()
    except Exception:  # noqa: BLE001 - bootstrap is best-effort for a diagnostic
        pass

    from cli.institutional_preflight import probe_institutional_access

    rpt = probe_institutional_access(doi=args.doi)
    for line in rpt.lines:
        print(f"  {line}")
    if rpt.probe_detail:
        print(f"  probe detail: {rpt.probe_detail}")

    if not rpt.ezproxy_configured:
        print(
            "\nEZproxy is NOT configured. Set GVF_EZPROXY_HOST (e.g. "
            "login.proxy.library.vanderbilt.edu) or GVF_EZPROXY_PREFIX, and log\n"
            "into the library proxy once in Chrome (cookie_loader reads "
            "proxy.library.vanderbilt.edu). Then re-run."
        )
        return 2

    if rpt.should_block:
        print(f"\nRESULT: institutional full-text access is NOT live — {rpt.reason}")
        return 1

    print(
        "\nRESULT: institutional full-text access is LIVE — Cloudflare cleared, "
        "real licensed full text returned."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
