"""Offline unit tests for the gvf-run institutional-access preflight + audit.

No network: the probe short-circuits to a block when EZproxy is unconfigured, and
the audit is pure filesystem. These lock in the two behaviors that matter:

  * the H1 fix — a login/SSO redirect (an expired session's tell) is detected as
    NOT-full-text, which the old scripts/check_ezproxy.py probe accepted; and
  * the audit threshold is high enough that a normal (~30-50% full-text) run is
    not falsely flagged degraded.
"""

from pathlib import Path

from cli.institutional_preflight import (
    _looks_cf,
    _looks_login,
    audit_source_integrity,
    format_block_message,
    probe_institutional_access,
)

_EZPROXY_ENV = (
    "GVF_EZPROXY_PREFIX",
    "GVF_EZPROXY_HOST",
    "PROXY_LOGIN_PREFIX",
    "PROXY_HOST",
)


def test_looks_login_detects_sso_and_proxy_login():
    # The exact false-pass the old probe missed: an expired session redirects here.
    assert _looks_login("https://login.microsoftonline.com/common/oauth2/authorize")
    assert _looks_login("https://login.proxy.library.vanderbilt.edu/login?url=x")
    assert _looks_login("https://shibboleth.example.edu/idp/profile")
    # A real Wiley full-text URL is not a login page.
    assert not _looks_login(
        "https://onlinelibrary.wiley.com/doi/full/10.1111/jce.14865"
    )


def test_looks_cf_detects_challenge_not_real_body():
    assert _looks_cf("Just a moment...", 200, {"server": "cloudflare"})
    assert _looks_cf("blocked", 403, {})
    assert not _looks_cf("A real article body with results and methods.", 200, {})


def test_probe_blocks_when_ezproxy_unconfigured(monkeypatch):
    for k in _EZPROXY_ENV:
        monkeypatch.delenv(k, raising=False)
    rpt = probe_institutional_access()
    assert rpt.should_block is True
    assert rpt.ezproxy_configured is False
    assert rpt.live_probe == "skipped"  # no network attempted
    assert "GVF PREFLIGHT" in format_block_message(rpt)


def _write(p: Path, first_line: str) -> None:
    p.write_text(first_line + "\nbody text\n", encoding="utf-8")


def test_audit_counts_and_high_threshold(tmp_path: Path):
    _write(tmp_path / "1_FULL_CONTEXT.md", "# ABSTRACT-ONLY FALLBACK")
    _write(tmp_path / "2_FULL_CONTEXT.md", "# A real paper")
    (tmp_path / "sub").mkdir()
    _write(tmp_path / "sub" / "3_FULL_CONTEXT.md", "# ABSTRACT-ONLY FALLBACK")

    rpt = audit_source_integrity(tmp_path)
    assert rpt.total == 3
    assert rpt.abstract_only == 2
    assert rpt.full_text == 1
    # 2/3 = 0.667 is normal-ish coverage; must NOT trip the default 0.85 threshold.
    assert rpt.degraded is False
    # A stricter threshold does trip.
    assert audit_source_integrity(tmp_path, threshold=0.5).degraded is True


def test_audit_flags_near_total_failure(tmp_path: Path):
    for i in range(9):
        _write(tmp_path / f"{i}_FULL_CONTEXT.md", "# ABSTRACT-ONLY FALLBACK")
    _write(tmp_path / "9_FULL_CONTEXT.md", "# A real paper")
    rpt = audit_source_integrity(tmp_path)  # 9/10 = 0.9 >= 0.85 default
    assert rpt.degraded is True


def test_audit_empty_dir_is_safe(tmp_path: Path):
    rpt = audit_source_integrity(tmp_path)
    assert rpt.total == 0 and rpt.degraded is False


def test_audit_dedupes_by_pmid_full_text_wins(tmp_path: Path):
    # The same PMID landing in two subdirs (a stub copy + a recovered full-text
    # copy) must count ONCE, as full text — not twice, which would skew the ratio.
    (tmp_path / "pmc_fulltext").mkdir()
    (tmp_path / "recovery").mkdir()
    _write(
        tmp_path / "pmc_fulltext" / "111_FULL_CONTEXT.md", "# ABSTRACT-ONLY FALLBACK"
    )
    _write(tmp_path / "recovery" / "111_FULL_CONTEXT.md", "# Recovered full text")
    _write(
        tmp_path / "pmc_fulltext" / "222_FULL_CONTEXT.md", "# ABSTRACT-ONLY FALLBACK"
    )
    rpt = audit_source_integrity(tmp_path)
    assert rpt.total == 2  # 111 deduped + 222
    assert rpt.full_text == 1  # 111: full text wins over its stub copy
    assert rpt.abstract_only == 1  # 222


def test_env_helpers_are_defensive(monkeypatch):
    # A malformed override must degrade to the default, never raise (this module is
    # imported inside a guarded run path).
    from cli.institutional_preflight import _float_env, _int_env

    monkeypatch.setenv("X_BAD", "notanumber")
    monkeypatch.setenv("X_EMPTY", "")
    assert _int_env("X_BAD", 45) == 45
    assert _int_env("X_EMPTY", 45) == 45
    assert _int_env("X_UNSET", 7) == 7
    assert _float_env("X_BAD", 0.85) == 0.85
    monkeypatch.setenv("X_GOOD", "12")
    assert _int_env("X_GOOD", 45) == 12


def test_audit_none_run_dir_is_safe():
    # run_dir is Optional[Path] at the call site; None must not raise (a bare
    # Path(None).rglob would TypeError, which the OSError handler wouldn't catch).
    rpt = audit_source_integrity(None)
    assert rpt.total == 0 and rpt.degraded is False
