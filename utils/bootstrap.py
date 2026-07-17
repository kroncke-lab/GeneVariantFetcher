"""Runtime bootstrap utilities for environment and logging initialization."""

from __future__ import annotations

import logging
import os
import threading
from pathlib import Path
from typing import Optional

from dotenv import load_dotenv

LLM_PROVIDER_KEY_ENV_VARS = (
    "OPENAI_API_KEY",
    "AZURE_AI_API_KEY",
    "ANTHROPIC_API_KEY",
)

_BOOTSTRAP_LOCK = threading.Lock()
_BOOTSTRAPPED = False
_logger = logging.getLogger(__name__)


def initialize_runtime(
    *,
    env_file: Optional[Path] = None,
    extra_env_var: str = "GVF_EXTRA_ENV_FILE",
) -> None:
    """Initialize process runtime exactly once.

    Loads `.env` from project root and optionally one extra env file from an
    environment variable. This avoids module-import side effects and hardcoded
    machine-specific paths.
    """
    global _BOOTSTRAPPED
    if _BOOTSTRAPPED:
        return

    with _BOOTSTRAP_LOCK:
        if _BOOTSTRAPPED:
            return

        project_root = Path(__file__).resolve().parent.parent
        primary_env = env_file or (project_root / ".env")
        if primary_env.exists():
            load_dotenv(primary_env, override=False)

        extra_env = os.getenv(extra_env_var)
        if extra_env:
            extra_path = Path(extra_env)
            if extra_path.exists():
                load_dotenv(extra_path, override=False)
            else:
                _logger.warning(
                    "Configured extra env file does not exist (%s=%s)",
                    extra_env_var,
                    extra_env,
                )

        _apply_network_hardening()

        _BOOTSTRAPPED = True


_TRUTHY = {"1", "true", "yes", "on"}


def _apply_network_hardening() -> None:
    """Opt-in resilience for flaky networks (default: no-op).

    On some hosts the IPv6 path to external services blackholes: TCP connections
    establish but reads hang forever, because most HTTP client calls set no read
    timeout. These env-gated knobs let an operator force IPv4 and/or install a
    process-wide default socket timeout so a stalled read fails fast and the
    caller's retry logic recovers instead of hanging indefinitely. Both are unset
    by default, so behavior is unchanged unless explicitly enabled.

    - ``GVF_FORCE_IPV4=1``       resolve hostnames to IPv4 (AF_INET) only.
    - ``GVF_SOCKET_TIMEOUT=<s>`` ``socket.setdefaulttimeout(<s>)`` for calls that
      do not set their own timeout. NOTE: httpx-based LLM clients set their own
      per-request timeouts and are unaffected by this default.
    """
    import socket

    if os.getenv("GVF_FORCE_IPV4", "").strip().lower() in _TRUTHY:
        _orig_getaddrinfo = socket.getaddrinfo

        def _ipv4_only(host, port, family=0, type=0, proto=0, flags=0):
            return _orig_getaddrinfo(host, port, socket.AF_INET, type, proto, flags)

        socket.getaddrinfo = _ipv4_only
        _logger.info("GVF_FORCE_IPV4 enabled: resolving hostnames to IPv4 only")

    raw_timeout = os.getenv("GVF_SOCKET_TIMEOUT", "").strip()
    if raw_timeout:
        try:
            timeout = float(raw_timeout)
        except ValueError:
            _logger.warning(
                "Ignoring invalid GVF_SOCKET_TIMEOUT=%r (not a number)", raw_timeout
            )
        else:
            if timeout > 0:
                socket.setdefaulttimeout(timeout)
                _logger.info(
                    "GVF_SOCKET_TIMEOUT enabled: default socket timeout %.0fs", timeout
                )


def llm_provider_key_status() -> dict[str, bool]:
    """Return which supported LLM provider keys are present in the process env."""

    return {key: bool(os.getenv(key, "").strip()) for key in LLM_PROVIDER_KEY_ENV_VARS}


def has_llm_provider_key() -> bool:
    """True when any supported LLM provider credential is available."""

    return any(llm_provider_key_status().values())
