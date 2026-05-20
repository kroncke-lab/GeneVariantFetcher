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

        _BOOTSTRAPPED = True


def llm_provider_key_status() -> dict[str, bool]:
    """Return which supported LLM provider keys are present in the process env."""

    return {key: bool(os.getenv(key, "").strip()) for key in LLM_PROVIDER_KEY_ENV_VARS}


def has_llm_provider_key() -> bool:
    """True when any supported LLM provider credential is available."""

    return any(llm_provider_key_status().values())
