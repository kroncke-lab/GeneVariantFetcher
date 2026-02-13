"""Runtime bootstrap utilities for environment and logging initialization."""

from __future__ import annotations

import logging
import os
import threading
from pathlib import Path
from typing import Optional

from dotenv import load_dotenv

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
