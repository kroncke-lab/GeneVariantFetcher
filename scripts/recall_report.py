#!/usr/bin/env python3
"""Compatibility wrapper for the canonical multi-gene recall runner.

The old version duplicated KCNH2-only recall logic and hardcoded a desktop gold
workbook path. Keep this entry point for existing shell history, but delegate to
``scripts/run_recall_suite.py`` so all recall reporting uses the same comparator
and the same normalized gold-standard package.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.run_recall_suite import main  # noqa: E402


if __name__ == "__main__":
    sys.exit(main())
