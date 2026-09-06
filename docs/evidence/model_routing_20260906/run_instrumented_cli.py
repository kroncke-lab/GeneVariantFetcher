"""Run the unchanged production CLI under the campaign-only budget wrapper."""

import runpy
import sys
from pathlib import Path

sys.path.insert(0, str(Path.cwd()))
from budget_guard import install

install()
runpy.run_module("cli", run_name="__main__")
