# Proposed Pipeline Package Refactor

## Proposed Directory Tree
```
GeneVariantFetcher/
├── README.md
├── pyproject.toml
├── requirements.txt
├── requirements_harvest.txt
├── config/
│   ├── __init__.py
│   └── settings.py
├── pipeline/
│   ├── __init__.py
│   ├── cli.py
│   ├── sourcing.py
│   ├── filters.py
│   ├── harvesting.py
│   ├── extraction.py
│   ├── aggregation.py
│   └── utils/
│       ├── __init__.py
│       ├── io.py
│       ├── parsing.py
│       └── rate_limits.py
├── docs/
│   └── index.md
├── tests/
│   ├── unit/
│   │   ├── test_filters.py
│   │   ├── test_parsers.py
│   │   └── test_utils.py
│   └── integration/
│       ├── test_pipeline_end_to_end.py
│       └── test_harvest_to_extract.py
└── scripts/
    └── dev_seeds.py
```

## CLI Entry Point (Typer skeleton)
```python
# pipeline/cli.py
import typer
from pipeline.sourcing import PaperSourcer
from pipeline.harvesting import PMCHarvester

app = typer.Typer(help="Gene Variant Fetcher pipeline CLI")


@app.command()
def source(query: str, limit: int = typer.Option(50, help="Max records to fetch")):
    """Run the sourcing stage and persist candidate papers."""
    sourcer = PaperSourcer()
    results = sourcer.run(query=query, limit=limit)
    typer.echo(f"Sourced {len(results)} papers for query '{query}'")


@app.command()
def harvest(input_path: str, output_path: str = typer.Option("data/harvests", help="Destination directory")):
    """Harvest full text for previously sourced papers."""
    harvester = PMCHarvester()
    count = harvester.run(input_path=input_path, output_path=output_path)
    typer.echo(f"Harvested {count} articles into {output_path}")


if __name__ == "__main__":
    app()
```

## Configuration Strategy (singleton loader)
```python
# config/settings.py
import os
from functools import lru_cache
from pathlib import Path
from dotenv import load_dotenv

_ENV_PATH = Path(__file__).resolve().parent.parent / ".env"
load_dotenv(_ENV_PATH)

class Settings:
    """Centralized application settings loaded from environment variables."""

    api_key: str
    rate_limit_per_minute: int

    def __init__(self):
        self.api_key = os.getenv("API_KEY", "")
        self.rate_limit_per_minute = int(os.getenv("RATE_LIMIT_PER_MINUTE", 60))


@lru_cache(maxsize=1)
def get_settings() -> Settings:
    """Return a cached Settings instance so configuration is loaded once."""
    return Settings()
```

## Migration Plan
1. Create the `pipeline/` package with submodules mirroring the architecture (sourcing, filtering, harvesting, extraction, aggregation) and move existing logic from flat modules into their respective files. Update `__init__.py` exports to preserve public interfaces.
2. Introduce `pipeline/cli.py` using Typer and replace `example_*.py` scripts with subcommands; add a console script entry in `pyproject.toml` pointing to `pipeline.cli:app`.
3. Add the `config/` package and `.env` loading via `config/settings.py`; replace direct environment lookups in code with `get_settings()` and update modules to accept settings via dependency injection where appropriate.
4. Refactor imports throughout the codebase to use the new package paths (e.g., `from pipeline.sourcing import PaperSourcer`), and provide backward-compatible shim imports in legacy modules if necessary during transition.
5. Establish `tests/unit/` and `tests/integration/` directories; move or create tests aligned with the new module locations, updating test imports accordingly.
6. Update `ARCHITECTURE.md` (or `docs/index.md`) to map each architectural stage to its concrete class and CLI command for discoverability, ensuring README points to the unified CLI.
