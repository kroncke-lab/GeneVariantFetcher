# Contributing to GeneVariantFetcher

## Getting Started

1. Fork and clone the repository
2. Create a virtual environment: `python -m venv venv && source venv/bin/activate`
3. Install in development mode: `pip install -e . && pip install -r gui/requirements.txt`
4. Install pre-commit hooks: `pip install pre-commit && pre-commit install`

## Code Style

- Python 3.11+
- Follow PEP 8 guidelines
- Use type hints for function signatures
- Use `logging.getLogger(__name__)` for logging

## Commit Messages

Use imperative mood:
- "Add feature" not "Added feature"
- "Fix bug" not "Fixed bug"

## Pull Requests

1. Create a feature branch from `main`
2. Make your changes with clear commit messages
3. Run tests: `pytest tests/`
4. Run linting: `ruff check .`
5. Submit PR with description of changes

## Testing

```bash
# Run all tests
pytest tests/

# Run specific test file
pytest tests/test_sqlite_migration.py -v

# Run with coverage
pytest tests/ --cov=pipeline --cov=harvesting
```

## Project Structure

| Directory | Purpose |
|-----------|---------|
| `pipeline/` | Core extraction logic (filters, extraction, aggregation) |
| `harvesting/` | Paper download and SQLite migration |
| `gene_literature/` | PubMed/PMC discovery and synonym finding |
| `gui/` | FastAPI web interface |
| `config/` | Pydantic settings |
| `utils/` | Shared utilities |

## Environment Variables

Copy from CLAUDE.md or create `.env` with required keys. Never commit actual API keys.
