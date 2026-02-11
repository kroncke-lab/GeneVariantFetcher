# Contributing to GeneVariantFetcher

Thank you for your interest in contributing! This guide covers the essentials.

## Development Setup

```bash
# Clone and enter repo
git clone https://github.com/your-org/GeneVariantFetcher.git
cd GeneVariantFetcher

# Create virtual environment (Python 3.11+ required)
python3.11 -m venv venv
source venv/bin/activate

# Install in development mode
pip install -e .
pip install -r gui/requirements.txt

# Install pre-commit hooks
pip install pre-commit && pre-commit install

# Copy environment template and add your API keys
cp .env.example .env
# Edit .env with your Elsevier, Springer, Wiley API keys
```

## Adding a New Gene

1. **Create gene configuration** in `config/genes/`:
   ```python
   # config/genes/YOUR_GENE.py
   GENE_CONFIG = {
       "symbol": "GENE",
       "aliases": ["ALIAS1", "ALIAS2"],
       "uniprot_id": "PXXXXX",
       "reference_transcript": "NM_XXXXXXX.X",
   }
   ```

2. **Add variant alias dictionary** (optional but recommended):
   ```
   resources/variant_aliases/GENE_aliases.tsv
   ```
   Format: `alias\tcanonical_variant` (one per line)

3. **Run discovery** to find papers:
   ```bash
   gvf extract GENE --max-pmids 500 --max-downloads 0
   ```

4. **Test extraction** on a few papers before full run:
   ```bash
   gvf extract GENE --limit 10 --scout-first
   ```

## Extending Extractors

The extraction pipeline lives in `pipeline/`:

| Module | Purpose |
|--------|---------|
| `pipeline/extraction.py` | Main LLM extraction + table regex |
| `utils/variant_scanner.py` | Comprehensive regex variant scanner |
| `utils/variant_normalizer.py` | Canonical normalizer (THE source of truth) |
| `pipeline/data_scout.py` | Pre-extraction quality filters + data zones |
| `pipeline/aggregation.py` | Post-extraction deduplication |

To add a new extraction pattern:
1. Add regex to `utils/variant_scanner.py` in the pattern definitions
2. Test against sample papers in `tests/fixtures/`
3. Ensure pattern output feeds into `variant_normalizer.normalize_variant()`

## Code Style

- **Python 3.11+** with type hints for all function signatures
- **PEP 8** compliance (enforced by ruff)
- **Logging**: Use `logging.getLogger(__name__)`, never `print()` in library code
- **Imports**: Standard library → third-party → local, alphabetized
- **Docstrings**: Google style for public functions

```python
def normalize_variant(raw: str, gene: str = "KCNH2") -> str | None:
    """Normalize a variant string to canonical form.
    
    Args:
        raw: Raw variant string (e.g., "p.A561V", "Ala561Val")
        gene: Gene symbol for context
        
    Returns:
        Canonical form (e.g., "p.Ala561Val") or None if unparseable
    """
```

## Running Tests

```bash
# Run all tests
pytest tests/ -v

# Run specific test file
pytest tests/test_normalizer.py -v

# Run with coverage
pytest tests/ --cov=pipeline --cov=harvesting --cov-report=term-missing

# Run linting
ruff check .
ruff format --check .
```

## Commit Messages

Use imperative mood:
- ✅ "Add fuzzy variant matching"
- ✅ "Fix Excel truncation at 100 rows"
- ❌ "Added feature" / "Fixed bug"

Format: `<type>: <description>`

Types: `feat`, `fix`, `docs`, `refactor`, `test`, `chore`

## Pull Request Process

1. Create a feature branch from `main`
2. Make focused changes with clear commits
3. Add/update tests for new functionality
4. Run full test suite and linting
5. Update CHANGELOG.md under "Unreleased" section
6. Submit PR with description of changes and any breaking changes

## Project Structure

```
GeneVariantFetcher/
├── cli/              # CLI commands (Typer)
├── config/           # Pydantic settings, gene configs
├── gene_literature/  # PubMed/PMC discovery
├── gui/              # FastAPI web interface
├── harvesting/       # Paper download, publisher APIs
├── pipeline/         # Core extraction logic
├── resources/        # Variant aliases, reference data
├── tests/            # Pytest test suite
└── utils/            # Shared utilities
```

## Questions?

Open an issue for bugs or feature requests. For architecture questions, check `CLAUDE.md` for detailed context.
