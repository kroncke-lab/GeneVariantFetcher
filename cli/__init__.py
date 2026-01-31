"""
CLI tools for GeneVariantFetcher.

This package contains standalone command-line utilities:
- automated_workflow: Full pipeline CLI
- fetch_manager: Semi-manual download helper for paywalled papers
- browser_fetch: Browser-automated download with Cloudflare/CAPTCHA handling
- scout: Standalone Data Scout for identifying high-value data zones
"""

from cli.automated_workflow import automated_variant_extraction_workflow
from cli.scout import run_scout

__all__ = ["automated_variant_extraction_workflow", "run_scout"]
