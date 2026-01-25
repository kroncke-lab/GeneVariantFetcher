"""
CLI tools for GeneVariantFetcher.

This package contains standalone command-line utilities:
- automated_workflow: Full pipeline CLI
- fetch_manager: Semi-manual download helper for paywalled papers
- browser_fetch: Browser-automated download with Cloudflare/CAPTCHA handling
"""

from cli.automated_workflow import automated_variant_extraction_workflow

__all__ = ["automated_variant_extraction_workflow"]
