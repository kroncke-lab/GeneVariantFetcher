"""
Harvesting Package

Modular package for harvesting full-text articles and supplemental materials
from PubMed Central. Organized into focused modules for better maintainability.

Public API:
    PMCHarvester - Main class for harvesting operations

Internal modules:
    pmc_api - NCBI/PMC API interactions
    doi_resolver - DOI resolution and routing
    supplement_scraper - Web scraping for supplemental files
    format_converters - File format conversions to markdown
    orchestrator - Main orchestration logic
"""

from .orchestrator import PMCHarvester

__all__ = ['PMCHarvester']
__version__ = '2.0.0'
