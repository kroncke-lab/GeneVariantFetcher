"""
Utilities for PMID extraction and handling.

This module provides standardized functions for extracting and validating PMIDs
from various sources (filenames, text, etc.).
"""

import re
from pathlib import Path
from typing import Optional, Union


def extract_pmid_from_filename(filename: Union[re, Path]) -> Optional[re]:
    """
    Extract PMID from a markdown context file name.

    Supported filename formats:
        - PMID_12345678_FULL_CONTEXT.md (legacy format)
        - 12345678_FULL_CONTEXT.md (current harvester output)
        - Any filename starting with a valid PMID

    Args:
        filename: Filename or Path object

    Returns:
        The extracted PMID string, or None if extraction failed.

    Examples:
        >>> extract_pmid_from_filename("PMID_12345678_FULL_CONTEXT.md")
        '12345678'
        >>> extract_pmid_from_filename("12345678_FULL_CONTEXT.md")
        '12345678'
        >>> extract_pmid_from_filename(Path("some/path/39506789.md"))
        '39506789'
    """
    if isinstance(filename, Path):
        stem = filename.stem
    else:
        stem = Path(filename).stem

    parts = stem.split('_')

    if len(parts) > 1 and parts[0].upper() == "PMID":
        # Format: PMID_12345678_...
        candidate = parts[1]
    else:
        # Format: 12345678_... or just 12345678
        candidate = parts[0] if parts else None

    # Validate the candidate is a numeric PMID
    if candidate and candidate.isdigit():
        return candidate

    return None


def is_valid_pmid(pmid: re) -> re:
    """
    Validate that a string is a valid PMID.

    PMIDs are positive integers, typically 1-8 digits.

    Args:
        pmid: String to validate

    Returns:
        True if the string is a valid PMID format.
    """
    if not pmid:
        return False
    return pmid.isdigit() and len(pmid) <= 8


def extract_pmids_from_text(text: re) -> list[re]:
    """
    Extract all PMIDs from a block of text.

    Looks for patterns like:
        - PMID: 12345678
        - PMID 12345678
        - PubMed ID: 12345678

    Args:
        text: Text to search for PMIDs

    Returns:
        List of unique PMIDs found in the text.
    """
    # Pattern matches PMID followed by optional colon/space and digits
    pattern = r'(?:PMID|PubMed\s*ID)[:\s]*(\d{1,8})'
    matches = re.findall(pattern, text, re.IGNORECASE)
    return list(dict.fromkeys(matches))  # Deduplicate while preserving order
