"""
Supplement Reference Parser Module

Parses main text for references to supplementary materials.
Used to:
1. Know what supplements to expect
2. Detect when supplements weren't downloaded
3. Extract explicit supplement URLs from text
"""

import re
from typing import Dict, Any, List, Optional


def parse_supplement_references(text: str) -> Dict[str, Any]:
    """
    Parse main text for references to supplementary materials.

    Args:
        text: Main article text (markdown or plain text)

    Returns:
        Dict with keys:
        - expected_tables: int (number of distinct supplement tables referenced)
        - expected_figures: int (number of distinct supplement figures referenced)
        - table_refs: List[str] (specific table references found)
        - figure_refs: List[str] (specific figure references found)
        - mentioned_variants: Optional[int] (variant count if explicitly mentioned near supplement ref)
        - supplement_urls: List[str] (explicit URLs found)
        - raw_refs: List[str] (all raw supplement references found)
    """
    result = {
        "expected_tables": 0,
        "expected_figures": 0,
        "table_refs": [],
        "figure_refs": [],
        "mentioned_variants": None,
        "supplement_urls": [],
        "raw_refs": [],
    }

    # Patterns for supplement table references
    table_patterns = [
        # "online suppl. table 1", "online supplementary table S1"
        r"online\s+suppl(?:ementary)?\.?\s+table\s+[S]?(\d+)",
        # "supplementary table 1", "supplemental table S1"
        r"supplement(?:ary|al)\s+table\s+[S]?(\d+)",
        # "Table S1", "Tables S1-S3"
        r"\bTable\s+S(\d+)",
        # "additional file 1: table"
        r"additional\s+file\s+(\d+).*?table",
        # "supporting information table S1"
        r"supporting\s+information\s+table\s+[S]?(\d+)",
        # "(suppl. table 1)"
        r"\(suppl\.?\s+table\s+[S]?(\d+)\)",
    ]

    # Patterns for supplement figure references
    figure_patterns = [
        r"online\s+suppl(?:ementary)?\.?\s+fig(?:ure)?\.?\s+[S]?(\d+)",
        r"supplement(?:ary|al)\s+fig(?:ure)?\.?\s+[S]?(\d+)",
        r"\bFig(?:ure)?\.?\s+S(\d+)",
        r"supporting\s+information\s+fig(?:ure)?\.?\s+[S]?(\d+)",
    ]

    # Extract table references
    table_numbers = set()
    for pattern in table_patterns:
        matches = re.findall(pattern, text, re.IGNORECASE)
        for match in matches:
            table_numbers.add(int(match))
            result["table_refs"].append(f"Table S{match}")

    result["expected_tables"] = len(table_numbers)
    result["table_refs"] = list(set(result["table_refs"]))  # Deduplicate

    # Extract figure references
    figure_numbers = set()
    for pattern in figure_patterns:
        matches = re.findall(pattern, text, re.IGNORECASE)
        for match in matches:
            figure_numbers.add(int(match))
            result["figure_refs"].append(f"Figure S{match}")

    result["expected_figures"] = len(figure_numbers)
    result["figure_refs"] = list(set(result["figure_refs"]))

    # Look for variant counts mentioned near supplement references
    # Pattern: "N variants/mutations in [supplement reference]"
    variant_count_patterns = [
        r"(\d+)\s+(?:KCNH2\s+)?(?:variants?|mutations?)\s+(?:in|are\s+listed\s+in|shown\s+in).*?(?:suppl|table\s+S)",
        r"(?:suppl|table\s+S).*?(?:lists?|contains?|shows?)\s+(\d+)\s+(?:variants?|mutations?)",
        r"(\d+)\s+(?:variants?|mutations?).*?\(.*?suppl",
    ]

    for pattern in variant_count_patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            try:
                count = int(match.group(1))
                if count > 0 and count < 10000:  # Sanity check
                    result["mentioned_variants"] = count
                    break
            except (ValueError, IndexError):
                continue

    # Extract explicit supplement URLs
    url_patterns = [
        r'https?://[^\s<>"]+(?:suppl|supplement)[^\s<>"]*',
        r'www\.[^\s<>"]+(?:suppl|supplement)[^\s<>"]*',
        r'(?:karger\.com|doi\.org)/[^\s<>"]*suppl[^\s<>"]*',
    ]

    for pattern in url_patterns:
        urls = re.findall(pattern, text, re.IGNORECASE)
        for url in urls:
            # Clean up URL
            url = url.rstrip(".,;:)")
            if not url.startswith("http"):
                url = "https://" + url
            if url not in result["supplement_urls"]:
                result["supplement_urls"].append(url)

    # Collect all raw references for debugging
    all_ref_patterns = [
        r"online\s+suppl(?:ementary)?\.?\s+(?:table|fig|material|info)[^.]*",
        r"supplement(?:ary|al)\s+(?:table|fig|material|data)[^.]*",
        r"supporting\s+information[^.]*",
        r"additional\s+file\s+\d+[^.]*",
    ]

    for pattern in all_ref_patterns:
        refs = re.findall(pattern, text, re.IGNORECASE)
        result["raw_refs"].extend(refs[:10])  # Limit to avoid huge lists

    result["raw_refs"] = list(set(result["raw_refs"]))[:20]  # Deduplicate and limit

    return result


def extract_supplement_urls_from_text(text: str) -> List[str]:
    """
    Extract all potential supplement URLs from article text.

    Args:
        text: Article text (markdown or plain text)

    Returns:
        List of supplement URLs found
    """
    urls = []

    patterns = [
        # Karger supplements
        r'(?:https?://)?(?:www\.)?karger\.com/doi/[^\s<>"]+',
        # General DOI-based supplements
        r'(?:https?://)?doi\.org/[^\s<>"]+suppl[^\s<>"]*',
        # Direct supplement links
        r'https?://[^\s<>"]+/suppl(?:ement)?(?:ary)?[^\s<>"]*',
        # ScienceDirect MMC
        r'https?://[^\s<>"]+/mmc\d+[^\s<>"]*',
    ]

    for pattern in patterns:
        matches = re.findall(pattern, text, re.IGNORECASE)
        for url in matches:
            url = url.rstrip(".,;:)")
            if not url.startswith("http"):
                url = "https://" + url
            if url not in urls:
                urls.append(url)

    return urls


def check_supplement_gap(
    text: str, downloaded_count: int, extracted_variant_count: int
) -> Dict[str, Any]:
    """
    Check if there's likely a gap between expected and actual supplements.

    Args:
        text: Main article text
        downloaded_count: Number of supplement files actually downloaded
        extracted_variant_count: Number of variants extracted

    Returns:
        Dict with gap analysis:
        - has_gap: bool
        - expected_tables: int
        - downloaded_count: int
        - mentioned_variants: Optional[int]
        - extracted_variants: int
        - warnings: List[str]
    """
    refs = parse_supplement_references(text)

    result = {
        "has_gap": False,
        "expected_tables": refs["expected_tables"],
        "downloaded_count": downloaded_count,
        "mentioned_variants": refs["mentioned_variants"],
        "extracted_variants": extracted_variant_count,
        "warnings": [],
    }

    # Check table count mismatch
    if refs["expected_tables"] > 0 and downloaded_count < refs["expected_tables"]:
        result["has_gap"] = True
        result["warnings"].append(
            f"Main text references {refs['expected_tables']} supplementary table(s) "
            f"but only {downloaded_count} were downloaded"
        )

    # Check variant count mismatch
    if refs["mentioned_variants"] is not None:
        if extracted_variant_count < refs["mentioned_variants"] * 0.5:
            result["has_gap"] = True
            result["warnings"].append(
                f"Paper mentions {refs['mentioned_variants']} variants but "
                f"extraction found only {extracted_variant_count} - likely missing supplements"
            )

    # Check for explicit URLs that might not have been followed
    if refs["supplement_urls"]:
        result["warnings"].append(
            f"Found {len(refs['supplement_urls'])} explicit supplement URL(s) in text: "
            f"{refs['supplement_urls'][:3]}"
        )

    return result


if __name__ == "__main__":
    # Test with example text
    test_text = """
    We identified 54 KCNH2 mutations (online suppl. table 1, 
    www.karger.com/doi/10.1159/000440608 for all online material).
    The mutation spectrum is shown in supplementary figure S1.
    See also Table S2 for detailed phenotype data.
    """

    result = parse_supplement_references(test_text)
    print("Parsed references:")
    for key, value in result.items():
        print(f"  {key}: {value}")

    # Test gap detection
    gap = check_supplement_gap(test_text, downloaded_count=0, extracted_variant_count=4)
    print("\nGap analysis:")
    for key, value in gap.items():
        print(f"  {key}: {value}")
