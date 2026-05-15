"""Shared fixtures for GVF recall tests."""

import os

import pytest
import pandas as pd
import re
import glob
from pathlib import Path
from collections import defaultdict

REPO_ROOT = Path(__file__).resolve().parents[2]

DEFAULT_GOLD_STANDARD_PATH = (
    REPO_ROOT
    / "comparison_results"
    / "KCNH2 HeterozygoteDatabase-4-clinical_w_paris_japan_mayo_and_italycohort.xls"
)

GOLD_STANDARD_PATH = Path(
    os.environ.get("GVF_GOLD_STANDARD_PATH", str(DEFAULT_GOLD_STANDARD_PATH))
).expanduser()
DOWNLOAD_DIR = Path(
    os.environ.get("GVF_RECALL_DOWNLOAD_DIR", str(REPO_ROOT / "results"))
).expanduser()
OUTPUT_DIR = Path(
    os.environ.get("GVF_RECALL_OUTPUT_DIR", str(REPO_ROOT / "results"))
).expanduser()


def normalize_variant(v):
    """Normalize variant nomenclature for comparison."""
    if pd.isna(v):
        return None
    v = str(v).upper().strip()
    v = v.replace("*", "X")
    v = v.replace("P.", "")
    # Remove p. prefix if present
    if v.startswith("P."):
        v = v[2:]
    # Handle frameshift variants
    if "FS" in v:
        v = re.sub(r"FS.*$", "FSX", v)
    return v


@pytest.fixture(scope="session")
def gold_standard():
    """Load and normalize gold standard data."""
    if not GOLD_STANDARD_PATH.exists():
        pytest.skip(
            "Gold standard workbook not found; set GVF_GOLD_STANDARD_PATH "
            "to run recall tests."
        )
    df = pd.read_excel(GOLD_STANDARD_PATH)

    # Normalize PMIDs
    df["PMID_int"] = df["PMID"].apply(lambda x: int(float(x)) if pd.notna(x) else None)

    # Normalize variants
    df["Variant_norm"] = df["Variant"].apply(normalize_variant)

    # Flag LQT-affected (LQT >= 1 means affected)
    df["is_LQT_affected"] = df["LQT"].apply(lambda x: x >= 1 if pd.notna(x) else False)

    pmids = set(df["PMID_int"].dropna().astype(int))
    variants = set(df["Variant_norm"].dropna())

    return {
        "df": df,
        "pmids": pmids,
        "variants": variants,
        "total_carriers": len(df),
        "carriers_with_pmid": len(df[df["PMID_int"].notna()]),
        "total_lqt_affected": int(df["is_LQT_affected"].sum()),
    }


@pytest.fixture(scope="session")
def downloaded_pmids():
    """Get set of downloaded PMIDs from PDF filenames."""
    if not DOWNLOAD_DIR.exists():
        pytest.skip(
            "Recall download directory not found; set GVF_RECALL_DOWNLOAD_DIR "
            "to run recall tests."
        )
    pdfs = glob.glob(str(DOWNLOAD_DIR / "*.pdf"))
    if not pdfs:
        pytest.skip(
            "No recall download PDFs found; set GVF_RECALL_DOWNLOAD_DIR "
            "to a populated verified download directory."
        )
    pmids = set()

    for pdf in pdfs:
        basename = Path(pdf).name
        # Extract PMID from filename (assumes PMID_*.pdf or PMID.pdf format)
        match = re.match(r"^(\d+)", basename)
        if match:
            pmids.add(int(match.group(1)))

    return pmids


@pytest.fixture(scope="session")
def downloaded_papers():
    """Get detailed info about downloaded papers."""
    if not DOWNLOAD_DIR.exists():
        pytest.skip(
            "Recall download directory not found; set GVF_RECALL_DOWNLOAD_DIR "
            "to run recall tests."
        )
    pdfs = glob.glob(str(DOWNLOAD_DIR / "*.pdf")) + glob.glob(
        str(DOWNLOAD_DIR / "*.txt")
    )
    if not pdfs:
        pytest.skip(
            "No recall download files found; set GVF_RECALL_DOWNLOAD_DIR "
            "to a populated verified download directory."
        )

    pmid_files = defaultdict(list)
    for f in pdfs:
        basename = Path(f).name
        match = re.match(r"^(\d+)", basename)
        if match:
            pmid = int(match.group(1))
            is_supp = "supp" in basename.lower() or "_s" in basename.lower()
            pmid_files[pmid].append(
                {
                    "file": basename,
                    "is_supplement": is_supp,
                    "is_fulltext": not is_supp,  # Simplified assumption
                }
            )
    if not pmid_files:
        pytest.skip(
            "No PMID-prefixed recall download files found; set "
            "GVF_RECALL_DOWNLOAD_DIR to a populated verified download directory."
        )

    papers = []
    for pmid, files in pmid_files.items():
        papers.append(
            {
                "pmid": pmid,
                "files": files,
                "has_fulltext": any(not f["is_supplement"] for f in files),
                "has_supplements": any(f["is_supplement"] for f in files),
            }
        )

    return papers


@pytest.fixture(scope="session")
def extraction_results():
    """Load current GVF extraction results."""
    if not OUTPUT_DIR.exists():
        pytest.skip(
            "Recall output directory not found; set GVF_RECALL_OUTPUT_DIR "
            "to run recall tests."
        )
    # Auto-detect latest extraction file
    candidates = sorted(OUTPUT_DIR.glob("kcnh2_variant_scan_*.tsv"))
    if not candidates:
        pytest.skip("No extraction results found")

    latest = candidates[-1]
    df = pd.read_csv(latest, sep="\t")

    # Normalize variants
    if "variant_raw" in df.columns:
        df["variant_norm"] = df["variant_raw"].apply(normalize_variant)

    return df


@pytest.fixture(scope="session")
def extracted_variants(extraction_results):
    """Get set of extracted variants."""
    if "variant_norm" not in extraction_results.columns:
        return set()
    return set(extraction_results["variant_norm"].dropna())
