"""Shared fixtures for GVF recall tests."""
import pytest
import pandas as pd
import re
import glob
from pathlib import Path
from collections import defaultdict

GOLD_STANDARD_PATH = Path(
    "/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher/comparison_results/"
    "KCNH2 HeterozygoteDatabase-4-clinical_w_paris_japan_mayo_and_italycohort.xls"
)

DOWNLOAD_DIR = Path("/mnt/temp2/kronckbm/gvf_output/verified_downloads_20260208")
OUTPUT_DIR = Path("/mnt/temp2/kronckbm/gvf_output")


def normalize_variant(v):
    """Normalize variant nomenclature for comparison."""
    if pd.isna(v):
        return None
    v = str(v).upper().strip()
    v = v.replace('*', 'X')
    v = v.replace('P.', '')
    # Remove p. prefix if present
    if v.startswith('P.'):
        v = v[2:]
    # Handle frameshift variants
    if 'FS' in v:
        v = re.sub(r'FS.*$', 'FSX', v)
    return v


@pytest.fixture(scope="session")
def gold_standard():
    """Load and normalize gold standard data."""
    df = pd.read_excel(GOLD_STANDARD_PATH)
    
    # Normalize PMIDs
    df['PMID_int'] = df['PMID'].apply(
        lambda x: int(float(x)) if pd.notna(x) else None
    )
    
    # Normalize variants
    df['Variant_norm'] = df['Variant'].apply(normalize_variant)
    
    # Flag LQT-affected (LQT >= 1 means affected)
    df['is_LQT_affected'] = df['LQT'].apply(lambda x: x >= 1 if pd.notna(x) else False)
    
    pmids = set(df['PMID_int'].dropna().astype(int))
    variants = set(df['Variant_norm'].dropna())
    
    return {
        'df': df,
        'pmids': pmids,
        'variants': variants,
        'total_carriers': len(df),
        'carriers_with_pmid': len(df[df['PMID_int'].notna()]),
        'total_lqt_affected': int(df['is_LQT_affected'].sum()),
    }


@pytest.fixture(scope="session")
def downloaded_pmids():
    """Get set of downloaded PMIDs from PDF filenames."""
    pdfs = glob.glob(str(DOWNLOAD_DIR / "*.pdf"))
    pmids = set()
    
    for pdf in pdfs:
        basename = Path(pdf).name
        # Extract PMID from filename (assumes PMID_*.pdf or PMID.pdf format)
        match = re.match(r'^(\d+)', basename)
        if match:
            pmids.add(int(match.group(1)))
    
    return pmids


@pytest.fixture(scope="session")
def downloaded_papers():
    """Get detailed info about downloaded papers."""
    pdfs = glob.glob(str(DOWNLOAD_DIR / "*.pdf")) + glob.glob(str(DOWNLOAD_DIR / "*.txt"))
    
    pmid_files = defaultdict(list)
    for f in pdfs:
        basename = Path(f).name
        match = re.match(r'^(\d+)', basename)
        if match:
            pmid = int(match.group(1))
            is_supp = 'supp' in basename.lower() or '_s' in basename.lower()
            pmid_files[pmid].append({
                'file': basename,
                'is_supplement': is_supp,
                'is_fulltext': not is_supp  # Simplified assumption
            })
    
    papers = []
    for pmid, files in pmid_files.items():
        papers.append({
            'pmid': pmid,
            'files': files,
            'has_fulltext': any(not f['is_supplement'] for f in files),
            'has_supplements': any(f['is_supplement'] for f in files),
        })
    
    return papers


@pytest.fixture(scope="session")
def extraction_results():
    """Load current GVF extraction results."""
    # Auto-detect latest extraction file
    candidates = sorted(OUTPUT_DIR.glob("kcnh2_variant_scan_*.tsv"))
    if not candidates:
        pytest.skip("No extraction results found")
    
    latest = candidates[-1]
    df = pd.read_csv(latest, sep='\t')
    
    # Normalize variants
    if 'variant_raw' in df.columns:
        df['variant_norm'] = df['variant_raw'].apply(normalize_variant)
    
    return df


@pytest.fixture(scope="session")
def extracted_variants(extraction_results):
    """Get set of extracted variants."""
    if 'variant_norm' not in extraction_results.columns:
        return set()
    return set(extraction_results['variant_norm'].dropna())
