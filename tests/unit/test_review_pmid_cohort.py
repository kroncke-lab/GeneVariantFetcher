"""Integrity checks for the committed initial paper-review cohort."""

from pathlib import Path

from utils.pmid_utils import is_valid_pmid


COHORT_DIR = (
    Path(__file__).resolve().parents[2]
    / "benchmarks"
    / "curated_extraction_eval"
    / "review_pmids_50"
)
EXPECTED_GENES = {
    "APOE",
    "BRCA1",
    "BRCA2",
    "KCNH2",
    "KCNQ1",
    "MYBPC3",
    "RYR2",
    "SCN5A",
}


def test_review_cohort_has_50_unique_valid_pmids_per_gene():
    cohort_files = sorted(COHORT_DIR.glob("*.txt"))

    assert {path.stem for path in cohort_files} == EXPECTED_GENES
    for path in cohort_files:
        pmids = [line.strip() for line in path.read_text().splitlines() if line.strip()]
        assert len(pmids) == 50, path
        assert len(set(pmids)) == 50, path
        assert all(is_valid_pmid(pmid) for pmid in pmids), path

    kcnh2_pmids = (COHORT_DIR / "KCNH2.txt").read_text().splitlines()
    assert "34546463" in kcnh2_pmids
    assert "PMC9522753" not in kcnh2_pmids
