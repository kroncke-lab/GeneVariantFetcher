"""Integrity checks for the committed initial paper-review cohort."""

from pathlib import Path

from utils.pmid_utils import is_valid_pmid


COHORT_DIR = (
    Path(__file__).resolve().parents[2]
    / "benchmarks"
    / "curated_extraction_eval"
    / "review_pmids_50"
)
ACTIVE_BRCA2_DIR = (
    Path(__file__).resolve().parents[2]
    / "benchmarks"
    / "curated_extraction_eval"
    / "review_pmids_20260811_brca2_provenance"
)
PRESENTATION_COHORT_DIR = (
    Path(__file__).resolve().parents[2]
    / "benchmarks"
    / "curated_extraction_eval"
    / "review_pmids_50_20260824_curated"
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


def test_active_brca2_review_scope_excludes_provenance_and_canine_papers():
    historical = set((COHORT_DIR / "BRCA2.txt").read_text().splitlines())
    active = [
        line.strip()
        for line in (ACTIVE_BRCA2_DIR / "BRCA2.txt").read_text().splitlines()
        if line.strip()
    ]
    excluded = {
        "10398279",
        "15365993",
        "18489799",
        "19944633",
        "21356067",
        "22655046",
        "25802882",
    }

    assert len(active) == 45
    assert len(set(active)) == 45
    assert all(is_valid_pmid(pmid) for pmid in active)
    assert set(active) == historical - excluded
    assert {"26833046", "26848529"} <= set(active)


def test_presentation_cohort_has_three_source_screened_50_paper_queues():
    cohort_files = sorted(PRESENTATION_COHORT_DIR.glob("*.txt"))
    assert {path.stem for path in cohort_files} == {"BRCA1", "BRCA2", "BMPR2"}

    for path in cohort_files:
        pmids = [line.strip() for line in path.read_text().splitlines() if line.strip()]
        assert len(pmids) == len(set(pmids)) == 50, path
        assert all(is_valid_pmid(pmid) for pmid in pmids), path

    excluded = {
        "16014699",
        "18779604",
        "19476645",
        "20807450",
        "22144684",
        "23519070",
        "25136594",
        "25543221",
        "25859162",
        "26109977",
        "27077124",
        "27767231",
        "29580235",
        "29936257",
        "30283497",
        "31757951",
        "32019284",
        "32068069",
        "32566972",
        "32766142",
        "34678164",
        "35229499",
        "37278169",
        "39367156",
        "39848555",
        "19941162",
        "14559878",
        "15059534",
        "23139147",
    }
    active = {
        pmid
        for path in cohort_files
        for pmid in path.read_text().splitlines()
        if pmid.strip()
    }
    assert excluded.isdisjoint(active)
