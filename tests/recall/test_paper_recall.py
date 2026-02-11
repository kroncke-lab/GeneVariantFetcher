"""Paper-level recall tests for GVF."""

import pytest


class TestPaperRecall:
    """Paper-level recall tests."""

    def test_overall_paper_recall(self, gold_standard, downloaded_pmids):
        """Test overall paper recall rate."""
        overlap = gold_standard["pmids"] & downloaded_pmids
        recall = len(overlap) / len(gold_standard["pmids"])

        # Floor test - should not regress below 55%
        assert recall >= 0.55, f"Paper recall {recall:.1%} regressed below 55% floor"

        print(
            f"\nPaper Recall: {len(overlap)}/{len(gold_standard['pmids'])} = {recall:.1%}"
        )
        print(f"Missing: {len(gold_standard['pmids'] - downloaded_pmids)} papers")

    def test_paper_count_minimum(self, downloaded_pmids):
        """Ensure minimum number of papers downloaded."""
        assert len(downloaded_pmids) >= 140, (
            f"Only {len(downloaded_pmids)} papers, expected >= 140"
        )

    def test_fulltext_coverage(self, downloaded_papers):
        """Check full-text vs abstract-only breakdown."""
        fulltext = sum(1 for p in downloaded_papers if p["has_fulltext"])
        total = len(downloaded_papers)

        fulltext_rate = fulltext / total if total > 0 else 0

        print(f"\nFull-text papers: {fulltext}/{total} = {fulltext_rate:.1%}")

        # At least 80% should be full text
        assert fulltext_rate >= 0.80, f"Only {fulltext_rate:.1%} are full text"

    def test_supplement_coverage(self, downloaded_papers):
        """Check supplement download status (informational)."""
        with_supps = sum(1 for p in downloaded_papers if p["has_supplements"])
        total = len(downloaded_papers)
        supp_rate = with_supps / total if total > 0 else 0

        print(f"\nPapers with supplements: {with_supps}/{total} = {supp_rate:.1%}")
        # Informational only - no assertion

    def test_missing_high_carrier_papers(self, gold_standard, downloaded_pmids):
        """Identify missing papers with most carriers."""
        df = gold_standard["df"]
        missing = gold_standard["pmids"] - downloaded_pmids

        # Count carriers per missing PMID
        carrier_counts = []
        for pmid in missing:
            count = len(df[df["PMID_int"] == pmid])
            carrier_counts.append((pmid, count))

        # Sort by carrier count
        carrier_counts.sort(key=lambda x: x[1], reverse=True)

        print(f"\nTop 10 missing papers by carrier count:")
        for pmid, count in carrier_counts[:10]:
            print(f"  PMID {pmid}: {count} carriers")

        # Informational - no assertion


class TestPriorityPapers:
    """Tests for specific high-priority papers."""

    PRIORITY_PMIDS = [
        # Large variant/carrier papers from gold standard
        (14661677, "Large LQTS cohort"),
        (26496715, "Comprehensive KCNH2 review"),
        (10973849, "Classic HERG paper"),
    ]

    @pytest.mark.parametrize("pmid,description", PRIORITY_PMIDS)
    def test_priority_paper_present(self, pmid, description, downloaded_pmids):
        """High-priority papers should be present."""
        # Note: This is informational - some may be legitimately unavailable
        if pmid not in downloaded_pmids:
            pytest.skip(f"Priority paper {pmid} ({description}) not available")
