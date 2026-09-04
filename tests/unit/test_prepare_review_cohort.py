from scripts.prepare_review_cohort import PaperCandidate, select_cohort


def _paper(
    index: int,
    *,
    complete: bool,
    source: str,
    decade: int,
    adjudicated: bool = False,
) -> PaperCandidate:
    year = decade + index % 10
    return PaperCandidate(
        pmid=str(10000000 + index),
        title=f"Paper {index}",
        journal="Journal",
        publication_date=str(year),
        year=year,
        decade=f"{decade}s",
        doi="",
        pmc_id="",
        source_origin=source,
        has_full_text=source != "abstract-only",
        variant_links=index + 1,
        unique_variants=index + 1,
        count_rows=1,
        carrier_rows=1,
        affected_rows=1,
        unaffected_rows=1 if complete else 0,
        complete_rows=1 if complete else 0,
        trusted_complete_rows=1 if complete else 0,
        total_carriers=index + 1,
        total_affected=index,
        total_unaffected=1 if complete else 0,
        extraction_selected=True,
        secondary_llm_adjudicated=adjudicated,
        claim_verifier_calls=1,
    )


def test_select_cohort_enforces_completeness_and_diversity():
    papers = [
        _paper(1, complete=True, source="pmc_xml", decade=2000),
        _paper(2, complete=True, source="pmc_xml", decade=2010),
        _paper(3, complete=True, source="elsevier-api", decade=2020),
        _paper(4, complete=True, source="abstract-only", decade=2010),
        _paper(5, complete=True, source="pmc_xml", decade=2020),
        _paper(6, complete=False, source="publisher-free", decade=2000),
        _paper(7, complete=False, source="wiley-api", decade=2010),
        _paper(
            8,
            complete=False,
            source="browser-html-generic",
            decade=2020,
            adjudicated=True,
        ),
    ]

    cohort = select_cohort(
        papers,
        size=7,
        complete_count=4,
        min_source_origins=4,
        min_decades=3,
    )

    assert len(cohort) == len({paper.pmid for paper in cohort}) == 7
    assert sum(paper.complete_for_quota for paper in cohort) >= 4
    assert len({paper.source_origin for paper in cohort}) >= 4
    assert {paper.decade for paper in cohort} == {"2000s", "2010s", "2020s"}
    assert any(paper.secondary_llm_adjudicated for paper in cohort)


def test_select_cohort_fails_when_complete_pool_is_too_small():
    papers = [
        _paper(1, complete=True, source="pmc_xml", decade=2000),
        _paper(2, complete=False, source="elsevier-api", decade=2010),
    ]

    try:
        select_cohort(
            papers,
            size=2,
            complete_count=2,
            min_source_origins=1,
            min_decades=1,
        )
    except ValueError as exc:
        assert "need 2 trusted complete papers, found 1" in str(exc)
    else:
        raise AssertionError(
            "selection should fail when the complete pool is undersized"
        )
