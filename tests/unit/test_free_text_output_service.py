"""Tests for shared free-text output writer helpers."""

import csv

from harvesting.free_text_output_service import (
    publisher_api_fallback_source,
    source_from_free_text_flags,
    write_free_text_output,
)


def test_source_from_free_text_flags():
    elsevier_source = source_from_free_text_flags(
        used_elsevier_api=True, used_wiley_api=False
    )
    assert elsevier_source.success_marker == "ELSEVIER_API"
    assert elsevier_source.status_source == "elsevier-api"
    assert elsevier_source.source_tag == "[via Elsevier API]"

    wiley_source = source_from_free_text_flags(
        used_elsevier_api=False, used_wiley_api=True
    )
    assert wiley_source.success_marker == "WILEY_API"
    assert wiley_source.status_source == "wiley-api"
    assert wiley_source.source_tag == "[via Wiley API]"

    publisher_source = source_from_free_text_flags(
        used_elsevier_api=False, used_wiley_api=False
    )
    assert publisher_source.success_marker == "PUBLISHER_FREE"
    assert publisher_source.status_source == "publisher-free"
    assert publisher_source.source_tag == "[from publisher]"


def test_write_free_text_output_writes_file_log_and_status(tmp_path):
    success_log = tmp_path / "success.csv"
    success_log.write_text("PMID,PMCID,Supplements_Downloaded\n", encoding="utf-8")

    # Content must be long enough (>1500 chars) and have paper indicators to pass validation
    main_markdown = """# MAIN TEXT

## Abstract

This is a test paper about KCNH2 variants and their role in Long QT syndrome. We identified 
the R534C mutation in 12 patients from our multicenter cohort study. The goal of this 
research was to characterize the clinical phenotype of carriers and determine penetrance.
Our findings have important implications for genetic counseling and risk stratification.

## Introduction

Long QT syndrome (LQTS) is a cardiac channelopathy characterized by prolonged ventricular 
repolarization. The HERG gene (KCNH2) encodes the alpha subunit of the rapid delayed 
rectifier potassium channel (IKr). Mutations in KCNH2 cause LQT2, which accounts for 
approximately 30-40% of all LQTS cases. Understanding genotype-phenotype correlations 
is essential for clinical management of affected individuals and their families.

## Methods

Genetic testing was performed using next-generation sequencing panels. Variants were 
classified according to ACMG guidelines. Clinical data including ECG parameters, symptoms, 
and family history were collected from all participants. Statistical analysis was 
performed using R software version 4.0.

## Results

Among 50 carriers of pathogenic variants, 30 (60%) experienced breakthrough cardiac events 
before age 40. The R534C variant was associated with a hazard ratio of 1.32 for cardiac 
events compared to other KCNH2 variants. QTc prolongation was observed in 85% of carriers.
Penetrance was estimated at 60% by age 40 and 75% by age 60.

## Discussion

Our findings support the pathogenicity of these KCNH2 variants and highlight the importance 
of early intervention in affected individuals. The high penetrance observed underscores 
the need for cascade genetic testing in families. Future studies should focus on 
identifying modifiers of disease expression.

## References

1. Test reference one
2. Test reference two
"""
    supplement_markdown = "\n\n## Supplement\n\nAdditional data here."
    expected_content = main_markdown + supplement_markdown

    status_calls = []
    source = publisher_api_fallback_source()
    output_file, unified_content = write_free_text_output(
        output_dir=tmp_path,
        success_log=success_log,
        pmid="123456",
        main_markdown=main_markdown,
        supplement_markdown=supplement_markdown,
        downloaded_count=2,
        source=source,
        write_pmid_status=lambda pmid, status, details: status_calls.append(
            (pmid, status, details)
        ),
    )

    assert output_file is not None, "output_file should not be None for valid content"
    assert output_file.exists()
    assert output_file.read_text(encoding="utf-8") == expected_content
    assert unified_content == expected_content

    # Also add a test for binary/invalid content rejection
    binary_content = "# MAIN TEXT\n\n\x00\x01\x02 binary garbage"
    bad_output, error_msg = write_free_text_output(
        output_dir=tmp_path,
        success_log=success_log,
        pmid="999999",
        main_markdown=binary_content,
        supplement_markdown="",
        downloaded_count=0,
        source=source,
    )
    assert bad_output is None, "Should return None for invalid content"
    assert "validation failed" in error_msg.lower()

    with open(success_log, "r", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    assert rows[-1] == ["123456", "publisher-api", "2"]

    assert len(status_calls) == 1
    pmid, status, details = status_calls[0]
    assert pmid == "123456"
    assert status == "extracted"
    assert details["variant_count"] == 0
    assert details["source"] == "publisher-api-fallback"
