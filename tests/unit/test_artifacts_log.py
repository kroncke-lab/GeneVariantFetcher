"""Tests for harvesting.artifacts_log."""

from __future__ import annotations

import json
from pathlib import Path

from harvesting.artifacts_log import (
    ArtifactsLog,
    FigureArtifact,
    SupplementArtifact,
)


def test_artifacts_log_writes_summary(tmp_path: Path):
    log = ArtifactsLog(
        pmid="12345",
        output_dir=tmp_path,
        pmcid="PMC999",
        doi="10.1/test",
        gene_symbol="KCNH2",
    )
    log.record_main_text(
        source="pmc_xml",
        chars=1234,
        figure_captions=2,
        table_captions=1,
        supplement_descriptions=3,
    )
    fig = tmp_path / "fig1.png"
    fig.write_bytes(b"\x89PNGfake")
    log.record_figure(
        ArtifactsLog.figure_artifact_from_path(
            fig,
            source="pmc_html",
            caption_label="Figure 1",
            caption_text="A pedigree carrying KCNH2 G604S",
        )
    )
    log.record_supplement_dict(
        filename="mmc1.xlsx",
        path=str(tmp_path / "mmc1.xlsx"),
        url="https://example.com/mmc1.xlsx",
        size_bytes=2048,
        source="elsevier_api",
        description="Supplementary Table S1: 86 missense variants",
        converted=True,
        converted_chars=8000,
        figures_extracted=0,
    )
    log.add_note("first run")

    out = log.save()
    assert out.exists()
    data = json.loads(out.read_text())
    assert data["pmid"] == "12345"
    assert data["pmcid"] == "PMC999"
    assert data["main_text"]["chars"] == 1234
    assert data["main_text"]["figure_captions_count"] == 2
    assert len(data["figures"]) == 1
    assert data["figures"][0]["caption_label"] == "Figure 1"
    assert data["figures"][0]["size_bytes"] == len(b"\x89PNGfake")
    assert len(data["supplements"]) == 1
    assert data["supplements"][0]["source"] == "elsevier_api"
    assert data["supplements"][0]["converted_chars"] == 8000
    assert data["summary"]["figure_count"] == 1
    assert data["summary"]["figures_with_captions"] == 1
    assert data["summary"]["supplement_count"] == 1
    assert data["summary"]["supplements_converted"] == 1
    assert data["summary"]["supplements_total_chars"] == 8000
    assert data["notes"] == ["first run"]


def test_supplement_artifact_from_missing_path_records_zero_size(tmp_path: Path):
    missing = tmp_path / "does-not-exist.zip"
    artifact = ArtifactsLog.supplement_artifact_from_path(
        filepath=missing,
        url="https://example.com/x.zip",
        source="scraper",
    )
    assert isinstance(artifact, SupplementArtifact)
    assert artifact.size_bytes == 0
    assert not artifact.converted


def test_figure_artifact_records_size(tmp_path: Path):
    p = tmp_path / "fig.png"
    p.write_bytes(b"x" * 17)
    fig = ArtifactsLog.figure_artifact_from_path(p, source="pmc_html")
    assert isinstance(fig, FigureArtifact)
    assert fig.size_bytes == 17
    assert fig.filename == "fig.png"


# ---------------------------------------------------------------------------
# supplement_links — the inventory of what the paper advertised
# ---------------------------------------------------------------------------


class _Supp:
    """Minimal stand-in for figure_extractor.SupplementDescription."""

    def __init__(self, label, href, title=None, media_type=None):
        self.label = label
        self.href = href
        self.title = title
        self.media_type = media_type


class _Captions:
    def __init__(self, supplements):
        self.supplements = supplements


def test_record_supplement_links_marks_unfetched_links(tmp_path: Path):
    log = ArtifactsLog(pmid="26484152", output_dir=tmp_path, pmcid="PMC4535901")
    log.record_supplement_dict(
        filename="mmc1.xls",
        path=str(tmp_path / "mmc1.xls"),
        url="https://pmc.ncbi.nlm.nih.gov/articles/PMC4535901/bin/mmc1.xls",
        converted=True,
        converted_chars=500,
    )
    log.record_supplement_links(
        _Captions([_Supp("Table S1", "mmc1.xls"), _Supp("Table S2", "mmc2.pdf")]),
        source="pmc",
    )

    data = json.loads(log.save().read_text())
    links = {link["label"]: link for link in data["supplement_links"]}
    assert links["Table S1"]["downloaded"] is True
    assert links["Table S2"]["downloaded"] is False
    # The recovery signal a later pass targets.
    assert data["summary"]["supplement_link_count"] == 2
    assert data["summary"]["supplement_links_unfetched"] == 1
    # Link rows must not inflate the processed-file count.
    assert data["summary"]["supplement_count"] == 1


def test_record_supplement_links_ignores_descriptions_without_an_href(tmp_path: Path):
    log = ArtifactsLog(pmid="1", output_dir=tmp_path)
    log.record_supplement_links(
        _Captions([_Supp("Supplementary Methods", None), _Supp("Table S1", "  ")]),
        source="publisher_html",
    )
    data = json.loads(log.save().read_text())
    assert data["supplement_links"] == []
    assert data["summary"]["supplement_links_unfetched"] == 0


def test_record_supplement_links_dedupes_repeated_hrefs(tmp_path: Path):
    log = ArtifactsLog(pmid="1", output_dir=tmp_path)
    caps = _Captions([_Supp("A", "mmc1.xls"), _Supp("B", "mmc1.xls")])
    log.record_supplement_links(caps, source="pmc")
    log.record_supplement_links(caps, source="pmc")
    data = json.loads(log.save().read_text())
    assert len(data["supplement_links"]) == 1


def test_older_manifests_still_report_zero_links(tmp_path: Path):
    log = ArtifactsLog(pmid="1", output_dir=tmp_path)
    data = json.loads(log.save().read_text())
    assert data["supplement_links"] == []
    assert data["summary"]["supplement_link_count"] == 0


def test_a_link_delivered_inside_an_archive_counts_as_downloaded(tmp_path: Path):
    """Europe PMC serves a whole article's supplements as one ZIP.

    The advertised per-file links then exist only inside the extracted
    directory. Matching on the archive's own name alone marked every one of
    them unfetched and would send a recovery pass after files already held.
    """
    log = ArtifactsLog(pmid="30309222", output_dir=tmp_path, pmcid="PMC6639209")
    extracted = tmp_path / "PMC6639209_supplements"
    log.record_supplement_dict(
        filename="PMC6639209_supplements.zip",
        path=str(tmp_path / "PMC6639209_supplements.zip"),
        url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6639209/supplementaryFiles",
        converted=True,
        converted_chars=4200,
        nested_files=[
            str(extracted / "crt-2018-312-suppl1.pdf"),
            str(extracted / "crt-2018-312-suppl2.pdf"),
        ],
    )
    log.record_supplement_links(
        _Captions(
            [
                _Supp(
                    "S1",
                    "https://pmc.ncbi.nlm.nih.gov/articles/PMC6639209/bin/crt-2018-312-suppl1.pdf",
                ),
                _Supp("S2", "crt-2018-312-suppl2.pdf"),
                _Supp("S3", "crt-2018-312-suppl3.pdf"),  # genuinely absent
            ]
        ),
        source="pmc",
    )

    data = json.loads(log.save().read_text())
    links = {link["label"]: link["downloaded"] for link in data["supplement_links"]}
    assert links == {"S1": True, "S2": True, "S3": False}
    # Only the genuinely absent one should drive a recovery pass.
    assert data["summary"]["supplement_links_unfetched"] == 1
