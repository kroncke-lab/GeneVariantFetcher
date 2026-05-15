"""Tests for FULL_CONTEXT.md assembly from body, captions, and supplements."""

from harvesting.figure_extractor import (
    CaptionExtractionResult,
    FigureCaption,
    SupplementDescription,
    TableCaption,
)
from harvesting.orchestrator import PMCHarvester


def test_unified_content_preserves_body_caption_supplement_order():
    harvester = PMCHarvester.__new__(PMCHarvester)
    captions = CaptionExtractionResult(
        figures=[
            FigureCaption(
                label="Figure 1",
                title="Pedigree of family carrying KCNH2 G604S",
                text="Filled symbols are affected; open symbols are unaffected.",
                image_url="https://example.org/fig1.png",
                figure_id="fig1",
            )
        ],
        tables=[
            TableCaption(
                label="Table 1",
                title="Variant carrier counts",
                text="Counts of affected and unaffected carriers.",
                table_id="t1",
            )
        ],
        supplements=[
            SupplementDescription(
                label="Supplementary Table S1",
                title="Extended KCNH2 variant list",
                text="Supplement contains p.Ala561Val and 86 other variants.",
                href="supp/S1.xlsx",
            )
        ],
    )

    unified = harvester._build_unified_content(
        "# MAIN TEXT\n\nBody with p.Arg176Trp.\n",
        captions,
        "\n\n# SUPPLEMENTAL FILE 1: S1.xlsx\n\np.Ala561Val appears here.\n",
    )

    assert unified.index("# MAIN TEXT") < unified.index("## FIGURE CAPTIONS")
    assert unified.index("## FIGURE CAPTIONS") < unified.index("## TABLE CAPTIONS")
    assert unified.index("## TABLE CAPTIONS") < unified.index(
        "## SUPPLEMENT DESCRIPTIONS"
    )
    assert unified.index("## SUPPLEMENT DESCRIPTIONS") < unified.index(
        "# SUPPLEMENTAL FILE 1"
    )
    assert "p.Arg176Trp" in unified
    assert "KCNH2 G604S" in unified
    assert "affected and unaffected carriers" in unified
    assert "p.Ala561Val" in unified
