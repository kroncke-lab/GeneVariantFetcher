"""
Pydantic models for the Genetic Data Scout module.

These models define the structure for identifying and tracking high-value
data zones in biomedical literature for variant extraction.
"""

from typing import Literal, Optional

from pydantic import BaseModel, Field


class DataZone(BaseModel):
    """
    Represents a high-value data zone identified in a scientific paper.

    Data zones are sections of text (tables or prose) that contain
    patient-level or variant-specific information worth extracting.
    """

    id: str = Field(
        description="Unique identifier for the zone (e.g., 'supp_table_1', 'results_para_3')"
    )
    type: Literal["TABLE", "TEXT"] = Field(
        description="Zone structure type: TABLE for tabular data, TEXT for prose/narratives"
    )
    keep: bool = Field(
        description="True if zone contains individual patient/variant data, False if aggregate only"
    )
    content: str = Field(
        description="Brief description of zone contents (e.g., 'Table 2: List of 15 variants')"
    )
    text_snippet: str = Field(
        description="First ~200 characters of the zone for preview", max_length=500
    )
    char_start: int = Field(description="Start character position in source text", ge=0)
    char_end: int = Field(description="End character position in source text", ge=0)
    relevance_score: float = Field(
        description="Confidence score for zone relevance (0.0 to 1.0)", ge=0.0, le=1.0
    )

    # Optional metadata
    gene_mentions: int = Field(
        default=0, description="Count of target gene mentions in this zone"
    )
    variant_mentions: int = Field(
        default=0,
        description="Count of variant notation patterns (c., p.) in this zone",
    )
    clinical_keywords: int = Field(
        default=0,
        description="Count of clinical keywords (patient, proband, affected, etc.)",
    )
    source_section: Optional[str] = Field(
        default=None,
        description="Parent section name if identifiable (e.g., 'Results', 'Supplementary')",
    )

    @property
    def char_range(self) -> tuple:
        """Return character range as tuple for convenience."""
        return (self.char_start, self.char_end)

    @property
    def length(self) -> int:
        """Return length of the zone in characters."""
        return self.char_end - self.char_start


class DataZoneReport(BaseModel):
    """
    Complete report from the Genetic Data Scout analysis.
    """

    pmid: Optional[str] = Field(
        default=None, description="PubMed ID of the analyzed paper"
    )
    gene_symbol: str = Field(description="Target gene symbol used for analysis")
    total_zones_found: int = Field(description="Total number of data zones identified")
    zones_kept: int = Field(description="Number of zones marked keep=True")
    zones_discarded: int = Field(description="Number of zones marked keep=False")
    total_chars_original: int = Field(
        description="Character count of original full text"
    )
    total_chars_condensed: int = Field(
        description="Character count of condensed (kept) zones"
    )
    compression_ratio: float = Field(
        description="Ratio of condensed to original size", ge=0.0, le=1.0
    )
    zones: list[DataZone] = Field(
        default_factory=list, description="List of all identified data zones"
    )

    @property
    def kept_zones(self) -> list[DataZone]:
        """Return only zones marked for keeping."""
        return [z for z in self.zones if z.keep]

    @property
    def discarded_zones(self) -> list[DataZone]:
        """Return only zones marked for discarding."""
        return [z for z in self.zones if not z.keep]
