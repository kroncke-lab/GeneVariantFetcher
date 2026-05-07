"""Per-PMID audit log of harvested figures and supplements.

Saves one JSON file per PMID at ``{output_dir}/{PMID}_artifacts.json`` with
sufficient detail to answer questions like:

  * Did we download all the supplement files this paper has?
  * Did the supplement get converted to text, and how many characters?
  * Were figure captions extracted? How many?
  * Where did each artifact come from (PMC API, Elsevier API, scraper)?

This is additive — it does not replace ``manifest.json`` (which tracks
per-PMID stage status across the whole batch) or
``successful_downloads.csv`` / ``paywalled_missing.csv``. It complements
them with file-level granularity for the figure/supplement audit.
"""

from __future__ import annotations

import datetime
import json
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


@dataclass
class FigureArtifact:
    filename: str
    path: str
    size_bytes: int = 0
    source_url: Optional[str] = None
    source: str = ""  # "pmc_html", "pmc_xml", "publisher_html", "supplement_pdf"
    caption_label: Optional[str] = None
    caption_title: Optional[str] = None
    caption_text: Optional[str] = None
    figure_id: Optional[str] = None


@dataclass
class SupplementArtifact:
    filename: str
    path: str
    url: str = ""
    size_bytes: int = 0
    source: str = ""  # "pmc_europepmc", "elsevier_api", "scraper", etc.
    description: Optional[str] = None
    converted: bool = False
    converted_chars: int = 0
    figures_extracted: int = 0
    nested_files: List[str] = field(default_factory=list)
    error: Optional[str] = None


@dataclass
class MainTextArtifact:
    source: str = ""  # "pmc_xml", "pmc_html", "elsevier_api", "wiley_api", "scraper"
    chars: int = 0
    figure_captions_count: int = 0
    table_captions_count: int = 0
    supplement_descriptions_count: int = 0


@dataclass
class ArtifactsManifest:
    pmid: str
    pmcid: Optional[str] = None
    doi: Optional[str] = None
    gene_symbol: Optional[str] = None
    timestamp: str = ""
    main_text: MainTextArtifact = field(default_factory=MainTextArtifact)
    figures: List[FigureArtifact] = field(default_factory=list)
    supplements: List[SupplementArtifact] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pmid": self.pmid,
            "pmcid": self.pmcid,
            "doi": self.doi,
            "gene_symbol": self.gene_symbol,
            "timestamp": self.timestamp,
            "main_text": asdict(self.main_text),
            "figures": [asdict(f) for f in self.figures],
            "supplements": [asdict(s) for s in self.supplements],
            "notes": list(self.notes),
            "summary": {
                "figure_count": len(self.figures),
                "figures_with_captions": sum(
                    1 for f in self.figures if (f.caption_text or f.caption_title)
                ),
                "supplement_count": len(self.supplements),
                "supplements_converted": sum(
                    1 for s in self.supplements if s.converted
                ),
                "supplements_total_chars": sum(
                    s.converted_chars for s in self.supplements
                ),
                "main_text_chars": self.main_text.chars,
            },
        }


class ArtifactsLog:
    """Builder for an ArtifactsManifest tied to one PMID."""

    def __init__(
        self,
        pmid: str,
        output_dir: Path,
        pmcid: Optional[str] = None,
        doi: Optional[str] = None,
        gene_symbol: Optional[str] = None,
    ):
        self.output_dir = Path(output_dir)
        self.manifest = ArtifactsManifest(
            pmid=pmid,
            pmcid=pmcid,
            doi=doi,
            gene_symbol=gene_symbol,
            timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        )

    # ------------------------------------------------------------------
    # Recording
    # ------------------------------------------------------------------

    def record_main_text(
        self,
        source: str,
        chars: int,
        figure_captions: int = 0,
        table_captions: int = 0,
        supplement_descriptions: int = 0,
    ) -> None:
        self.manifest.main_text = MainTextArtifact(
            source=source,
            chars=chars,
            figure_captions_count=figure_captions,
            table_captions_count=table_captions,
            supplement_descriptions_count=supplement_descriptions,
        )

    def record_figure(self, artifact: FigureArtifact) -> None:
        self.manifest.figures.append(artifact)

    def record_figure_dict(self, **kwargs: Any) -> None:
        self.manifest.figures.append(FigureArtifact(**kwargs))

    def record_supplement(self, artifact: SupplementArtifact) -> None:
        self.manifest.supplements.append(artifact)

    def record_supplement_dict(self, **kwargs: Any) -> None:
        self.manifest.supplements.append(SupplementArtifact(**kwargs))

    def add_note(self, note: str) -> None:
        self.manifest.notes.append(note)

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def path(self) -> Path:
        return self.output_dir / f"{self.manifest.pmid}_artifacts.json"

    def save(self) -> Path:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        out_path = self.path()
        try:
            out_path.write_text(
                json.dumps(self.manifest.to_dict(), indent=2),
                encoding="utf-8",
            )
            logger.info(
                "Wrote artifacts log for PMID %s: %d figures, %d supplements",
                self.manifest.pmid,
                len(self.manifest.figures),
                len(self.manifest.supplements),
            )
        except Exception as exc:
            logger.warning(
                "Failed to write artifacts log for PMID %s: %s",
                self.manifest.pmid,
                exc,
            )
        return out_path

    # ------------------------------------------------------------------
    # Convenience factories
    # ------------------------------------------------------------------

    @classmethod
    def figure_artifact_from_path(
        cls,
        filepath: Path,
        source: str,
        source_url: Optional[str] = None,
        caption_label: Optional[str] = None,
        caption_title: Optional[str] = None,
        caption_text: Optional[str] = None,
        figure_id: Optional[str] = None,
    ) -> FigureArtifact:
        try:
            size = filepath.stat().st_size if filepath.exists() else 0
        except OSError:
            size = 0
        return FigureArtifact(
            filename=filepath.name,
            path=str(filepath),
            size_bytes=size,
            source_url=source_url,
            source=source,
            caption_label=caption_label,
            caption_title=caption_title,
            caption_text=caption_text,
            figure_id=figure_id,
        )

    @classmethod
    def supplement_artifact_from_path(
        cls,
        filepath: Path,
        url: str,
        source: str,
        description: Optional[str] = None,
        converted_chars: int = 0,
        figures_extracted: int = 0,
        nested_files: Optional[List[str]] = None,
        error: Optional[str] = None,
    ) -> SupplementArtifact:
        try:
            size = filepath.stat().st_size if filepath.exists() else 0
        except OSError:
            size = 0
        return SupplementArtifact(
            filename=filepath.name,
            path=str(filepath),
            url=url,
            size_bytes=size,
            source=source,
            description=description,
            converted=converted_chars > 0,
            converted_chars=converted_chars,
            figures_extracted=figures_extracted,
            nested_files=list(nested_files or []),
            error=error,
        )
