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
from urllib.parse import urlparse

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
class SupplementLinkArtifact:
    """A supplement link the paper advertised, recorded whether or not we got it.

    Distinct from :class:`SupplementArtifact`, which describes a file we
    actually processed. These rows are the *inventory*: what the markup said
    existed. Kept structured (rather than only rendered into ``FULL_CONTEXT.md``
    as a ``_link_:`` line) so a recovery pass can resolve and fetch them
    afterwards — see ``harvesting.supplement_link_resolver``.
    """

    label: str
    href: str
    title: Optional[str] = None
    media_type: Optional[str] = None
    source: str = ""  # "pmc_xml", "pmc_html", "publisher_html", "browser_html"
    downloaded: bool = False


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
    supplement_links: List[SupplementLinkArtifact] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        downloaded_names = {
            Path(s.path).name.lower() for s in self.supplements if s.converted or s.path
        }
        links_unfetched = sum(
            1
            for link in self.supplement_links
            if not link.downloaded
            and Path(urlparse(link.href).path or link.href).name.lower()
            not in downloaded_names
        )
        return {
            "pmid": self.pmid,
            "pmcid": self.pmcid,
            "doi": self.doi,
            "gene_symbol": self.gene_symbol,
            "timestamp": self.timestamp,
            "main_text": asdict(self.main_text),
            "figures": [asdict(f) for f in self.figures],
            "supplements": [asdict(s) for s in self.supplements],
            "supplement_links": [asdict(link) for link in self.supplement_links],
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
                # The recovery signal: the paper advertised links we never
                # fetched. Non-zero here is what a source-recovery pass targets.
                "supplement_link_count": len(self.supplement_links),
                "supplement_links_unfetched": links_unfetched,
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

    def record_supplement_links(
        self,
        captions: Any,
        source: str,
    ) -> None:
        """Record the supplement links a paper advertised, from a caption result.

        Takes a ``figure_extractor.CaptionExtractionResult``. Only entries with
        an href are recorded — a description with no link is already covered by
        ``main_text.supplement_descriptions_count``. Marks each link
        ``downloaded`` when a supplement of the same filename was processed, so
        a recovery pass can target the genuine gaps.
        """
        fetched = {
            Path(urlparse(s.url or "").path or s.url or "").name.lower()
            for s in self.manifest.supplements
            if s.url
        } | {Path(s.path).name.lower() for s in self.manifest.supplements if s.path}
        fetched.discard("")

        seen = {link.href for link in self.manifest.supplement_links}
        for supp in getattr(captions, "supplements", []) or []:
            href = (getattr(supp, "href", None) or "").strip()
            if not href or href in seen:
                continue
            seen.add(href)
            name = Path(urlparse(href).path or href).name.lower()
            self.manifest.supplement_links.append(
                SupplementLinkArtifact(
                    label=getattr(supp, "label", "") or "",
                    href=href,
                    title=getattr(supp, "title", None) or None,
                    media_type=getattr(supp, "media_type", None),
                    source=source,
                    downloaded=bool(name) and name in fetched,
                )
            )

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
