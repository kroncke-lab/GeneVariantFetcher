"""Figure caption and supplement-description extraction.

Variant names frequently appear in figure legends and supplement descriptions
("Pedigree of family carrying KCNH2 G604S", "Table S1: 86 KCNH2 missense
variants identified in screening cohort"). The standard XML/HTML→markdown
converters in `format_converters.py` walk `<sec>`/`<p>` elements only and
silently drop these regions. This module extracts them so they get included
in `FULL_CONTEXT.md` for both the regex pre-scanner and the LLM extractor.

Sources handled:
  - PMC JATS XML: <fig>, <table-wrap>, <supplementary-material>
  - PMC / publisher HTML: <figure>, <figcaption>, supplementary-material divs

The output is markdown with a stable `## FIGURE CAPTIONS` /
`## SUPPLEMENT DESCRIPTIONS` heading so downstream consumers can locate it.
"""

from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import List, Optional

from bs4 import BeautifulSoup, Tag


@dataclass
class FigureCaption:
    """Caption metadata for a single figure."""

    label: str  # "Figure 1", "Fig. 2A", etc. (best-effort)
    title: str  # First sentence / bold heading of the caption
    text: str  # Full caption body
    image_url: Optional[str] = None  # Resolved if present in source
    figure_id: Optional[str] = None  # XML/HTML id attribute, if any


@dataclass
class SupplementDescription:
    """Description text associated with a supplementary material reference."""

    label: str  # e.g. "Supplementary Table S1"
    title: str  # First sentence of the description
    text: str  # Full description body
    href: Optional[str] = None  # Source URL or xlink:href
    media_type: Optional[str] = None  # e.g. "application/vnd.ms-excel"


@dataclass
class TableCaption:
    """Caption metadata for an inline (non-supplementary) table."""

    label: str
    title: str
    text: str
    table_id: Optional[str] = None
    image_url: Optional[str] = None


@dataclass
class CaptionExtractionResult:
    """All caption-like content extracted from a single source document."""

    figures: List[FigureCaption] = field(default_factory=list)
    tables: List[TableCaption] = field(default_factory=list)
    supplements: List[SupplementDescription] = field(default_factory=list)

    def is_empty(self) -> bool:
        return not (self.figures or self.tables or self.supplements)

    def to_dict(self) -> dict:
        return {
            "figures": [asdict(f) for f in self.figures],
            "tables": [asdict(t) for t in self.tables],
            "supplements": [asdict(s) for s in self.supplements],
        }


# ---------------------------------------------------------------------------
# JATS XML extraction (PMC)
# ---------------------------------------------------------------------------

# JATS uses optional namespaces; strip them when matching tag names.
_NS_RE = re.compile(r"^\{[^}]+\}")


def _local(tag: str) -> str:
    return _NS_RE.sub("", tag)


def _itertext_clean(elem: Optional[ET.Element]) -> str:
    if elem is None:
        return ""
    text = "".join(elem.itertext())
    text = re.sub(r"\s+", " ", text).strip()
    return text


def _first_sentence(text: str, max_chars: int = 200) -> str:
    if not text:
        return ""
    if len(text) <= max_chars:
        return text
    # Sentence-end on . ! ? followed by space + capital, but be lenient.
    match = re.search(r"(.{20,200}?[.!?])\s", text[: max_chars + 50])
    if match:
        return match.group(1).strip()
    return text[:max_chars].rstrip() + "…"


def _figure_label_from_jats(fig: ET.Element, fallback_index: int) -> str:
    label_elem = None
    for child in fig.iter():
        if _local(child.tag) == "label":
            label_elem = child
            break
    label = _itertext_clean(label_elem)
    if label:
        return label
    return f"Figure {fallback_index}"


def _resolve_jats_graphic_href(fig: ET.Element) -> Optional[str]:
    for child in fig.iter():
        if _local(child.tag) != "graphic":
            continue
        for attr_name, attr_val in child.attrib.items():
            if attr_name.endswith("href"):
                return attr_val
    return None


def _caption_text_from_jats(fig: ET.Element) -> tuple[str, str]:
    """Return (title, full_body) for a JATS <fig>/<table-wrap>/<supp> caption."""
    caption_elem = None
    for child in fig:
        if _local(child.tag) == "caption":
            caption_elem = child
            break

    title_text = ""
    body_parts: List[str] = []

    if caption_elem is not None:
        for sub in caption_elem:
            if _local(sub.tag) == "title":
                title_text = _itertext_clean(sub)
            else:
                txt = _itertext_clean(sub)
                if txt:
                    body_parts.append(txt)
        # Some JATS files put text directly in <caption>
        direct = (caption_elem.text or "").strip()
        if direct:
            body_parts.insert(0, re.sub(r"\s+", " ", direct))

    # Some figures use <p> or <alt-text> outside <caption>.
    for child in fig:
        local = _local(child.tag)
        if local in {"caption", "label", "graphic"}:
            continue
        if local in {"p", "alt-text"}:
            txt = _itertext_clean(child)
            if txt:
                body_parts.append(txt)

    full_body = " ".join(part for part in body_parts if part).strip()
    if not title_text:
        title_text = _first_sentence(full_body)
    return title_text, full_body


def extract_from_jats_xml(xml_content: str) -> CaptionExtractionResult:
    """Extract figure, table, and supplement captions from PMC JATS XML."""
    result = CaptionExtractionResult()
    if not xml_content or not xml_content.strip():
        return result

    try:
        root = ET.fromstring(xml_content)
    except ET.ParseError:
        return result

    fig_idx = 0
    for elem in root.iter():
        local = _local(elem.tag)
        if local == "fig":
            fig_idx += 1
            label = _figure_label_from_jats(elem, fig_idx)
            title, body = _caption_text_from_jats(elem)
            href = _resolve_jats_graphic_href(elem)
            fig_id = elem.attrib.get("id")
            if title or body:
                result.figures.append(
                    FigureCaption(
                        label=label,
                        title=title,
                        text=body,
                        image_url=href,
                        figure_id=fig_id,
                    )
                )

    tbl_idx = 0
    for elem in root.iter():
        if _local(elem.tag) != "table-wrap":
            continue
        tbl_idx += 1
        label_elem = None
        for child in elem:
            if _local(child.tag) == "label":
                label_elem = child
                break
        label = _itertext_clean(label_elem) or f"Table {tbl_idx}"
        title, body = _caption_text_from_jats(elem)
        if title or body:
            result.tables.append(
                TableCaption(
                    label=label,
                    title=title,
                    text=body,
                    table_id=elem.attrib.get("id"),
                )
            )

    supp_idx = 0
    for elem in root.iter():
        if _local(elem.tag) != "supplementary-material":
            continue
        supp_idx += 1
        label_elem = None
        media_type = None
        href = None
        for child in elem.iter():
            local = _local(child.tag)
            if label_elem is None and local == "label":
                label_elem = child
            if local == "media" or local == "supplementary-material":
                for attr_name, attr_val in child.attrib.items():
                    if attr_name.endswith("href"):
                        href = href or attr_val
                    if attr_name.endswith("mimetype") or attr_name.endswith(
                        "mime-subtype"
                    ):
                        media_type = (
                            (media_type or "") + ("/" if media_type else "") + attr_val
                        )
        label = _itertext_clean(label_elem) or f"Supplementary Material {supp_idx}"
        title, body = _caption_text_from_jats(elem)
        if title or body or href:
            result.supplements.append(
                SupplementDescription(
                    label=label,
                    title=title,
                    text=body,
                    href=href,
                    media_type=media_type,
                )
            )

    return result


# ---------------------------------------------------------------------------
# HTML extraction (PMC HTML pages and publisher HTML)
# ---------------------------------------------------------------------------


_FIGURE_LABEL_RE = re.compile(r"^\s*(figure|fig\.?|f)\s*\d+[a-z]?", re.IGNORECASE)
_TABLE_LABEL_RE = re.compile(r"^\s*table\s*\d+[a-z]?", re.IGNORECASE)
_SUPP_LABEL_RE = re.compile(
    r"^\s*(supplement(ary|al)?|appendix|supporting)\s",
    re.IGNORECASE,
)


def _html_text(elem: Optional[Tag]) -> str:
    if elem is None:
        return ""
    text = elem.get_text(" ", strip=True)
    return re.sub(r"\s+", " ", text).strip()


def _resolve_image_url(img: Tag) -> Optional[str]:
    for attr in ("src", "data-src", "data-lazy-src", "data-original"):
        val = img.get(attr)
        if val:
            return val.strip()
    return None


def _figcaption_text(fig: Tag) -> tuple[str, str]:
    figcap = fig.find("figcaption")
    if figcap is not None:
        # Title: first <strong>/<b>/<h*>, otherwise first sentence
        title_node = figcap.find(["strong", "b", "h1", "h2", "h3", "h4"])
        title = _html_text(title_node) if title_node else ""
        body = _html_text(figcap)
        if not title:
            title = _first_sentence(body)
        return title, body
    # Some publisher HTML uses adjacent .caption div
    caption_div = fig.find(class_=re.compile(r"caption|legend", re.IGNORECASE))
    if caption_div is not None:
        body = _html_text(caption_div)
        return _first_sentence(body), body
    return "", ""


def _extract_figures_from_html(soup: BeautifulSoup) -> List[FigureCaption]:
    out: List[FigureCaption] = []
    seen_ids: set[str] = set()

    # Pass 1: <figure> elements (preferred — semantic).
    for idx, fig in enumerate(soup.find_all("figure"), 1):
        title, body = _figcaption_text(fig)
        if not (title or body):
            continue
        img = fig.find("img")
        href = _resolve_image_url(img) if img else None
        fid = fig.get("id")
        label = ""
        # Try to find a label inside (e.g. "Figure 1.")
        label_node = fig.find(class_=re.compile(r"label", re.IGNORECASE))
        if label_node:
            label = _html_text(label_node)
        if not label:
            # First word/two of caption often is "Figure N."
            head = (title or body).split(".", 1)[0]
            if _FIGURE_LABEL_RE.match(head):
                label = head.strip()
        if not label:
            label = f"Figure {idx}"
        if fid and fid in seen_ids:
            continue
        if fid:
            seen_ids.add(fid)
        out.append(
            FigureCaption(
                label=label,
                title=title,
                text=body,
                image_url=href,
                figure_id=fid,
            )
        )

    # Pass 2: divs that *look* like figures but aren't <figure>.
    if not out:
        for idx, div in enumerate(
            soup.find_all(
                "div",
                class_=re.compile(r"(^|[\s_-])fig(\b|ure)", re.IGNORECASE),
            ),
            1,
        ):
            cap = div.find(class_=re.compile(r"caption|legend", re.IGNORECASE))
            if cap is None:
                continue
            body = _html_text(cap)
            if not body:
                continue
            img = div.find("img")
            href = _resolve_image_url(img) if img else None
            fid = div.get("id")
            head = body.split(".", 1)[0]
            label = head.strip() if _FIGURE_LABEL_RE.match(head) else f"Figure {idx}"
            out.append(
                FigureCaption(
                    label=label,
                    title=_first_sentence(body),
                    text=body,
                    image_url=href,
                    figure_id=fid,
                )
            )

    return out


def _extract_table_captions_from_html(soup: BeautifulSoup) -> List[TableCaption]:
    out: List[TableCaption] = []
    for idx, tbl in enumerate(soup.find_all("table"), 1):
        cap_elem = tbl.find("caption")
        body = _html_text(cap_elem) if cap_elem else ""
        if not body:
            # Some publishers put captions in a sibling div labeled .table-caption
            parent = tbl.find_parent(
                class_=re.compile(r"table.*wrap|table.*caption", re.IGNORECASE)
            )
            if parent is not None:
                cap_div = parent.find(
                    class_=re.compile(r"caption|legend|title", re.IGNORECASE)
                )
                body = _html_text(cap_div)
        if not body:
            continue
        head = body.split(".", 1)[0]
        label = head.strip() if _TABLE_LABEL_RE.match(head) else f"Table {idx}"
        out.append(
            TableCaption(
                label=label,
                title=_first_sentence(body),
                text=body,
                table_id=tbl.get("id"),
            )
        )

    # PMC's current HTML often renders tables as image-only ``section.tw``
    # blocks with a caption div and an <img> rather than as text tables.
    # Those image URLs are valuable input for the optional figure OCR pass.
    seen_keys = {((tbl.table_id or "").lower(), tbl.text[:120].lower()) for tbl in out}
    image_table_sections = soup.select(
        "section.tw, section.table-wrap, div.table-wrap, div[class*='table-wrap']"
    )
    for idx, section in enumerate(image_table_sections, len(out) + 1):
        cap_elem = (
            section.find(class_=re.compile(r"caption", re.IGNORECASE))
            or section.find("caption")
            or section.find(["figcaption"])
        )
        body = _html_text(cap_elem) if cap_elem else ""
        if not body:
            continue
        key = ((section.get("id") or "").lower(), body[:120].lower())
        if key in seen_keys:
            continue
        seen_keys.add(key)

        label_node = cap_elem.find(["strong", "b"]) if cap_elem else None
        label = _html_text(label_node) if label_node else ""
        if not label:
            head = body.split(".", 1)[0]
            label = head.strip() if _TABLE_LABEL_RE.match(head) else f"Table {idx}"

        img = section.find("img")
        href = _resolve_image_url(img) if img else None
        out.append(
            TableCaption(
                label=label,
                title=_first_sentence(body),
                text=body,
                table_id=section.get("id"),
                image_url=href,
            )
        )
    return out


def _extract_supplement_descriptions_from_html(
    soup: BeautifulSoup,
) -> List[SupplementDescription]:
    out: List[SupplementDescription] = []

    # Look for sections that look like supplementary-materials.
    candidates: List[Tag] = []
    selectors = [
        "section[id*='supplement']",
        "section[class*='supplement']",
        "div[id*='supplement']",
        "div[class*='supplement']",
        "div#mc-supplementary-materials",
        "section#supplementary-information",
    ]
    seen: set[int] = set()
    for selector in selectors:
        for elem in soup.select(selector):
            if id(elem) in seen:
                continue
            seen.add(id(elem))
            candidates.append(elem)

    for section in candidates:
        # Each supplement entry is typically a <p> or <div> with a label
        # ("Supplementary Table 1", "File S1") followed by a description
        # and a download link.
        for entry in section.find_all(["p", "li", "div"], recursive=True):
            text = _html_text(entry)
            if not text:
                continue
            if not _SUPP_LABEL_RE.match(text) and "supplement" not in text.lower():
                continue
            # Extract a label from leading bold/strong if present
            label_node = entry.find(["strong", "b"])
            label = _html_text(label_node) if label_node else text.split(":", 1)[0]
            label = label.strip().rstrip(".:")
            if len(label) > 80:
                label = label[:80] + "…"
            link = entry.find("a", href=True)
            href = link["href"].strip() if link else None
            out.append(
                SupplementDescription(
                    label=label or "Supplement",
                    title=_first_sentence(text),
                    text=text,
                    href=href,
                )
            )

    # Deduplicate by (label, title).
    deduped: List[SupplementDescription] = []
    seen_keys: set[tuple[str, str]] = set()
    for item in out:
        key = (item.label.lower(), item.title.lower())
        if key in seen_keys:
            continue
        seen_keys.add(key)
        deduped.append(item)
    return deduped


def extract_from_html(html_content: str) -> CaptionExtractionResult:
    """Extract figure / table / supplement captions from a publisher HTML page."""
    result = CaptionExtractionResult()
    if not html_content or not html_content.strip():
        return result

    soup = BeautifulSoup(html_content, "html.parser")
    result.figures = _extract_figures_from_html(soup)
    result.tables = _extract_table_captions_from_html(soup)
    result.supplements = _extract_supplement_descriptions_from_html(soup)
    return result


# ---------------------------------------------------------------------------
# Markdown rendering
# ---------------------------------------------------------------------------


def _format_figure_caption_md(fig: FigureCaption) -> str:
    parts = [f"### {fig.label}"]
    if fig.title and fig.title != fig.label:
        parts.append(f"**{fig.title}**")
    if fig.text and fig.text not in (fig.title, fig.label):
        parts.append(fig.text)
    if fig.image_url:
        parts.append(f"_image_: {fig.image_url}")
    return "\n\n".join(parts)


def _format_table_caption_md(tbl: TableCaption) -> str:
    parts = [f"### {tbl.label}"]
    if tbl.title and tbl.title != tbl.label:
        parts.append(f"**{tbl.title}**")
    if tbl.text and tbl.text not in (tbl.title, tbl.label):
        parts.append(tbl.text)
    if tbl.image_url:
        parts.append(f"_image_: {tbl.image_url}")
    return "\n\n".join(parts)


def _format_supplement_desc_md(supp: SupplementDescription) -> str:
    parts = [f"### {supp.label}"]
    if supp.title and supp.title != supp.label:
        parts.append(f"**{supp.title}**")
    if supp.text and supp.text not in (supp.title, supp.label):
        parts.append(supp.text)
    if supp.href:
        parts.append(f"_link_: {supp.href}")
    return "\n\n".join(parts)


def render_captions_markdown(result: CaptionExtractionResult) -> str:
    """Render a CaptionExtractionResult as a markdown block.

    Empty sections are omitted. Returns "" when nothing was extracted.
    """
    if result.is_empty():
        return ""

    sections: List[str] = []

    if result.figures:
        sections.append(
            "## FIGURE CAPTIONS\n\n"
            + "\n\n".join(_format_figure_caption_md(f) for f in result.figures)
        )

    if result.tables:
        sections.append(
            "## TABLE CAPTIONS\n\n"
            + "\n\n".join(_format_table_caption_md(t) for t in result.tables)
        )

    if result.supplements:
        sections.append(
            "## SUPPLEMENT DESCRIPTIONS\n\n"
            + "\n\n".join(_format_supplement_desc_md(s) for s in result.supplements)
        )

    return "\n\n" + "\n\n".join(sections) + "\n\n"


def merge_results(*results: CaptionExtractionResult) -> CaptionExtractionResult:
    """Merge multiple extraction results, deduplicating by content fingerprint."""
    merged = CaptionExtractionResult()
    seen_fig: set[str] = set()
    seen_tbl: set[str] = set()
    seen_supp: set[str] = set()

    for r in results:
        for f in r.figures:
            key = (f.figure_id or "") + "|" + f.text[:120].lower()
            if key in seen_fig:
                continue
            seen_fig.add(key)
            merged.figures.append(f)
        for t in r.tables:
            key = (t.table_id or "") + "|" + t.text[:120].lower()
            if key in seen_tbl:
                continue
            seen_tbl.add(key)
            merged.tables.append(t)
        for s in r.supplements:
            key = s.label.lower() + "|" + s.text[:120].lower()
            if key in seen_supp:
                continue
            seen_supp.add(key)
            merged.supplements.append(s)

    return merged


def save_captions_json(result: CaptionExtractionResult, dest: Path) -> None:
    """Persist captions as JSON next to the figures directory."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(result.to_dict(), indent=2), encoding="utf-8")
