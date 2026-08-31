"""
Genetic Data Scout module for identifying high-value data zones in biomedical literature.

This module analyzes full-text markdown (main text + supplements) to identify sections
containing patient-level data or variant lists for a target gene. It outputs a JSON
list of "Data Zones" and can generate a condensed markdown file with only high-value content.

Role: Genetic Data Scout
Task: Analyze text (or supplement chunk). Identify ANY section that contains specific
      patient-level data or variant lists for the target gene.
"""

import json
import logging
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from config.constants import SCOUT_CLINICAL_KEYWORDS, SCOUT_PREFER_CLEANED_ABOVE_CHARS
from pipeline.source_quality import demote_empty_full_context
from utils.scout_models import DataZone, DataZoneReport

logger = logging.getLogger(__name__)

# Hard bounds so a single pathological file cannot wedge a whole run (RYR2
# PMID 32508047: a 14.5 MB binary-garbage supplement fold produced ~18K fake
# pipe tables and the per-candidate bookkeeping went quadratic — 2h45m+ at
# 100% CPU). Every regex call gets bounded input, every per-candidate lookup
# gets a bounded window, and scan() enforces a wall-clock budget between them.
#
# Longest regex input taken from one physical line. Real table rows top out
# around a few KB; anything longer is a converter artifact, and one match in
# the first 20K chars is enough to classify the line.
_MAX_LINE_SCAN_CHARS = 20_000
# Longest regex input for signal *counting* (paragraph triggers, zone
# scoring). Zones keep their full extent; only the heuristic counts are
# computed on a bounded prefix.
_MAX_SIGNAL_SCAN_CHARS = 200_000
# _find_parent_section only accepts a header ending within 10_000 chars of
# the zone, so searching a slightly larger window is semantics-preserving
# (the slack covers the header match's own length).
_PARENT_SECTION_DISTANCE_CHARS = 10_000
_PARENT_SECTION_WINDOW_CHARS = _PARENT_SECTION_DISTANCE_CHARS + 100
# Caption/header walk-backs above a detected table. Real captions sit within
# a few lines; without a cap, adjacent pipe-line runs make each separator row
# re-walk the whole region above it (quadratic).
_CAPTION_LOOKBACK_LINES = 200
_PSEUDO_HEADER_LOOKBACK_LINES = 50
# How often the line loops poll the scan deadline.
_DEADLINE_CHECK_EVERY_LINES = 256


class ScanBudgetExceeded(Exception):
    """Raised internally when a scan exceeds its wall-clock budget."""


def select_scout_source_path(
    full_context_path: Path,
    prefer_cleaned_above_chars: int = SCOUT_PREFER_CLEANED_ABOVE_CHARS,
) -> Path:
    """Return the best text source for scouting a FULL_CONTEXT artifact.

    The raw FULL_CONTEXT file remains the canonical audit artifact, but the
    scout stage can safely use an existing CLEANED sibling when the raw file is
    very large. Preprocessing is deterministic and non-destructive.
    """

    # An empty/truncated FULL_CONTEXT would trivially pass the size check below
    # and scout nothing while a populated CLEANED sits beside it.
    full_context_path = demote_empty_full_context(full_context_path)

    if full_context_path.stat().st_size <= prefer_cleaned_above_chars:
        return full_context_path

    cleaned_path = full_context_path.with_name(
        full_context_path.name.replace("_FULL_CONTEXT.md", "_CLEANED.md")
    )
    if cleaned_path.exists() and cleaned_path.stat().st_size > 0:
        logger.warning(
            "Using cleaned scout source for oversized context %s "
            "(raw %.1f MB, cleaned %.1f MB)",
            full_context_path.name,
            full_context_path.stat().st_size / 1_000_000,
            cleaned_path.stat().st_size / 1_000_000,
        )
        return cleaned_path

    logger.warning(
        "Oversized context %s has no cleaned sibling; scouting raw %.1f MB file",
        full_context_path.name,
        full_context_path.stat().st_size / 1_000_000,
    )
    return full_context_path


@dataclass
class ZoneCandidate:
    """Internal representation of a potential data zone before scoring."""

    start: int
    end: int
    zone_type: str  # "TABLE" or "TEXT"
    raw_text: str
    source_section: Optional[str] = None
    table_caption: Optional[str] = None


class GeneticDataScout:
    """
    Identifies and prioritizes high-value data zones in scientific papers.

    Scans for:
    - Section headers (Results, Supplementary, Tables)
    - Table captions and tabular structure
    - Clinical descriptors (patient, case, carrier, affected)
    - Variant nomenclature density (c., p. notation)
    """

    # Patterns for detecting individual-level data (keep=True signals)
    #
    # Repeat-bounding rule for every pattern class below: a repeat whose
    # failure forces the engine to re-count at each start position (an
    # unbounded token repeat followed by a required literal, or chained
    # `\s*`/optional-group runs) must carry an explicit upper bound. On
    # normal text the bounds are far above anything real; on pathological
    # content (multi-MB uppercase/digit runs, huge whitespace pads from
    # converter output) they turn quadratic scans into linear ones.
    INDIVIDUAL_PATTERNS = [
        r"\bpatient[- ]?\d+",  # Patient 1, Patient-1, Patient1
        r"\bcase[- ]?\d+",  # Case 1, Case-62, Case1
        r"\bproband\b",  # proband
        r"\bindex\s*case\b",  # index case
        r"\b[IVX]{1,7}-\d{1,4}\b",  # Pedigree notation: II-1, III-2, IV-1
        r"\b[Pp]\d+\b",  # P1, P2 (patient IDs)
        r"\bsubject\s*\d+",  # Subject 1
        r"\bfamily\s*[A-Z]\b",  # Family A, Family B
        r"\bindividual\s*\d+",  # Individual 1
        r"\d+-year-old\s*(male|female|man|woman|boy|girl)",  # Age/sex pattern
        # Age at onset: 45 (bounded gaps: three chained \s* around optionals
        # backtrack quadratically over long whitespace pads)
        r"\bage\s{0,10}(at\s{0,10})?(onset|diagnosis|evaluation)\s{0,10}[:=]?\s{0,10}\d+",
        r"\bpresented\s+with\b",  # "presented with"
        r"\bdiagnosed\s+with\b",  # "diagnosed with"
        r"\bcarrier\s+(status|of)\b",  # carrier status, carrier of
        r"\baffected\s+(status|individual|member|sibling|parent)",  # affected individual
        r"\bunaffected\s+(individual|member|sibling|parent|carrier)",  # unaffected carrier
    ]

    # Patterns for detecting aggregate data (keep=False signals)
    AGGREGATE_PATTERNS = [
        r"\b\d+%\s+of\s+(patients|carriers|individuals|cases)",  # 20% of patients
        r"\bon\s+average\b",  # on average
        r"\boverall\b",  # overall
        r"\bmean\s{0,10}(age|value)?\s{0,10}[:=]?\s{0,10}\d+",  # mean age, mean value
        r"\bmedian\s{0,10}(age|value)?\s{0,10}[:=]?\s{0,10}\d+",  # median age
        r"\bstatistical(ly)?\s*(analysis|significant)",  # statistical analysis
        r"\bp\s*[<>=]\s*0\.\d+",  # p-value notation
        r"\bpreviously\s+reported\b",  # previously reported
        r"\bas\s+described\s+(previously|elsewhere|in)\b",  # as described previously
        r"\breview\s+of\s+(the\s+)?literature\b",  # review of literature
        r"\bmeta-analysis\b",  # meta-analysis
        r"\bpooled\s+(analysis|data)\b",  # pooled analysis
    ]

    # Patterns for variant nomenclature (expanded for better detection)
    VARIANT_PATTERNS = [
        r"c\.\d+[ACGT]>[ACGT]",  # cDNA: c.1234G>A
        r"c\.\d+[_\+\-]\d*[delinsdup]+",  # cDNA: c.1234_1235del
        r"c\.\d+[\+\-]\d+[ACGT]>[ACGT]",  # Intronic: c.1234+5G>A
        r"p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}",  # Protein: p.Arg412His
        r"p\.[A-Z][a-z]{2}\d+\*",  # Nonsense: p.Arg412*
        r"p\.[A-Z][a-z]{2}\d+fs",  # Frameshift: p.Arg412fs
        r"p\.\([A-Z][a-z]{2}\d+[A-Z][a-z]{2}\)",  # Protein with parens: p.(Arg412His)
        r"[A-Z]\d+[A-Z]",  # Single-letter: R412H, G628S
        r"[A-Z][a-z]{2}\d+[A-Z][a-z]{2}",  # Three-letter without p.: Arg412His
        # Nucleotide position: 1234 G>A (bounded: unbounded \d+/\s* re-count
        # long digit runs at every start position when no >[ACGT] follows)
        r"\d{1,10}\s{0,5}[ACGT]\s{0,5}>\s{0,5}[ACGT]",
        r"del[A-Z]{2,}",  # Deletion: delAG, delCT
        r"ins[A-Z]{2,}",  # Insertion: insAG
        r"IVS\d+[+-]\d+[ACGT]>[ACGT]",  # Intronic legacy: IVS2+1G>A
        # Gene-mutation: KCNH2-T895M, SCN5A-G1743R (bounded gene token: an
        # unbounded [A-Z0-9]{2,} re-counts multi-MB uppercase/digit runs —
        # DNA/protein sequence — at every start position, O(n^2))
        r"[A-Z][A-Z0-9]{2,11}-[A-Z]\d{1,5}[A-Z]",
    ]

    # Clinical keywords for relevance scoring (from centralized constants)
    CLINICAL_KEYWORDS = SCOUT_CLINICAL_KEYWORDS

    # Section headers to identify source context
    SECTION_PATTERNS = {
        "results": r"#{1,4}\s*results?\b",
        "methods": r"#{1,4}\s*methods?\b",
        "discussion": r"#{1,4}\s*discussion\b",
        "supplementary": r"#{1,4}\s*(supplementary|supplemental|supporting)\b",
        "tables": r"#{1,4}\s*table\s*\d*",
        "case_report": r"#{1,4}\s*case\s*(report|description|study)",
        "clinical": r"#{1,4}\s*clinical\s*(data|characteristics|features)",
        "patients": r"#{1,4}\s*patients?\s*(and\s*methods)?",
    }

    # Table caption patterns (high-value signals)
    TABLE_CAPTION_PATTERNS = [
        r"table\s*\d+[.:]\s*",
        r"supplementary\s*table\s*\w*[.:]\s*",
        r"supplemental\s*table\s*\w*[.:]\s*",
        r"s\d+\s*table[.:]\s*",
    ]

    def __init__(
        self,
        gene_symbol: str,
        min_relevance_score: float = 0.1,
        max_zones: int = 30,
        scan_budget_seconds: Optional[float] = None,
    ):
        """
        Initialize the Genetic Data Scout.

        Args:
            gene_symbol: Target gene symbol to search for (e.g., "KCNQ1", "BRCA1")
            min_relevance_score: Minimum score (0.0-1.0) to include a zone
            max_zones: Maximum number of zones to return
            scan_budget_seconds: Wall-clock budget for a single scan() call.
                None resolves from settings (SCOUT_SCAN_BUDGET_SECONDS);
                0 or negative disables the budget. When the budget runs out,
                scan() returns an empty report with ``scan_truncated=True``
                instead of wedging the run — callers must then fall back to
                the un-condensed source rather than write a partial
                DATA_ZONES file that silently hides the unscanned tail.
        """
        self.gene_symbol = gene_symbol.upper()
        self.scan_budget_seconds = scan_budget_seconds
        self._deadline: Optional[float] = None
        self.gene_pattern = re.compile(
            rf"\b{re.escape(self.gene_symbol)}\b", re.IGNORECASE
        )
        self.min_relevance_score = min_relevance_score
        self.max_zones = max_zones

        # Compile patterns for efficiency
        self._individual_patterns = [
            re.compile(p, re.IGNORECASE) for p in self.INDIVIDUAL_PATTERNS
        ]
        self._aggregate_patterns = [
            re.compile(p, re.IGNORECASE) for p in self.AGGREGATE_PATTERNS
        ]
        self._variant_patterns = [
            re.compile(p, re.IGNORECASE) for p in self.VARIANT_PATTERNS
        ]
        self._section_patterns = {
            k: re.compile(v, re.IGNORECASE) for k, v in self.SECTION_PATTERNS.items()
        }
        self._table_caption_patterns = [
            re.compile(p, re.IGNORECASE) for p in self.TABLE_CAPTION_PATTERNS
        ]

    def scan(self, text: str, pmid: Optional[str] = None) -> DataZoneReport:
        """
        Scan text for high-value data zones.

        Args:
            text: Full markdown text (main text + supplements)
            pmid: Optional PubMed ID for the report

        Returns:
            DataZoneReport with all identified zones and metadata
        """
        logger.info(
            f"Scanning for data zones in {len(text)} characters for gene {self.gene_symbol}"
        )

        budget = self.scan_budget_seconds
        if budget is None:
            budget = self._settings_scan_budget_seconds()
        started = time.monotonic()
        self._deadline = started + budget if budget and budget > 0 else None

        try:
            # Detect all candidate zones
            candidates = []
            candidates.extend(self._detect_tables(text))
            candidates.extend(self._detect_clinical_text(text))

            # Merge overlapping candidates
            candidates = self._merge_overlapping(candidates)

            # Score and convert to DataZone objects
            zones = []
            for candidate in candidates:
                self._check_deadline()
                zone = self._score_and_create_zone(candidate, text, len(zones))
                # Always keep TABLE zones (they often contain variant data)
                # For TEXT zones, apply relevance threshold
                if (
                    candidate.zone_type == "TABLE"
                    or zone.relevance_score >= self.min_relevance_score
                ):
                    zones.append(zone)
        except ScanBudgetExceeded:
            elapsed = time.monotonic() - started
            logger.warning(
                "Data scout scan budget exhausted for PMID %s after %.1fs "
                "(%d chars, budget %.0fs); returning no zones so extraction "
                "falls back to the un-condensed source",
                pmid,
                elapsed,
                len(text),
                budget or 0.0,
            )
            return DataZoneReport(
                pmid=pmid,
                gene_symbol=self.gene_symbol,
                total_zones_found=0,
                zones_kept=0,
                zones_discarded=0,
                total_chars_original=len(text),
                total_chars_condensed=0,
                compression_ratio=0.0,
                zones=[],
                scan_truncated=True,
            )
        finally:
            self._deadline = None

        # Sort by relevance score (descending) and limit
        zones.sort(key=lambda z: z.relevance_score, reverse=True)
        zones = zones[: self.max_zones]

        # Re-sort by position for output
        zones.sort(key=lambda z: z.char_start)

        # Calculate report metrics
        kept_zones = [z for z in zones if z.keep]
        condensed_chars = sum(z.length for z in kept_zones)

        report = DataZoneReport(
            pmid=pmid,
            gene_symbol=self.gene_symbol,
            total_zones_found=len(zones),
            zones_kept=len(kept_zones),
            zones_discarded=len(zones) - len(kept_zones),
            total_chars_original=len(text),
            total_chars_condensed=condensed_chars,
            compression_ratio=condensed_chars / len(text) if len(text) > 0 else 0.0,
            zones=zones,
        )

        logger.info(
            f"Found {report.total_zones_found} zones, keeping {report.zones_kept} "
            f"({report.compression_ratio:.1%} of original)"
        )

        return report

    @staticmethod
    def _settings_scan_budget_seconds() -> float:
        """Resolve the default scan budget from settings (0 = unlimited)."""
        try:
            from config.settings import get_settings

            settings = get_settings()
            return float(settings.scout_scan_budget_seconds) if settings else 0.0
        except Exception:  # noqa: BLE001 - budget is a guard, never a crash
            return 0.0

    def _check_deadline(self) -> None:
        """Raise ScanBudgetExceeded once the scan's wall-clock budget is spent."""
        if self._deadline is not None and time.monotonic() > self._deadline:
            raise ScanBudgetExceeded

    @staticmethod
    def _line_offsets(lines: list[str]) -> list[int]:
        """Prefix sums: offsets[i] is the char offset where line i starts.

        Replaces the per-candidate ``sum(len(lines[j]) + 1 for j in range(k))``
        bookkeeping, which re-walked every preceding line for every detected
        table (quadratic once a pathological file yields thousands of them).
        Matches that sum exactly, including the trailing +1 for the final
        line's absent newline.
        """
        offsets = [0] * (len(lines) + 1)
        total = 0
        for i, line in enumerate(lines):
            total += len(line) + 1
            offsets[i + 1] = total
        return offsets

    def _detect_tables(self, text: str) -> list[ZoneCandidate]:
        """
        Detect markdown tables and their surrounding context.

        Looks for:
        - Markdown table markers (|---|)
        - Table captions (Table 1:, Supplementary Table S1:)
        - Table content rows
        - Pseudo-tabular data (whitespace-aligned or TSV-like from PDF extraction)
        """
        candidates = []
        lines = text.split("\n")
        offsets = self._line_offsets(lines)

        # First pass: detect markdown tables
        candidates.extend(self._detect_markdown_tables(lines, text, offsets))

        # Second pass: detect pseudo-tabular data (common in PDF extractions)
        candidates.extend(self._detect_pseudo_tables(lines, text, offsets))

        return candidates

    def _detect_pseudo_tables(
        self, lines: list[str], text: str, offsets: list[int]
    ) -> list[ZoneCandidate]:
        """
        Detect pseudo-tabular data that isn't in markdown format.

        This catches variant tables extracted from PDFs that appear as:
        - Tab-separated values
        - Whitespace-aligned columns
        - Lines with multiple variant-like patterns
        """
        candidates = []
        i = 0

        while i < len(lines):
            if i % _DEADLINE_CHECK_EVERY_LINES == 0:
                self._check_deadline()
            # Bounded regex input: one match in the head of an absurdly long
            # "line" (a converter artifact, not a table row) is enough.
            line = lines[i][:_MAX_LINE_SCAN_CHARS]

            # Check for lines with tab-separated or multi-column variant data
            has_tabs = "\t" in line
            has_multiple_spaces = "  " in line  # Double-space column separator
            variant_matches = sum(1 for p in self._variant_patterns if p.search(line))

            # Detect start of pseudo-table: line with tabs/spaces and variant patterns
            if (has_tabs or has_multiple_spaces) and variant_matches >= 1:
                table_start = i

                # Look backward for table header/caption (bounded: without a
                # cap, keyword-dense regions make every detection re-walk the
                # whole region above it)
                while table_start > max(0, i - _PSEUDO_HEADER_LOOKBACK_LINES):
                    prev_line = lines[table_start - 1].strip().lower()
                    if not prev_line:
                        break
                    # Table headers often contain these keywords
                    if any(
                        kw in prev_line
                        for kw in [
                            "mutation",
                            "variant",
                            "nucleotide",
                            "amino acid",
                            "patient",
                            "phenotype",
                            "table",
                        ]
                    ):
                        table_start -= 1
                    else:
                        break

                # Look forward for table end
                table_end = i + 1
                while table_end < len(lines):
                    if table_end % _DEADLINE_CHECK_EVERY_LINES == 0:
                        self._check_deadline()
                    next_line = lines[table_end][:_MAX_LINE_SCAN_CHARS]
                    next_has_structure = "\t" in next_line or "  " in next_line
                    next_variant_matches = sum(
                        1 for p in self._variant_patterns if p.search(next_line)
                    )
                    # Continue if line has similar structure or is short (table row)
                    if (
                        next_has_structure
                        or next_variant_matches >= 1
                        or (len(next_line.strip()) < 100 and next_line.strip())
                    ):
                        table_end += 1
                    else:
                        break

                # Only create zone if we found a substantial table (3+ rows)
                if table_end - table_start >= 3:
                    char_start = offsets[table_start]
                    char_end = offsets[table_end]
                    raw_text = "\n".join(lines[table_start:table_end])

                    candidates.append(
                        ZoneCandidate(
                            start=char_start,
                            end=char_end,
                            zone_type="TABLE",
                            raw_text=raw_text,
                            table_caption=lines[table_start].strip()
                            if table_start > 0
                            else None,
                            source_section=self._find_parent_section(text, char_start),
                        )
                    )

                    i = table_end
                    continue

            i += 1

        return candidates

    def _detect_markdown_tables(
        self, lines: list[str], text: str, offsets: list[int]
    ) -> list[ZoneCandidate]:
        """Detect standard markdown tables with |---| separators."""
        candidates = []

        i = 0
        while i < len(lines):
            if i % _DEADLINE_CHECK_EVERY_LINES == 0:
                self._check_deadline()
            line = lines[i][:_MAX_LINE_SCAN_CHARS]

            # Check for table separator row (|---|---|)
            if re.match(r"\s*\|[-:| ]+\|", line):
                # Found a table - find its boundaries
                table_start_line = i

                # Look backward for table caption and header (bounded: in a
                # pipe-dense region — e.g. binary garbage folded as text —
                # every separator row would otherwise re-walk all the pipe
                # lines above it, going quadratic)
                caption = None
                while table_start_line > max(0, i - _CAPTION_LOOKBACK_LINES):
                    prev_line = lines[table_start_line - 1].strip()
                    if not prev_line:
                        break
                    if prev_line.startswith("|") or any(
                        p.search(prev_line[:_MAX_LINE_SCAN_CHARS])
                        for p in self._table_caption_patterns
                    ):
                        table_start_line -= 1
                        if any(
                            p.search(prev_line[:_MAX_LINE_SCAN_CHARS])
                            for p in self._table_caption_patterns
                        ):
                            caption = prev_line
                    else:
                        break

                # Look forward for table end
                table_end_line = i + 1
                while table_end_line < len(lines):
                    if table_end_line % _DEADLINE_CHECK_EVERY_LINES == 0:
                        self._check_deadline()
                    next_line = lines[table_end_line].strip()
                    if next_line.startswith("|"):
                        table_end_line += 1
                    elif not next_line:
                        # Empty line - check if there's a note/footnote
                        if table_end_line + 1 < len(lines):
                            following = lines[table_end_line + 1].strip()
                            if following.startswith("*") or following.startswith("†"):
                                table_end_line += 2
                            else:
                                break
                        else:
                            break
                    else:
                        break

                # Calculate character positions
                char_start = offsets[table_start_line]
                char_end = offsets[table_end_line]

                raw_text = "\n".join(lines[table_start_line:table_end_line])

                candidates.append(
                    ZoneCandidate(
                        start=char_start,
                        end=char_end,
                        zone_type="TABLE",
                        raw_text=raw_text,
                        table_caption=caption,
                        source_section=self._find_parent_section(text, char_start),
                    )
                )

                i = table_end_line
            else:
                i += 1

        return candidates

    def _detect_clinical_text(self, text: str) -> list[ZoneCandidate]:
        """
        Detect text sections containing clinical/patient data.

        Looks for:
        - Case report paragraphs
        - Patient descriptions with IDs
        - Results sections with individual data
        """
        candidates = []

        # Split by paragraphs (double newline or section headers)
        paragraph_pattern = re.compile(r"\n\s*\n|\n(?=#{1,4}\s)")
        paragraphs = paragraph_pattern.split(text)

        char_pos = 0
        for para in paragraphs:
            self._check_deadline()
            if not para.strip():
                char_pos += len(para)
                continue

            # Bounded regex input: signals from the head of a giant paragraph
            # (converter artifact) are representative enough to trigger a zone.
            scan_para = para[:_MAX_SIGNAL_SCAN_CHARS]

            # Check if paragraph has individual-level data signals
            individual_matches = sum(
                1 for p in self._individual_patterns if p.search(scan_para)
            )

            # Also check for variant mentions
            variant_matches = sum(
                1 for p in self._variant_patterns if p.search(scan_para)
            )

            # Check for gene mentions
            gene_matches = len(self.gene_pattern.findall(scan_para))

            # Check for clinical outcome keywords (penetrance, affected, carrier, etc.)
            para_lower = scan_para.lower()
            clinical_keyword_count = sum(
                1 for kw in self.CLINICAL_KEYWORDS if kw in para_lower
            )

            # If paragraph has strong signals, include it
            # RELAXED: Also include paragraphs with multiple clinical keywords
            # (catches penetrance prose like "50 carriers, 30 affected")
            has_clinical_outcome_signal = clinical_keyword_count >= 3
            if (
                individual_matches >= 2
                or (
                    individual_matches >= 1
                    and (variant_matches >= 1 or gene_matches >= 1)
                )
                or has_clinical_outcome_signal
            ):
                # Extend to include context (preceding header if any)
                start_pos = char_pos
                end_pos = char_pos + len(para)

                # Look for section header before this paragraph
                header_match = re.search(
                    r"(#{1,4}\s+[^\n]+\n)", text[max(0, start_pos - 200) : start_pos]
                )
                if header_match:
                    header_offset = text[max(0, start_pos - 200) : start_pos].rfind(
                        header_match.group(1)
                    )
                    if header_offset >= 0:
                        start_pos = max(0, start_pos - 200) + header_offset

                candidates.append(
                    ZoneCandidate(
                        start=start_pos,
                        end=end_pos,
                        zone_type="TEXT",
                        raw_text=text[start_pos:end_pos],
                        source_section=self._find_parent_section(text, start_pos),
                    )
                )

            char_pos += len(para)
            # Account for the delimiter we split on (pattern.match with a pos
            # argument — slicing text[char_pos:] here copied the remaining
            # megabytes once per paragraph)
            if char_pos < len(text):
                match = paragraph_pattern.match(text, char_pos)
                if match:
                    char_pos += len(match.group(0))

        return candidates

    def _merge_overlapping(
        self, candidates: list[ZoneCandidate]
    ) -> list[ZoneCandidate]:
        """Merge overlapping zone candidates."""
        if not candidates:
            return []

        # Sort by start position
        sorted_candidates = sorted(candidates, key=lambda c: c.start)

        merged = [sorted_candidates[0]]
        for current in sorted_candidates[1:]:
            prev = merged[-1]

            # Check for overlap (with some tolerance for nearby zones)
            if current.start <= prev.end + 100:
                # Merge: extend the previous zone
                prev.end = max(prev.end, current.end)
                prev.raw_text = (
                    prev.raw_text
                )  # Keep original text, will be recalculated
                # Prefer TABLE type if either is a table
                if current.zone_type == "TABLE":
                    prev.zone_type = "TABLE"
                    prev.table_caption = current.table_caption or prev.table_caption
            else:
                merged.append(current)

        return merged

    def _score_and_create_zone(
        self, candidate: ZoneCandidate, full_text: str, index: int
    ) -> DataZone:
        """
        Score a zone candidate and create a DataZone object.

        Scoring factors:
        - Gene mentions (high weight)
        - Variant nomenclature (high weight)
        - Individual-level patterns (high weight for keep)
        - Aggregate patterns (negative weight for keep)
        - Clinical keywords (medium weight)
        - Table with clinical data (bonus)
        - Supplementary section (bonus)
        """
        text = full_text[candidate.start : candidate.end]

        # Count signals on a bounded prefix: the zone keeps its full extent,
        # but heuristic counting over a multi-megabyte merged zone re-scans
        # it once per pattern (~35 passes) for scores that saturate anyway.
        scan_text = text[:_MAX_SIGNAL_SCAN_CHARS]
        gene_mentions = len(self.gene_pattern.findall(scan_text))
        variant_mentions = sum(
            len(p.findall(scan_text)) for p in self._variant_patterns
        )
        individual_signals = sum(
            len(p.findall(scan_text)) for p in self._individual_patterns
        )
        aggregate_signals = sum(
            len(p.findall(scan_text)) for p in self._aggregate_patterns
        )
        clinical_keywords = sum(
            scan_text.lower().count(kw) for kw in self.CLINICAL_KEYWORDS
        )

        # Calculate relevance score (0.0 to 1.0)
        score = 0.0

        # Gene presence is critical
        if gene_mentions > 0:
            score += 0.3

        # Variant nomenclature is strong signal
        if variant_mentions > 0:
            score += min(0.25, variant_mentions * 0.05)

        # Individual-level data patterns
        score += min(0.25, individual_signals * 0.05)

        # Clinical keyword density
        score += min(0.1, clinical_keywords * 0.01)

        # Table type bonus
        if candidate.zone_type == "TABLE":
            score += 0.1
            # Clinical table caption bonus
            if candidate.table_caption:
                caption_lower = candidate.table_caption.lower()
                if any(
                    kw in caption_lower
                    for kw in [
                        "clinical",
                        "patient",
                        "carrier",
                        "variant",
                        "genotype",
                        "phenotype",
                    ]
                ):
                    score += 0.1

        # Supplementary section bonus (often richest data source)
        if (
            candidate.source_section
            and "supplementary" in candidate.source_section.lower()
        ):
            score += 0.1

        # Normalize to 0.0-1.0
        score = min(1.0, score)

        # Determine keep/discard based on content signals
        # IMPORTANT: Be permissive - variant tables are often in supplements without prose
        keep = True

        # Only discard if aggregate data dominates AND no variant mentions
        # RELAXED: Also keep if section has clinical outcome keywords (affected, penetrance)
        # because aggregate penetrance data ("50 carriers, 30 affected") is valuable
        # even without explicit variant nomenclature
        has_clinical_outcome_keywords = (
            clinical_keywords >= 2
        )  # At least 2 clinical terms
        if aggregate_signals > individual_signals * 2 and variant_mentions == 0:
            if not has_clinical_outcome_keywords:
                keep = False
            # If clinical keywords present, keep the section but lower score
            else:
                score *= 0.8  # Slight penalty but don't discard

        # Require SOME signal: gene OR variant OR individual mentions OR clinical keywords
        # (relaxed from requiring all three - clinical outcome sections are valuable)
        if individual_signals == 0 and gene_mentions == 0 and variant_mentions == 0:
            if clinical_keywords < 2:
                keep = False

        # CRITICAL: Tables with variant mentions should ALWAYS be kept
        # (they likely contain the data we're looking for)
        if candidate.zone_type == "TABLE" and variant_mentions >= 2:
            keep = True
            score = max(score, 0.5)  # Ensure minimum score for variant tables

        # Methods section is usually not useful
        if candidate.source_section and "methods" in candidate.source_section.lower():
            if individual_signals < 2 and variant_mentions < 3:  # Keep if has variants
                keep = False
                score *= 0.5

        # Discussion section usually summarizes (but keep if has variant data)
        if (
            candidate.source_section
            and "discussion" in candidate.source_section.lower()
        ):
            if individual_signals < 3 and variant_mentions < 2:
                keep = False
                score *= 0.7

        # Generate zone ID
        zone_id = self._generate_zone_id(candidate, index)

        # Generate content description
        content = self._generate_content_description(
            candidate, gene_mentions, variant_mentions
        )

        # Create snippet
        snippet = text[:200].replace("\n", " ").strip()
        if len(text) > 200:
            snippet += "..."

        return DataZone(
            id=zone_id,
            type=candidate.zone_type,
            keep=keep,
            content=content,
            text_snippet=snippet,
            char_start=candidate.start,
            char_end=candidate.end,
            relevance_score=round(score, 2),
            gene_mentions=gene_mentions,
            variant_mentions=variant_mentions,
            clinical_keywords=clinical_keywords,
            source_section=candidate.source_section,
        )

    def _find_parent_section(self, text: str, position: int) -> Optional[str]:
        """Find the parent section header for a given position.

        Only a header ending within ``_PARENT_SECTION_DISTANCE_CHARS`` of the
        position counts, so searching just that window (plus slack for the
        header match's own length) is semantics-preserving. Scanning the whole
        ``text[:position]`` prefix here — 8 patterns, once per candidate —
        was the dominant cost that wedged the scout on pathological folds.
        """
        window_start = max(0, position - _PARENT_SECTION_WINDOW_CHARS)
        search_text = text[window_start:position]

        for section_name, pattern in self._section_patterns.items():
            matches = list(pattern.finditer(search_text))
            if matches:
                # Get the last (nearest) match
                last_match = matches[-1]
                # Only consider if it's within reasonable distance (10000 chars)
                if position - (window_start + last_match.end()) < (
                    _PARENT_SECTION_DISTANCE_CHARS
                ):
                    return section_name

        return None

    def _generate_zone_id(self, candidate: ZoneCandidate, index: int) -> str:
        """Generate a unique ID for the zone."""
        if candidate.zone_type == "TABLE":
            # Try to extract table number from caption
            if candidate.table_caption:
                match = re.search(
                    r"(?:supplementary\s*)?table\s*(\w+)",
                    candidate.table_caption,
                    re.IGNORECASE,
                )
                if match:
                    table_num = match.group(1).lower()
                    if (
                        "supplementary" in candidate.table_caption.lower()
                        or "s" in table_num
                    ):
                        return f"supp_table_{table_num}"
                    return f"table_{table_num}"
            return f"table_{index + 1}"
        else:
            section = candidate.source_section or "text"
            return f"{section}_para_{index + 1}"

    def _generate_content_description(
        self, candidate: ZoneCandidate, gene_mentions: int, variant_mentions: int
    ) -> str:
        """Generate a brief description of the zone content."""
        parts = []

        if candidate.zone_type == "TABLE":
            if candidate.table_caption:
                parts.append(candidate.table_caption[:100])
            else:
                parts.append("Table")
        else:
            if candidate.source_section:
                parts.append(candidate.source_section.replace("_", " ").title())
            else:
                parts.append("Text section")

        details = []
        if gene_mentions > 0:
            details.append(f"{gene_mentions} {self.gene_symbol} mentions")
        if variant_mentions > 0:
            details.append(f"{variant_mentions} variant notations")

        if details:
            parts.append(f"({', '.join(details)})")

        return " ".join(parts)

    def to_json(self, report: DataZoneReport, indent: int = 2) -> str:
        """
        Convert DataZoneReport to JSON string.

        Returns the simplified JSON format as specified:
        [
          { "id": "...", "type": "TABLE", "keep": true, "text_snippet": "..." },
          ...
        ]
        """
        simplified = []
        for zone in report.zones:
            simplified.append(
                {
                    "id": zone.id,
                    "type": zone.type,
                    "keep": zone.keep,
                    "content": zone.content,
                    "text_snippet": zone.text_snippet,
                    "relevance_score": zone.relevance_score,
                }
            )
        return json.dumps(simplified, indent=indent)

    def to_full_json(self, report: DataZoneReport, indent: int = 2) -> str:
        """Convert DataZoneReport to full JSON with all metadata."""
        return report.model_dump_json(indent=indent)

    def format_markdown(self, report: DataZoneReport, full_text: str) -> str:
        """
        Generate condensed markdown with only high-value zones.

        Args:
            report: DataZoneReport from scan()
            full_text: Original full text

        Returns:
            Condensed markdown string with zone headers and content
        """
        lines = [
            f"# Data Zones for PMID {report.pmid or 'Unknown'}",
            f"**Gene:** {report.gene_symbol}",
            f"**Zones identified:** {report.total_zones_found} ({report.zones_kept} kept, {report.zones_discarded} discarded)",
            f"**Compression:** {report.compression_ratio:.1%} of original text",
            "",
            "---",
            "",
        ]

        kept_zones = report.kept_zones
        if not kept_zones:
            lines.append("*No high-value data zones identified.*")
            return "\n".join(lines)

        for zone in kept_zones:
            zone_text = full_text[zone.char_start : zone.char_end]

            lines.extend(
                [
                    f"## [{zone.type}] {zone.id}",
                    f"**Source:** {zone.content}",
                    f"**Relevance:** {zone.relevance_score:.2f}",
                    "",
                    zone_text.strip(),
                    "",
                    "---",
                    "",
                ]
            )

        return "\n".join(lines)
