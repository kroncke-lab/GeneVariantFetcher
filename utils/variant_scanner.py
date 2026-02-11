"""
Comprehensive Variant Scanner for GeneVariantFetcher

Scans ALL text (not just markdown tables) for genetic variants using
regex-based pattern matching. Designed to catch variants that appear in:
- Narrative text ("the R534C mutation was found...")
- Figure captions
- Methods sections
- Non-table data dumps
- Any other text context

This scanner supplements LLM extraction by:
1. Pre-extracting variant hints to pass to the LLM prompt
2. Adding low-confidence scanner results to extracted variants
3. Ensuring high recall for variant detection

CREATED: 2026-02-10
"""

import logging
import re
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)

# Import from variant_normalizer for AA codes and normalization
try:
    from utils.variant_normalizer import (
        AA_MAP,
        AA_MAP_REVERSE,
        PROTEIN_LENGTHS,
        VariantNormalizer,
        get_variant_type,
        normalize_variant,
    )
except ImportError:
    # Fallback for standalone testing
    from variant_normalizer import (
        AA_MAP,
        AA_MAP_REVERSE,
        PROTEIN_LENGTHS,
        VariantNormalizer,
        get_variant_type,
        normalize_variant,
    )


@dataclass
class ScannedVariant:
    """A variant found by the scanner."""

    raw_text: str  # Original matched text
    normalized: str  # Normalized form (e.g., A561V)
    variant_type: str  # missense, frameshift, nonsense, etc.
    notation_type: str  # protein, cdna, splice
    position: Optional[int]  # Amino acid or nucleotide position
    context: str  # Surrounding text (for debugging)
    confidence: float  # 0.0-1.0 based on pattern quality
    source: str  # Where it was found (narrative, table, etc.)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "raw_text": self.raw_text,
            "normalized": self.normalized,
            "variant_type": self.variant_type,
            "notation_type": self.notation_type,
            "position": self.position,
            "confidence": self.confidence,
            "source": self.source,
        }


@dataclass
class ScanResult:
    """Result of scanning a document."""

    variants: List[ScannedVariant] = field(default_factory=list)
    unique_normalized: Set[str] = field(default_factory=set)
    stats: Dict[str, int] = field(default_factory=dict)

    def get_hints_for_prompt(self, max_hints: int = 50) -> str:
        """Format variants as hints for LLM prompt."""
        if not self.variants:
            return ""

        # Deduplicate and sort by confidence
        seen = set()
        unique_variants = []
        for v in sorted(self.variants, key=lambda x: -x.confidence):
            if v.normalized not in seen:
                seen.add(v.normalized)
                unique_variants.append(v)

        unique_variants = unique_variants[:max_hints]

        lines = [
            "\n\n--- PRE-SCANNED VARIANT HINTS (regex detection) ---",
            f"Pattern matching found {len(unique_variants)} potential {self.stats.get('gene', 'target gene')} variants.",
            "Verify these and extract full clinical details:",
            "",
        ]

        for i, v in enumerate(unique_variants, 1):
            conf_label = (
                "HIGH"
                if v.confidence >= 0.8
                else "MED"
                if v.confidence >= 0.5
                else "LOW"
            )
            lines.append(f"  {i}. {v.normalized} [{v.variant_type}] ({conf_label})")

        lines.append("\n--- END PRE-SCANNED HINTS ---\n")
        return "\n".join(lines)

    def to_variant_dicts(self, gene_symbol: str) -> List[Dict[str, Any]]:
        """Convert to variant dicts compatible with extraction pipeline."""
        seen = set()
        results = []

        for v in sorted(self.variants, key=lambda x: -x.confidence):
            if v.normalized in seen:
                continue
            seen.add(v.normalized)

            variant_dict = {
                "gene_symbol": gene_symbol,
                "protein_notation": v.normalized
                if v.notation_type == "protein"
                else None,
                "cdna_notation": v.normalized if v.notation_type == "cdna" else None,
                "clinical_significance": "unknown",
                "evidence_level": "scanner",
                "source_location": f"Text scan ({v.source})",
                "additional_notes": f"Auto-detected via pattern scanning (confidence: {v.confidence:.2f})",
                "patients": {},
                "penetrance_data": {},
                "individual_records": [],
                "functional_data": {"summary": "", "assays": []},
                "key_quotes": [],
                "_scanner_confidence": v.confidence,
                "_scanner_raw": v.raw_text,
            }
            results.append(variant_dict)

        return results


class VariantScanner:
    """
    Comprehensive regex-based variant scanner.

    Finds variants in all text, not just structured tables.
    """

    # ==========================================================================
    # PROTEIN VARIANT PATTERNS
    # ==========================================================================

    # Full HGVS protein notation: p.Arg534Cys, p.Ala561Val, p.Leu987fs, etc.
    PROTEIN_HGVS_FULL = re.compile(
        r"\bp\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|fs\*?\d*|del|dup|ins|Ter|\*)"
        r"(?:\*?\d*)?",  # Optional frameshift extension
        re.IGNORECASE,
    )

    # Short HGVS with p.: p.R534C, p.A561V, p.L987fs
    PROTEIN_HGVS_SHORT = re.compile(
        r"\bp\.([A-Z])(\d+)([A-Z]|fs[X\*]?\d*|del|dup|ins|\*|X)", re.IGNORECASE
    )

    # Three-letter AA without p. prefix: Arg534Cys, Ala561Val, Leu987fs
    PROTEIN_THREE_LETTER = re.compile(
        r"\b([A-Z][a-z]{2})(\d{2,4})([A-Z][a-z]{2}|fs\*?\d*|del|dup|ins|Ter|\*)"
        r"(?:[X\*]?\d*)?",
        re.IGNORECASE,
    )

    # Single-letter variants: R534C, A561V, L987fs, W1001X
    # Must have 2-4 digit position to avoid false positives
    PROTEIN_SINGLE_LETTER = re.compile(
        r"\b([A-Z])(\d{2,4})([A-Z]|fs[X\*]?\d*|del|dup|ins|\*|X)\b", re.IGNORECASE
    )

    # Concatenated gene+variant: HERGG604S, KCNH2A561T, hERGT613M, Kv11.1R534C
    # Gene name prefix is consumed but not captured as part of the variant
    CONCATENATED_GENE_VARIANT = re.compile(
        r"\b(?:HERG|hERG|KCNH2|kcnh2|Kv11\.1)"
        r"[-_]?"  # Optional separator
        r"([A-Z])"  # Ref AA
        r"(\d{2,4})"  # Position
        r"([A-Z]|fs[X\*]?\d*|del|dup|ins|\*|X)"  # Alt AA or special
        r"\b",
        re.IGNORECASE,
    )

    # Frameshift variants with various notations
    FRAMESHIFT_PATTERNS = re.compile(
        r"\b(?:p\.)?"
        r"([A-Z][a-z]{2}|[A-Z])"  # Ref AA (3-letter or 1-letter)
        r"(\d{2,4})"  # Position
        r"(?:[A-Z][a-z]{2})?"  # Optional second AA (for Profs type)
        r"(fs"  # Start of frameshift
        r"(?:[X\*]?\d*|Ter\d*)?)"  # Optional extension: fsX, fs*10, fsTer10
        r"\b",
        re.IGNORECASE,
    )

    # Nonsense/stop variants: W1001X, W1001*, p.Trp1001Ter, R864stop
    NONSENSE_PATTERNS = re.compile(
        r"\b(?:p\.)?"
        r"([A-Z][a-z]{2}|[A-Z])"  # Ref AA
        r"(\d{2,4})"  # Position
        r"(\*|X|Ter|stop|sp)"  # Stop indicator
        r"\b",
        re.IGNORECASE,
    )

    # Deletion variants: L552del, p.Leu552del, del552, del I642-V644
    DELETION_PATTERNS = re.compile(
        r"\b(?:p\.)?"
        r"(?:([A-Z][a-z]{2}|[A-Z])(\d{2,4})del"  # Leu552del or L552del
        r"|del\s*([A-Z])?(\d{2,4})"  # del552 or del L552
        r"(?:[_\-]([A-Z])?(\d{2,4}))?)"  # Optional range: del I642-V644
        r"\b",
        re.IGNORECASE,
    )

    # Duplication variants
    DUPLICATION_PATTERNS = re.compile(
        r"\b(?:p\.)?"
        r"([A-Z][a-z]{2}|[A-Z])?"
        r"(\d{2,4})"
        r"(?:_([A-Z][a-z]{2}|[A-Z])?(\d{2,4}))?"
        r"dup"
        r"\b",
        re.IGNORECASE,
    )

    # Insertion variants
    INSERTION_PATTERNS = re.compile(
        r"\b(?:p\.)?"
        r"([A-Z][a-z]{2}|[A-Z])?"
        r"(\d{2,4})"
        r"(?:_([A-Z][a-z]{2}|[A-Z])?(\d{2,4}))?"
        r"ins"
        r"([A-Z][a-z]{2}|[A-Z])*"  # Inserted residues
        r"\b",
        re.IGNORECASE,
    )

    # ==========================================================================
    # cDNA VARIANT PATTERNS
    # ==========================================================================

    # Standard cDNA substitution: c.1600C>T, c.526C>T
    # Note: using case-sensitive for nucleotides to reduce false positives
    CDNA_SUBSTITUTION = re.compile(r"c\.(\d+)([ACGTacgt])>([ACGTacgt])")

    # cDNA with position only (no c. prefix): 1600C>T
    CDNA_NO_PREFIX = re.compile(r"\b(\d{3,5})([ACGT])>([ACGT])\b", re.IGNORECASE)

    # Intronic cDNA: c.1234+1G>A, c.1234-2A>G
    CDNA_INTRONIC = re.compile(
        r"\bc\.(\d+)([\+\-]\d+)([ACGT])>([ACGT])\b", re.IGNORECASE
    )

    # cDNA deletion/duplication: c.1234del, c.1234_1235delAG, c.1234dup
    CDNA_INDEL = re.compile(
        r"\bc\.(\d+)(?:_(\d+))?(del|dup|ins)([ACGT]*)\b", re.IGNORECASE
    )

    # ==========================================================================
    # SPLICE VARIANT PATTERNS (IVS notation)
    # ==========================================================================

    # IVS notation: IVS9+1G>A, IVS5-2A>G
    IVS_SPLICE = re.compile(r"\bIVS(\d+)([\+\-]\d+)([ACGT])>([ACGT])\b", re.IGNORECASE)

    # IVS indel: IVS10+1del
    IVS_INDEL = re.compile(
        r"\bIVS(\d+)([\+\-]\d+)(del|dup|ins)([ACGT]*)\b", re.IGNORECASE
    )

    # ==========================================================================
    # NARRATIVE CONTEXT PATTERNS
    # ==========================================================================

    # Patterns that indicate a variant mention in narrative text
    NARRATIVE_CONTEXTS = [
        # "the R534C mutation", "a R534C variant"
        re.compile(
            r"\b(?:the|a|an)\s+([A-Z]\d{2,4}[A-Z])\s+(?:mutation|variant|substitution)",
            re.IGNORECASE,
        ),
        # "carrying the p.Arg534Cys variant"
        re.compile(
            r"carrying\s+(?:the\s+)?(\bp\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2})\b",
            re.IGNORECASE,
        ),
        # "mutation R534C was found"
        re.compile(r"\bmutation\s+([A-Z]\d{2,4}[A-Z])\s+(?:was|is|has)", re.IGNORECASE),
        # "identified the R534C"
        re.compile(r"\bidentified\s+(?:the\s+)?([A-Z]\d{2,4}[A-Z])\b", re.IGNORECASE),
        # "c.1600C>T mutation"
        re.compile(
            r"\b(c\.\d+[ACGT]>[ACGT])\s+(?:mutation|variant|change)", re.IGNORECASE
        ),
        # "mutation at position 534"
        re.compile(
            r"\bmutation\s+(?:at\s+)?(?:position\s+)?(\d{2,4})\b", re.IGNORECASE
        ),
    ]

    # ==========================================================================
    # FALSE POSITIVE FILTERS
    # ==========================================================================

    # Patterns that look like variants but aren't
    FALSE_POSITIVE_PATTERNS = [
        re.compile(r"^[STFR]\d{1}$"),  # S1, T2, F1 (figure/table refs)
        re.compile(r"^Table\s*\d", re.I),  # Table 1, Table 2
        re.compile(r"^Fig", re.I),  # Figure refs
        re.compile(r"^[A-Z]\d{1,2}$"),  # Short codes like A1, B12
        re.compile(r"^\d+\s*°?[CF]$"),  # Temperatures: 37C, 37°C
        re.compile(r"^\d+\s*[mkμn]?[gGlLmM]"),  # Concentrations: 5mg, 10mL
        re.compile(r"^IC\d+$"),  # IC50 etc.
        re.compile(r"^EC\d+$"),  # EC50 etc.
        re.compile(r"^LD\d+$"),  # LD50 etc.
        re.compile(r"^K\d+$"),  # K+ channel nomenclature
        re.compile(r"^HEK\d+"),  # Cell lines
        re.compile(r"^CHO\s"),  # Cell lines
        re.compile(r"^\d+[xX]$"),  # Magnification: 10x, 100X
    ]

    # ==========================================================================
    # NON-TARGET GENE HOTSPOTS (from variant_normalizer)
    # ==========================================================================

    NON_TARGET_HOTSPOTS = {
        # TP53 (most common cancer gene mutations)
        "R175H",
        "R248W",
        "R248Q",
        "R249S",
        "R273H",
        "R273C",
        "R282W",
        "G245S",
        "G245D",
        "Y220C",
        "C176F",
        "C176Y",
        "C242F",
        "C242S",
        "H179R",
        "H179Y",
        "M237I",
        "S241F",
        "S241C",
        "C277F",
        "Y163C",
        # KRAS
        "G12D",
        "G12V",
        "G12C",
        "G12R",
        "G12A",
        "G12S",
        "G13D",
        "Q61H",
        "Q61L",
        "Q61R",
        "Q61G",
        # BRAF
        "V600E",
        "V600K",
        # PIK3CA
        "E545K",
        "H1047R",
        "H1047L",
    }

    def __init__(self, gene_symbol: str = "KCNH2"):
        """
        Initialize the variant scanner.

        Args:
            gene_symbol: Target gene for position validation
        """
        self.gene_symbol = gene_symbol.upper()
        self.protein_length = PROTEIN_LENGTHS.get(self.gene_symbol, 9999)
        self.normalizer = VariantNormalizer(self.gene_symbol)

    def scan(self, text: str, source: str = "full_text") -> ScanResult:
        """
        Scan text for all variants.

        Args:
            text: Full text to scan
            source: Source identifier (for logging)

        Returns:
            ScanResult with all found variants
        """
        result = ScanResult()
        result.stats["gene"] = self.gene_symbol
        result.stats["text_length"] = len(text)
        result.stats["source"] = source

        # Pre-process: normalize Unicode arrows (→ to >) so cDNA patterns match
        # Papers often use Unicode arrows: "1810G→A" instead of "1810G>A"
        text = text.replace("\u2192", ">").replace("\u2190", "<").replace("\u21d2", ">")

        # Track what we've found to avoid duplicates
        seen_normalized: Set[str] = set()

        # Run all pattern matchers
        protein_variants = self._scan_protein_variants(text)
        cdna_variants = self._scan_cdna_variants(text)
        splice_variants = self._scan_splice_variants(text)
        narrative_variants = self._scan_narrative_variants(text)

        # Combine all findings
        all_candidates = (
            protein_variants + cdna_variants + splice_variants + narrative_variants
        )

        # Filter and deduplicate
        filtered_count = 0
        for v in all_candidates:
            # Skip false positives
            if self._is_false_positive(v.raw_text):
                filtered_count += 1
                continue

            # Skip non-target hotspots
            if v.normalized.upper() in self.NON_TARGET_HOTSPOTS:
                filtered_count += 1
                logger.debug(f"Scanner: filtered non-target hotspot {v.raw_text}")
                continue

            # Skip invalid positions (only for protein variants)
            # cDNA positions are nucleotide positions which can be much larger
            if (
                v.notation_type == "protein"
                and v.position
                and v.position > self.protein_length
            ):
                filtered_count += 1
                logger.debug(
                    f"Scanner: filtered invalid position {v.raw_text} (pos={v.position})"
                )
                continue

            # Skip duplicates
            if v.normalized in seen_normalized:
                continue

            seen_normalized.add(v.normalized)
            result.variants.append(v)
            result.unique_normalized.add(v.normalized)

        # Update stats
        result.stats["total_candidates"] = len(all_candidates)
        result.stats["filtered"] = filtered_count
        result.stats["unique_variants"] = len(result.unique_normalized)
        result.stats["protein_variants"] = sum(
            1 for v in result.variants if v.notation_type == "protein"
        )
        result.stats["cdna_variants"] = sum(
            1 for v in result.variants if v.notation_type == "cdna"
        )
        result.stats["splice_variants"] = sum(
            1 for v in result.variants if v.notation_type == "splice"
        )

        logger.info(
            f"Variant scanner found {len(result.variants)} unique variants in {source} "
            f"(candidates: {len(all_candidates)}, filtered: {filtered_count})"
        )

        return result

    def _scan_protein_variants(self, text: str) -> List[ScannedVariant]:
        """Scan for protein variants."""
        variants = []

        # Full HGVS: p.Arg534Cys
        for m in self.PROTEIN_HGVS_FULL.finditer(text):
            ref, pos, alt = m.group(1), m.group(2), m.group(3)
            raw = m.group(0)
            position = int(pos)

            # Normalize
            normalized = self.normalizer.normalize_to_single_letter(raw)
            if not normalized:
                normalized = f"{ref[0].upper()}{pos}{alt[0].upper() if alt else 'X'}"

            variants.append(
                ScannedVariant(
                    raw_text=raw,
                    normalized=normalized,
                    variant_type=self._classify_variant_type(alt),
                    notation_type="protein",
                    position=position,
                    context=self._get_context(text, m.start(), m.end()),
                    confidence=0.95,  # High confidence for full HGVS
                    source="protein_hgvs_full",
                )
            )

        # Short HGVS: p.R534C
        for m in self.PROTEIN_HGVS_SHORT.finditer(text):
            ref, pos, alt = m.group(1), m.group(2), m.group(3)
            raw = m.group(0)
            position = int(pos)

            normalized = self.normalizer.normalize_to_single_letter(raw)
            if not normalized:
                normalized = f"{ref.upper()}{pos}{alt.upper()}"

            variants.append(
                ScannedVariant(
                    raw_text=raw,
                    normalized=normalized,
                    variant_type=self._classify_variant_type(alt),
                    notation_type="protein",
                    position=position,
                    context=self._get_context(text, m.start(), m.end()),
                    confidence=0.90,
                    source="protein_hgvs_short",
                )
            )

        # Three-letter without prefix: Arg534Cys
        for m in self.PROTEIN_THREE_LETTER.finditer(text):
            ref, pos, alt = m.group(1), m.group(2), m.group(3)
            raw = m.group(0)

            # Must be valid amino acids
            ref_single = AA_MAP_REVERSE.get(ref.capitalize())
            if not ref_single:
                continue

            position = int(pos)
            normalized = self.normalizer.normalize_to_single_letter(raw)
            if not normalized:
                alt_single = AA_MAP_REVERSE.get(
                    alt.capitalize(), alt[0].upper() if alt else "X"
                )
                normalized = f"{ref_single}{pos}{alt_single}"

            variants.append(
                ScannedVariant(
                    raw_text=raw,
                    normalized=normalized,
                    variant_type=self._classify_variant_type(alt),
                    notation_type="protein",
                    position=position,
                    context=self._get_context(text, m.start(), m.end()),
                    confidence=0.85,
                    source="protein_three_letter",
                )
            )

        # Single-letter: R534C, A561V
        for m in self.PROTEIN_SINGLE_LETTER.finditer(text):
            ref, pos, alt = m.group(1), m.group(2), m.group(3)
            raw = m.group(0)
            position = int(pos)

            # Basic validation
            if ref.upper() not in AA_MAP:
                continue
            if len(alt) == 1 and alt.upper() not in AA_MAP and alt.upper() not in "X*":
                continue

            normalized = f"{ref.upper()}{pos}{alt.upper()}"

            # Lower confidence for standalone single-letter (more false positives)
            confidence = 0.70
            # But higher if in a variant-like context
            context = self._get_context(text, m.start(), m.end())
            if any(
                word in context.lower()
                for word in ["mutation", "variant", "substitution", "carrier"]
            ):
                confidence = 0.85

            variants.append(
                ScannedVariant(
                    raw_text=raw,
                    normalized=normalized,
                    variant_type=self._classify_variant_type(alt),
                    notation_type="protein",
                    position=position,
                    context=context,
                    confidence=confidence,
                    source="protein_single_letter",
                )
            )

        # Frameshift patterns
        for m in self.FRAMESHIFT_PATTERNS.finditer(text):
            ref, pos, fs_suffix = m.group(1), m.group(2), m.group(3)
            raw = m.group(0)
            position = int(pos)

            # Convert to single letter if needed
            if len(ref) == 3:
                ref_single = AA_MAP_REVERSE.get(ref.capitalize(), ref[0].upper())
            else:
                ref_single = ref.upper()

            normalized = f"{ref_single}{pos}fsX"

            variants.append(
                ScannedVariant(
                    raw_text=raw,
                    normalized=normalized,
                    variant_type="frameshift",
                    notation_type="protein",
                    position=position,
                    context=self._get_context(text, m.start(), m.end()),
                    confidence=0.90,
                    source="frameshift",
                )
            )

        # Nonsense patterns
        for m in self.NONSENSE_PATTERNS.finditer(text):
            ref, pos, stop = m.group(1), m.group(2), m.group(3)
            raw = m.group(0)
            position = int(pos)

            if len(ref) == 3:
                ref_single = AA_MAP_REVERSE.get(ref.capitalize(), ref[0].upper())
            else:
                ref_single = ref.upper()

            normalized = f"{ref_single}{pos}X"

            variants.append(
                ScannedVariant(
                    raw_text=raw,
                    normalized=normalized,
                    variant_type="nonsense",
                    notation_type="protein",
                    position=position,
                    context=self._get_context(text, m.start(), m.end()),
                    confidence=0.90,
                    source="nonsense",
                )
            )

        # Deletion patterns
        for m in self.DELETION_PATTERNS.finditer(text):
            raw = m.group(0)
            groups = m.groups()

            # Parse different deletion formats
            if groups[0]:  # Leu552del or L552del
                ref, pos = groups[0], groups[1]
            elif groups[3]:  # del552 or del L552
                ref = groups[2] or "?"
                pos = groups[3]
            else:
                continue

            position = int(pos)

            if len(ref) == 3:
                ref_single = AA_MAP_REVERSE.get(ref.capitalize(), ref[0].upper())
            elif len(ref) == 1:
                ref_single = ref.upper()
            else:
                ref_single = "?"

            normalized = f"{ref_single}{pos}del"

            variants.append(
                ScannedVariant(
                    raw_text=raw,
                    normalized=normalized,
                    variant_type="deletion",
                    notation_type="protein",
                    position=position,
                    context=self._get_context(text, m.start(), m.end()),
                    confidence=0.85,
                    source="deletion",
                )
            )

        # Concatenated gene+variant: HERGG604S, KCNH2A561T, hERGT613M
        for m in self.CONCATENATED_GENE_VARIANT.finditer(text):
            ref, pos, alt = m.group(1), m.group(2), m.group(3)
            raw = m.group(0)
            position = int(pos)

            normalized = f"{ref.upper()}{pos}{alt.upper()}"

            variants.append(
                ScannedVariant(
                    raw_text=raw,
                    normalized=normalized,
                    variant_type=self._classify_variant_type(alt),
                    notation_type="protein",
                    position=position,
                    context=self._get_context(text, m.start(), m.end()),
                    confidence=0.90,  # High confidence - gene name prefix confirms target gene
                    source="concatenated_gene_variant",
                )
            )

        return variants

    def _scan_cdna_variants(self, text: str) -> List[ScannedVariant]:
        """Scan for cDNA variants."""
        variants = []

        # Standard c.1600C>T
        for m in self.CDNA_SUBSTITUTION.finditer(text):
            pos, ref, alt = m.group(1), m.group(2), m.group(3)
            raw = m.group(0)

            normalized = f"c.{pos}{ref.upper()}>{alt.upper()}"

            variants.append(
                ScannedVariant(
                    raw_text=raw,
                    normalized=normalized,
                    variant_type="substitution",
                    notation_type="cdna",
                    position=int(pos),
                    context=self._get_context(text, m.start(), m.end()),
                    confidence=0.95,
                    source="cdna_substitution",
                )
            )

        # Intronic: c.1234+1G>A
        for m in self.CDNA_INTRONIC.finditer(text):
            pos, offset, ref, alt = m.group(1), m.group(2), m.group(3), m.group(4)
            raw = m.group(0)

            normalized = f"c.{pos}{offset}{ref.upper()}>{alt.upper()}"

            variants.append(
                ScannedVariant(
                    raw_text=raw,
                    normalized=normalized,
                    variant_type="splice",
                    notation_type="cdna",
                    position=int(pos),
                    context=self._get_context(text, m.start(), m.end()),
                    confidence=0.95,
                    source="cdna_intronic",
                )
            )

        # cDNA indels
        for m in self.CDNA_INDEL.finditer(text):
            pos1 = m.group(1)
            pos2 = m.group(2)  # May be None
            indel_type = m.group(3).lower()
            bases = m.group(4) or ""

            raw = m.group(0)

            if pos2:
                normalized = f"c.{pos1}_{pos2}{indel_type}{bases.upper()}"
            else:
                normalized = f"c.{pos1}{indel_type}{bases.upper()}"

            variants.append(
                ScannedVariant(
                    raw_text=raw,
                    normalized=normalized,
                    variant_type=indel_type,
                    notation_type="cdna",
                    position=int(pos1),
                    context=self._get_context(text, m.start(), m.end()),
                    confidence=0.90,
                    source="cdna_indel",
                )
            )

        return variants

    def _scan_splice_variants(self, text: str) -> List[ScannedVariant]:
        """Scan for splice/IVS variants."""
        variants = []

        # IVS notation
        for m in self.IVS_SPLICE.finditer(text):
            ivs_num, offset, ref, alt = m.group(1), m.group(2), m.group(3), m.group(4)
            raw = m.group(0)

            normalized = f"IVS{ivs_num}{offset}{ref.upper()}>{alt.upper()}"

            variants.append(
                ScannedVariant(
                    raw_text=raw,
                    normalized=normalized,
                    variant_type="splice",
                    notation_type="splice",
                    position=int(ivs_num),  # IVS number as position
                    context=self._get_context(text, m.start(), m.end()),
                    confidence=0.95,
                    source="ivs_splice",
                )
            )

        # IVS indels
        for m in self.IVS_INDEL.finditer(text):
            ivs_num, offset, indel_type, bases = (
                m.group(1),
                m.group(2),
                m.group(3),
                m.group(4) or "",
            )
            raw = m.group(0)

            normalized = f"IVS{ivs_num}{offset}{indel_type.lower()}{bases.upper()}"

            variants.append(
                ScannedVariant(
                    raw_text=raw,
                    normalized=normalized,
                    variant_type="splice",
                    notation_type="splice",
                    position=int(ivs_num),
                    context=self._get_context(text, m.start(), m.end()),
                    confidence=0.90,
                    source="ivs_indel",
                )
            )

        return variants

    def _scan_narrative_variants(self, text: str) -> List[ScannedVariant]:
        """Scan for variants mentioned in narrative context."""
        variants = []

        for pattern in self.NARRATIVE_CONTEXTS:
            for m in pattern.finditer(text):
                raw = m.group(1)

                # Try to normalize
                normalized = normalize_variant(raw, self.gene_symbol)
                if not normalized or normalized == raw.upper():
                    # Couldn't normalize well, try as-is
                    normalized = raw.upper()

                # Extract position
                pos_match = re.search(r"(\d{2,4})", raw)
                position = int(pos_match.group(1)) if pos_match else None

                # Determine notation type
                notation_type = "protein"
                if raw.lower().startswith("c."):
                    notation_type = "cdna"
                elif raw.lower().startswith("ivs"):
                    notation_type = "splice"

                variants.append(
                    ScannedVariant(
                        raw_text=raw,
                        normalized=normalized,
                        variant_type=get_variant_type(normalized),
                        notation_type=notation_type,
                        position=position,
                        context=self._get_context(text, m.start(), m.end(), window=100),
                        confidence=0.75,  # Lower confidence for narrative mentions
                        source="narrative",
                    )
                )

        return variants

    def _classify_variant_type(self, suffix: str) -> str:
        """Classify variant type from the suffix."""
        if not suffix:
            return "unknown"

        suffix_lower = suffix.lower()

        if "fs" in suffix_lower:
            return "frameshift"
        if suffix_lower in ("*", "x", "ter", "stop", "sp"):
            return "nonsense"
        if "del" in suffix_lower:
            return "deletion"
        if "dup" in suffix_lower:
            return "duplication"
        if "ins" in suffix_lower:
            return "insertion"
        if len(suffix) == 1 or len(suffix) == 3:
            return "missense"

        return "unknown"

    def _is_false_positive(self, text: str) -> bool:
        """Check if a match is a false positive."""
        for pattern in self.FALSE_POSITIVE_PATTERNS:
            if pattern.match(text):
                return True

        # Additional heuristics
        if len(text) < 3:
            return True

        return False

    def _get_context(self, text: str, start: int, end: int, window: int = 50) -> str:
        """Get surrounding context for a match."""
        ctx_start = max(0, start - window)
        ctx_end = min(len(text), end + window)
        return text[ctx_start:ctx_end]


def scan_document_for_variants(
    text: str, gene_symbol: str = "KCNH2", source: str = "full_text"
) -> ScanResult:
    """
    Convenience function to scan a document for variants.

    Args:
        text: Document text to scan
        gene_symbol: Target gene for position validation
        source: Source identifier

    Returns:
        ScanResult with found variants
    """
    scanner = VariantScanner(gene_symbol)
    return scanner.scan(text, source)


def merge_scanner_results(
    extracted_data: Dict[str, Any],
    scan_result: ScanResult,
    gene_symbol: str,
    min_confidence: float = 0.5,
) -> Dict[str, Any]:
    """
    Merge scanner-found variants with LLM-extracted variants.

    Only adds variants that weren't already extracted by the LLM.
    Scanner variants are added with lower evidence level.

    Args:
        extracted_data: LLM extraction result dict
        scan_result: Scanner result
        gene_symbol: Target gene symbol
        min_confidence: Minimum scanner confidence to include

    Returns:
        Updated extracted_data with merged variants
    """
    existing_variants = extracted_data.get("variants", [])

    # Build set of existing variant keys (normalized)
    existing_keys = set()
    for v in existing_variants:
        protein = normalize_variant(v.get("protein_notation", "") or "", gene_symbol)
        cdna = normalize_variant(v.get("cdna_notation", "") or "", gene_symbol)
        if protein:
            existing_keys.add(protein.upper())
        if cdna:
            existing_keys.add(cdna.upper())

    # Add scanner variants not already present
    added_count = 0
    for sv in scan_result.variants:
        if sv.confidence < min_confidence:
            continue

        if sv.normalized.upper() in existing_keys:
            continue

        # Create variant dict
        new_variant = {
            "gene_symbol": gene_symbol,
            "protein_notation": sv.normalized
            if sv.notation_type == "protein"
            else None,
            "cdna_notation": sv.normalized
            if sv.notation_type in ("cdna", "splice")
            else None,
            "clinical_significance": "unknown",
            "evidence_level": "scanner",
            "source_location": f"Text scan ({sv.source})",
            "additional_notes": f"Auto-detected via pattern scanning (confidence: {sv.confidence:.2f}, raw: '{sv.raw_text}')",
            "patients": {},
            "penetrance_data": {},
            "individual_records": [],
            "functional_data": {"summary": "", "assays": []},
            "key_quotes": [],
        }

        existing_variants.append(new_variant)
        existing_keys.add(sv.normalized.upper())
        added_count += 1

    if added_count > 0:
        logger.info(f"Scanner merge added {added_count} variants not found by LLM")

        # Update metadata
        if "extraction_metadata" in extracted_data:
            extracted_data["extraction_metadata"]["total_variants_found"] = len(
                existing_variants
            )
            extracted_data["extraction_metadata"]["scanner_added"] = added_count

    extracted_data["variants"] = existing_variants
    return extracted_data


# ==========================================================================
# TESTING
# ==========================================================================

if __name__ == "__main__":
    # Basic test with various variant formats
    test_text = """
    The patient carried the p.Arg534Cys mutation, also known as R534C. 
    This variant c.1600C>T was previously reported in LQT2 families.
    We also identified A561V, Gly628Ser, and the frameshift mutation L987fsX.
    The IVS9+1G>A splice variant was found in 3 families.
    The intronic variant c.2398+1G>A causes splicing defects.
    The W1001X nonsense mutation leads to truncation.
    In Table 2, we list the mutations: T613M, N470D, and the deletion p.Leu552del.
    Also found: c.526C>T (R176W), p.Thr613Met, and c.1234del.
    """

    result = scan_document_for_variants(test_text, "KCNH2")

    print(f"\nFound {len(result.variants)} variants:")
    for v in result.variants:
        print(
            f"  {v.normalized:15} [{v.variant_type:12}] conf={v.confidence:.2f} ({v.source})"
        )

    print(f"\nStats: {result.stats}")
    print(f"\nHints for prompt:\n{result.get_hints_for_prompt()}")
