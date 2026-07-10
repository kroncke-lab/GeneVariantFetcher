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
    from utils.gene_metadata import (
        gene_alias_regex,
        get_gene_aliases,
        known_gene_aliases,
    )
    from utils.variant_normalizer import (
        AA_MAP,
        AA_MAP_REVERSE,
        NON_TARGET_HOTSPOT_GENES,
        PROTEIN_LENGTHS,
        VariantNormalizer,
        get_variant_type,
        normalize_variant,
    )
except ImportError:
    # Fallback for standalone testing
    gene_alias_regex = None
    get_gene_aliases = None
    known_gene_aliases = None
    from variant_normalizer import (
        AA_MAP,
        AA_MAP_REVERSE,
        NON_TARGET_HOTSPOT_GENES,
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
    notation_type: str  # protein, cdna, splice, structural
    position: Optional[int]  # Amino acid or nucleotide position
    context: str  # Surrounding text (for debugging)
    confidence: float  # 0.0-1.0 based on pattern quality
    source: str  # Where it was found (narrative, table, etc.)
    variant_class: Optional[str] = None  # closed taxonomy when known
    structural_description: Optional[str] = None  # free-text structural event

    def to_dict(self) -> Dict[str, Any]:
        return {
            "raw_text": self.raw_text,
            "normalized": self.normalized,
            "variant_type": self.variant_type,
            "notation_type": self.notation_type,
            "position": self.position,
            "confidence": self.confidence,
            "source": self.source,
            "variant_class": self.variant_class,
            "structural_description": self.structural_description,
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

            context_quote = " ".join((v.context or v.raw_text or "").split())
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
                "key_quotes": [context_quote] if context_quote else [],
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

    GENE_CONTEXT_TERMS = {gene: {gene} for gene in PROTEIN_LENGTHS}
    GENE_CONTEXT_TERMS.setdefault("KCNH2", set()).update({"HERG", "hERG", "Kv11.1"})
    GENE_CONTEXT_TERMS.setdefault("KCNQ1", set()).update({"KvLQT1", "Kv7.1"})
    GENE_CONTEXT_TERMS.setdefault("SCN5A", set()).update({"Nav1.5", "Na v 1.5"})

    @classmethod
    def _known_context_genes(cls) -> set[str]:
        genes = set(PROTEIN_LENGTHS) | set(cls.GENE_CONTEXT_TERMS)
        if known_gene_aliases is not None:
            try:
                genes.update(known_gene_aliases(include_query_aliases=True))
            except Exception:
                pass
        return genes

    @classmethod
    def _gene_context_terms(cls, gene_symbol: str) -> set[str]:
        gene = gene_symbol.upper()
        terms = set(cls.GENE_CONTEXT_TERMS.get(gene, {gene_symbol}))
        if get_gene_aliases is not None:
            try:
                terms.update(get_gene_aliases(gene, include_query_aliases=True))
            except Exception:
                pass
        return terms

    # ==========================================================================
    # PROTEIN VARIANT PATTERNS
    # ==========================================================================

    # Full HGVS protein notation: p.Arg534Cys, p.Ala561Val, p.Leu987fs, etc.
    PROTEIN_HGVS_FULL = re.compile(
        r"\bp\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|fs\*?\d*|del|dup|ins|Ter|\*)"
        r"(?:\*?\d*)?",  # Optional frameshift extension
        re.IGNORECASE,
    )

    # Parenthesized HGVS protein notation: p.(Arg176Trp), p.(Ser906Leu),
    # p.(Ala1198fs). HGVS uses parentheses to mark *predicted* protein
    # consequences (vs. experimentally confirmed) and this form is pervasive
    # in supplement tables generated by VEP / Annovar / classical clinical
    # genetics pipelines. PROTEIN_HGVS_FULL above cannot match this — its
    # leading "\bp\.[A-Z][a-z]{2}" requires a letter immediately after
    # "p." and "(" breaks that. Spot-check on the KCNH2 corpus
    # (results/KCNH2/20260506_102238) found 937 occurrences across 492
    # unique variants invisible to the scanner.
    PROTEIN_HGVS_PAREN = re.compile(
        r"\bp\.\(([A-Z][a-z]{2})(\d+)"
        r"([A-Z][a-z]{2}|fs\*?\d*|del|dup|ins|Ter|\*|=)"
        r"(?:\*?\d*)?\)",
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
    RANGE_DELETION_PATTERNS = re.compile(
        r"\b(?:p\.)?"
        r"([A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWY])"
        r"(\d{2,4})"
        r"[_-]"
        r"([A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWY])"
        r"(\d{2,4})"
        r"del[A-Z]*\b",
        re.IGNORECASE,
    )
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

    # cDNA delins: c.123_456delXXXinsYYY or c.123delAinsT
    CDNA_DELINS = re.compile(
        r"\bc\.(\d+)(?:_(\d+))?del([ACGT]*)ins([ACGT]+)\b", re.IGNORECASE
    )

    # ==========================================================================
    # STRUCTURAL / EXON-LEVEL / PREFIXLESS BIC PATTERNS (tightly gated)
    # ==========================================================================

    # Prefixless BIC-style indel: 185delAG, 5382insC (need digit length ≥ 3 + bases)
    BIC_PREFIXLESS_INDEL = re.compile(
        r"\b(\d{3,5})(del|dup|ins)([ACGT]{1,20})\b", re.IGNORECASE
    )

    # Exon-level deletion/duplication: deletion of exons 3-5, del exons 3–5
    EXON_EVENT = re.compile(
        r"\b(?:(?:large\s+)?(?:deletion|duplication)\s+of\s+)?"
        r"(?:del(?:etion)?|dup(?:lication)?)\s*"
        r"(?:of\s+)?"
        r"exons?\s*"
        r"(\d{1,3})"
        r"(?:\s*[-–—to]+\s*(\d{1,3}))?"
        r"\b",
        re.IGNORECASE,
    )

    # Whole-gene / large CNV phrasing (require gene-ish context nearby is handled later)
    WHOLE_GENE_DEL = re.compile(
        r"\b(?:whole[\s-]?gene|entire\s+gene)\s+(?:deletion|del)\b"
        r"|\b(?:deletion|del)\s+of\s+(?:the\s+)?(?:entire|whole)\s+gene\b",
        re.IGNORECASE,
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
        structural_variants = self._scan_structural_variants(text)

        # Combine all findings
        all_candidates = (
            protein_variants
            + cdna_variants
            + splice_variants
            + narrative_variants
            + structural_variants
        )

        # Filter and deduplicate
        filtered_count = 0
        for v in all_candidates:
            # Skip false positives
            if self._is_false_positive(v.raw_text):
                filtered_count += 1
                continue

            if self._has_conflicting_gene_context(v.context, v.raw_text):
                filtered_count += 1
                logger.debug(
                    "Scanner: filtered %s due to conflicting gene context",
                    v.raw_text,
                )
                continue

            # Skip common comparator hotspots only when they are not the target gene.
            hotspot_genes = NON_TARGET_HOTSPOT_GENES.get(v.normalized.upper())
            if hotspot_genes and self.gene_symbol not in hotspot_genes:
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

    def _context_mentions_gene(self, context: str, gene_symbol: str) -> bool:
        if gene_alias_regex is not None:
            try:
                return bool(
                    gene_alias_regex(gene_symbol, include_query_aliases=True).search(
                        context
                    )
                )
            except Exception:
                pass
        terms = self._gene_context_terms(gene_symbol)
        for term in terms:
            pattern = rf"(?<![A-Za-z0-9]){re.escape(term)}(?![A-Za-z0-9])"
            if re.search(pattern, context, re.IGNORECASE):
                return True
        return False

    def _gene_assigned_to_variant(self, context: str, raw_text: str) -> Optional[str]:
        """Infer the closest gene label attached to this variant mention."""
        if not context or not raw_text:
            return None

        raw_match = re.search(re.escape(raw_text), context, re.IGNORECASE)
        if not raw_match:
            return None

        variant_start, variant_end = raw_match.span()
        mentions: List[tuple[str, int, int]] = []
        for gene in self._known_context_genes():
            terms = self._gene_context_terms(gene)
            for term in terms:
                pattern = rf"(?<![A-Za-z0-9]){re.escape(term)}(?![A-Za-z0-9])"
                for match in re.finditer(pattern, context, re.IGNORECASE):
                    mentions.append((gene.upper(), match.start(), match.end()))

        def same_clause_or_row(separator: str, *, has_table: bool) -> bool:
            if len(separator) > 180:
                return False
            if has_table and "|" in separator:
                return True
            return not re.search(r";|\.\s", separator)

        before: List[tuple[int, str]] = []
        for gene, start, end in mentions:
            if end <= variant_start:
                separator = context[end:variant_start]
                if same_clause_or_row(separator, has_table="|" in context):
                    before.append((variant_start - end, gene))
        if before:
            before.sort(key=lambda item: item[0])
            return before[0][1]

        after: List[tuple[int, str]] = []
        for gene, start, _end in mentions:
            if start >= variant_end:
                separator = context[variant_end:start]
                if len(separator) <= 80 and not re.search(r";|\.\s", separator):
                    after.append((start - variant_end, gene))
        if after:
            after.sort(key=lambda item: item[0])
            return after[0][1]

        return None

    def _has_conflicting_gene_context(self, context: str, raw_text: str = "") -> bool:
        """Return true when local context names another gene but not target."""
        assigned_gene = self._gene_assigned_to_variant(context, raw_text)
        if assigned_gene:
            return assigned_gene != self.gene_symbol

        if self._context_mentions_gene(context, self.gene_symbol):
            return False
        return any(
            self._context_mentions_gene(context, other_gene)
            for other_gene in self._known_context_genes()
            if other_gene != self.gene_symbol
        )

    def _scan_protein_variants(self, text: str) -> List[ScannedVariant]:
        """Scan for protein variants."""
        variants = []

        # Parenthesised HGVS first — most specific. p.(Arg176Trp), p.(Ser906Leu).
        # Downstream `seen_normalized` dedup in scan() collapses the inner
        # 3-letter match emitted by PROTEIN_THREE_LETTER, so emitting both is safe.
        for m in self.PROTEIN_HGVS_PAREN.finditer(text):
            ref, pos, alt = m.group(1), m.group(2), m.group(3)
            raw = m.group(0)
            position = int(pos)

            # Same fallback path as the three-letter scanner below: map the
            # captured 3-letter codes through AA_MAP_REVERSE rather than taking
            # the first character (which silently produces e.g. A176T from
            # "Arg176Trp" because Arg/Trp both map to wrong single letters).
            normalized = self.normalizer.normalize_to_single_letter(raw)
            if not normalized:
                ref_single = AA_MAP_REVERSE.get(ref.capitalize())
                if not ref_single:
                    continue
                if alt and alt.lower() in {
                    "fs",
                    "del",
                    "dup",
                    "ins",
                    "ter",
                    "*",
                    "=",
                }:
                    alt_single = alt
                else:
                    alt_single = AA_MAP_REVERSE.get(alt.capitalize(), "X")
                normalized = f"{ref_single}{pos}{alt_single}"

            variants.append(
                ScannedVariant(
                    raw_text=raw,
                    normalized=normalized,
                    variant_type=self._classify_variant_type(alt),
                    notation_type="protein",
                    position=position,
                    context=self._get_context(text, m.start(), m.end()),
                    confidence=0.95,
                    source="protein_hgvs_paren",
                )
            )

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

        # Range deletion patterns: p.Lys1505_Gln1507del, K1505_Q1507del
        for m in self.RANGE_DELETION_PATTERNS.finditer(text):
            ref1, pos1, ref2, pos2 = m.group(1), m.group(2), m.group(3), m.group(4)
            raw = m.group(0)

            ref1_single = (
                AA_MAP_REVERSE.get(ref1.capitalize(), ref1[0].upper())
                if len(ref1) == 3
                else ref1.upper()
            )
            ref2_single = (
                AA_MAP_REVERSE.get(ref2.capitalize(), ref2[0].upper())
                if len(ref2) == 3
                else ref2.upper()
            )
            if not ref1_single or not ref2_single:
                continue

            variants.append(
                ScannedVariant(
                    raw_text=raw,
                    normalized=f"{ref1_single}{pos1}_{ref2_single}{pos2}del",
                    variant_type="deletion",
                    notation_type="protein",
                    position=int(pos1),
                    context=self._get_context(text, m.start(), m.end()),
                    confidence=0.90,
                    source="range_deletion",
                )
            )

        # Single-residue deletion patterns
        for m in self.DELETION_PATTERNS.finditer(text):
            if m.start() > 0 and text[m.start() - 1] in "_-":
                continue

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

        # cDNA delins
        for m in self.CDNA_DELINS.finditer(text):
            pos1, pos2, deleted, inserted = (
                m.group(1),
                m.group(2),
                m.group(3) or "",
                m.group(4),
            )
            raw = m.group(0)
            if pos2:
                normalized = f"c.{pos1}_{pos2}del{deleted.upper()}ins{inserted.upper()}"
            else:
                normalized = f"c.{pos1}del{deleted.upper()}ins{inserted.upper()}"
            variants.append(
                ScannedVariant(
                    raw_text=raw,
                    normalized=normalized,
                    variant_type="delins",
                    notation_type="cdna",
                    position=int(pos1),
                    context=self._get_context(text, m.start(), m.end()),
                    confidence=0.92,
                    source="cdna_delins",
                    variant_class="complex",
                )
            )

        return variants

    def _scan_structural_variants(self, text: str) -> List[ScannedVariant]:
        """Scan for exon-level, large del/dup, and prefixless BIC-style indels."""
        variants: List[ScannedVariant] = []

        for m in self.BIC_PREFIXLESS_INDEL.finditer(text):
            pos, op, bases = m.group(1), m.group(2).lower(), m.group(3).upper()
            ctx = self._get_context(text, m.start(), m.end(), window=80)
            if self.gene_symbol and self.gene_symbol.upper() not in ctx.upper():
                aliases: set[str] = set()
                try:
                    if get_gene_aliases:
                        aliases = {
                            a.upper()
                            for a in get_gene_aliases(
                                self.gene_symbol, include_query_aliases=True
                            )
                        }
                except Exception:
                    aliases = set()
                if not any(a in ctx.upper() for a in aliases):
                    continue
            normalized = f"c.{pos}{op}{bases}"
            variants.append(
                ScannedVariant(
                    raw_text=m.group(0),
                    normalized=normalized,
                    variant_type=op,
                    notation_type="cdna",
                    position=int(pos),
                    context=ctx,
                    confidence=0.75,
                    source="bic_prefixless",
                    variant_class=(
                        "frameshift" if op in {"del", "ins"} else "inframe_indel"
                    ),
                )
            )

        for m in self.EXON_EVENT.finditer(text):
            exon1, exon2 = m.group(1), m.group(2)
            raw = m.group(0)
            op = "deletion" if re.search(r"del", raw, re.I) else "duplication"
            if exon2:
                desc = f"{op} of exons {exon1}-{exon2}"
                key = f"{'del' if op == 'deletion' else 'dup'}:exon{exon1}-{exon2}"
            else:
                desc = f"{op} of exon {exon1}"
                key = f"{'del' if op == 'deletion' else 'dup'}:exon{exon1}"
            vclass = "exon_deletion" if op == "deletion" else "exon_duplication"
            variants.append(
                ScannedVariant(
                    raw_text=raw,
                    normalized=key,
                    variant_type=op,
                    notation_type="structural",
                    position=int(exon1),
                    context=self._get_context(text, m.start(), m.end()),
                    confidence=0.85,
                    source="exon_event",
                    variant_class=vclass,
                    structural_description=desc,
                )
            )

        for m in self.WHOLE_GENE_DEL.finditer(text):
            raw = m.group(0)
            variants.append(
                ScannedVariant(
                    raw_text=raw,
                    normalized="del:wholegene",
                    variant_type="deletion",
                    notation_type="structural",
                    position=None,
                    context=self._get_context(text, m.start(), m.end()),
                    confidence=0.80,
                    source="whole_gene_del",
                    variant_class="large_deletion",
                    structural_description="whole-gene deletion",
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
        line_start = text.rfind("\n", 0, start) + 1
        line_end = text.find("\n", end)
        if line_end == -1:
            line_end = len(text)
        if 0 <= line_start <= start and line_end - line_start <= 500:
            line_context = text[line_start:line_end]
            if any(
                self._context_mentions_gene(line_context, gene)
                for gene in self._known_context_genes()
            ):
                return line_context

            # Legacy Word/table converters can wrap a single table row across
            # physical lines. Pull a few previous lines into the context so
            # gene-column labels such as "SCN5A" or "KCNH2" stay attached to
            # their wrapped variant cells.
            prev_lines: List[str] = []
            cursor = line_start
            for _ in range(3):
                if cursor <= 0:
                    break
                prev_end = cursor - 1
                prev_start = text.rfind("\n", 0, prev_end) + 1
                prev_line = text[prev_start:prev_end]
                candidate_lines = [prev_line] + prev_lines + [line_context]
                candidate = "\n".join(candidate_lines)
                if len(candidate) > 800:
                    break
                prev_lines.insert(0, prev_line)
                if any(
                    self._context_mentions_gene(candidate, gene)
                    for gene in self._known_context_genes()
                ):
                    return candidate
                cursor = prev_start

            if prev_lines:
                return "\n".join(prev_lines + [line_context])
            return line_context

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
        context_quote = " ".join((sv.context or sv.raw_text or "").split())
        new_variant = {
            "gene_symbol": gene_symbol,
            "protein_notation": sv.normalized
            if sv.notation_type == "protein"
            else None,
            "cdna_notation": sv.normalized
            if sv.notation_type in ("cdna", "splice")
            else None,
            "variant_class": getattr(sv, "variant_class", None),
            "structural_description": getattr(sv, "structural_description", None),
            "clinical_significance": "unknown",
            "evidence_level": "scanner",
            "source_location": f"Text scan ({sv.source})",
            "additional_notes": f"Auto-detected via pattern scanning (confidence: {sv.confidence:.2f}, raw: '{sv.raw_text}')",
            "patients": {},
            "penetrance_data": {},
            "individual_records": [],
            "functional_data": {"summary": "", "assays": []},
            "key_quotes": [context_quote] if context_quote else [],
        }
        if sv.notation_type == "structural" and sv.normalized:
            # structural key lives in structural_description; keep class
            if not new_variant["structural_description"]:
                new_variant["structural_description"] = sv.normalized

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
