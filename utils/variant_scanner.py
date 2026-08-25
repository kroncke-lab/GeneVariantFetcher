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
        structural_variant_identity,
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
        structural_variant_identity,
    )


@dataclass
class ScannedVariant:
    """A variant found by the scanner."""

    raw_text: str  # Original matched text
    normalized: str  # Normalized form (e.g., A561V)
    variant_type: str  # missense, frameshift, nonsense, etc.
    notation_type: str  # protein, cdna, splice, legacy, structural
    position: Optional[int]  # Amino acid or nucleotide position
    context: str  # Surrounding text (for debugging)
    confidence: float  # 0.0-1.0 based on pattern quality
    source: str  # Where it was found (narrative, table, etc.)
    variant_class: Optional[str] = None  # closed taxonomy when known
    structural_description: Optional[str] = None  # free-text structural event
    start: Optional[int] = None  # first exact mention offset when available
    end: Optional[int] = None

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
            "start": self.start,
            "end": self.end,
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
                "cdna_notation": v.normalized
                if v.notation_type in {"cdna", "splice"}
                else None,
                "legacy_notation": v.normalized
                if v.notation_type == "legacy"
                else None,
                "variant_class": v.variant_class,
                "structural_description": v.structural_description,
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
            if v.notation_type == "legacy":
                variant_dict["source_notation"] = v.raw_text
            if (
                v.notation_type == "structural"
                and not variant_dict["structural_description"]
            ):
                variant_dict["structural_description"] = v.normalized
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
        r"(?:\*?\d*)?(?!\w)",  # Optional extension + complete-token boundary
        re.IGNORECASE,
    )

    # Parenthesized HGVS protein notation: p.(Arg176Trp), p.(Ser906Leu),
    # p.(Ala1198fs). HGVS uses parentheses to mark *predicted* protein
    # consequences (vs. experimentally confirmed) and this form is pervasive
    # in supplement tables generated by VEP / Annovar / classical clinical
    # genetics pipelines. PROTEIN_HGVS_FULL above cannot match this — its
    # leading "\bp\.[A-Z][a-z]{2}" requires a letter immediately after
    # "p." and "(" breaks that. A spot-check on the consolidated KCNH2 corpus
    # found 937 occurrences across 492
    # unique variants invisible to the scanner.
    PROTEIN_HGVS_PAREN = re.compile(
        r"\bp\.\(([A-Z][a-z]{2})(\d+)"
        r"([A-Z][a-z]{2}|fs\*?\d*|del|dup|ins|Ter|\*|=)"
        r"(?:\*?\d*)?\)",
        re.IGNORECASE,
    )

    # Short HGVS with p.: p.R534C, p.A561V, p.L987fs
    PROTEIN_HGVS_SHORT = re.compile(
        r"\bp\.([A-Z])(\d+)([A-Z]|fs[X\*]?\d*|del|dup|ins|\*|X)(?!\w)",
        re.IGNORECASE,
    )

    # Three-letter AA without p. prefix: Arg534Cys, Ala561Val, Leu987fs
    PROTEIN_THREE_LETTER = re.compile(
        r"\b([A-Z][a-z]{2})(\d{2,4})([A-Z][a-z]{2}|fs\*?\d*|del|dup|ins|Ter|\*)"
        r"(?:[X\*]?\d*)?(?!\w)",
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

    # Duplication variants. Literature sometimes restates the duplicated run
    # after "dup" ("p.R360_Q361dupQKQR"); allow the trailing residue run so
    # the token still matches (normalize_duplication drops the redundant run).
    DUPLICATION_PATTERNS = re.compile(
        r"\b(?:p\.)?"
        r"([A-Z][a-z]{2}|[A-Z])?"
        r"(\d{2,4})"
        r"(?:_([A-Z][a-z]{2}|[A-Z])?(\d{2,4}))?"
        r"dup"
        r"(?:(?:[A-Z][a-z]{2})+|[ACDEFGHIKLMNPQRSTVWY]+)?"
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

    # Prefixless simple indel: BRCA BIC labels such as 185delAG/4184del4, or
    # omitted-c. HGVS in other gene literature. Bases must be UPPERCASE (BIC
    # convention); the (?-i:) scope keeps [ACGT] case-sensitive even under the
    # global IGNORECASE, so "120delta" (del + "ta") no longer matches.
    BIC_PREFIXLESS_INDEL = re.compile(
        r"(?<![A-Za-z0-9_.])(\d{3,5})(del|dup|ins)"
        r"((?-i:[ACGT]{1,20})|\d{1,3})\b",
        re.IGNORECASE,
    )

    # Older cardiac literature also places the affected base before the
    # operation (``2768Cdel``). This is a source-coordinate allele, not a
    # protein call. It is admitted only from a positive, non-negated target-gene
    # finding and later must still pass current-study attribution in
    # ``merge_scanner_results``.
    PREFIXLESS_POSITION_BASE_INDEL = re.compile(
        r"(?<![A-Za-z0-9_.])(\d{3,6})([ACGT])(del|dup)\b",
        re.IGNORECASE,
    )

    # Negation preceding a free-text structural event ("no deletion of exon 5",
    # "ruled out a whole-gene deletion"). Sentence-break aware so it only looks
    # within the current clause.
    NEGATION_BEFORE_RE = re.compile(
        r"(?:\bno\b|\bnot\b|\bnegative\s+for\b|\bwithout\b|\babsen(?:ce|t)\b|"
        r"\bexclud(?:e[ds]?|ing)\b|\brule[ds]?\s+out\b|\bno\s+evidence\b)"
        r"[^.;:!?]{0,40}$",
        re.IGNORECASE,
    )

    # Negation can also follow the structural noun phrase: "deletion of exons
    # 3-5 was not detected". Keep the look-ahead clause-bounded so a negative
    # statement later in the paragraph cannot suppress a real finding.
    NEGATION_AFTER_RE = re.compile(
        r"^[^.;:!?\n]{0,60}\b(?:not|never)\s+(?:be\s+|been\s+)?"
        r"(?:detected|identified|found|observed|seen|confirmed|demonstrated|present)\b"
        r"|^[^.;:!?\n]{0,60}\b(?:was|were|is|are)\s+"
        r"(?:absent|excluded|ruled\s+out)\b",
        re.IGNORECASE,
    )

    # Generic structural phrases are useful only when the surrounding statement
    # says they were a positive finding. This deliberately excludes methods and
    # background prose such as "MLPA detects exon deletions" unless the paper
    # also reports that one was found/observed in the target gene.
    STRUCTURAL_POSITIVE_FINDING_RE = re.compile(
        r"\b(?:identified|detected|found|observed|reported|confirmed|revealed|"
        r"showed|demonstrated|harbou?rs?|harbou?red|carries|carried|segregates|"
        r"segregated)\b"
        r"|\bwe\s+(?:identify|detect|observe|report|confirm|demonstrate)\b"
        r"|\b(?:is|was|were)\s+present\b|\bpresent\s+in\b",
        re.IGNORECASE,
    )

    # Exon-level deletion/duplication is matched with the shared
    # utils.structural_alleles.EXON_EVENT_RE, which supports operation-first
    # ("deletion of exons 3-5") and exons-first ("RyR2 exon 3 deletion")
    # phrasing alike; see _scan_structural_variants.

    # Spelled-out nucleotide events: "insertion of a guanine at position
    # 3108", "deletion of 4 bp at nucleotide position 1261". Older papers
    # narrate single-base events instead of printing HGVS.
    PROSE_NT_INDEL = re.compile(
        r"\b(?P<op>insertion|deletion|duplication)\s+of\s+"
        r"(?:(?:a|an|one|single)\s+)*"
        r"(?:(?P<base>adenine|guanine|cytosine|thymine)|"
        r"(?P<bp>\d{1,4})\s*(?:bp|base\s*pairs?|nucleotides?|bases?))\s+"
        r"at\s+(?:nucleotide\s+)?position\s+(?P<pos>\d{1,6})\b",
        re.IGNORECASE,
    )

    # Companion clause naming the protein-level consequence of the event
    # above ("... truncation of the protein at amino acid position 1036");
    # anchors a residue-less frameshift at that codon.
    PROSE_TRUNCATION = re.compile(
        r"\btruncat(?:ed|es|ion|ing)\b[^.;!?]{0,160}?"
        r"\bat\s+(?:amino[\s-]*acid\s+)?position\s+(?P<aa>\d{1,5})\b",
        re.IGNORECASE,
    )

    #: Explicitly non-coding context. "position 1261" inside a promoter,
    #: genomic or vector statement is not a cDNA coordinate, and emitting
    #: ``c.`` there asserts a coding allele the paper never stated.
    PROSE_NONCODING_CONTEXT = re.compile(
        r"\bpromoter\b|\bgenomic\s+(?:dna|sequence|position|coordinate)|"
        r"\bg\.\d|\bchr(?:omosome)?\s*\d|\b(?:5|3)['\u2032]\s*UTR\b|"
        r"\bconstruct\b|\bvector\b|\bplasmid\b|\bcloning\s+site\b|"
        r"\bprimer\b",
        re.IGNORECASE,
    )

    #: A bare length ("4 bp at position 1261") names no sequence, so the
    #: coordinate system has to be stated somewhere in the sentence. A named
    #: base ("a guanine at position 3108") is nucleotide-level by
    #: construction and needs no further witness.
    PROSE_CDNA_CONTEXT = re.compile(
        r"\bc\.\d|\bcDNA\b|\bcoding\s+(?:sequence|region|dna)\b|"
        r"\bnucleotide\b|\bexon\s+\d|\btranscript\b|\bNM_\d",
        re.IGNORECASE,
    )

    #: A negated or merely predicted consequence is not an observation.
    PROSE_NEGATION = re.compile(
        r"\b(?:no|not|never|without|absence\s+of|failed\s+to|neither|"
        r"did\s+not|does\s+not|was\s+not|were\s+not)\b",
        re.IGNORECASE,
    )

    PROSE_BASE_LETTERS = {
        "adenine": "A",
        "guanine": "G",
        "cytosine": "C",
        "thymine": "T",
    }

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
        # Dense tables can repeat the same literal notation thousands of times.
        # Document attribution scans the full source for every occurrence, so
        # recomputing it for each candidate creates an avoidable
        # candidates-times-document-length pass. Attribution depends only on
        # this document, target gene, and literal token; cache it for this scan.
        document_attribution_cache: dict[str, set[str]] = {}

        # Run all pattern matchers
        protein_variants = self._scan_protein_variants(text)
        cdna_variants = self._scan_cdna_variants(text)
        splice_variants = self._scan_splice_variants(text)
        narrative_variants = self._scan_narrative_variants(text)
        structural_variants = self._scan_structural_variants(text)
        prose_indel_variants = self._scan_prose_indel_events(text)

        # Combine all findings
        all_candidates = (
            protein_variants
            + cdna_variants
            + splice_variants
            + narrative_variants
            + structural_variants
            + prose_indel_variants
        )

        # Filter and deduplicate
        filtered_count = 0
        for v in all_candidates:
            # Once one representation of a normalized identity has survived
            # every guard, later representations cannot change the emitted
            # result. Put this before even the local gene-context work: dense
            # tables can repeat an accepted identity thousands of times. If
            # the first representation is rejected it is never added to
            # ``seen_normalized``, so a later target-gene occurrence still gets
            # the full evaluation.
            if v.normalized in seen_normalized:
                continue

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

            attribution_key = v.raw_text.casefold()
            assigned_genes = document_attribution_cache.get(attribution_key)
            if assigned_genes is None:
                assigned_genes = self._document_gene_attribution(text, v.raw_text)
                document_attribution_cache[attribution_key] = assigned_genes
            if assigned_genes and self.gene_symbol.upper() not in assigned_genes:
                filtered_count += 1
                logger.debug(
                    "Scanner: filtered %s; document attributes it to %s",
                    v.raw_text,
                    ", ".join(sorted(assigned_genes)),
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

            # Matchers already know these offsets but the original scanner object
            # discarded them. Persist a conservative first exact occurrence so
            # downstream attribution can recover document position without
            # trusting the shortened debug context.
            if v.start is None and v.raw_text:
                mention = re.search(re.escape(v.raw_text), text, re.IGNORECASE)
                if mention:
                    v.start, v.end = mention.span()

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

        # A gene token joined to the variant with at most one punctuation
        # character and no whitespace ("Q2958R-RyR2", "KCNE1/G38S",
        # "R14C(KCNQ1)") labels that specific mention; it outranks any
        # detached clause-mate in either direction. The nearest-preceding
        # preference below stays as the tiebreak for detached mentions only.
        attached: List[tuple[int, str]] = []
        for gene, start, end in mentions:
            if end <= variant_start:
                separator = context[end:variant_start]
            elif start >= variant_end:
                separator = context[variant_end:start]
            else:
                continue
            if len(separator) <= 1 and not re.search(r"[\w\s]", separator):
                attached.append((len(separator), gene))
        if attached:
            attached.sort(key=lambda item: (item[0], item[1]))
            return attached[0][1]

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

    # A ±50-character window is enough to see "KCNQ1 L187P" but not
    # "...nucleotide 560 of KCNQ1, which results in a substitution of amino
    # acid residue leucine by proline (L187P)". The document-level pass below
    # re-reads every occurrence of the token with a window wide enough to hold
    # a full clause.
    DOCUMENT_ATTRIBUTION_WINDOW = 240

    def _document_gene_attribution(self, text: str, raw_text: str) -> set[str]:
        """Genes this document attaches to ``raw_text``, across all mentions.

        The per-hit filter only sees one narrow window, and papers that report
        two genes routinely name the gene once per sentence rather than next to
        every repetition of the variant. Reading every occurrence recovers the
        attribution that a single window misses.
        """
        if not text or not raw_text:
            return set()
        pattern = rf"(?<![A-Za-z0-9]){re.escape(raw_text)}(?![A-Za-z0-9])"
        assigned: set[str] = set()
        for match in re.finditer(pattern, text, re.IGNORECASE):
            start = max(0, match.start() - self.DOCUMENT_ATTRIBUTION_WINDOW)
            end = min(len(text), match.end() + self.DOCUMENT_ATTRIBUTION_WINDOW)
            window = text[start:end]
            gene = self._gene_assigned_to_variant(window, raw_text)
            if not gene:
                continue
            # ``_gene_assigned_to_variant`` prefers the nearest mention
            # *before* the token and only falls back to the one after it, so a
            # list like "P297S KCNH2 and P1177L SCN5A" assigns P1177L to
            # KCNH2. When a window names another gene but also names the target
            # gene, that single reading is not trustworthy: the window
            # abstains. It must not vote for the target either, or one
            # ambiguous sentence would rescue a token every other sentence in
            # the paper attributes elsewhere.
            if gene != self.gene_symbol.upper() and self._context_mentions_gene(
                window, self.gene_symbol
            ):
                continue
            assigned.add(gene)
        return assigned

    def _document_assigns_variant_elsewhere(self, text: str, raw_text: str) -> bool:
        """True when the whole document only ever attaches this token to
        another gene.

        Fail-open by construction: a single mention that names the target gene,
        or a document that names no gene at all, keeps the candidate. It fires
        only on the unambiguous case -- every attributed mention belongs to a
        different gene -- which is how ``L187P in KCNQ1`` stops being scanned
        as a KCNH2 variant on a paper that sequenced both.
        """
        assigned = self._document_gene_attribution(text, raw_text)
        if not assigned:
            return False
        return self.gene_symbol.upper() not in assigned

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

    def _structural_statement_context(
        self, text: str, start: int, end: int, window: int = 240
    ) -> str:
        """Return the clause-sized statement containing a structural match."""
        lower = max(0, start - window)
        upper = min(len(text), end + window)
        before = text[lower:start]
        after = text[end:upper]
        # A bare period is not a safe boundary because HGVS c./p. prefixes are
        # common immediately before the event. Only a period followed by space
        # (or the end of the window) closes the statement.
        boundary = re.compile(r"\n|[;?!]|\.(?=\s|$)")
        left_matches = list(boundary.finditer(before))
        left = lower + left_matches[-1].end() if left_matches else lower
        right_match = boundary.search(after)
        right = end + right_match.start() if right_match else upper
        return text[left:right].strip()

    def _structural_match_is_negated(
        self, text: str, start: int, end: int, window: int = 80
    ) -> bool:
        """Return true for negation immediately before or after an event."""
        before = text[max(0, start - window) : start]
        after = text[end : min(len(text), end + window)]
        return bool(
            self.NEGATION_BEFORE_RE.search(before)
            or self.NEGATION_AFTER_RE.search(after)
        )

    def _positive_structural_context(
        self, text: str, start: int, end: int
    ) -> Optional[str]:
        """Return evidence context only for target-gene positive findings."""
        context = self._structural_statement_context(text, start, end)
        if not self._context_mentions_gene(context, self.gene_symbol):
            return None
        if self._structural_match_is_negated(text, start, end):
            return None
        if not self.STRUCTURAL_POSITIVE_FINDING_RE.search(context):
            return None
        return context

    def _scan_structural_variants(self, text: str) -> List[ScannedVariant]:
        """Scan for exon-level, large del/dup, delta-peptide, and prefixless
        BIC-style indels."""
        variants: List[ScannedVariant] = []

        from utils.legacy_notation import gene_supports_legacy_notation
        from utils.structural_alleles import (
            DELTA_RE,
            EXON_EVENT_RE,
            parse_delta_peptide,
            parse_exon_event,
        )

        for m in self.BIC_PREFIXLESS_INDEL.finditer(text):
            pos, op, payload = m.group(1), m.group(2).lower(), m.group(3).upper()
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
            source_notation = f"{pos}{op}{payload}"
            is_legacy = gene_supports_legacy_notation(self.gene_symbol)
            # A numeric payload (for example BRCA2 BIC ``4184del4``) reports
            # the number of affected bases, not literal HGVS sequence. Outside
            # genes with this known legacy system, prefixing it with ``c.``
            # fabricates an allele. Base payloads such as ``185delAG`` remain a
            # valid omitted-prefix cDNA hint for other genes.
            if payload.isdigit() and not is_legacy:
                continue
            normalized = source_notation if is_legacy else f"c.{source_notation}"
            variants.append(
                ScannedVariant(
                    raw_text=m.group(0),
                    normalized=normalized,
                    variant_type=op,
                    notation_type="legacy" if is_legacy else "cdna",
                    position=int(pos),
                    context=ctx,
                    confidence=0.75,
                    source="bic_prefixless" if is_legacy else "prefixless_cdna_indel",
                    variant_class=(
                        "frameshift" if op in {"del", "ins"} else "inframe_indel"
                    ),
                )
            )

        for m in self.PREFIXLESS_POSITION_BASE_INDEL.finditer(text):
            pos, base, op = m.group(1), m.group(2).upper(), m.group(3).lower()
            ctx = self._positive_structural_context(text, m.start(), m.end())
            if ctx is None:
                continue
            # ``reported`` is part of the general positive-finding vocabulary,
            # but a prior-literature statement must not become a HIGH prompt
            # hint before the downstream merge has a chance to reject it.
            if re.search(
                r"\b(?:previously|prior(?:ly)?|earlier)\s+(?:been\s+)?"
                r"(?:reported|described|published|identified|found|observed)\b",
                ctx,
                re.IGNORECASE,
            ):
                continue
            variants.append(
                ScannedVariant(
                    raw_text=m.group(0),
                    normalized=f"c.{pos}{op}{base}",
                    variant_type=op,
                    notation_type="cdna",
                    position=int(pos),
                    context=ctx,
                    confidence=0.80,
                    source="prefixless_position_base_indel",
                    variant_class="frameshift",
                )
            )

        for m in EXON_EVENT_RE.finditer(text):
            event = parse_exon_event(m.group(0))
            if not event:
                continue
            context = self._positive_structural_context(text, m.start(), m.end())
            if context is None:
                continue
            op_key, first, last = event
            op = "deletion" if op_key == "del" else "duplication"
            if last != first:
                desc = f"{op} of exons {first}-{last}"
                key = f"{op_key}:exon{first}-{last}"
            else:
                desc = f"{op} of exon {first}"
                key = f"{op_key}:exon{first}"
            vclass = "exon_deletion" if op_key == "del" else "exon_duplication"
            variants.append(
                ScannedVariant(
                    raw_text=m.group(0),
                    normalized=key,
                    variant_type=op,
                    notation_type="structural",
                    position=first,
                    context=context,
                    confidence=0.85,
                    source="exon_event",
                    variant_class=vclass,
                    structural_description=desc,
                )
            )

        # Delta-peptide alleles ("ΔKPQ", "delta KPQ"): record only the raw
        # token; expand_structural_keys resolves it to the protein range
        # downstream, so no residue arithmetic belongs here.
        for m in DELTA_RE.finditer(text):
            peptide = parse_delta_peptide(m.group(0))
            if not peptide:
                continue
            context = self._positive_structural_context(text, m.start(), m.end())
            if context is None:
                continue
            raw = m.group(0)
            variants.append(
                ScannedVariant(
                    raw_text=raw,
                    normalized=f"DELTA{peptide}",
                    variant_type="deletion",
                    notation_type="structural",
                    position=None,
                    context=context,
                    confidence=0.85,
                    source="delta_notation",
                    variant_class="inframe_indel",
                    structural_description=raw,
                )
            )

        for m in self.WHOLE_GENE_DEL.finditer(text):
            context = self._positive_structural_context(text, m.start(), m.end())
            if context is None:
                continue
            raw = m.group(0)
            variants.append(
                ScannedVariant(
                    raw_text=raw,
                    normalized="del:wholegene",
                    variant_type="deletion",
                    notation_type="structural",
                    position=None,
                    context=context,
                    confidence=0.80,
                    source="whole_gene_del",
                    variant_class="large_deletion",
                    structural_description="whole-gene deletion",
                )
            )

        return variants

    def _scan_prose_indel_events(self, text: str) -> List[ScannedVariant]:
        """Scan spelled-out nucleotide indels and their truncation clauses."""
        variants: List[ScannedVariant] = []

        for m in self.PROSE_NT_INDEL.finditer(text):
            statement = self._structural_statement_context(text, m.start(), m.end())
            if self.PROSE_NONCODING_CONTEXT.search(statement):
                continue
            phrase_at = statement.find(m.group(0))
            preamble = statement[:phrase_at] if phrase_at > 0 else statement
            if self.PROSE_NEGATION.search(preamble):
                continue
            op = m.group("op").lower()[:3]
            pos = int(m.group("pos"))
            base = m.group("base")
            if not base and not self.PROSE_CDNA_CONTEXT.search(statement):
                # Bare-length form with no stated coordinate system.
                continue
            if base:
                normalized = f"c.{pos}{op}{self.PROSE_BASE_LETTERS[base.lower()]}"
            else:
                length = int(m.group("bp"))
                # An insertion length names no reference span; a fabricated
                # range would assert an allele the paper never printed.
                if op == "ins":
                    continue
                if length <= 1:
                    normalized = f"c.{pos}{op}"
                else:
                    normalized = f"c.{pos}_{pos + length - 1}{op}"
            variants.append(
                ScannedVariant(
                    raw_text=m.group(0),
                    normalized=normalized,
                    variant_type=op,
                    notation_type="cdna",
                    position=pos,
                    context=self._get_context(text, m.start(), m.end(), window=100),
                    confidence=0.70,
                    source="prose_nt_indel",
                )
            )

            # A same-sentence truncation clause anchors the protein-level
            # consequence at a codon even without a reference residue.
            truncation = self.PROSE_TRUNCATION.search(statement)
            if truncation and not self.PROSE_NEGATION.search(
                statement[: truncation.start()]
            ):
                aa_pos = int(truncation.group("aa"))
                variants.append(
                    ScannedVariant(
                        raw_text=truncation.group(0),
                        normalized=f"{aa_pos}fs",
                        variant_type="frameshift",
                        notation_type="protein",
                        position=aa_pos,
                        context=statement[:300],
                        confidence=0.60,
                        source="prose_truncation",
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
                    variant_class="splice",
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
                    variant_class="splice",
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
    document_text: str = "",
) -> Dict[str, Any]:
    """
    Merge scanner-found variants with LLM-extracted variants.

    Only adds variants that weren't already extracted by an authoritative layer
    and whose source text independently supports both the target gene and a
    current-study observation. Regex confidence is identity confidence, not
    gene/study attribution; silence therefore fails closed.

    Args:
        extracted_data: LLM extraction result dict
        scan_result: Scanner result
        gene_symbol: Target gene symbol
        min_confidence: Minimum scanner confidence to include
        document_text: Full source text used for gene/study attribution

    Returns:
        Updated extracted_data with merged variants
    """
    existing_variants = extracted_data.get("variants", [])
    from utils.legacy_notation import preserve_source_only_legacy_identity

    for variant in existing_variants:
        if isinstance(variant, dict):
            preserve_source_only_legacy_identity(variant, gene_symbol)
    target_gene = gene_symbol.upper()
    scanner = VariantScanner(target_gene)
    source_text = document_text or "\n".join(
        v.context for v in scan_result.variants if v.context
    )

    variant_token_re = re.compile(
        r"(?:[cp]\.\(?[A-Za-z*][A-Za-z]{0,2}\d{1,6}[^\s,;)]*"
        r"|\bIVS\d+[^\s,;)]*"
        r"|\b[A-Z][a-z]{0,2}\d{2,5}(?:[A-Z*]|fs|del|dup|ins)[^\s,;)]*)",
        re.IGNORECASE,
    )
    study_subject = (
        r"(?:we|our\s+(?:study|cohort|series)|(?:the\s+)?(?:current|present)\s+study|"
        r"patients?|probands?|carriers?|cases?|subjects?|participants?|individuals?|"
        r"index\s+(?:case|patient)|famil(?:y|ies))"
    )
    observation = (
        r"(?:identif(?:y|ied)|detect(?:ed)?|found|carri(?:ed|es)|harbou?r(?:ed|s)?|"
        r"observ(?:ed|es)|present(?:ed)?\s+with|was\s+positive\s+for)"
    )
    positive_study_re = re.compile(
        rf"(?:{study_subject}.{{0,180}}{observation}|"
        rf"{observation}.{{0,180}}{study_subject})",
        re.IGNORECASE | re.DOTALL,
    )
    background_re = re.compile(
        r"\b(?:previous(?:ly)?|prior|earlier)\s+(?:been\s+)?"
        r"(?:reported|described|published|identified|found|observed)\b|"
        r"\b(?:literature|review|overview|catalog(?:ue)?|compilation|"
        r"founder\s+(?:panel|mutation)|testing\s+panel|screened\s+for)\b|"
        r"\b(?:reported|described|published)\s+by\b",
        re.IGNORECASE,
    )
    reference_re = re.compile(
        r"\bet\s+al\b|\b(?:references|bibliography)\b",
        re.IGNORECASE,
    )
    negative_finding_re = re.compile(
        r"\b(?:could\s+not|did\s+not|was\s+not|were\s+not|not)\s+"
        r"(?:confirm(?:ed)?|replicat(?:ed|e)|validat(?:ed|e))\b|"
        r"\b(?:sequencing|DNA|PCR|technical)\s+artifacts?\b|"
        r"\bfalse[- ]positive\b",
        re.IGNORECASE,
    )
    # External-evidence signals inside a table row: a citation ("Smith et al.
    # 2010"), an rs id, or a population/clinical repository name marks the row
    # as a compilation of other studies' findings, never a current-study
    # observation. A BARE year is deliberately not a signal here: clinical
    # rows carry years of diagnosis, onset and follow-up, so vetoing on one
    # discards exactly the current-study rows this lane exists to admit. A
    # year only counts when it sits in a citation shape.
    row_external_signal_re = re.compile(
        r"\bet\s+al\b|\brs\d{3,}\b|"
        r"\((?:[^()]{0,40}?\b(?:19|20)\d{2}[a-z]?)\)|"
        r"\b(?:gnomad|exac|clinvar|dbsnp|hgmd|topmed|bravo|esp6500|lovd|bic|"
        r"1000\s*genomes|uk\s*biobank)\b",
        re.IGNORECASE,
    )
    # A caption or heading that scopes its table to previously published or
    # database-sourced variants. Applied to the nearest preceding caption
    # because a physical row cannot show it.
    compilation_scope_re = re.compile(
        r"\b(?:previously|prior(?:ly)?|already)\s+(?:reported|published|"
        r"described|identified)\b|"
        r"\b(?:catalog(?:ue)?|compilation|compendium|summary|review|list(?:ing)?)"
        r"\s+of\b|"
        r"\breported\s+(?:in\s+the\s+)?literature\b|"
        r"\bliterature[- ]reported\b|\bpublished\s+(?:variants?|mutations?)\b",
        re.IGNORECASE,
    )

    document_genes = {
        gene
        for gene in scanner._known_context_genes()
        if scanner._context_mentions_gene(source_text, gene)
    }
    single_gene_scope = document_genes == {target_gene}
    target_terms = {
        re.sub(r"[^A-Z0-9]", "", term.upper())
        for term in scanner._gene_context_terms(target_gene)
    }
    generic_gene_stopwords = {
        "AIC",
        "ANOVA",
        "BIC",
        "CDNA",
        "CI",
        "DNA",
        "EDTA",
        "HGVS",
        "HR",
        "MRI",
        "NGS",
        "OR",
        "PCR",
        "PMID",
        "RNA",
        "SNP",
        "USA",
    }

    def scanner_key(sv: ScannedVariant) -> str:
        if sv.notation_type == "structural":
            structural = structural_variant_identity(
                sv.structural_description or sv.normalized
            )
            return f"STRUCT:{structural}"
        if sv.notation_type == "legacy":
            from utils.legacy_notation import normalize_legacy_notation

            legacy = normalize_legacy_notation(sv.normalized or sv.raw_text)
            return f"LEGACY:{legacy or sv.raw_text}".upper()
        normalized = normalize_variant(sv.normalized or sv.raw_text, target_gene)
        return (normalized or sv.normalized or sv.raw_text).upper()

    def preceding_scope(line_start: int) -> str:
        """Nearest caption/heading above a row, or "" when none is close.

        A row's own text cannot show that its table is a catalogue of prior
        publications; the caption above it can. Bounded to a few lines so an
        unrelated paragraph far above never scopes a table.
        """

        if line_start <= 0:
            return ""
        head = source_text[:line_start].rsplit("\n", 9)[:-1]
        for candidate in reversed(head):
            text = candidate.strip()
            if not text:
                continue
            if re.match(
                r"^[#*_>\s]*(?:e?table|figure|fig\.?|supplement(?:al|ary)?)\b",
                text,
                re.IGNORECASE,
            ):
                return text[:400]
            if text.startswith("#"):
                return text[:400]
        return ""

    def mention_contexts(sv: ScannedVariant) -> List[Tuple[str, str, str]]:
        """Return (physical line, evidence window, table scope) per mention."""
        if not source_text or not sv.raw_text:
            context = sv.context or ""
            return [(context, context, "")] if context else []
        pattern = re.compile(re.escape(sv.raw_text), re.IGNORECASE)
        contexts: List[Tuple[str, str, str]] = []
        for match in list(pattern.finditer(source_text))[:100]:
            line_start = source_text.rfind("\n", 0, match.start()) + 1
            line_end = source_text.find("\n", match.end())
            if line_end == -1:
                line_end = len(source_text)
            line = source_text[line_start:line_end]
            window_start = max(line_start, match.start() - 320)
            window_end = min(line_end, match.end() + 320)
            contexts.append(
                (
                    line,
                    source_text[window_start:window_end],
                    preceding_scope(line_start),
                )
            )
        if not contexts and sv.context:
            contexts.append((sv.context, sv.context, ""))
        return contexts

    def nearby_generic_gene(
        window: str, raw_text: str
    ) -> Optional[Tuple[str, int, int]]:
        """Nearest gene-like symbol (with span) absent from built-in metadata.

        This is an abstention aid, not a new attribution oracle. A nearby
        all-caps biomedical symbol such as PALB2 may veto a BRCA2 merge; it can
        never positively establish the requested gene.
        """
        raw_match = re.search(re.escape(raw_text), window, re.IGNORECASE)
        if not raw_match:
            return None
        candidates: List[Tuple[int, str, int, int]] = []
        for match in re.finditer(
            r"(?<![A-Za-z0-9])[A-Z][A-Z0-9-]{1,11}(?![A-Za-z0-9])", window
        ):
            if match.end() > raw_match.start() and match.start() < raw_match.end():
                continue
            label = re.sub(r"[^A-Z0-9]", "", match.group(0).upper())
            if label in generic_gene_stopwords or label.isdigit():
                continue
            # Nearby variant tokens are not gene symbols. Without this guard,
            # a valid list such as ``RYR2 H469Y, L2299F`` makes each allele
            # falsely veto the other as an unknown gene label.
            if re.fullmatch(r"[A-Z]\d{1,6}[A-Z*]", label):
                continue
            distance = min(
                abs(raw_match.start() - match.end()),
                abs(match.start() - raw_match.end()),
            )
            if distance <= 100:
                candidates.append((distance, label, match.start(), match.end()))
        if not candidates:
            return None
        candidates.sort(key=lambda item: item[0])
        return candidates[0][1], candidates[0][2], candidates[0][3]

    def generic_intercepts_target(
        window: str, raw_text: str, generic_span: Tuple[int, int]
    ) -> bool:
        """True when the unknown symbol sits between the nearest target-gene
        mention and the variant token.

        An affirmative target assignment outranks an unknown-symbol hunch (a
        phenotype cell such as "Yes, VT"), but an unknown symbol that
        intercepts the assignment path is the label that actually owns the
        token -- compilation rows list several genes' variants on one line.
        Fail closed: no locatable variant or target mention keeps the veto.
        """
        raw_match = re.search(re.escape(raw_text), window, re.IGNORECASE)
        if not raw_match:
            return True
        nearest: Optional[Tuple[int, int, int]] = None
        for term in scanner._gene_context_terms(target_gene):
            pattern = rf"(?<![A-Za-z0-9]){re.escape(term)}(?![A-Za-z0-9])"
            for match in re.finditer(pattern, window, re.IGNORECASE):
                distance = min(
                    abs(raw_match.start() - match.end()),
                    abs(match.start() - raw_match.end()),
                )
                if nearest is None or distance < nearest[0]:
                    nearest = (distance, match.start(), match.end())
        if nearest is None:
            return True
        _, target_start, target_end = nearest
        if target_end <= raw_match.start():
            region = (target_end, raw_match.start())
        elif raw_match.end() <= target_start:
            region = (raw_match.end(), target_start)
        else:
            return False
        generic_start, generic_end = generic_span
        return generic_start < region[1] and generic_end > region[0]

    def is_table_like(line: str) -> bool:
        token_count = len(variant_token_re.findall(line))
        return bool(
            line.count("|") >= 2
            or token_count >= 3
            or (
                re.search(
                    r"\b(?:e?table|supplement(?:al|ary)?\s+table)\s*\d*",
                    line,
                    re.IGNORECASE,
                )
                and token_count
            )
            or (len(line) > 500 and token_count >= 2)
            or (
                len(re.findall(r"\bet\s+al\b|\b(?:19|20)\d{2}\b", line, re.I)) >= 2
                and token_count >= 2
            )
            or bool(re.search(r"-{4,}\d+-{1,}", line))
        )

    def is_row_shaped(line: str, window: str) -> bool:
        """Table rows never contain study prose; recognize them structurally.

        Deliberately stricter than ``is_table_like``, which also treats bare
        variant-token density as table-ish. That is right for flagging a
        mention as table-adjacent, but it is not evidence of a row: a prose
        sentence listing three variants ("mutations R176W, G584S and N470D
        have been reported") would otherwise enter the row lane and skip the
        study-prose requirement entirely.
        """

        if line.count("|") >= 2 or re.search(r"-{4,}\d+-{1,}", line):
            return True
        # PDF/Word extraction can flatten a genotype table to one space per
        # cell. Require the full row fingerprint: numeric row id, a known gene,
        # a genomic coordinate, ref/alt base cells, frequency, and an HGVS
        # protein token. This is deliberately much narrower than merely seeing
        # several variants or numbers in prose.
        single_space_genotype_row = bool(
            re.match(r"^\s*\d+\*?(?:\s+[A-Z])?\s+", line)
            and len(line.split()) >= 11
            and any(
                scanner._context_mentions_gene(line, gene)
                for gene in scanner._known_context_genes()
            )
            and re.search(r"\b\d{5,}\b", line)
            and re.search(
                r"\s[ACGT-]+\s+[ACGT-]+\s+"
                r"(?:0(?:\.\d+)?|1(?:\.0+)?)\s+p\.",
                line,
                re.IGNORECASE,
            )
            and not re.search(
                r"\b(?:we|patients?|probands?|subjects?|participants?|cohort|study|"
                r"mutations?|variants?|identified|detected|found|observed|reported|"
                r"had|was|were)\b",
                line,
                re.IGNORECASE,
            )
        )
        if single_space_genotype_row:
            return True
        # Wrapped/converted rows keep cell shape without pipes: no sentence
        # punctuation and several short delimited fields.
        if re.search(r"[.;!?](?:\s|$)", window):
            return False
        cells = [cell.strip() for cell in re.split(r"\||\t|\s{2,}", window)]
        cells = [cell for cell in cells if cell]
        return len(cells) >= 3 and all(len(cell) <= 40 for cell in cells)

    def attribution_decision(sv: ScannedVariant) -> Tuple[bool, str, str]:
        """Fail-closed decision; any independently supported mention may pass."""
        saw_table = False
        saw_reference = False
        saw_background = False
        saw_wrong_gene = False
        saw_gene_support = False
        best_quote = sv.context or sv.raw_text

        # A wrong-gene prefix embedded in the token itself ("KCNH2_R176W")
        # dooms every mention; no window can overrule it.
        embedded_gene = re.match(
            r"^([A-Z][A-Z0-9]{1,11})[_-]", sv.raw_text, re.IGNORECASE
        )
        if embedded_gene and embedded_gene.group(1).upper() not in target_terms:
            return False, "wrong_gene", " ".join(best_quote.split())[:240]

        for line, window, scope in mention_contexts(sv):
            quote = " ".join(window.split())[:240]
            best_quote = quote or best_quote
            if is_row_shaped(line, window):
                saw_table = True
                # The caption scopes the whole table; a compilation caption
                # disqualifies every row under it even though no row can
                # carry that wording itself.
                if scope and (
                    compilation_scope_re.search(scope)
                    or reference_re.search(scope)
                    or background_re.search(scope)
                ):
                    saw_reference = True
                    continue
                # Table rows by construction never contain study prose, so
                # the prose lane below can never admit them. A row is
                # attributable on its own when it affirms the target gene and
                # carries no external-evidence signal marking it as a
                # compilation of other studies' or repositories' findings.
                if (
                    row_external_signal_re.search(window)
                    or reference_re.search(window)
                    or negative_finding_re.search(window)
                    or background_re.search(window)
                ):
                    continue
                assigned_gene = scanner._gene_assigned_to_variant(window, sv.raw_text)
                if assigned_gene != target_gene:
                    if assigned_gene:
                        saw_wrong_gene = True
                    continue
                generic = nearby_generic_gene(window, sv.raw_text)
                if (
                    generic
                    and generic[0] not in target_terms
                    and generic_intercepts_target(window, sv.raw_text, generic[1:])
                ):
                    saw_wrong_gene = True
                    continue
                return True, "table_row_target_gene", quote
            if reference_re.search(window):
                saw_reference = True
                continue
            if negative_finding_re.search(window):
                saw_background = True
                continue

            # "not previously reported" supports novelty; remove that phrase
            # before applying the prior-literature veto.
            veto_text = re.sub(
                r"\bnot\s+(?:been\s+)?previously\s+"
                r"(?:reported|described|published|identified|found|observed)\b",
                "",
                window,
                flags=re.IGNORECASE,
            )
            if background_re.search(veto_text):
                saw_background = True
                continue

            assigned_gene = scanner._gene_assigned_to_variant(window, sv.raw_text)
            generic = nearby_generic_gene(window, sv.raw_text)
            if generic and generic[0] not in target_terms:
                # An affirmative target assignment outranks an unknown-symbol
                # hunch unless the symbol intercepts the assignment path.
                if assigned_gene != target_gene or generic_intercepts_target(
                    window, sv.raw_text, generic[1:]
                ):
                    saw_wrong_gene = True
                    continue
            window_genes = {
                gene
                for gene in scanner._known_context_genes()
                if scanner._context_mentions_gene(window, gene)
            }
            gene_supported = assigned_gene == target_gene or (
                target_gene in window_genes and len(window_genes) == 1
            )
            if assigned_gene and assigned_gene != target_gene:
                saw_wrong_gene = True
                continue
            if not gene_supported and single_gene_scope:
                gene_supported = True
            if not gene_supported:
                continue
            saw_gene_support = True

            if positive_study_re.search(window):
                return True, "current_study_target_gene", quote

        if saw_wrong_gene and not saw_gene_support:
            reason = "wrong_gene"
        elif saw_table:
            reason = "table_like"
        elif saw_reference:
            reason = "reference_list"
        elif saw_background:
            reason = "background_mention"
        elif saw_gene_support:
            reason = "study_unattributed"
        else:
            reason = "gene_unattributed"
        return False, reason, " ".join(best_quote.split())[:240]

    from utils.legacy_notation import normalize_legacy_notation

    # Build set of existing variant keys (normalized). Structural-only rows need
    # their description in this identity set; otherwise every scanner replay
    # appends a duplicate because both point-notation fields are empty.
    existing_keys = set()
    for v in existing_variants:
        protein = normalize_variant(v.get("protein_notation", "") or "", gene_symbol)
        cdna = normalize_variant(v.get("cdna_notation", "") or "", gene_symbol)
        if protein:
            existing_keys.add(protein.upper())
        if cdna:
            existing_keys.add(cdna.upper())
        structural = structural_variant_identity(v.get("structural_description"))
        if structural:
            existing_keys.add(f"STRUCT:{structural}")
        legacy = str(v.get("legacy_notation") or "").strip()
        if legacy:
            existing_keys.add(f"LEGACY:{legacy}".upper())
        # Deterministic table lanes still prefix a bare legacy label with "c."
        # when they cannot see the verbatim cell. Such a row and a scanner
        # legacy identity are the same mutation, so register the legacy form of
        # any cDNA that is exactly "c." + a strict legacy label. Without this
        # the two lanes emit a duplicate pair: the table row carries the counts
        # and the legacy row carries none.
        if cdna:
            bare = normalize_legacy_notation(
                cdna[2:] if cdna[:2].upper() == "C." else ""
            )
            from utils.legacy_notation import gene_supports_legacy_notation

            if bare and gene_supports_legacy_notation(gene_symbol):
                existing_keys.add(f"LEGACY:{bare}".upper())

    # Add scanner variants not already present
    added_count = 0
    skipped: List[Dict[str, Any]] = []
    for sv in scan_result.variants:
        if sv.confidence < min_confidence:
            skipped.append(
                {
                    "variant": sv.normalized,
                    "reason": "below_confidence",
                    "quote": " ".join((sv.context or sv.raw_text).split())[:240],
                }
            )
            continue

        candidate_key = scanner_key(sv)

        equivalent_existing = candidate_key in existing_keys
        if not equivalent_existing and re.match(
            r"^C\.\d+(?:_\d+)?(?:DEL|DUP|INS)$", candidate_key
        ):
            equivalent_existing = any(
                existing_key.startswith(candidate_key) for existing_key in existing_keys
            )

        if equivalent_existing:
            skipped.append(
                {
                    "variant": sv.normalized,
                    "reason": "already_extracted",
                    "quote": " ".join((sv.context or sv.raw_text).split())[:240],
                }
            )
            continue

        allowed, reason, quote = attribution_decision(sv)
        if not allowed:
            skipped.append({"variant": sv.normalized, "reason": reason, "quote": quote})
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
            "legacy_notation": sv.normalized if sv.notation_type == "legacy" else None,
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
        if sv.notation_type == "legacy":
            new_variant["source_notation"] = sv.raw_text
        if sv.notation_type == "structural" and sv.normalized:
            # structural key lives in structural_description; keep class
            if not new_variant["structural_description"]:
                new_variant["structural_description"] = sv.normalized

        existing_variants.append(new_variant)
        existing_keys.add(candidate_key)
        added_count += 1

    if added_count > 0:
        logger.info(f"Scanner merge added {added_count} variants not found by LLM")

        # Update metadata
        if "extraction_metadata" in extracted_data:
            extracted_data["extraction_metadata"]["total_variants_found"] = len(
                existing_variants
            )
            extracted_data["extraction_metadata"]["scanner_added"] = added_count

    metadata = extracted_data.setdefault("extraction_metadata", {})
    metadata["scanner_merge"] = {
        "policy": "target_gene_current_study_v1",
        "candidate_count": len(scan_result.variants),
        "hint_count": min(len(scan_result.variants), 50),
        "added": added_count,
        "skipped": skipped,
    }
    metadata["scanner_added"] = added_count

    extracted_data["variants"] = existing_variants
    return extracted_data


# ==========================================================================
# TESTING
# ==========================================================================
