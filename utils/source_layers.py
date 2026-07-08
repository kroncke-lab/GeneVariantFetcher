"""Source-layer labeling and cheap notation-quality gates."""

from __future__ import annotations

import re
from typing import Iterable

SOURCE_LAYERS = {
    "llm_text",
    "llm_table",
    "regex_text",
    "regex_table",
    "clinvar",
    "pubtator",
    "figure",
    "mixed",
}

JUNK_REJECT_LAYERS = {"figure", "regex_table"}

KNOWN_GENE_SYMBOLS = {
    "ABCC9",
    "ACTC1",
    "AKAP9",
    "ANK2",
    "APC",
    "APOB",
    "BRAF",
    "BRCA1",
    "BRCA2",
    "CACNA1C",
    "CACNA2D1",
    "CACNB2",
    "CALM1",
    "CALM2",
    "CALM3",
    "CASQ2",
    "CAV3",
    "DES",
    "DSC2",
    "DSG2",
    "DSP",
    "EGFR",
    "FLNC",
    "HCN4",
    "JUP",
    "KCNE1",
    "KCNE2",
    "KCNH2",
    "KCNJ2",
    "KCNJ5",
    "KCNQ1",
    "KRAS",
    "LDLR",
    "LMNA",
    "MLH1",
    "MSH2",
    "MSH6",
    "MYBPC3",
    "MYH7",
    "MYL2",
    "MYL3",
    "PCSK9",
    "PIK3CA",
    "PKP2",
    "PLN",
    "PMS2",
    "PTEN",
    "RYR2",
    "SCN1B",
    "SCN2B",
    "SCN3B",
    "SCN4B",
    "SCN5A",
    "SNTA1",
    "TECRL",
    "TNNC1",
    "TNNI3",
    "TNNT2",
    "TP53",
    "TPM1",
    "TRDN",
    "TTN",
}

_REGEX_TABLE_MARKERS = (
    "regex extraction",
    "table regex",
    "deterministic",
    "fixed-width",
    "router-first",
    "vertical gene table",
    "lqt mutation table",
)
_LLM_TEXT_MARKERS = (
    "functional study",
    "discussion",
    "abstract",
    "results",
    "introduction",
    "additional file",
    "supplementary material",
    "supplement",
)
_RESIDUE_PROSE_RE = re.compile(r"\bresidues?\b", re.IGNORECASE)


def normalize_source_layer(layer: object) -> str | None:
    """Normalize a persisted or inferred source-layer label."""

    text = str(layer or "").strip().lower()
    if not text:
        return None
    if text == "manual_or_legacy":
        return "llm_text"
    if text in SOURCE_LAYERS:
        return text
    return None


def source_layer_tokens(layer: object) -> set[str]:
    """Split a layer string such as ``figure,regex_table`` into labels."""

    tokens: set[str] = set()
    for part in re.split(r"[,;|]", str(layer or "")):
        normalized = normalize_source_layer(part)
        if normalized:
            tokens.add(normalized)
    return tokens


def infer_source_layer_from_text(
    *,
    source_location: object = None,
    additional_notes: object = None,
    extraction_source: object = None,
    source_layer: object = None,
) -> str:
    """Infer the stable source layer for a variant-paper link."""

    explicit = normalize_source_layer(source_layer)
    if explicit:
        return explicit

    text = " ".join(
        str(value or "")
        for value in (source_location, additional_notes, extraction_source)
    ).lower()
    if "clinvar" in text:
        return "clinvar"
    if "pubtator" in text:
        return "pubtator"
    if "figure" in text:
        return "figure"
    if any(marker in text for marker in _REGEX_TABLE_MARKERS):
        return "regex_table"
    if "scanner" in text or "text scan" in text:
        return "regex_text"
    if "table" in text:
        return "llm_table"
    if any(marker in text for marker in _LLM_TEXT_MARKERS):
        return "llm_text"
    return "llm_text"


def source_layer_sql_case(
    source_location_expr: str,
    source_layer_expr: str | None = None,
    *,
    additional_notes_expr: str | None = None,
) -> str:
    """Return a SQLite CASE expression mirroring ``infer_source_layer_from_text``."""

    explicit = (
        f"NULLIF(LOWER(TRIM({source_layer_expr})), '')" if source_layer_expr else "NULL"
    )
    loc = f"LOWER(COALESCE({source_location_expr}, ''))"
    notes = (
        f"LOWER(COALESCE({additional_notes_expr}, ''))"
        if additional_notes_expr
        else "''"
    )
    text = f"({loc} || ' ' || {notes})"
    return f"""
        COALESCE(
            CASE
                WHEN {explicit} = 'manual_or_legacy' THEN 'llm_text'
                WHEN {explicit} IN (
                    'llm_text', 'llm_table', 'regex_text', 'regex_table',
                    'clinvar', 'pubtator', 'figure', 'mixed'
                ) THEN {explicit}
                ELSE NULL
            END,
            CASE
                WHEN {text} LIKE '%clinvar%' THEN 'clinvar'
                WHEN {text} LIKE '%pubtator%' THEN 'pubtator'
                WHEN {text} LIKE '%figure%' THEN 'figure'
                WHEN {text} LIKE '%regex%' OR {text} LIKE '%deterministic%'
                     OR {text} LIKE '%fixed-width%' OR {text} LIKE '%router%'
                     OR {text} LIKE '%vertical gene table%'
                     OR {text} LIKE '%lqt mutation table%'
                     THEN 'regex_table'
                WHEN {text} LIKE '%scanner%' OR {text} LIKE '%text scan%'
                     THEN 'regex_text'
                WHEN {text} LIKE '%table%' THEN 'llm_table'
                WHEN {text} LIKE '%functional study%' OR {text} LIKE '%discussion%'
                     OR {text} LIKE '%abstract%' OR {text} LIKE '%results%'
                     OR {text} LIKE '%introduction%'
                     OR {text} LIKE '%additional file%'
                     OR {text} LIKE '%supplementary material%'
                     OR {text} LIKE '%supplement%'
                     THEN 'llm_text'
                ELSE 'llm_text'
            END
        )
    """


def _clean_gene_candidate(value: object) -> str:
    return re.sub(r"\s+", "", str(value or "").strip()).upper()


def _clean_protein_short_candidate(value: object) -> str:
    text = re.sub(r"\s+", "", str(value or "").strip())
    if text.lower().startswith("p."):
        text = text[2:]
    return text


def _known_genes_for(row_gene_symbol: object) -> set[str]:
    genes = set(KNOWN_GENE_SYMBOLS)
    row_gene = _clean_gene_candidate(row_gene_symbol)
    if row_gene:
        genes.add(row_gene)
    return genes


def junk_notation_reason(
    *,
    source_layer: object,
    protein_notation: object = None,
    cdna_notation: object = None,
    variant: object = None,
    gene_symbol: object = None,
) -> str | None:
    """Return a stable reject reason for known junk in noisy source layers."""

    if not (source_layer_tokens(source_layer) & JUNK_REJECT_LAYERS):
        return None

    genes = _known_genes_for(gene_symbol)
    for value in (protein_notation, cdna_notation, variant):
        if _clean_gene_candidate(value) in genes:
            return "gene_symbol_as_variant"

    protein = _clean_protein_short_candidate(protein_notation)
    if protein and len(protein) <= 2:
        return "short_protein_notation"

    for value in (protein_notation, cdna_notation, variant):
        if _RESIDUE_PROSE_RE.search(str(value or "")):
            return "residue_prose"

    return None


def count_reasons(reasons: Iterable[str]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for reason in reasons:
        counts[reason] = counts.get(reason, 0) + 1
    return counts
