"""Target-gene specificity checks for extracted variant rows.

Association-heavy genes can pull in variants from nearby disease literature
even when the row is saved under the target gene. These checks are deliberately
conservative: they annotate risk by default and do not delete data unless a
caller opts into a stricter policy later.
"""

from __future__ import annotations

import os
import re
from typing import Any, Optional

from utils.gene_metadata import (
    gene_alias_regex,
    get_gene_aliases,
    get_gene_metadata,
    lookup_variantfeatures_residue,
    normalize_gene_symbol,
)

APOE_ISOFORM_RE = re.compile(
    r"(?:APOE\s*)?(?:[εe]\s*[234]|epsilon\s*[234]|E[234])\b|"
    r"rs429358|rs7412|c\.334T>C|c\.388T>C|p\.Cys112Arg|p\.Arg158Cys|"
    r"Cys112Arg|Arg158Cys|C112R|R158C",
    re.IGNORECASE,
)

PROTEIN_POS_RE = re.compile(
    r"(?:p\.)?(?:[A-Z][a-z]{2}|[A-Z])(\d+)(?:[A-Z][a-z]{2}|[A-Z]|\*|Ter|fs)",
    re.IGNORECASE,
)

PROTEIN_REF_POS_RE = re.compile(
    r"(?:p\.)?([A-Z][a-z]{2}|[A-Z])(\d+)(?:[A-Z][a-z]{2}|[A-Z]|\*|Ter|fs)",
    re.IGNORECASE,
)

CONTEXT_SCAN_MAX_SOURCE_CHARS = int(
    os.environ.get("TARGET_GENE_SPECIFICITY_CONTEXT_MAX_SOURCE_CHARS", "2000000")
)
CONTEXT_SCAN_MAX_VARIANTS = int(
    os.environ.get("TARGET_GENE_SPECIFICITY_CONTEXT_MAX_VARIANTS", "2000")
)
VARIANTFEATURES_MAX_VARIANTS = int(
    os.environ.get("TARGET_GENE_SPECIFICITY_VARIANTFEATURES_MAX_VARIANTS", "2000")
)

AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "TER": "*",
    "STOP": "*",
}


def _alias_pattern(gene_symbol: str) -> re.Pattern[str]:
    return gene_alias_regex(gene_symbol, include_query_aliases=True)


def _compact_gene_name(value: str) -> str:
    return re.sub(r"[^A-Z0-9]+", "", value.upper())


def _variant_gene_matches_target(variant_gene: str, gene_symbol: str) -> bool:
    target_aliases = get_gene_aliases(gene_symbol, include_query_aliases=True)
    normalized_aliases = {_compact_gene_name(alias) for alias in target_aliases}
    normalized_aliases.add(_compact_gene_name(gene_symbol))
    return _compact_gene_name(variant_gene) in normalized_aliases


def _variant_notations(variant: dict[str, Any]) -> list[str]:
    values: list[str] = []
    for key in (
        "cdna_notation",
        "protein_notation",
        "genomic_position",
        "variant",
        "hgvs",
    ):
        value = variant.get(key)
        if value and str(value).strip().lower() not in {"none", "n/a", "unknown"}:
            values.append(str(value).strip())
    return values


def _protein_positions(notations: list[str]) -> list[int]:
    positions: list[int] = []
    for notation in notations:
        for match in PROTEIN_POS_RE.finditer(notation):
            try:
                positions.append(int(match.group(1)))
            except ValueError:
                continue
    return positions


def _cdna_notation(notations: list[str]) -> str | None:
    for notation in notations:
        if notation.lower().startswith(("c.", "nm_")):
            return notation
    return None


def _aa_to_one_letter(value: str) -> str | None:
    token = (value or "").strip()
    if len(token) == 1 and token.isalpha():
        return token.upper()
    return AA3_TO_1.get(token.upper())


def _protein_reference_positions(notations: list[str]) -> list[tuple[str, int, str]]:
    parsed: list[tuple[str, int, str]] = []
    for notation in notations:
        for match in PROTEIN_REF_POS_RE.finditer(notation):
            ref = _aa_to_one_letter(match.group(1))
            if not ref:
                continue
            try:
                parsed.append((ref, int(match.group(2)), notation))
            except ValueError:
                continue
    return parsed


def _context_has_gene(
    source_text: str, notation: str, gene_pattern: re.Pattern[str]
) -> bool:
    if not source_text or not notation:
        return False
    escaped = re.escape(notation)
    for match in re.finditer(escaped, source_text, flags=re.IGNORECASE):
        start = max(0, match.start() - 300)
        end = min(len(source_text), match.end() + 300)
        if gene_pattern.search(source_text[start:end]):
            return True
    return False


def assess_variant_specificity(
    variant: dict[str, Any],
    *,
    gene_symbol: str,
    source_text: str = "",
    check_variantfeatures: bool = True,
) -> dict[str, Any]:
    """Return a specificity annotation for one extracted variant row."""

    gene = normalize_gene_symbol(gene_symbol)
    variant_gene = str(variant.get("gene") or variant.get("gene_symbol") or "").strip()
    notations = _variant_notations(variant)
    reasons: list[str] = []
    status = "supported"
    risk = "low"

    if variant_gene and not _variant_gene_matches_target(variant_gene, gene):
        return {
            "status": "off_target_gene",
            "risk": "high",
            "reasons": [
                f"variant gene {variant_gene} does not match target gene {gene}"
            ],
        }

    if not notations:
        status = "weak"
        risk = "medium"
        reasons.append("no usable variant notation")

    if gene == "APOE" and any(APOE_ISOFORM_RE.search(n) for n in notations):
        reasons.append("recognized APOE isoform/defining allele notation")
        return {"status": status, "risk": risk, "reasons": reasons}

    gene_metadata = get_gene_metadata(gene)
    length = gene_metadata.protein_length
    protein_positions = _protein_positions(notations)
    if length and any(pos > length for pos in protein_positions):
        status = "off_target_risk"
        risk = "high"
        reasons.append(f"protein position exceeds known {gene} length ({length} aa)")

    variantfeatures_checks: list[dict[str, Any]] = []
    cdna = _cdna_notation(notations)
    if check_variantfeatures and status != "off_target_risk":
        for ref, position, protein_notation in _protein_reference_positions(notations):
            support = lookup_variantfeatures_residue(
                gene,
                position=position,
                protein_notation=protein_notation,
                cdna_notation=cdna,
            )
            if not support:
                continue
            expected = set(support.reference_residues)
            check = {
                "position": position,
                "extracted_reference": ref,
                "reference_residues": list(support.reference_residues),
                "alternate_residues": list(support.alternate_residues),
                "transcripts": list(support.transcripts[:10]),
                "matched_hgvs_p": support.matched_hgvs_p,
                "matched_hgvs_c": support.matched_hgvs_c,
            }
            variantfeatures_checks.append(check)
            if expected and ref not in expected:
                if status == "supported":
                    status = "weak"
                    risk = "medium"
                reasons.append(
                    f"VariantFeatures reference residue mismatch at {gene} position "
                    f"{position}: extracted {ref}, expected {'/'.join(sorted(expected))}"
                )

    protein_only = bool(protein_positions) and not any(
        n.lower().startswith("c.")
        or n.lower().startswith("nm_")
        or n.lower().startswith("rs")
        for n in notations
    )
    if protein_only and source_text:
        gene_pattern = _alias_pattern(gene)
        if not any(
            _context_has_gene(source_text, notation, gene_pattern)
            for notation in notations
        ):
            if status == "supported":
                status = "weak"
                risk = "medium"
            reasons.append("protein-only notation lacks nearby target-gene context")

    if not reasons:
        reasons.append("target gene and notation passed specificity checks")
    result: dict[str, Any] = {"status": status, "risk": risk, "reasons": reasons}
    if variantfeatures_checks:
        result["variantfeatures"] = variantfeatures_checks
    if length:
        result["protein_length"] = length
    if gene_metadata.canonical_transcript:
        result["canonical_transcript"] = gene_metadata.canonical_transcript
    return result


def apply_target_gene_specificity(
    extracted_data: Optional[dict[str, Any]],
    *,
    gene_symbol: str,
    source_text: str = "",
    policy: Optional[str] = None,
) -> dict[str, int]:
    """Annotate extracted variants with target-gene specificity metadata.

    policy values:
    - off/none: do nothing
    - flag: annotate only (default)
    - clear: drop high-risk off-target rows after annotating counts in metadata
    """

    selected_policy = (
        policy or os.getenv("TARGET_GENE_SPECIFICITY_POLICY") or "flag"
    ).lower()
    if selected_policy in {"off", "none"}:
        return {"checked": 0, "weak": 0, "off_target_risk": 0, "cleared": 0}
    if not isinstance(extracted_data, dict):
        return {"checked": 0, "weak": 0, "off_target_risk": 0, "cleared": 0}
    variants = extracted_data.get("variants")
    if not isinstance(variants, list):
        return {"checked": 0, "weak": 0, "off_target_risk": 0, "cleared": 0}

    kept: list[dict[str, Any]] = []
    stats = {"checked": 0, "weak": 0, "off_target_risk": 0, "cleared": 0}
    source_for_context = source_text
    context_scan_skipped = False
    if source_text and (
        len(source_text) > CONTEXT_SCAN_MAX_SOURCE_CHARS
        or len(variants) > CONTEXT_SCAN_MAX_VARIANTS
    ):
        source_for_context = ""
        context_scan_skipped = True
    variantfeatures_skipped = len(variants) > VARIANTFEATURES_MAX_VARIANTS
    for variant in variants:
        if not isinstance(variant, dict):
            kept.append(variant)
            continue
        assessment = assess_variant_specificity(
            variant,
            gene_symbol=gene_symbol,
            source_text=source_for_context,
            check_variantfeatures=not variantfeatures_skipped,
        )
        variant["target_gene_specificity"] = assessment
        stats["checked"] += 1
        if assessment["status"] == "weak":
            stats["weak"] += 1
        if assessment["status"] in {"off_target_risk", "off_target_gene"}:
            stats["off_target_risk"] += 1
            if selected_policy == "clear":
                stats["cleared"] += 1
                continue
        kept.append(variant)

    if selected_policy == "clear":
        extracted_data["variants"] = kept
    if stats["checked"]:
        metadata = extracted_data.setdefault("extraction_metadata", {})
        metadata["target_gene_specificity"] = {
            "policy": selected_policy,
            **stats,
        }
        if context_scan_skipped:
            metadata["target_gene_specificity_context_scan"] = {
                "skipped": True,
                "reason": "source_or_variant_count_exceeds_context_scan_limit",
                "source_chars": len(source_text),
                "variant_count": len(variants),
                "max_source_chars": CONTEXT_SCAN_MAX_SOURCE_CHARS,
                "max_variants": CONTEXT_SCAN_MAX_VARIANTS,
            }
        if variantfeatures_skipped:
            metadata["target_gene_specificity_variantfeatures"] = {
                "skipped": True,
                "reason": "variant_count_exceeds_variantfeatures_lookup_limit",
                "variant_count": len(variants),
                "max_variants": VARIANTFEATURES_MAX_VARIANTS,
            }
    return stats
