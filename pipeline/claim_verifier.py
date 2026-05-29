"""Small-context verification for variant carrier-count claims.

The verifier is intentionally narrower than extraction: it receives one
variant/count claim plus the exact nearby evidence lines and decides whether
each structured field is supported, ambiguous, unsupported, or source-blocked.
"""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass
from typing import Any

from utils.llm_utils import BaseLLMCaller, clamp_max_tokens

SUPPORTED_VERDICTS = {"directly_supported", "inferred_supported"}
UNTRUSTED_VERDICTS = {"ambiguous", "unsupported", "source_missing"}
FIELD_NAMES = ("variant", "total_carriers", "affected", "unaffected")


@dataclass
class VariantClaimCard:
    gene: str
    disease: str | None
    pmid: str
    title: str | None
    variant: str
    extracted: dict[str, int | None]
    evidence: str
    source_location: str | None = None

    def to_prompt_json(self) -> str:
        return json.dumps(asdict(self), ensure_ascii=False, indent=2)


def coerce_int(value: Any) -> int | None:
    if value is None or value == "":
        return None
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return int(value)
    match = re.search(r"-?\d+", str(value))
    return int(match.group()) if match else None


def variant_display(variant: dict[str, Any]) -> str:
    for key in ("protein_notation", "cdna_notation", "genomic_position", "variant"):
        value = variant.get(key)
        if value is not None and str(value).strip():
            return str(value).strip()
    return ""


def extracted_counts(variant: dict[str, Any]) -> dict[str, int | None]:
    patients = variant.get("patients") or {}
    pdata = variant.get("penetrance_data") or {}
    return {
        "total_carriers": coerce_int(
            pdata.get("total_carriers_observed", patients.get("count"))
        ),
        "affected": coerce_int(pdata.get("affected_count")),
        "unaffected": coerce_int(pdata.get("unaffected_count")),
    }


def _variant_terms(variant: str) -> set[str]:
    terms = {variant.strip()} if variant and variant.strip() else set()
    for term in list(terms):
        if term.startswith("p."):
            terms.add(term[2:])
        if len(term) >= 3:
            terms.add(term.replace("p.", "").replace(" ", ""))
    return {term for term in terms if term}


def _is_markdown_table_line(line: str) -> bool:
    stripped = line.strip()
    return stripped.startswith("|") and stripped.endswith("|")


def _is_markdown_separator(line: str) -> bool:
    stripped = line.strip()
    return bool(re.fullmatch(r"\|[\s:.\-|\+]+\|", stripped)) and "---" in stripped


def _table_context_indices(lines: list[str], idx: int, max_scan: int = 60) -> set[int]:
    """Return nearby header/title indices for a markdown table row."""
    if not _is_markdown_table_line(lines[idx]):
        return set()
    start = max(0, idx - max_scan)
    for pos in range(idx - 1, start - 1, -1):
        if _is_markdown_separator(lines[pos]) and pos > 0:
            header = pos - 1
            indices = {header, pos}
            for title_idx in range(max(0, header - 4), header):
                if lines[title_idx].strip() and not _is_markdown_table_line(
                    lines[title_idx]
                ):
                    indices.add(title_idx)
            return indices
        if lines[pos].strip() and not _is_markdown_table_line(lines[pos]):
            break
    return set()


def build_evidence_snippet(
    *,
    source_text: str,
    gene: str | None,
    variant: str,
    counts: dict[str, int | None],
    max_chars: int = 4000,
    window: int = 2,
) -> str:
    """Select nearby lines for one variant/count claim."""
    lines = source_text.splitlines()
    gene_l = (gene or "").lower()
    variant_terms = _variant_terms(variant)
    count_values = {str(value) for value in counts.values() if value is not None}
    count_words = (
        "carrier",
        "carriers",
        "affected",
        "unaffected",
        "asymptomatic",
        "symptomatic",
        "patient",
        "patients",
        "proband",
        "probands",
        "control",
        "controls",
        "table",
        "mutation",
        "variant",
    )

    scored: list[tuple[int, int]] = []
    for idx, line in enumerate(lines):
        compact = line.replace(" ", "")
        lower = line.lower()
        score = 0
        if gene_l and gene_l in lower:
            score += 1
        if any(
            term.lower() in lower or term.replace(" ", "") in compact
            for term in variant_terms
        ):
            score += 8
        if any(
            value and re.search(rf"(?<!\d){re.escape(value)}(?!\d)", line)
            for value in count_values
        ):
            score += 2
        if any(word in lower for word in count_words):
            score += 1
        if score:
            scored.append((score, idx))

    selected: list[str] = []
    seen: set[int] = set()
    for _score, idx in sorted(scored, key=lambda item: (-item[0], item[1])):
        context_indices = set(
            range(max(0, idx - window), min(len(lines), idx + window + 1))
        )
        context_indices.update(_table_context_indices(lines, idx))
        for nearby in sorted(context_indices):
            if nearby in seen:
                continue
            seen.add(nearby)
            selected.append(f"L{nearby + 1}: {lines[nearby][:700]}")
        if sum(len(item) + 1 for item in selected) >= max_chars:
            break
    return "\n".join(selected)[:max_chars]


def build_claim_card(
    *,
    source_text: str,
    gene: str,
    disease: str | None,
    pmid: str,
    title: str | None,
    variant: dict[str, Any],
    max_evidence_chars: int = 4000,
) -> VariantClaimCard | None:
    display = variant_display(variant)
    if not display:
        return None
    counts = extracted_counts(variant)
    evidence = build_evidence_snippet(
        source_text=source_text,
        gene=gene,
        variant=display,
        counts=counts,
        max_chars=max_evidence_chars,
    )
    return VariantClaimCard(
        gene=gene,
        disease=disease,
        pmid=pmid,
        title=title,
        variant=display,
        extracted=counts,
        evidence=evidence,
        source_location=variant.get("source_location"),
    )


def build_claim_verification_prompt(card: VariantClaimCard) -> str:
    return f"""You are verifying one biomedical variant/count claim.

Do not do broad extraction from the paper. Verify and, when the local evidence
clearly supports it, correct this one claim. Field verdict labels refer to the
final corrected value you return, not necessarily the original extracted value:
- directly_supported: the exact table row/sentence supports the value.
- inferred_supported: the value follows from a very constrained clinical case inference.
- ambiguous: evidence mentions the field but the semantics are unclear.
- unsupported: evidence does not support the extracted value.
- source_missing: the evidence says the relevant table/supplement is missing or absent.

Important rules:
- Count variant-positive people, not everyone enrolled, sampled, sequenced, or
  screened in the cohort. If a paper says 39 cases and 46 controls were
  sampled, but only 17 cases and 9 controls were heterozygous for the variant,
  then total_carriers=26, affected=17, unaffected=9 for that variant.
- Do not copy study-wide, family-set, domain, or mutation-class totals onto a variant.
- Do not copy aggregate carrier counts across several variants/families onto
  one variant unless a row, pedigree, or sentence gives that variant-specific
  count.
- If the original extracted value is wrong but the evidence supports a concrete
  replacement, mark that field directly_supported or inferred_supported and put
  the replacement integer in corrected_values.
- Control counts are unaffected only when the evidence says those controls carry the variant.
- Assay replicate/cell counts are not carrier counts.
- Do not interpret unlabeled table numbers as carrier counts. First use the
  table header/title: prediction scores, classifications, exon/domain columns,
  allele frequencies, and functional scores are not carrier counts.
- Affected means carriers with the target disease/phenotype, not a subset with
  syncope, symptoms, sudden death, arrhythmic events, or severe presentation
  unless that symptom/event is the target phenotype.
- If evidence says N target-disease patients/probands carry the variant, then
  affected=N even if only a smaller subset is symptomatic.
- Example: "112 LQTS patients, of which 18 were symptomatic" means
  total_carriers=112 and affected=112 for Long QT syndrome; 18 is a symptom or
  event subset, not affected_count.
- Example: "16 carriers in an unselected population cohort" means
  total_carriers=16, and unaffected=16 if no target-disease diagnosis is
  reported for those carriers.
- If evidence says N carriers come from an unselected population, control, or
  blood-donor cohort and does not report target-disease diagnosis in those
  carriers, then unaffected=N for the target disease.
- Background/prior-literature, population-database, or referenced-study counts
  should not be treated as the present paper's primary cohort unless the claim
  explicitly describes that cohort.
- If counts are in a pedigree, figure, or table that is referenced but missing
  from the evidence packet, return null for the count fields rather than using
  a nearby family/haplotype aggregate.
- Return null for corrected count values unless the evidence supports a concrete integer.
- The variant can be supported even when total/affected/unaffected are unsupported.

Claim card:
{card.to_prompt_json()}

Return strict JSON only:
{{
  "verdict": "directly_supported|inferred_supported|ambiguous|unsupported|source_missing",
  "field_verdicts": {{
    "variant": "directly_supported|inferred_supported|ambiguous|unsupported|source_missing",
    "total_carriers": "directly_supported|inferred_supported|ambiguous|unsupported|source_missing",
    "affected": "directly_supported|inferred_supported|ambiguous|unsupported|source_missing",
    "unaffected": "directly_supported|inferred_supported|ambiguous|unsupported|source_missing"
  }},
  "corrected_values": {{
    "total_carriers": "integer or null",
    "affected": "integer or null",
    "unaffected": "integer or null"
  }},
  "reason": "brief evidence-based explanation",
  "evidence_quote": "short quote or line reference supporting the verdict"
}}
"""


class VariantClaimVerifier(BaseLLMCaller):
    """LLM wrapper for one-card field-level verification."""

    def __init__(
        self,
        model: str,
        temperature: float = 0.0,
        max_tokens: int = 2500,
        reasoning_effort: str | None = None,
    ):
        super().__init__(
            model=model,
            temperature=temperature,
            max_tokens=clamp_max_tokens(model, max_tokens),
            reasoning_effort=reasoning_effort,
        )

    def verify(self, card: VariantClaimCard) -> dict[str, Any]:
        raw = self.call_llm_json(build_claim_verification_prompt(card))
        return normalize_verification(raw, card=card)


def _disease_aliases(disease: str | None) -> set[str]:
    text = (disease or "").lower()
    aliases = {part.strip() for part in re.split(r"[,;/]", text) if part.strip()}
    if "long qt" in text:
        aliases.update({"lqts", "long qt syndrome"})
    if "short qt" in text:
        aliases.update({"sqts", "short qt syndrome"})
    if "brugada" in text:
        aliases.update({"brugada", "brugada syndrome"})
    if "catecholaminergic polymorphic ventricular tachycardia" in text:
        aliases.update({"cpvt", "cpvt1"})
    return {alias for alias in aliases if len(alias) >= 3}


def _mentions_target_disease_group(
    evidence: str, disease: str | None, total: int
) -> bool:
    evidence_l = evidence.lower()
    aliases = _disease_aliases(disease)
    group_words = r"(?:patients|probands|cases|individuals|members|subjects)"
    for alias in aliases:
        alias_re = re.escape(alias)
        patterns = (
            rf"\b{total}\s+{alias_re}\s+{group_words}\b",
            rf"\b{total}\s+clinically\s+affected\s+{group_words}\b",
            rf"\b{alias_re}\s+{group_words}\b.{{0,80}}\b{total}\b",
        )
        if any(re.search(pattern, evidence_l) for pattern in patterns):
            return True
    return False


def _apply_consistency_guards(
    normalized: dict[str, Any], card: VariantClaimCard | None
) -> dict[str, Any]:
    if card is None or not card.disease:
        return normalized

    corrected = normalized["corrected_values"]
    field_verdicts = normalized["field_verdicts"]
    total = corrected.get("total_carriers") or card.extracted.get("total_carriers")
    affected = corrected.get("affected") or card.extracted.get("affected")
    original_affected = card.extracted.get("affected")
    if total is None or affected is None:
        return normalized
    has_smaller_symptom_subset = affected < total or (
        corrected.get("affected") == total
        and original_affected is not None
        and original_affected < total
    )
    if not has_smaller_symptom_subset:
        return normalized

    symptom_subset_words = (
        "symptomatic",
        "symptoms",
        "syncope",
        "syncopal",
        "sudden death",
        "arrhythmic event",
        "ventricular fibrillation",
    )
    evidence_l = card.evidence.lower()
    if not any(word in evidence_l for word in symptom_subset_words):
        return normalized
    if not _mentions_target_disease_group(card.evidence, card.disease, total):
        return normalized

    corrected["affected"] = total
    field_verdicts["affected"] = "directly_supported"
    if field_verdicts.get("unaffected") not in SUPPORTED_VERDICTS:
        corrected["unaffected"] = 0
        field_verdicts["unaffected"] = "inferred_supported"
    normalized["reason"] = (
        normalized["reason"]
        + " Consistency guard: target-disease patient/proband count overrides "
        "a smaller symptom/event subset for affected_count."
    ).strip()
    return normalized


def _apply_count_identity_guard(normalized: dict[str, Any]) -> dict[str, Any]:
    corrected = normalized["corrected_values"]
    field_verdicts = normalized["field_verdicts"]

    def supported(field: str) -> bool:
        return (
            corrected.get(field) is not None
            and field_verdicts.get(field) in SUPPORTED_VERDICTS
        )

    total = corrected.get("total_carriers")
    affected = corrected.get("affected")
    unaffected = corrected.get("unaffected")
    if (
        supported("total_carriers")
        and supported("affected")
        and not supported("unaffected")
    ):
        inferred = total - affected
        if inferred >= 0:
            corrected["unaffected"] = inferred
            field_verdicts["unaffected"] = "inferred_supported"
    if (
        supported("total_carriers")
        and supported("unaffected")
        and not supported("affected")
    ):
        inferred = total - unaffected
        if inferred >= 0:
            corrected["affected"] = inferred
            field_verdicts["affected"] = "inferred_supported"
    if (
        supported("affected")
        and supported("unaffected")
        and not supported("total_carriers")
    ):
        corrected["total_carriers"] = affected + unaffected
        field_verdicts["total_carriers"] = "inferred_supported"
    return normalized


def normalize_verification(
    raw: dict[str, Any], card: VariantClaimCard | None = None
) -> dict[str, Any]:
    field_verdicts = raw.get("field_verdicts") or {}
    normalized_fields = {}
    for field in FIELD_NAMES:
        verdict = str(field_verdicts.get(field) or raw.get("verdict") or "ambiguous")
        verdict = verdict.strip().lower()
        normalized_fields[field] = (
            verdict
            if verdict in SUPPORTED_VERDICTS | UNTRUSTED_VERDICTS
            else "ambiguous"
        )

    corrected = raw.get("corrected_values") or {}
    normalized = {
        "verdict": str(raw.get("verdict") or "ambiguous").strip().lower(),
        "field_verdicts": normalized_fields,
        "corrected_values": {
            "total_carriers": coerce_int(corrected.get("total_carriers")),
            "affected": coerce_int(corrected.get("affected")),
            "unaffected": coerce_int(corrected.get("unaffected")),
        },
        "reason": str(raw.get("reason") or ""),
        "evidence_quote": str(raw.get("evidence_quote") or ""),
    }
    if normalized["verdict"] not in SUPPORTED_VERDICTS | UNTRUSTED_VERDICTS:
        normalized["verdict"] = "ambiguous"
    if normalized["verdict"] in SUPPORTED_VERDICTS:
        for field in ("total_carriers", "affected", "unaffected"):
            if (
                normalized["corrected_values"].get(field) is not None
                and normalized["field_verdicts"].get(field) not in SUPPORTED_VERDICTS
            ):
                normalized["field_verdicts"][field] = "inferred_supported"
    normalized = _apply_consistency_guards(normalized, card)
    return _apply_count_identity_guard(normalized)


def apply_verification_to_variant(
    variant: dict[str, Any], verification: dict[str, Any]
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Return a copy with unsupported count fields nulled/corrected."""
    updated = dict(variant)
    patients = dict(updated.get("patients") or {})
    pdata = dict(updated.get("penetrance_data") or {})
    field_verdicts = verification.get("field_verdicts") or {}
    corrected = verification.get("corrected_values") or {}
    changes: dict[str, Any] = {}

    field_map = {
        "total_carriers": ("total_carriers_observed", "count"),
        "affected": ("affected_count", None),
        "unaffected": ("unaffected_count", None),
    }
    for field, (pdata_key, patient_key) in field_map.items():
        verdict = field_verdicts.get(field, "ambiguous")
        old = coerce_int(
            pdata.get(pdata_key, patients.get(patient_key) if patient_key else None)
        )
        if verdict in SUPPORTED_VERDICTS:
            new = corrected.get(field)
            if new is not None and new != old:
                pdata[pdata_key] = new
                if patient_key:
                    patients[patient_key] = new
                changes[field] = {"old": old, "new": new, "verdict": verdict}
        else:
            if old is not None:
                pdata[pdata_key] = None
                if patient_key:
                    patients[patient_key] = None
                changes[field] = {"old": old, "new": None, "verdict": verdict}

    updated["patients"] = patients
    updated["penetrance_data"] = pdata
    updated["claim_verification"] = verification
    return updated, changes
