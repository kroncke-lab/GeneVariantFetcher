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

from utils.llm_trace import attempt_link_summary, llm_attempt_ledger, llm_trace_scope
from utils.llm_utils import (
    BaseLLMCaller,
    clamp_max_tokens,
    normalize_reasoning_effort,
)

SUPPORTED_VERDICTS = {"directly_supported", "inferred_supported"}
UNTRUSTED_VERDICTS = {"ambiguous", "unsupported", "source_missing"}
FIELD_NAMES = ("variant", "total_carriers", "affected", "unaffected")
CLAIM_VERIFIER_DEFAULT_MAX_TOKENS = 2_500
CLAIM_VERIFIER_XHIGH_MAX_TOKENS = 64_000
AA3_TO_1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
}
AA1_TO_3 = {one: three for three, one in AA3_TO_1.items()}


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
    for key in (
        "protein_notation",
        "cdna_notation",
        "legacy_notation",
        "genomic_position",
        "variant",
    ):
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
        three_letter = re.fullmatch(
            r"(?:p\.)?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter)",
            term,
            flags=re.IGNORECASE,
        )
        if three_letter:
            ref, position, alt = three_letter.groups()
            ref = ref.title()
            alt = alt.title()
            ref_one = AA3_TO_1.get(ref)
            alt_one = "*" if alt == "Ter" else AA3_TO_1.get(alt)
            if ref_one and alt_one:
                one_letter_aliases = {
                    f"{ref_one}{position}{alt_one}",
                    f"p.{ref_one}{position}{alt_one}",
                }
                if alt == "Ter":
                    one_letter_aliases.update(
                        {f"{ref_one}{position}X", f"p.{ref_one}{position}X"}
                    )
                terms.update(one_letter_aliases)
        one_letter = re.fullmatch(r"(?:p\.)?([A-Z])(\d+)([A-Z*X])", term)
        if one_letter:
            ref, position, alt = one_letter.groups()
            ref_three = AA1_TO_3.get(ref)
            alt_three = "Ter" if alt in {"*", "X"} else AA1_TO_3.get(alt)
            if ref_three and alt_three:
                terms.update(
                    {
                        f"{ref_three}{position}{alt_three}",
                        f"p.{ref_three}{position}{alt_three}",
                    }
                )
                if alt in {"*", "X"}:
                    terms.update(
                        {
                            f"{ref}{position}*",
                            f"p.{ref}{position}*",
                            f"{ref}{position}X",
                            f"p.{ref}{position}X",
                        }
                    )
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


def _centered_line_excerpt(
    line: str, *, preferred_terms: set[str], max_chars: int = 700
) -> str:
    """Keep the variant mention when source conversion creates a giant paragraph."""
    if len(line) <= max_chars:
        return line
    lower = line.lower()
    positions = [
        lower.find(term.lower())
        for term in sorted(preferred_terms, key=len, reverse=True)
        if term and lower.find(term.lower()) >= 0
    ]
    if not positions:
        return line[:max_chars]
    target = min(positions)
    start = max(0, target - max_chars // 3)
    end = min(len(line), start + max_chars)
    start = max(0, end - max_chars)
    excerpt = line[start:end]
    if start:
        excerpt = "…" + excerpt[1:]
    if end < len(line):
        excerpt = excerpt[:-1] + "…"
    return excerpt


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
            excerpt = (
                _centered_line_excerpt(
                    lines[nearby], preferred_terms=variant_terms, max_chars=700
                )
                if nearby == idx
                else lines[nearby][:700]
            )
            selected.append(f"L{nearby + 1}: {excerpt}")
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


# Byte-stable across cards: rules and the output schema live in the system
# turn so the provider can cache the prefix, and the schema now precedes the
# card instead of trailing it. Only the card itself varies per call.
CLAIM_VERIFICATION_SYSTEM_PROMPT = """You are verifying one biomedical variant/count claim.

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
- Affected means carriers the evidence explicitly assigns the target
  disease/phenotype. Enrollment as a patient/proband or a smaller subset with
  syncope, symptoms, sudden death, arrhythmic events, or severe presentation is
  not enough to manufacture a complete affected/unaffected split.
- Example: "112 LQTS patients, of which 18 were symptomatic" supports
  total_carriers=112, but the affected/unaffected split is null unless the
  source explicitly defines which status is being counted.
- An unselected population, control, or blood-donor cohort is not automatically
  unaffected. Absence of a reported diagnosis is not evidence of absence; set
  unaffected only when the source explicitly calls those carriers unaffected,
  healthy, asymptomatic, or clinically negative for the target phenotype.
- Background/prior-literature, population-database, or referenced-study counts
  should not be treated as the present paper's primary cohort unless the claim
  explicitly describes that cohort.
- If counts are in a pedigree, figure, or table that is referenced but missing
  from the evidence packet, return null for the count fields rather than using
  a nearby family/haplotype aggregate.
- Return null for corrected count values unless the evidence supports a concrete integer.
- The variant can be supported even when total/affected/unaffected are unsupported.

Return strict JSON only:
{
  "verdict": "directly_supported|inferred_supported|ambiguous|unsupported|source_missing",
  "field_verdicts": {
    "variant": "directly_supported|inferred_supported|ambiguous|unsupported|source_missing",
    "total_carriers": "directly_supported|inferred_supported|ambiguous|unsupported|source_missing",
    "affected": "directly_supported|inferred_supported|ambiguous|unsupported|source_missing",
    "unaffected": "directly_supported|inferred_supported|ambiguous|unsupported|source_missing"
  },
  "corrected_values": {
    "total_carriers": "integer or null",
    "affected": "integer or null",
    "unaffected": "integer or null"
  },
  "reason": "brief evidence-based explanation",
  "evidence_quote": "short quote or line reference supporting the verdict"
}
"""


def build_claim_verification_prompt(card: VariantClaimCard) -> str:
    return f"""Claim card:
{card.to_prompt_json()}

Verify this one claim now. Return ONLY the strict JSON object specified in the
system message.
"""


class VariantClaimVerifier(BaseLLMCaller):
    """LLM wrapper for one-card field-level verification."""

    def __init__(
        self,
        model: str,
        temperature: float = 0.0,
        max_tokens: int = CLAIM_VERIFIER_DEFAULT_MAX_TOKENS,
        reasoning_effort: str | None = None,
    ):
        # GPT-5.6 xhigh spends output-budget tokens on hidden reasoning before it
        # emits the small JSON verdict.  A 2.5k cap therefore returns an empty or
        # truncated response even for a compact claim card.  Preserve the cheap
        # default for ordinary models/efforts, but give xhigh the same validated
        # reasoning headroom as the final paper check and count recovery.
        requested_max_tokens = max_tokens
        if (
            "gpt-5.6" in model.lower()
            and normalize_reasoning_effort(reasoning_effort) == "xhigh"
        ):
            requested_max_tokens = max(
                requested_max_tokens, CLAIM_VERIFIER_XHIGH_MAX_TOKENS
            )
        super().__init__(
            model=model,
            temperature=temperature,
            max_tokens=clamp_max_tokens(model, requested_max_tokens),
            reasoning_effort=reasoning_effort,
        )

    def verify(self, card: VariantClaimCard) -> dict[str, Any]:
        with (
            llm_trace_scope(
                gene=card.gene,
                pmid=card.pmid,
                stage="variant_claim_verification",
                component=self.__class__.__name__,
                variant=card.variant,
            ),
            llm_attempt_ledger(),
        ):
            raw = self.call_llm_json(
                build_claim_verification_prompt(card),
                system_message=CLAIM_VERIFICATION_SYSTEM_PROMPT,
            )
            normalized = normalize_verification(raw, card=card)
            self.record_llm_decision(
                "variant_claim_verification_decision",
                {**normalized, **attempt_link_summary()},
            )
            return normalized


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

    corrected["affected"] = None
    field_verdicts["affected"] = "ambiguous"
    corrected["unaffected"] = None
    field_verdicts["unaffected"] = "ambiguous"
    normalized["reason"] = (
        normalized["reason"]
        + " Consistency guard: target-disease enrollment and a smaller "
        "symptom/event subset do not prove an affected/unaffected partition; "
        "the split was cleared."
    ).strip()
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
    for field in ("total_carriers", "affected", "unaffected"):
        if normalized["field_verdicts"].get(field) in UNTRUSTED_VERDICTS:
            normalized["corrected_values"][field] = None
    # Phenotype partitions must be directly supported. Inference from cohort
    # labels, missing diagnosis text, or arithmetic is specifically forbidden by
    # the extraction contract.
    for field in ("affected", "unaffected"):
        if normalized["field_verdicts"].get(field) == "inferred_supported":
            normalized["field_verdicts"][field] = "ambiguous"
            normalized["corrected_values"][field] = None
    normalized = _apply_consistency_guards(normalized, card)
    # Never complete a partition arithmetically. This applies equally to
    # production claim cards and older card-free audit callers.
    return normalized


def apply_verification_to_variant(
    variant: dict[str, Any], verification: dict[str, Any]
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Return a copy with unsupported carrier fields nulled/corrected.

    A typed non-carrier observation (for example, ``family_count``) is not a
    claim about variant-positive individuals. Preserve that source value and
    let the trust gate mask it from carrier-facing consumers. Otherwise claim
    verification can erase an exact table cell merely because it correctly
    rejects the cell as an individual carrier count.
    """
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
    provenance = updated.get("count_provenance")
    provenance = provenance if isinstance(provenance, dict) else {}
    provenance_keys = {
        "total_carriers": "carriers_count_type",
        "affected": "affected_count_type",
        "unaffected": "unaffected_count_type",
    }
    for field, (pdata_key, patient_key) in field_map.items():
        verdict = field_verdicts.get(field, "ambiguous")
        old = coerce_int(
            pdata.get(pdata_key, patients.get(patient_key) if patient_key else None)
        )
        count_type = str(provenance.get(provenance_keys[field]) or "").strip().lower()
        # Trust Gate tg5 has one deliberate raw-observation exception: a
        # source-reported number of families remains auditable in the carrier
        # slot and is field-masked from trusted carrier consumers. Other typed
        # values are not covered by that exception; verification must still
        # clear a case ID, cohort denominator, or mis-typed phenotype count.
        if field == "total_carriers" and count_type == "family_count":
            continue
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
