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

from pipeline.count_provenance import (
    PATIENT_ROW_PHENOTYPE_SOURCE,
    SOURCE_BOUND_PHENOTYPE_SOURCE,
)
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
    derivation: dict[str, Any] | None = None
    # High-confidence phenotype terms taken from the paper title itself.  This
    # is deliberately separate from ``disease``: the latter is a gene-level
    # candidate vocabulary assembled from ClinGen and manual aliases, and may
    # contain several syndromes that are not the phenotype studied by this
    # particular paper.
    paper_target_phenotypes: list[str] | None = None
    paper_title_scope: str | None = None

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
            pdata.get(
                "total_carriers_observed",
                patients.get("count", variant.get("carriers")),
            )
        ),
        "affected": coerce_int(pdata.get("affected_count", variant.get("affected"))),
        "unaffected": coerce_int(
            pdata.get("unaffected_count", variant.get("unaffected"))
        ),
        "uncertain": coerce_int(pdata.get("uncertain_count", variant.get("uncertain"))),
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


def _identity_aliases(variant: dict[str, Any]) -> set[str]:
    """Return source-facing aliases without inventing an identity mapping."""
    aliases: set[str] = set()
    for key in (
        "protein_notation",
        "cdna_notation",
        "legacy_notation",
        "source_notation",
        "genomic_position",
        "variant",
    ):
        value = str(variant.get(key) or "").strip()
        if not value:
            continue
        aliases.add(value)
        aliases.update(_variant_terms(value))
        for part in re.split(r"[;|,]", value):
            part = part.strip(" ()[]")
            if 4 <= len(part) <= 100:
                aliases.add(part)
                aliases.update(_variant_terms(part))
    for fact in variant.get("fact_provenance") or []:
        if not isinstance(fact, dict) or fact.get("fact_type") != "variant_identity":
            continue
        value = str(fact.get("fact_value") or "").strip()
        if value:
            aliases.add(value)
            aliases.update(_variant_terms(value))
            for part in re.split(r"[;|,]", value):
                part = part.strip(" ()[]")
                if 4 <= len(part) <= 100:
                    aliases.add(part)
    return {alias for alias in aliases if len(alias) >= 4}


def _source_fold(value: str) -> str:
    """Fold presentation-only differences for quote/source comparison."""
    folded = str(value or "").translate(
        str.maketrans({"\u2018": "'", "\u2019": "'", "\u201c": '"', "\u201d": '"'})
    )
    return re.sub(r"\s+", " ", folded).strip().casefold()


def _quote_is_in_source(quote: str, source_text: str) -> bool:
    quote_folded = _source_fold(quote).strip("\"' ")
    return len(quote_folded) >= 16 and quote_folded in _source_fold(source_text)


def _quote_has_alias(quote: str, aliases: set[str]) -> bool:
    quote_folded = _source_fold(quote).replace(" ", "")
    return any(
        _source_fold(alias).replace(" ", "") in quote_folded
        for alias in aliases
        if len(_source_fold(alias).replace(" ", "")) >= 4
    )


def _quote_has_count(quote: str, value: int) -> bool:
    if re.search(rf"(?<!\d){value}(?!\d)", quote):
        return True
    if value != 1:
        return False
    return bool(
        re.search(
            r"\b(?:a|an|one)\s+(?:patient|proband|carrier|individual|subject|child)\b",
            quote,
            flags=re.IGNORECASE,
        )
    )


def _quote_mentions_target_disease(quote: str, disease: str | None) -> bool:
    folded = _source_fold(quote)
    return any(alias in folded for alias in _disease_aliases(disease))


_PHENOTYPE_EVENT_WORDS = re.compile(
    r"\b(?:symptoms?|symptomatic|asymptomatic|syncope|syncopal|cardiac events?|"
    r"arrhythmic events?|sudden cardiac death|aborted cardiac arrest|"
    r"ventricular fibrillation|torsade(?:s)? de pointes?|icd shocks?)\b",
    re.IGNORECASE,
)
_NONEXACT_COUNT_WORDS = re.compile(
    r"(?:\b(?:about|approximately|roughly|nearly|at least|at most|more than|"
    r"fewer than|up to)\b|[~<>]=?)",
    re.IGNORECASE,
)
_UNASSESSED_RESIDUE_WORDS = re.compile(
    r"\b(?:unknown|unassessed|not assessed|not phenotyped|without phenotype data|"
    r"lost to follow-up|remaining carriers? (?:were )?not evaluated)\b",
    re.IGNORECASE,
)
_ASSESSED_SUBSET_CARRIER_WORDS = re.compile(
    r"\b(?:known|evaluated|phenotyped|genotyped|followed|available)\b.{0,80}\bcarriers?\b|"
    r"\bcarriers?\s+(?:with|who had)\s+(?:available\s+)?(?:ecg|phenotype|follow-up|data)\b",
    re.IGNORECASE,
)
_NONPRIMARY_COUNT_WORDS = re.compile(
    r"\b(?:previously reported|prior (?:study|report|publication)|"
    r"reported elsewhere|adapted from)\b",
    re.IGNORECASE,
)
_SURVIVAL_ESTIMATE_WORDS = re.compile(
    r"\b(?:kaplan[–-]meier|age-adjusted|time-to-event|survival estimate)\b",
    re.IGNORECASE,
)


def _source_sentences(source_text: str) -> list[str]:
    """Return source-verbatim prose units without treating whole paragraphs as claims."""
    sentences: list[str] = []
    for raw_line in source_text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("|") or re.fullmatch(r"[-:| ]+", line):
            continue
        parts = re.split(r"(?<=[.!?])\s+(?=[\"'([]*[A-Z0-9])", line)
        sentences.extend(part.strip() for part in parts if part.strip())
    return sentences


# Conservative clinical phenotype vocabulary used only to identify terms in a
# paper title.  It never creates a count.  Its purpose is to keep a verifier
# from treating a negative test for one gene-associated syndrome as a negative
# finding for a different phenotype explicitly studied by the paper.
_PAPER_TARGET_PHENOTYPES: tuple[tuple[str, re.Pattern[str]], ...] = (
    (
        "catecholaminergic polymorphic ventricular tachycardia",
        re.compile(
            r"\bcatecholaminergic polymorphic ventricular tachycardia\b|\bCPVT\b",
            re.I,
        ),
    ),
    ("long QT syndrome", re.compile(r"\blong QT syndrome\b|\bLQTS?\d*\b", re.I)),
    ("short QT syndrome", re.compile(r"\bshort QT syndrome\b|\bSQTS?\b", re.I)),
    ("Brugada syndrome", re.compile(r"\bBrugada syndrome\b", re.I)),
    ("arrhythmic storm", re.compile(r"\barrhythmic storm\b", re.I)),
    (
        "ventricular fibrillation",
        re.compile(r"\bventricular fibrillation\b|\bVF\b", re.I),
    ),
    (
        "ventricular tachycardia",
        re.compile(r"\bventricular tachycardia\b|\bVT\b", re.I),
    ),
    ("torsades de pointes", re.compile(r"\btorsade(?:s)? de pointes\b|\bTdP\b", re.I)),
    ("atrial fibrillation", re.compile(r"\batrial fibrillation\b|\bAF\b", re.I)),
    ("cardiac conduction disease", re.compile(r"\bcardiac conduction disease\b", re.I)),
    ("sick sinus syndrome", re.compile(r"\bsick sinus syndrome\b", re.I)),
    ("dilated cardiomyopathy", re.compile(r"\bdilated cardiomyopathy\b", re.I)),
    (
        "arrhythmogenic right ventricular cardiomyopathy",
        re.compile(
            r"\barrhythmogenic right ventricular cardiomyopathy\b|\bARVC\b", re.I
        ),
    ),
)


def _paper_title_from_source(source_text: str, title: str | None) -> str | None:
    """Return a real metadata title or the first non-structural Markdown heading."""
    title_text = str(title or "").strip()
    if title_text and not re.fullmatch(r"Paper\s+\d+", title_text, re.I):
        return title_text
    for raw_line in source_text.splitlines()[:120]:
        match = re.match(r"^#{1,3}\s+(.+?)\s*$", raw_line.strip())
        if not match:
            continue
        heading = match.group(1).strip()
        if heading.casefold() in {
            "main text",
            "abstract",
            "introduction",
            "background",
            "methods",
            "results",
        }:
            continue
        return heading
    return title_text or None


def _paper_target_phenotypes(title_scope: str | None) -> list[str]:
    """Extract title-bound phenotype names without inferring any patient count."""
    text = str(title_scope or "")
    return [name for name, pattern in _PAPER_TARGET_PHENOTYPES if pattern.search(text)]


def _number_is_exact(sentence: str, value: int) -> bool:
    """Require an unhedged integer occurrence in the source sentence."""
    for match in re.finditer(rf"(?<!\d){value}(?!\d)", sentence):
        if re.match(r"\s*%", sentence[match.end() :]):
            continue
        nearby = sentence[max(0, match.start() - 24) : match.end() + 12]
        if not _NONEXACT_COUNT_WORDS.search(nearby):
            return True
    return False


def _direct_carrier_binding(sentence: str, value: int) -> bool:
    number = rf"(?<!\d){value}(?!\d)"
    carrier = r"(?:variant[- ]positive\s+|mutation[- ]positive\s+)?(?:carriers?|heterozygotes?)"
    return bool(
        re.search(rf"{number}.{{0,48}}\b{carrier}\b", sentence, re.I)
        or re.search(rf"\b{carrier}\b.{{0,48}}{number}", sentence, re.I)
    )


def _direct_disease_patient_binding(
    sentence: str, value: int, disease: str | None
) -> bool:
    if not _quote_mentions_target_disease(sentence, disease):
        return False
    number = rf"(?<!\d){value}(?!\d)"
    people = r"(?:patients?|cases?|probands?|individuals?|subjects?)"
    aliases = sorted(_disease_aliases(disease), key=len, reverse=True)
    return any(
        re.search(
            rf"{number}.{{0,24}}\b{re.escape(alias)}\b.{{0,24}}\b{people}\b",
            sentence,
            re.I,
        )
        or re.search(
            rf"{number}.{{0,24}}\b{people}\b.{{0,40}}\b(?:with|diagnosed with|affected by)\s+{re.escape(alias)}\b",
            sentence,
            re.I,
        )
        for alias in aliases
    )


def _direct_target_manifestation(sentence: str, disease: str | None) -> bool:
    """Recognize a narrow disease-defining measurement stated at 100% penetrance."""

    aliases = _disease_aliases(disease)
    long_qt_context = any(
        alias == "lqts" or alias.startswith("lqt") or "long qt syndrome" in alias
        for alias in aliases
    )
    if not long_qt_context:
        return False
    folded = _source_fold(sentence)
    return bool(
        re.search(r"\b(?:prolonged\s+qtc|qtc\s+prolongation)\b", folded)
        and re.search(r"\ball\s+(?:manifested|had|showed|exhibited)\b", folded)
    )


def _population_anchors(sentence: str) -> set[str]:
    folded = _source_fold(sentence)
    anchors: set[str] = set()
    if re.search(r"\b(?:the |this |current |present )?study population\b", folded):
        anchors.add("study population")
    if re.search(r"\b(?:the |this |current |present )?study cohort\b", folded):
        anchors.add("study cohort")
    if re.search(r"\b(?:the |this |current |present )?patient series\b", folded):
        anchors.add("patient series")
    return anchors


def _same_closed_population(
    carrier_sentence: str, disease_sentence: str, value: int
) -> bool:
    if carrier_sentence == disease_sentence:
        return True
    shared_anchors = _population_anchors(carrier_sentence) & _population_anchors(
        disease_sentence
    )
    if shared_anchors:
        return True
    return bool(
        re.search(
            rf"\b(?:all|these|those|the same)\s+{value}\b", disease_sentence, re.I
        )
    )


def _scope_denial(sentence: str) -> bool:
    return bool(
        _NONPRIMARY_COUNT_WORDS.search(sentence)
        or _SURVIVAL_ESTIMATE_WORDS.search(sentence)
        or _UNASSESSED_RESIDUE_WORDS.search(sentence)
        or re.search(
            r"\b\d+\s+(?:families|kindreds|alleles|chromosomes|assays?|"
            r"replicates?|cells?|samples?|participants screened|subjects screened)\b",
            sentence,
            re.I,
        )
    )


def detect_source_bound_phenotype_claims(
    variant: dict[str, Any], source_text: str, disease: str | None
) -> list[dict[str, Any]]:
    """Detect tightly closed phenotype claims independently of emitted nulls.

    Detection is source-only and gold-free.  It deliberately returns complete
    candidate bundles so the write path can reject the whole promotion when an
    existing integer conflicts with any member of the bundle.
    """
    if not source_text or not isinstance(variant, dict) or not disease:
        return []
    aliases = _identity_aliases(variant)
    if not aliases or not _disease_aliases(disease):
        return []
    sentences = _source_sentences(source_text)
    candidates: list[dict[str, Any]] = []

    # A. The exact same closed set is explicitly called both N carriers of the
    # variant and N target-disease patients.  Symptom/event subsets do not alter
    # disease status and are never used as the affected predicate.
    total = extracted_counts(variant).get("total_carriers")
    if total is not None and total >= 1:
        carrier_sentences = [
            sentence
            for sentence in sentences
            if _quote_has_alias(sentence, aliases)
            and _number_is_exact(sentence, total)
            and _direct_carrier_binding(sentence, total)
            and not _scope_denial(sentence)
        ]
        disease_sentences = [
            sentence
            for sentence in sentences
            if _number_is_exact(sentence, total)
            and _direct_disease_patient_binding(sentence, total, disease)
            and not _scope_denial(sentence)
        ]
        coreferences = {
            (carrier_sentence, disease_sentence)
            for carrier_sentence in carrier_sentences
            for disease_sentence in disease_sentences
            if _same_closed_population(carrier_sentence, disease_sentence, total)
        }
        if coreferences:
            # A paper may restate the same carrier population in Methods and
            # Results.  Repeated support for the same N is not ambiguity; pick
            # the shortest exact pair for the audit quote.
            carrier_sentence, disease_sentence = min(
                coreferences, key=lambda pair: (len(pair[0]) + len(pair[1]), pair)
            )
            candidates.append(
                {
                    "method": "identity_disease_coreference_equal_n",
                    "values": {"total_carriers": total, "affected": total},
                    "quote": (
                        carrier_sentence
                        if carrier_sentence == disease_sentence
                        else f"{carrier_sentence} || {disease_sentence}"
                    ),
                    "quotes": [carrier_sentence, disease_sentence],
                }
            )
        manifestation_sentences: list[tuple[str, str]] = []
        for index, sentence in enumerate(sentences):
            if sentence not in carrier_sentences or not _direct_target_manifestation(
                sentence, disease
            ):
                continue
            support_window = " ".join(sentences[index : index + 3])
            if not re.search(
                r"\bpenetrance\s+(?:in\s+this\s+family\s+was\s+therefore\s+"
                r"assumed\s+to\s+be\s+|of\s+)?100\s*%",
                _source_fold(support_window),
            ):
                continue
            manifestation_sentences.append((sentence, support_window))
        if manifestation_sentences:
            sentence, support_window = min(
                manifestation_sentences, key=lambda item: (len(item[1]), item)
            )
            candidates.append(
                {
                    "method": "identity_target_manifestation_equal_n",
                    "values": {"total_carriers": total, "affected": total},
                    "quote": sentence,
                    "quotes": [sentence, support_window],
                }
            )

    # B. The source itself defines an X-of-Y target-disease partition of this
    # exact variant's carriers.  This is the only permitted complement path.
    partitions: dict[tuple[int, int], str] = {}
    for sentence in sentences:
        sentence_l = _source_fold(sentence)
        if (
            not _quote_has_alias(sentence, aliases)
            or not _quote_mentions_target_disease(sentence, disease)
            or not re.search(r"\bcarriers?\b", sentence_l)
            or not re.search(r"\b(?:affected|diagnosed|phenotype)\b", sentence_l)
            or _PHENOTYPE_EVENT_WORDS.search(sentence)
            or _ASSESSED_SUBSET_CARRIER_WORDS.search(sentence)
            or _scope_denial(sentence)
        ):
            continue
        for match in re.finditer(
            r"\b(\d+)\s+(?:out\s+of|of(?:\s+the)?)\s+(\d+)\b", sentence_l
        ):
            affected, partition_total = int(match.group(1)), int(match.group(2))
            if not (0 <= affected <= partition_total) or partition_total < 1:
                continue
            if not (
                _number_is_exact(sentence, affected)
                and _number_is_exact(sentence, partition_total)
            ):
                continue
            tail = sentence_l[match.end() :]
            if not re.search(r"\bcarriers?\b", tail) or not re.search(
                r"\b(?:were|are|was|is)\s+(?:clinically\s+)?"
                r"(?:affected\s+(?:by|with)|diagnosed\s+with)\b",
                tail,
            ):
                continue
            partitions[(affected, partition_total)] = sentence
        for match in re.finditer(
            r"\b(\d+)(?:\s+[a-z-]+){0,4}\s+"
            r"carrier(?:s|\s+(?:subjects?|patients?|individuals?))?\b"
            r".{0,220}?\b(\d+)\s+of\s+whom\b",
            sentence_l,
        ):
            partition_total, affected = int(match.group(1)), int(match.group(2))
            if not (0 <= affected <= partition_total) or partition_total < 1:
                continue
            known_total = extracted_counts(variant).get("total_carriers")
            if known_total is not None and partition_total != known_total:
                continue
            if not (
                _number_is_exact(sentence, affected)
                and _number_is_exact(sentence, partition_total)
            ):
                continue
            tail = sentence_l[match.end() : match.end() + 220]
            if not re.search(
                r"\b(?:affected|diagnosed|phenotype|phenotype-positive)\b", tail
            ):
                continue
            partitions[(affected, partition_total)] = sentence
    if len(partitions) == 1:
        (affected, partition_total), sentence = next(iter(partitions.items()))
        candidates.append(
            {
                "method": "closed_variant_disease_partition",
                "values": {
                    "total_carriers": partition_total,
                    "affected": affected,
                    "unaffected": partition_total - affected,
                },
                "quote": sentence,
                "quotes": [sentence],
            }
        )
    return candidates


def _write_promoted_count(
    variant: dict[str, Any], field: str, value: int | None
) -> None:
    key = {
        "total_carriers": "total_carriers_observed",
        "affected": "affected_count",
        "unaffected": "unaffected_count",
    }[field]
    pdata = variant.get("penetrance_data")
    if isinstance(pdata, dict):
        pdata[key] = value
        return
    flat_key = {
        "total_carriers": "carriers",
        "affected": "affected",
        "unaffected": "unaffected",
    }[field]
    variant[flat_key] = value


def promote_source_bound_phenotype_counts(
    variant: dict[str, Any], source_text: str, disease: str | None
) -> dict[str, dict[str, Any]]:
    """Verify or fill an unambiguous source-closed claim.

    The returned mapping is also consumed by the phenotype copy guard.  A
    matching, already-populated value therefore must remain visible as a
    source verification even when there is no null to promote.  Otherwise a
    valid ``N carriers == N disease patients`` claim can be detected here and
    then immediately erased by the downstream anti-copy rule.

    One overwrite is permitted: when the source independently identifies the
    same closed population as N target-disease patients, an emitted smaller
    symptomatic/asymptomatic partition is not an affected/unaffected disease
    partition.  Replace the symptom count with N affected and clear its
    complement rather than publishing symptom absence as disease absence.
    """
    candidates = detect_source_bound_phenotype_claims(variant, source_text, disease)
    if not candidates:
        return {}
    counts = extracted_counts(variant)
    uncertain = counts.get("uncertain")
    if uncertain not in {None, 0}:
        return {}

    # Prefer the explicit X-of-Y partition when it agrees with any equal-N
    # candidate.  Conflicting candidate bundles fail closed.
    unique_values = {
        tuple(sorted(candidate["values"].items())) for candidate in candidates
    }
    if len(unique_values) > 1:
        compatible = []
        for candidate in candidates:
            values = candidate["values"]
            if all(
                other["values"].get(field, value) == value
                for other in candidates
                for field, value in values.items()
            ):
                compatible.append(candidate)
        if not compatible:
            return {}
        candidate = max(compatible, key=lambda item: len(item["values"]))
    else:
        candidate = max(candidates, key=lambda item: len(item["values"]))

    desired = candidate["values"]
    corrected: dict[str, dict[str, int | None]] = {}
    conflicts = {
        field: (counts.get(field), value)
        for field, value in desired.items()
        if counts.get(field) is not None and counts[field] != value
    }
    if conflicts:
        emitted_affected = counts.get("affected")
        emitted_unaffected = counts.get("unaffected")
        emitted_total = counts.get("total_carriers")
        closed_symptom_split = (
            candidate.get("method")
            in {
                "identity_disease_coreference_equal_n",
                "identity_target_manifestation_equal_n",
            }
            and set(conflicts) == {"affected"}
            and emitted_total == desired.get("total_carriers")
            and isinstance(emitted_affected, int)
            and isinstance(emitted_unaffected, int)
            and 0 <= emitted_affected < emitted_total
            and emitted_affected + emitted_unaffected == emitted_total
            and any(
                _number_is_exact(sentence, emitted_affected)
                and _number_is_exact(sentence, emitted_unaffected)
                and _PHENOTYPE_EVENT_WORDS.search(sentence)
                for sentence in _source_sentences(source_text)
            )
        )
        if not closed_symptom_split:
            return {}
        _write_promoted_count(variant, "affected", emitted_total)
        _write_promoted_count(variant, "unaffected", None)
        corrected = {
            "affected": {"old": emitted_affected, "new": emitted_total},
            "unaffected": {"old": emitted_unaffected, "new": None},
        }
        counts = extracted_counts(variant)
        variant["source_bound_phenotype_correction"] = {
            "method": candidate["method"],
            "fields_corrected": corrected,
            "quote": candidate["quote"],
            "policy": "replace_closed_symptom_split_with_disease_status",
        }

    filled: list[str] = []
    for field, value in desired.items():
        # A zero complement is not a reported observation.  Never manufacture
        # it; a literal zero must arrive with its own source evidence.
        if counts.get(field) is None and not (field == "unaffected" and value == 0):
            _write_promoted_count(variant, field, value)
            filled.append(field)
    method = str(candidate["method"])
    if method in {
        "closed_variant_disease_partition",
        "identity_disease_coreference_equal_n",
        "identity_target_manifestation_equal_n",
    }:
        provenance = variant.setdefault("count_provenance", {})
        if isinstance(provenance, dict):
            for field in ("affected", "unaffected"):
                if field in desired and desired[field] > 0:
                    provenance[f"{field}_count_type"] = "closed_variant_partition"
                    provenance[f"{field}_source"] = SOURCE_BOUND_PHENOTYPE_SOURCE
                    provenance[f"{field}_column_label"] = str(candidate["quote"])[:1000]

    if filled:
        audit = {
            "method": method,
            "fields_filled": filled,
            "values": {field: desired[field] for field in filled},
            "quote": candidate["quote"],
            "policy": "fill_null_only",
        }
        variant["source_bound_phenotype_promotion"] = audit
    return {
        field: {
            "value": value,
            "method": method,
            "quote": candidate["quote"],
            "promoted": field in filled,
            **({"corrected": True} if field in corrected else {}),
        }
        for field, value in desired.items()
        if not (field == "unaffected" and value == 0)
    }


def source_verified_claims(
    variant: dict[str, Any], source_text: str, disease: str | None
) -> dict[str, dict[str, Any]]:
    """Find exact, identity-bound count claims in the harvested source.

    This is a keep/verification rule, not a broad count parser. It validates
    values already emitted by extraction. The sole arithmetic exception is an
    already-emitted remainder from an explicit, closed variant/disease split.
    """
    if not source_text or not isinstance(variant, dict):
        return {}
    counts = extracted_counts(variant)
    aliases = _identity_aliases(variant)
    if not aliases:
        return {}
    claims: dict[str, dict[str, Any]] = {}
    fact_types = {
        "total_carriers_observed": "total_carriers",
        "carrier_count": "total_carriers",
        "patient_count": "total_carriers",
        "affected_count": "affected",
        "unaffected_count": "unaffected",
    }
    verified_facts: list[tuple[str, int, str]] = []
    for fact in variant.get("fact_provenance") or []:
        if not isinstance(fact, dict):
            continue
        field = fact_types.get(str(fact.get("fact_type") or ""))
        value = coerce_int(fact.get("fact_value"))
        quote = str(fact.get("evidence_quote") or "").strip()
        if (
            field is None
            or value is None
            or counts.get(field) != value
            or not _quote_has_count(quote, value)
            or not _quote_is_in_source(quote, source_text)
        ):
            continue
        verified_facts.append((field, value, quote))
        has_identity = _quote_has_alias(quote, aliases)
        quote_l = _source_fold(quote)
        if (
            field == "total_carriers"
            and has_identity
            and re.search(
                r"\b(?:carrier|carriers|heterozygote|heterozygotes|mutation-positive|variant-positive|patients?)\b",
                quote_l,
            )
        ):
            claims[field] = {
                "value": value,
                "method": "verbatim_identity_bound_fact",
                "quote": quote,
            }
        elif (
            field == "affected"
            and has_identity
            and _quote_mentions_target_disease(quote, disease)
            and re.search(r"\b(?:affected|patients?|cases?|diagnosed)\b", quote_l)
        ):
            claims[field] = {
                "value": value,
                "method": "verbatim_identity_bound_fact",
                "quote": quote,
            }
        elif (
            field == "unaffected"
            and has_identity
            and re.search(
                r"\b(?:unaffected|asymptomatic|healthy|clinically negative|phenotype-negative)\b",
                quote_l,
            )
        ):
            claims[field] = {
                "value": value,
                "method": "verbatim_identity_bound_fact",
                "quote": quote,
            }

    # Same-population co-reference: one sentence identifies N carriers of this
    # variant and another identifies that exact population as N target-disease
    # patients. A disease-cohort sentence alone is not enough.
    total_facts = [
        (value, quote)
        for field, value, quote in verified_facts
        if field == "total_carriers"
        and _quote_has_alias(quote, aliases)
        and re.search(r"\bcarriers?\b", quote, flags=re.IGNORECASE)
    ]
    affected_facts = [
        (value, quote)
        for field, value, quote in verified_facts
        if field == "affected"
        and _quote_mentions_target_disease(quote, disease)
        and re.search(r"\b(?:patients?|cases?)\b", quote, flags=re.IGNORECASE)
    ]
    subset_words = re.compile(
        r"\b(?:of whom|including|among|subset|identified in)\b", re.I
    )
    population_words = re.compile(r"\b(?:current\s+)?study population\b", re.I)
    for total_value, total_quote in total_facts:
        for affected_value, affected_quote in affected_facts:
            if (
                total_value == affected_value == counts.get("affected")
                and population_words.search(total_quote)
                and population_words.search(affected_quote)
                and not subset_words.search(total_quote)
                and not subset_words.search(affected_quote)
            ):
                claims["affected"] = {
                    "value": affected_value,
                    "method": "same_population_identity_disease_coreference",
                    "quote": f"{total_quote} || {affected_quote}",
                }

    # Explicit closed target-disease partition. It validates an already
    # extracted remainder and never fills a missing unaffected value.
    for line in source_text.splitlines():
        line_l = _source_fold(line)
        if (
            not _quote_has_alias(line, aliases)
            or not _quote_mentions_target_disease(line, disease)
            or not re.search(r"\bcarriers?\b", line_l)
            or not re.search(r"\baffected\b", line_l)
        ):
            continue
        for match in re.finditer(r"\b(\d+)\s+out\s+of\s+(\d+)\b", line_l):
            affected, total = (int(match.group(1)), int(match.group(2)))
            if affected > total:
                continue
            quote = line.strip()
            if counts.get("total_carriers") == total:
                claims["total_carriers"] = {
                    "value": total,
                    "method": "closed_variant_disease_partition",
                    "quote": quote,
                }
            if counts.get("affected") == affected:
                claims["affected"] = {
                    "value": affected,
                    "method": "closed_variant_disease_partition",
                    "quote": quote,
                }
            if counts.get("unaffected") == total - affected:
                claims["unaffected"] = {
                    "value": total - affected,
                    "method": "closed_variant_disease_partition",
                    "quote": quote,
                }
    return claims


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
    variant_aliases: set[str] | None = None,
    max_chars: int = 4000,
    window: int = 2,
) -> str:
    """Select nearby lines for one variant/count claim."""
    lines = source_text.splitlines()
    gene_l = (gene or "").lower()
    variant_terms = _variant_terms(variant) | set(variant_aliases or set())
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
    count_provenance = variant.get("count_provenance")
    phenotype_derivation = variant.get("phenotype_derivation")
    derivation = None
    if isinstance(count_provenance, dict) or isinstance(phenotype_derivation, dict):
        derivation = {
            "count_provenance": (
                count_provenance if isinstance(count_provenance, dict) else {}
            ),
            "phenotype_derivation": (
                phenotype_derivation if isinstance(phenotype_derivation, dict) else None
            ),
        }
    aliases = _identity_aliases(variant)
    paper_title_scope = _paper_title_from_source(source_text, title)
    paper_target_phenotypes = _paper_target_phenotypes(paper_title_scope)
    verified_claims = source_verified_claims(variant, source_text, disease)
    evidence = build_evidence_snippet(
        source_text=source_text,
        gene=gene,
        variant=display,
        counts=counts,
        variant_aliases=aliases,
        max_chars=max_evidence_chars,
    )
    verified_lines: list[str] = []
    seen_quotes: set[str] = set()
    for field, claim in verified_claims.items():
        quote = str(claim.get("quote") or "").strip()
        if quote and quote not in seen_quotes:
            seen_quotes.add(quote)
            verified_lines.append(f"SOURCE-VERIFIED {field}: {quote}")
    if verified_lines:
        evidence = ("\n".join(verified_lines) + "\n" + evidence)[:max_evidence_chars]
        if derivation is None:
            derivation = {}
        derivation["source_verified_claims"] = verified_claims
    return VariantClaimCard(
        gene=gene,
        disease=disease,
        pmid=pmid,
        title=title,
        variant=display,
        extracted=counts,
        evidence=evidence,
        source_location=variant.get("source_location"),
        derivation=derivation,
        paper_target_phenotypes=paper_target_phenotypes,
        paper_title_scope=paper_title_scope,
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
- The claim card's `disease` field is gene-level candidate vocabulary. It can
  contain several syndromes and is not automatically the phenotype studied by
  this paper. First bind the paper-defined target from `paper_title_scope`,
  `paper_target_phenotypes`, and the source evidence. An explicit paper title,
  objective, result, or conclusion associating the variant with phenotype T
  takes precedence over unrelated gene-level aliases.
- A negative test for syndrome A cannot make a carrier unaffected for a
  different paper-defined phenotype B. For example, a carrier explicitly
  assigned ischemia-induced arrhythmic storm remains affected for that target
  even if challenge tests for Brugada and long QT syndromes are negative.
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
- Keep diagnosis separate from symptom severity. A specifically identified
  variant carrier whom the source explicitly calls a target-disease patient is
  affected for that disease even when described as asymptomatic or without
  disease-associated symptoms. Do not turn that wording into unaffected. This
  does not promote an enrollment total: the evidence must bind the diagnosis
  to this variant carrier or to the exact same closed carrier population.
- Event totals are not phenotype totals. A patient having six VT/VF episodes is
  one affected carrier when the paper assigns that carrier the target; six is
  never the affected count.
- One narrow inferred-supported path exists for a claim card whose structured
  derivation uses count_type=derived_from_patient_rows. Accept it only when a
  complete one-variant patient table, closed operational rule, explicit
  affected/unaffected/uncertain row subtotals, predicate tallies, and any
  non-overlapping add-on carriers reconcile exactly to total_carriers. Missing
  cells belong in uncertain; never derive unaffected by subtraction. Treat the
  card's structured audit as the arithmetic evidence and the source snippet as
  support for table/variant scope; do not pretend the snippet contains all rows.
- Example: "112 LQTS patients, of which 18 were symptomatic" supports
  total_carriers=112, but the affected/unaffected split is null unless the
  source explicitly defines which status is being counted.
- An unselected population, control, or blood-donor cohort is not automatically
  unaffected. Absence of a reported diagnosis is not evidence of absence; set
  unaffected only when the source explicitly calls those carriers unaffected,
  healthy, asymptomatic, or clinically negative for the target phenotype.
- A literal zero is a counted clinical claim, not an abstention. Never introduce
  affected=0 or unaffected=0 as the arithmetic complement of another field.
  Zero requires an explicit phenotype column, a source-stated closed X-of-Y
  partition, or a complete audited patient-row table. If the carrier denominator
  is unknown, use null rather than zero.
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
    # Callers may supply either the expanded disease name or only the acronym.
    # Treating ``LQTS`` as narrower than ``long QT syndrome`` caused source-
    # verified LQT1/LQT2/LQT3 cohorts to be missed in otherwise identical runs.
    if "long qt" in text or re.search(r"\blqts?\d*\b", text):
        aliases.update({"lqt", "lqts", "lqt1", "lqt2", "lqt3", "long qt syndrome"})
    if "short qt" in text or re.search(r"\bsqts?\d*\b", text):
        aliases.update({"sqts", "short qt syndrome"})
    if "brugada" in text:
        aliases.update({"brugada", "brugada syndrome"})
    if "catecholaminergic polymorphic ventricular tachycardia" in text or re.search(
        r"\bcpvt\d*\b", text
    ):
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
    if _valid_patient_row_derivation(card):
        return normalized
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


_NEGATIVE_FINDING_WORDS = re.compile(
    r"\b(?:negative|not detected|no evidence of|without|absence of|excluded)\b",
    re.I,
)


def _term_is_negated(text: str, term: str) -> bool:
    """Return whether a phenotype term is locally governed by a negative."""
    folded = _source_fold(text)
    term_re = re.escape(_source_fold(term))
    return bool(
        re.search(
            rf"(?:negative|not detected|no evidence of|without|absence of|excluded)"
            rf"[^.!?]{{0,100}}\b{term_re}\b",
            folded,
        )
        or re.search(
            rf"\b{term_re}\b[^.!?]{{0,60}}(?:was|were|tested)?\s*negative\b",
            folded,
        )
    )


def _target_is_variant_bound(card: VariantClaimCard, target: str) -> bool:
    """Require affirmative, local source support tying the paper target to the claim."""
    evidence = _source_fold(card.evidence)
    target_folded = _source_fold(target)
    if target_folded not in evidence or _term_is_negated(evidence, target):
        return False
    aliases = _variant_terms(card.variant)
    anchors = [
        *(_source_fold(alias) for alias in aliases),
        "mutation carrier",
        "variant carrier",
        "the carrier",
        "the patient",
    ]
    for match in re.finditer(re.escape(target_folded), evidence):
        window = evidence[max(0, match.start() - 280) : match.end() + 280]
        if any(anchor and anchor in window for anchor in anchors):
            return True
        if re.search(r"\b(?:mutation|variant)\b.{0,100}\bassociated with\b", window):
            return True
    return False


def _apply_paper_target_guard(
    normalized: dict[str, Any], card: VariantClaimCard | None
) -> dict[str, Any]:
    """Reject an all-affected -> all-unaffected flip based on an off-target negative.

    This is intentionally a refusal guard, not a new extraction lane. It only
    preserves a coherent original closed partition when the paper title names
    a phenotype, the local variant evidence affirmatively binds that phenotype,
    and the verifier's reversal is explicitly based on a negative finding that
    does *not* negate the title-bound phenotype.
    """
    if card is None or not card.paper_target_phenotypes:
        return normalized
    total = card.extracted.get("total_carriers")
    original_affected = card.extracted.get("affected")
    original_unaffected = card.extracted.get("unaffected")
    corrected = normalized.get("corrected_values") or {}
    if (
        total is None
        or total < 1
        or original_affected != total
        or original_unaffected != 0
        or corrected.get("affected") != 0
        or corrected.get("unaffected") != total
    ):
        return normalized

    verifier_basis = " ".join(
        [
            str(normalized.get("reason") or ""),
            str(normalized.get("evidence_quote") or ""),
        ]
    )
    if not _NEGATIVE_FINDING_WORDS.search(verifier_basis):
        return normalized

    supported_targets = [
        target
        for target in card.paper_target_phenotypes
        if _target_is_variant_bound(card, target)
        and not _term_is_negated(verifier_basis, target)
    ]
    if not supported_targets:
        return normalized

    corrected["affected"] = original_affected
    # The off-target negative cannot reverse the paper phenotype, but it also
    # cannot establish that zero carriers were unaffected for that phenotype.
    # Preserve the supported affected value and abstain on the complement.
    corrected["unaffected"] = None
    normalized["field_verdicts"]["affected"] = "directly_supported"
    normalized["field_verdicts"]["unaffected"] = "ambiguous"
    normalized["paper_target_scope_overrides"] = ["affected"]
    normalized["paper_target_scope"] = {
        "title": card.paper_title_scope,
        "phenotypes": supported_targets,
        "reason": "off_target_negative_cannot_reverse_bound_paper_phenotype",
    }
    normalized["reason"] = (
        str(normalized.get("reason") or "")
        + " Paper-target guard: a negative finding for another syndrome cannot "
        + "reverse the source-bound title phenotype."
    ).strip()
    return normalized


def _zero_has_closed_source(
    card: VariantClaimCard | None,
    field: str,
    verified_claims: dict[str, Any],
    *,
    valid_row_derivation: bool,
) -> bool:
    """Return whether a literal zero has a non-arithmetic source contract."""
    if card is None:
        return False
    claim = verified_claims.get(field) if isinstance(verified_claims, dict) else None
    if isinstance(claim, dict) and coerce_int(claim.get("value")) == 0:
        return True
    if valid_row_derivation:
        return True
    return False


def _apply_zero_claim_guards(
    normalized: dict[str, Any],
    card: VariantClaimCard | None,
    verified_claims: dict[str, Any],
    *,
    valid_row_derivation: bool,
) -> dict[str, Any]:
    """Refuse open-set or verifier-invented literal zero complements."""
    if card is None:
        return normalized
    corrected = normalized["corrected_values"]
    verdicts = normalized["field_verdicts"]
    reasons: list[str] = []
    total = corrected.get("total_carriers")
    affected = corrected.get("affected")
    verifier_basis = " ".join(
        [
            str(normalized.get("reason") or ""),
            str(normalized.get("evidence_quote") or ""),
        ]
    )
    for field in ("affected", "unaffected"):
        if corrected.get(field) != 0 or _zero_has_closed_source(
            card,
            field,
            verified_claims,
            valid_row_derivation=valid_row_derivation,
        ):
            continue
        original = card.extracted.get(field)
        introduced = original != 0
        open_set = total is None or (field == "unaffected" and affected is None)
        same_target_negative = (
            field == "affected"
            and total is not None
            and corrected.get("unaffected") == total
            and any(
                _term_is_negated(verifier_basis, target)
                for target in (card.paper_target_phenotypes or [])
            )
        )
        if same_target_negative:
            continue
        if not (introduced or open_set):
            continue
        corrected[field] = None
        verdicts[field] = "ambiguous"
        reasons.append(
            f"{field}=0 refused: "
            + (
                "the verifier introduced an unsourced zero complement"
                if introduced
                else "the assessed carrier set was not closed"
            )
        )
    if reasons:
        normalized["reason"] = (
            str(normalized.get("reason") or "")
            + " Zero-claim guard: "
            + "; ".join(reasons)
            + "."
        ).strip()
    return normalized


def _valid_patient_row_derivation(card: VariantClaimCard | None) -> bool:
    """Return whether a claim carries a closed, arithmetically valid row audit.

    This deliberately validates only the exceptional structured path. Ordinary
    prose inference still fails closed in :func:`normalize_verification`.
    """
    if card is None or not isinstance(card.derivation, dict):
        return False
    provenance = card.derivation.get("count_provenance")
    audit = card.derivation.get("phenotype_derivation")
    if not isinstance(provenance, dict) or not isinstance(audit, dict):
        return False
    derived_type = "derived_from_patient_rows"
    if any(
        str(provenance.get(key) or "").strip().lower() != derived_type
        for key in ("affected_count_type", "unaffected_count_type")
    ):
        return False
    if any(
        str(provenance.get(key) or "").strip().lower() != PATIENT_ROW_PHENOTYPE_SOURCE
        for key in ("affected_source", "unaffected_source")
    ):
        return False
    if any(
        not str(provenance.get(key) or "").strip()
        for key in ("affected_column_label", "unaffected_column_label")
    ):
        return False
    if str(audit.get("method") or "").strip().lower() != derived_type:
        return False
    if audit.get("complete_table") is not True:
        return False
    if not str(audit.get("source_table") or "").strip():
        return False
    if not str(audit.get("operational_rule") or "").strip():
        return False

    integer_keys = (
        "table_total",
        "table_affected",
        "table_unaffected",
        "table_uncertain",
        "additional_carriers",
        "additional_affected",
        "additional_unaffected",
        "additional_uncertain",
    )
    values: dict[str, int] = {}
    for key in integer_keys:
        value = audit.get(key)
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            return False
        values[key] = value
    tallies = audit.get("predicate_tallies")
    if not isinstance(tallies, dict) or not tallies:
        return False
    if any(
        isinstance(value, bool) or not isinstance(value, int) or value < 0
        for value in tallies.values()
    ):
        return False

    if (
        values["table_affected"]
        + values["table_unaffected"]
        + values["table_uncertain"]
        != values["table_total"]
    ):
        return False
    if (
        values["additional_affected"]
        + values["additional_unaffected"]
        + values["additional_uncertain"]
        != values["additional_carriers"]
    ):
        return False

    expected = {
        "total_carriers": values["table_total"] + values["additional_carriers"],
        "affected": values["table_affected"] + values["additional_affected"],
        "unaffected": values["table_unaffected"] + values["additional_unaffected"],
        "uncertain": values["table_uncertain"] + values["additional_uncertain"],
    }
    return all(card.extracted.get(field) == value for field, value in expected.items())


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
    # Phenotype partitions must be directly supported except for the one
    # structured, arithmetically checked patient-row derivation path. Inference
    # from cohort labels, missing diagnosis text, or unstructured arithmetic
    # remains forbidden by the extraction contract.
    valid_row_derivation = _valid_patient_row_derivation(card)
    for field in ("affected", "unaffected"):
        if (
            normalized["field_verdicts"].get(field) == "inferred_supported"
            and not valid_row_derivation
        ):
            normalized["field_verdicts"][field] = "ambiguous"
            normalized["corrected_values"][field] = None
        if (
            valid_row_derivation
            and normalized["corrected_values"].get(field) is not None
            and normalized["corrected_values"][field] != card.extracted.get(field)
        ):
            normalized["field_verdicts"][field] = "ambiguous"
            normalized["corrected_values"][field] = None
    normalized = _apply_consistency_guards(normalized, card)
    # The verifier cannot erase a field already rebound to an exact
    # identity/role/count quote in the harvested source.
    verified_claims = (
        card.derivation.get("source_verified_claims", {})
        if card is not None and isinstance(card.derivation, dict)
        else {}
    )
    overrides: list[str] = []
    for field in ("total_carriers", "affected", "unaffected"):
        claim = (
            verified_claims.get(field) if isinstance(verified_claims, dict) else None
        )
        if not isinstance(claim, dict):
            continue
        value = coerce_int(claim.get("value"))
        if value is None or value != card.extracted.get(field):
            continue
        normalized["corrected_values"][field] = value
        normalized["field_verdicts"][field] = (
            "inferred_supported"
            if field == "unaffected"
            and claim.get("method") == "closed_variant_disease_partition"
            else "directly_supported"
        )
        overrides.append(field)
    if overrides:
        normalized["source_verified_overrides"] = overrides
    normalized = _apply_paper_target_guard(normalized, card)
    normalized = _apply_zero_claim_guards(
        normalized,
        card,
        verified_claims if isinstance(verified_claims, dict) else {},
        valid_row_derivation=valid_row_derivation,
    )
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
