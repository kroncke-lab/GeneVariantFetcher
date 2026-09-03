"""Audited phenotype counts from complete, variant-specific patient tables.

This is deliberately a narrow exception to the usual source-bound count rule.
It never treats an omitted diagnosis as unaffected and never completes a
partition by subtraction.  A table is eligible only when it enumerates a
variant-scoped carrier cohort with unique contiguous patient IDs and the prose
independently confirms both the table population and carrier total.
"""

from __future__ import annotations

import copy
import re
from typing import Any

from pipeline.count_provenance import PATIENT_ROW_PHENOTYPE_SOURCE


_MISSING = {
    "",
    "-",
    "--",
    "—",
    "na",
    "n/a",
    "nan",
    "not available",
    "not assessed",
    "not applicable",
    "not done",
    "not performed",
    "not recorded",
    "unknown",
}
_NEGATIVE = {
    "no",
    "none",
    "negative",
    "normal",
    "absent",
    "asymptomatic",
    "unaffected",
    "0",
}
_POSITIVE_RE = re.compile(
    r"\b(?:yes|positive|present|abnormal|affected|symptomatic|syncope|"
    r"dizz(?:y|iness)|palpitations?|seizures?|cardiac arrest|sudden (?:cardiac )?death|"
    r"arrhythmias?|tachycardia|fibrillation|ectopy)\b",
    re.IGNORECASE,
)
_SYMPTOM_HEADER_RE = re.compile(
    r"\b(?:previous |prior |clinical )?symptoms?\b|"
    r"^symptomatic(?:\s+at\s+.+)?$|^cpvt[ -]?syncope$",
    re.IGNORECASE,
)
_OBJECTIVE_HEADER_RE = re.compile(
    r"\b(?:ventricular\s+arrhythmias?|arrhythmias?|va|"
    r"aborted\s+cardiac\s+arrest|aca|sudden\s+cardiac\s+death|scd)\b",
    re.IGNORECASE,
)
_PATIENT_HEADER_RE = re.compile(
    r"^(?:patient|subject|participant|individual)(?:\s+id)?$", re.I
)
_ELIGIBLE_CAPTION_RE = re.compile(
    r"\b(?:phenotyp\w*|clinical)\b.*\b(?:mutation|variant)[ -]?positive\b|"
    r"\b(?:mutation|variant)[ -]?positive\b.*\b(?:phenotyp\w*|clinical)\b",
    re.IGNORECASE,
)
_DIVIDER_RE = re.compile(r"^:?-{3,}:?$")
PROTOCOL_VERSION = PATIENT_ROW_PHENOTYPE_SOURCE
_AA3_TO_1 = {
    "ala": "A",
    "arg": "R",
    "asn": "N",
    "asp": "D",
    "cys": "C",
    "gln": "Q",
    "glu": "E",
    "gly": "G",
    "his": "H",
    "ile": "I",
    "leu": "L",
    "lys": "K",
    "met": "M",
    "phe": "F",
    "pro": "P",
    "ser": "S",
    "thr": "T",
    "trp": "W",
    "tyr": "Y",
    "val": "V",
}


def _clean_cell(value: str) -> str:
    value = re.sub(r"[*_`]", "", value)
    return re.sub(r"\s+", " ", value).strip()


def _cells(line: str) -> list[str]:
    return [_clean_cell(cell) for cell in line.strip().strip("|").split("|")]


def _is_table_line(line: str) -> bool:
    stripped = line.strip()
    return stripped.startswith("|") and stripped.endswith("|")


def _is_divider_row(cells: list[str]) -> bool:
    return bool(cells) and all(
        not cell or _DIVIDER_RE.fullmatch(cell) for cell in cells
    )


def _caption_before(lines: list[str], start: int) -> str:
    for index in range(start - 1, max(-1, start - 6), -1):
        value = _clean_cell(lines[index])
        if not value:
            continue
        if _is_table_line(lines[index]):
            break
        return value.lstrip("# ")
    return ""


def _markdown_tables(source_text: str) -> list[dict[str, Any]]:
    lines = source_text.splitlines()
    tables: list[dict[str, Any]] = []
    index = 0
    while index < len(lines):
        if not _is_table_line(lines[index]):
            index += 1
            continue
        start = index
        block: list[list[str]] = []
        while index < len(lines) and _is_table_line(lines[index]):
            block.append(_cells(lines[index]))
            index += 1
        if len(block) < 3:
            continue

        header_index: int | None = None
        # Ordinary markdown puts the divider after the header. Some PMC
        # supplements first emit an empty spanning row, then a divider, then
        # the actual bold header; accept both layouts.
        for row_index, row in enumerate(block):
            if any(_PATIENT_HEADER_RE.fullmatch(cell) for cell in row):
                header_index = row_index
                break
        if header_index is None:
            continue
        header = block[header_index]
        data_rows = [
            row
            for row in block[header_index + 1 :]
            if not _is_divider_row(row) and any(cell for cell in row)
        ]
        tables.append(
            {
                "caption": _caption_before(lines, start),
                "header": header,
                "rows": data_rows,
            }
        )
    return tables


def _cell_status(value: str) -> str:
    normalized = _clean_cell(value).lower().rstrip(".")
    if normalized in _MISSING:
        return "missing"
    if normalized in _NEGATIVE or re.fullmatch(r"no(?:\s+symptoms?)?", normalized):
        return "negative"
    if _POSITIVE_RE.search(normalized):
        return "positive"
    return "unmapped"


def _slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def _selected_status_columns(
    header: list[str], rows: list[list[str]]
) -> list[tuple[int, str]]:
    symptom = [
        (index, label)
        for index, label in enumerate(header)
        if _SYMPTOM_HEADER_RE.search(label)
    ]
    objective = [
        (index, label)
        for index, label in enumerate(header)
        if _OBJECTIVE_HEADER_RE.search(label)
    ]
    # The exception is intentionally limited to tables containing both a
    # historical symptom field and an objective disease-phenotype field.
    candidates = symptom + [item for item in objective if item not in symptom]
    if not symptom or not objective:
        return []

    selected: list[tuple[int, str]] = []
    for column_index, label in candidates:
        statuses = [
            _cell_status(row[column_index] if column_index < len(row) else "")
            for row in rows
        ]
        if "unmapped" in statuses:
            continue
        if "positive" not in statuses or "negative" not in statuses:
            continue
        selected.append((column_index, label))
    if not any(item in symptom for item in selected) or not any(
        item in objective for item in selected
    ):
        return []
    return selected


def _table_population_confirmed(source_text: str, table_total: int) -> bool:
    """Require an independent statement describing the table population."""
    patterns = (
        rf"\b{table_total}\s+living\b.{{0,120}}\b(?:mutation|variant)[ -]?positive\b",
        rf"\b{table_total}\s+living\b.{{0,120}}\b(?:mutation|variant)-positive\b",
        rf"\bidentified\b.{{0,160}}\b{table_total}\s+living\b.{{0,160}}\b(?:mutation|variant)[ -]?positive\b",
        rf"\b{table_total}\s+(?:subjects|patients|individuals|carriers)\s+"
        r"(?:were\s+)?(?:eligible|included|analysed|analyzed|enrolled|participated)\b",
    )
    return any(re.search(pattern, source_text, re.I | re.S) for pattern in patterns)


def _variant_terms(variant: dict[str, Any]) -> list[str]:
    terms: list[str] = []
    for key in ("source_notation", "protein_notation", "cdna_notation"):
        value = str(variant.get(key) or "").strip()
        if value:
            terms.extend((value, value.removeprefix("p.")))
    gene = str(variant.get("gene_symbol") or "").strip()
    if gene:
        terms.append(gene)
    return sorted(set(terms), key=len, reverse=True)


def _normalized_identity(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def _identity_aliases(variant: dict[str, Any]) -> set[str]:
    aliases: set[str] = set()
    for term in _variant_terms(variant):
        normalized = _normalized_identity(term)
        if normalized:
            aliases.add(normalized)
            aliases.add(normalized.removeprefix("p"))
        for match in re.finditer(
            r"(?<![A-Za-z0-9])p?\.?([A-Z][a-z]{2}|[A-Z])(\d+)"
            r"([A-Z][a-z]{2}|[A-Z])(?=$|[^A-Za-z0-9])",
            term,
        ):
            ref, position, alt = match.groups()
            ref_one = _AA3_TO_1.get(ref.lower(), ref.upper())
            alt_one = _AA3_TO_1.get(alt.lower(), alt.upper())
            short = f"{ref_one}{position}{alt_one}".lower()
            aliases.update((short, "p" + short))
    return {alias for alias in aliases if len(alias) >= 4}


def _text_identifies_variant(text: str, variant: dict[str, Any]) -> bool:
    normalized_text = _normalized_identity(text)
    return any(identity in normalized_text for identity in _identity_aliases(variant))


def _carrier_total_from_source(
    source_text: str, variant: dict[str, Any], *, minimum: int
) -> tuple[int | None, str]:
    """Recover one exact carrier total from direct carrier language.

    This deliberately ignores screening denominators.  In multi-variant papers,
    the matched sentence must also identify the target variant.
    """
    patterns = (
        re.compile(
            r"\b(?:a\s+)?total\s+of\s+(\d{1,5})\s+(?:living\s+)?"
            r"(?:genotyped\s+)?(?:mutation|variant)[ -]?positive\s+"
            r"(?:individuals|subjects|patients|carriers)\b",
            re.I,
        ),
        re.compile(
            r"\b(?:a\s+)?total\s+of\s+(\d{1,5})\s+(?:living\s+)?"
            r"(?:genotyped\s+)?(?:mutation|variant)\s+carriers\b",
            re.I,
        ),
        re.compile(
            r"\b(\d{1,5})\s+(?:living\s+)?(?:genotyped\s+)?"
            r"(?:mutation|variant)[ -]?positive\s+"
            r"(?:individuals|subjects|patients|carriers)\b",
            re.I,
        ),
        re.compile(
            r"\b(\d{1,5})\s+(?:were\s+)?carriers\s+of\s+(?:the\s+)?"
            r"(?:mutation|variant)\b",
            re.I,
        ),
    )
    candidates: dict[int, str] = {}
    for pattern in patterns:
        for match in pattern.finditer(source_text):
            count = int(match.group(1))
            if count < minimum:
                continue
            start = max(0, match.start() - 180)
            end = min(len(source_text), match.end() + 180)
            window = re.sub(r"\s+", " ", source_text[start:end]).strip()
            if _text_identifies_variant(window, variant):
                candidates[count] = window
    if len(candidates) != 1:
        return None, ""
    return next(iter(candidates.items()))


def _table_is_variant_scoped(
    table: dict[str, Any],
    variant: dict[str, Any],
    *,
    paper_title: str,
    paper_variant_count: int,
) -> bool:
    caption = str(table.get("caption") or "")
    if _text_identifies_variant(caption, variant):
        return True
    if not _ELIGIBLE_CAPTION_RE.search(caption):
        # Spreadsheet supplements often expose only ``Sheet: Data``.  That is
        # usable only when the paper itself is unambiguously about one variant.
        return paper_variant_count == 1 and _text_identifies_variant(
            paper_title, variant
        )
    # A generic "mutation-positive phenotype" caption can bind implicitly only
    # when there is exactly one extracted target variant in the paper.
    return paper_variant_count == 1


def _fatal_additional_carriers(
    source_text: str, additional: int, variant: dict[str, Any]
) -> tuple[int, str]:
    """Return an exact additional affected count only for genotyped fatal cases."""
    if additional <= 0:
        return 0, ""
    fatal_re = re.compile(
        rf"\b{additional}\s+(?:genotyped\s+)?(?:SCD|sudden cardiac death)\s+"
        r"(?:cases|victims|individuals|subjects)\b",
        re.IGNORECASE,
    )
    terms = _variant_terms(variant)
    for match in fatal_re.finditer(source_text):
        start = max(0, match.start() - 260)
        end = min(len(source_text), match.end() + 260)
        window = re.sub(r"\s+", " ", source_text[start:end]).strip()
        linked = any(term.lower() in window.lower() for term in terms)
        genotype_language = re.search(
            r"\b(?:mutation|variant|genetic analysis|identified|positive|had)\b",
            window,
            re.IGNORECASE,
        )
        if linked and genotype_language:
            return additional, window[:700]
    return 0, ""


def _derive_table_partition(
    table: dict[str, Any],
    source_text: str,
    variant: dict[str, Any],
    *,
    paper_title: str,
    paper_variant_count: int,
) -> dict[str, Any] | None:
    caption = str(table.get("caption") or "")
    if not _table_is_variant_scoped(
        table,
        variant,
        paper_title=paper_title,
        paper_variant_count=paper_variant_count,
    ):
        return None
    header = table["header"]
    rows = table["rows"]
    try:
        patient_index = next(
            index
            for index, label in enumerate(header)
            if _PATIENT_HEADER_RE.fullmatch(label)
        )
    except StopIteration:
        return None

    patient_ids: list[int] = []
    for row in rows:
        if patient_index >= len(row) or not re.fullmatch(r"\d+", row[patient_index]):
            return None
        patient_ids.append(int(row[patient_index]))
    if not patient_ids or sorted(patient_ids) != list(range(1, len(patient_ids) + 1)):
        return None
    table_total = len(rows)
    if not _table_population_confirmed(source_text, table_total):
        return None

    selected = _selected_status_columns(header, rows)
    if not selected:
        return None
    affected = unaffected = uncertain = 0
    predicate_tallies = {f"{_slug(label)}_positive": 0 for _, label in selected}
    for row in rows:
        statuses: list[str] = []
        for column_index, label in selected:
            status = _cell_status(row[column_index] if column_index < len(row) else "")
            statuses.append(status)
            if status == "positive":
                predicate_tallies[f"{_slug(label)}_positive"] += 1
        if "positive" in statuses:
            affected += 1
        elif statuses and all(status == "negative" for status in statuses):
            unaffected += 1
        else:
            uncertain += 1
    predicate_tallies.update(
        {
            "any_selected_positive": affected,
            "all_selected_negative": unaffected,
            "missing_without_positive": uncertain,
        }
    )
    return {
        "caption": caption,
        "selected_labels": [label for _, label in selected],
        "table_total": table_total,
        "table_affected": affected,
        "table_unaffected": unaffected,
        "table_uncertain": uncertain,
        "predicate_tallies": predicate_tallies,
    }


def derive_patient_row_phenotype_counts(
    extracted_data: dict[str, Any], source_text: str
) -> dict[str, Any]:
    """Populate an audited phenotype partition when the strict table gate passes."""
    variants = extracted_data.get("variants") or []
    if not isinstance(variants, list):
        return extracted_data
    result = copy.deepcopy(extracted_data)
    tables = _markdown_tables(source_text)
    paper_title = str((result.get("paper_metadata") or {}).get("title") or "")
    outcomes: list[dict[str, Any]] = []
    newly_applied = 0
    derived = 0
    for variant_index, target in enumerate(result.get("variants") or []):
        if not isinstance(target, dict):
            outcomes.append(
                {"variant_index": variant_index, "status": "invalid_variant"}
            )
            continue
        identity = next(
            (
                str(target.get(key))
                for key in ("source_notation", "protein_notation", "cdna_notation")
                if target.get(key)
            ),
            f"variant_{variant_index + 1}",
        )
        penetrance = target.get("penetrance_data") or {}
        if (
            penetrance.get("affected_count") is not None
            or penetrance.get("unaffected_count") is not None
        ):
            existing_derivation = target.get("phenotype_derivation") or {}
            if (
                existing_derivation.get("method") == "derived_from_patient_rows"
                and existing_derivation.get("complete_table") is True
            ):
                derived += 1
                outcomes.append({"variant": identity, "status": "already_derived"})
            else:
                outcomes.append(
                    {"variant": identity, "status": "phenotype_already_populated"}
                )
            continue

        partitions = [
            partition
            for table in tables
            if (
                partition := _derive_table_partition(
                    table,
                    source_text,
                    target,
                    paper_title=paper_title,
                    paper_variant_count=len(variants),
                )
            )
            is not None
        ]
        if len(partitions) != 1:
            outcomes.append(
                {
                    "variant": identity,
                    "status": (
                        "no_eligible_complete_patient_table"
                        if not partitions
                        else "ambiguous_multiple_patient_tables"
                    ),
                    "eligible_table_count": len(partitions),
                }
            )
            continue
        partition = partitions[0]
        table_total = partition["table_total"]
        total = penetrance.get("total_carriers_observed")
        total_quote = ""
        if isinstance(total, bool) or not isinstance(total, int) or total <= 0:
            total, total_quote = _carrier_total_from_source(
                source_text, target, minimum=table_total
            )
        if total is None or total < table_total:
            outcomes.append(
                {
                    "variant": identity,
                    "status": "missing_or_inconsistent_carrier_total",
                    "table_total": table_total,
                    "carrier_total": total,
                }
            )
            continue

        additional = total - table_total
        additional_affected, additional_quote = _fatal_additional_carriers(
            source_text, additional, target
        )
        additional_uncertain = additional - additional_affected
        target_penetrance = target.setdefault("penetrance_data", {})
        target_penetrance["total_carriers_observed"] = total
        target_penetrance["affected_count"] = (
            partition["table_affected"] + additional_affected
        )
        target_penetrance["unaffected_count"] = partition["table_unaffected"]
        target_penetrance["uncertain_count"] = (
            partition["table_uncertain"] + additional_uncertain
        )

        labels = " | ".join(partition["selected_labels"])
        provenance = target.setdefault("count_provenance", {})
        provenance.update(
            {
                "carriers_column_label": provenance.get("carriers_column_label"),
                "carriers_count_type": "per_variant_carrier",
                "affected_column_label": labels,
                "affected_count_type": "derived_from_patient_rows",
                "affected_source": PROTOCOL_VERSION,
                "unaffected_column_label": labels,
                "unaffected_count_type": "derived_from_patient_rows",
                "unaffected_source": PROTOCOL_VERSION,
            }
        )
        rule = (
            "Affected when any selected phenotype cell is positive; unaffected only "
            "when every selected phenotype cell is explicitly negative; missing or "
            "unmapped cells without another positive are uncertain. Exact genotyped "
            "fatal cases outside the enumerated table are affected."
        )
        target["phenotype_derivation"] = {
            "protocol_version": PROTOCOL_VERSION,
            "method": "derived_from_patient_rows",
            "source_table": partition["caption"],
            "operational_rule": rule,
            "complete_table": True,
            "table_total": table_total,
            "table_affected": partition["table_affected"],
            "table_unaffected": partition["table_unaffected"],
            "table_uncertain": partition["table_uncertain"],
            "additional_carriers": additional,
            "additional_affected": additional_affected,
            "additional_unaffected": 0,
            "additional_uncertain": additional_uncertain,
            "predicate_tallies": partition["predicate_tallies"],
        }
        facts = target.setdefault("fact_provenance", [])
        if not isinstance(facts, list):
            facts = []
            target["fact_provenance"] = facts
        audit_quote = (
            f"{partition['caption']}; {table_total} contiguous patient rows; {rule}"
        )
        if total_quote:
            audit_quote += f" Carrier-total evidence: {total_quote}"
        if additional_quote:
            audit_quote += f" Additional-case evidence: {additional_quote}"
        for fact_type, fact_value in (
            ("affected_count", target_penetrance["affected_count"]),
            ("unaffected_count", target_penetrance["unaffected_count"]),
        ):
            facts.append(
                {
                    "fact_type": fact_type,
                    "fact_value": str(fact_value),
                    "source_location": partition["caption"],
                    "source_table": partition["caption"],
                    "source_row": "complete patient rows 1-" + str(table_total),
                    "source_column": labels,
                    "evidence_quote": audit_quote[:2000],
                }
            )
        newly_applied += 1
        derived += 1
        outcomes.append(
            {
                "variant": identity,
                "status": "applied",
                "source_table": partition["caption"],
                "table_total": table_total,
                "carrier_total": total,
            }
        )

    metadata = result.setdefault("extraction_metadata", {})
    metadata["patient_row_phenotype_derivation"] = {
        "protocol_version": PROTOCOL_VERSION,
        "attempted": True,
        "applied": bool(derived),
        "paper_variant_count": len(variants),
        "parsed_patient_table_count": len(tables),
        "applied_variant_count": derived,
        "newly_applied_variant_count": newly_applied,
        "outcomes": outcomes,
    }
    return result
