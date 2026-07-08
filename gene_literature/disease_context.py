"""Gene-disease context assembly for cold-start literature runs."""

from __future__ import annotations

import csv
import io
import json
import logging
import os
import re
import urllib.error
import urllib.request
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable

logger = logging.getLogger(__name__)

CLINGEN_GENE_VALIDITY_CSV_URL = (
    "https://search.clinicalgenome.org/kb/gene-validity/download"
)

SUPPORTIVE_CLASSIFICATIONS = {
    "Definitive",
    "Strong",
    "Moderate",
    "Limited",
}

# Curated aliases for common user-facing disease phrases that do not map cleanly
# to ClinGen's preferred disease labels.
MANUAL_DISEASE_ALIASES: dict[str, dict[str, list[str]]] = {
    "APOE": {
        "DEFAULT": [
            "Alzheimer disease",
            "Alzheimer's disease",
            "late-onset Alzheimer disease",
            "late onset Alzheimer disease",
            "LOAD",
            "Alzheimer dementia",
        ],
        "ALZHEIMER": [
            "Alzheimer disease",
            "Alzheimer's disease",
            "late-onset Alzheimer disease",
            "late onset Alzheimer disease",
            "LOAD",
            "Alzheimer dementia",
        ],
    },
    "BRCA1": {
        "DEFAULT": [
            "BRCA1-related cancer predisposition",
            "hereditary breast and ovarian cancer",
            "hereditary breast and ovarian cancer syndrome",
            "HBOC",
            "breast cancer",
            "ovarian cancer",
            "pancreatic cancer",
            "prostate cancer",
        ],
        "HEREDITARY BREAST AND OVARIAN CANCER": [
            "BRCA1-related cancer predisposition",
            "hereditary breast and ovarian cancer",
            "hereditary breast and ovarian cancer syndrome",
            "HBOC",
            "breast cancer",
            "ovarian cancer",
        ],
    },
    "MYBPC3": {
        "DEFAULT": [
            "hypertrophic cardiomyopathy",
            "familial hypertrophic cardiomyopathy",
            "hypertrophic obstructive cardiomyopathy",
            "cardiomyopathy, hypertrophic",
            "HCM",
            "HOCM",
        ],
        "HYPERTROPHIC CARDIOMYOPATHY": [
            "hypertrophic cardiomyopathy",
            "familial hypertrophic cardiomyopathy",
            "hypertrophic obstructive cardiomyopathy",
            "cardiomyopathy, hypertrophic",
            "HCM",
            "HOCM",
        ],
    },
}


@dataclass
class ClinGenDiseaseCuration:
    gene_symbol: str
    disease_label: str
    disease_id: str
    mode_of_inheritance: str
    classification: str
    classification_date: str
    gcep: str
    online_report: str


@dataclass
class GeneDiseaseContext:
    gene_symbol: str
    requested_disease: str | None
    disease_terms: list[str] = field(default_factory=list)
    prompt_disease: str | None = None
    clinigen_gene_validity_url: str = CLINGEN_GENE_VALIDITY_CSV_URL
    clinigen_curations: list[ClinGenDiseaseCuration] = field(default_factory=list)
    selected_clinigen_curations: list[ClinGenDiseaseCuration] = field(
        default_factory=list
    )
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2) + "\n", encoding="utf-8")


def build_gene_disease_context(
    gene_symbol: str,
    disease: str | None,
    *,
    include_all_clinigen_phenotypes: bool | None = None,
    timeout_s: int = 30,
) -> GeneDiseaseContext:
    """Build disease aliases and ClinGen gene-validity context for a run."""

    gene = gene_symbol.upper()
    requested = disease.strip() if disease and disease.strip() else None
    if include_all_clinigen_phenotypes is None:
        include_all_clinigen_phenotypes = os.environ.get(
            "GVF_INCLUDE_ALL_CLINGEN_PHENOTYPES", ""
        ).strip().lower() in {"1", "true", "yes", "on", "all"}
    warnings: list[str] = []

    curations: list[ClinGenDiseaseCuration] = []
    try:
        curations = fetch_clinigen_gene_validity_curations(gene, timeout_s=timeout_s)
    except Exception as exc:
        warnings.append(f"ClinGen gene-validity fetch failed: {exc}")
        logger.warning("ClinGen gene-validity fetch failed for %s: %s", gene, exc)

    terms: list[str] = []
    if requested:
        terms.append(requested)

    terms.extend(_manual_aliases_for(gene, requested))

    selected_curations = _select_curations_for_requested_disease(
        curations, requested, terms
    )

    # If no disease was supplied, collect every ClinGen disease label for the
    # gene into the run context. Some labels may later prove weak, disputed, or
    # clinically irrelevant, but cold-start discovery should preserve the full
    # phenotype surface before downstream adjudication narrows it.
    if include_all_clinigen_phenotypes or not requested:
        selected_curations = _dedupe_curations([*selected_curations, *curations])
        terms.extend(c.disease_label for c in selected_curations)
    else:
        # With an explicit disease, unrelated ClinGen rows stay recorded but do
        # not broaden the query/filter context.
        terms.extend(
            c.disease_label
            for c in selected_curations
            if c.classification in SUPPORTIVE_CLASSIFICATIONS
        )

    disease_terms = _dedupe_terms(terms)
    prompt_disease = _build_prompt_disease(disease_terms)

    return GeneDiseaseContext(
        gene_symbol=gene,
        requested_disease=requested,
        disease_terms=disease_terms,
        prompt_disease=prompt_disease,
        clinigen_curations=curations,
        selected_clinigen_curations=selected_curations,
        warnings=warnings,
    )


def fetch_clinigen_gene_validity_curations(
    gene_symbol: str,
    *,
    timeout_s: int = 30,
) -> list[ClinGenDiseaseCuration]:
    """Fetch ClinGen Gene-Disease Validity rows for one gene."""

    with urllib.request.urlopen(CLINGEN_GENE_VALIDITY_CSV_URL, timeout=timeout_s) as r:
        text = r.read().decode("utf-8-sig", errors="replace")
    return parse_clinigen_gene_validity_csv(text, gene_symbol)


def parse_clinigen_gene_validity_csv(
    csv_text: str, gene_symbol: str
) -> list[ClinGenDiseaseCuration]:
    """Parse ClinGen's decorated gene-validity CSV export."""

    gene = gene_symbol.upper()
    lines = csv_text.splitlines()
    header_idx = None
    for idx, line in enumerate(lines):
        if "GENE SYMBOL" in line and "DISEASE LABEL" in line:
            header_idx = idx
            break
    if header_idx is None:
        raise ValueError("ClinGen gene-validity CSV header not found")

    data_lines = [lines[header_idx]]
    for line in lines[header_idx + 1 :]:
        if not line or set(line.replace('"', "").replace(",", "")) <= {"+"}:
            continue
        data_lines.append(line)

    rows: list[ClinGenDiseaseCuration] = []
    reader = csv.DictReader(io.StringIO("\n".join(data_lines)))
    for row in reader:
        row_gene = (row.get("GENE SYMBOL") or "").strip().upper()
        if row_gene != gene:
            continue
        rows.append(
            ClinGenDiseaseCuration(
                gene_symbol=row_gene,
                disease_label=(row.get("DISEASE LABEL") or "").strip(),
                disease_id=(row.get("DISEASE ID (MONDO)") or "").strip(),
                mode_of_inheritance=(row.get("MOI") or "").strip(),
                classification=(row.get("CLASSIFICATION") or "").strip(),
                classification_date=(row.get("CLASSIFICATION DATE") or "").strip(),
                gcep=(row.get("GCEP") or "").strip(),
                online_report=(row.get("ONLINE REPORT") or "").strip(),
            )
        )
    return rows


def _manual_aliases_for(gene: str, requested: str | None) -> list[str]:
    gene_aliases = MANUAL_DISEASE_ALIASES.get(gene, {})
    aliases: list[str] = []
    if requested:
        req_norm = _norm(requested)
        for key, values in gene_aliases.items():
            if key == "DEFAULT":
                continue
            key_norm = _norm(key)
            if key_norm in req_norm or req_norm in key_norm:
                aliases.extend(values)
    aliases.extend(gene_aliases.get("DEFAULT", []))
    return aliases


def _select_curations_for_requested_disease(
    curations: Iterable[ClinGenDiseaseCuration],
    requested: str | None,
    terms: Iterable[str],
) -> list[ClinGenDiseaseCuration]:
    if not requested and not list(terms):
        return []
    norm_terms = [_norm(t) for t in terms if t]
    selected: list[ClinGenDiseaseCuration] = []
    for curation in curations:
        disease_norm = _norm(curation.disease_label)
        if any(_terms_match(disease_norm, term_norm) for term_norm in norm_terms):
            selected.append(curation)
    return selected


def _dedupe_curations(
    curations: Iterable[ClinGenDiseaseCuration],
) -> list[ClinGenDiseaseCuration]:
    seen: set[tuple[str, str, str, str]] = set()
    deduped: list[ClinGenDiseaseCuration] = []
    for curation in curations:
        key = (
            curation.gene_symbol,
            _norm(curation.disease_label),
            curation.disease_id,
            curation.classification,
        )
        if key in seen:
            continue
        seen.add(key)
        deduped.append(curation)
    return deduped


def _terms_match(left: str, right: str) -> bool:
    if not left or not right:
        return False
    if left in right or right in left:
        return True
    left_tokens = set(left.split())
    right_tokens = set(right.split())
    informative = {
        token
        for token in left_tokens & right_tokens
        if len(token) >= 5
        and token not in {"disease", "syndrome", "cancer", "cardiomyopathy", "disorder"}
    }
    return bool(informative)


def _build_prompt_disease(terms: list[str]) -> str | None:
    if not terms:
        return None
    if len(terms) == 1:
        return terms[0]
    return f"{terms[0]} (aliases: {', '.join(terms[1:])})"


def _dedupe_terms(terms: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for term in terms:
        cleaned = " ".join(str(term).split())
        if not cleaned:
            continue
        key = _norm(cleaned)
        if key in seen:
            continue
        seen.add(key)
        out.append(cleaned)
    return out


def _norm(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", text.lower()).strip()
