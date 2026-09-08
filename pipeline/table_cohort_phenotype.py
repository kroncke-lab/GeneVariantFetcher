"""Project a code-classified table cohort onto affected / unaffected counts.

The deterministic table parsers can read a per-variant people count without
assigning a phenotype. A case-series caption or explicit case count can supply
that context. A patient/proband noun by itself cannot: mutation catalogues may
include several phenotypes or observations compiled from other publications.

This module closes that gap without a model call and without arithmetic on
other counts. It classifies the SOURCE TABLE's cohort from what the paper
itself printed -- the caption (looked up in the source text when the parser
stored only ``Table N``) and the count-column header -- and only then projects
the row's people count onto the phenotype class that cohort defines:

* ``case``    -> ``affected = N``; ``unaffected`` stays NULL because the paper
  never assessed anyone else, and a zero it did not count is not emitted.
* ``control`` -> ``unaffected = N`` and, when the count is people rather than
  alleles, ``affected = 0``: the cohort is defined as unaffected, so the
  partition is closed by the study design rather than derived.

Tables that signal a different population refuse: family / relative /
segregation tables, mutation-carrier characteristic tables, phenotype-split
columns, population or biobank screening, allele and frequency columns,
functional catalogues, per-person clinical rows, and captions that mix cases
with controls. Paper-level ascertainment (title / abstract) is a separate
tier that is OFF by default: reviewers judged that projecting the abstract's
cohort onto a table that names no cohort is the manufactured partition the
extraction contract forbids.

Every derived value is distinguishable from a literal phenotype column. It
carries ``count_type`` ``case`` / ``control`` / ``unaffected_control``, the
code-owned source stamp ``table_cohort_phenotype_v1`` (scrubbed from model
output before this runs), a ``phenotype_derivation`` audit block quoting the
caption and column, and ``fact_provenance`` rows.

Public API::

    derive_table_cohort_phenotype_counts(extracted_data, source_text, ...)
    classify_table_cohort(caption, count_label, headers, ...) -> TableCohort
"""

from __future__ import annotations

import copy
import re
from dataclasses import dataclass, field
from typing import Any, Iterable, Optional

from pipeline.count_provenance import TABLE_COHORT_PHENOTYPE_SOURCE
from pipeline.phenotype_count_guard import phenotype_fields_to_clear

PROTOCOL_VERSION = TABLE_COHORT_PHENOTYPE_SOURCE
METHOD = "derived_from_table_cohort"
IMPLICIT_ROW_CARRIER_LABEL = "implicit one carrier per clinical row"

# Tiers, strongest evidence first. Only ``paper_ascertainment`` reads outside
# the table itself, and it is disabled unless the setting turns it on.
TIER_COLUMN = "count_column_names_cases"
TIER_CAPTION = "caption_names_disease_cases"
TIER_CAPTION_DISEASE = "caption_names_disease_column_counts_people"
TIER_PAPER = "paper_ascertainment"
# The paper's own title defines a disease proband/patient cohort ("Probands
# With Brugada Syndrome"). Stronger than a sentence anywhere in the head of the
# text, and the only paper-level evidence that may certify one-row-per-person
# mutation catalogues.
TIER_TITLE = "title_ascertainment"
TABLE_LOCAL_TIERS = frozenset({TIER_COLUMN, TIER_CAPTION, TIER_CAPTION_DISEASE})

_ROUTER_LOCATION_SUFFIX_RE = re.compile(
    r",\s*row\s+\d+\s*\(router\+deterministic\)\s*$", re.IGNORECASE
)
_ROUTER_TABLE_ID_RE = re.compile(r"^Table\s+(T\d+)$", re.IGNORECASE)
_TABLE_LABEL_RE = re.compile(
    r"^(?:#+\s*)?(?:\*\*)?\s*(?:e|supplementa\w+\s+|supp\.?\s+|online\s+)?table\s+"
    r"[a-z]?\d+[a-z]?(?:\s*[\.:]|\s|$)",
    re.IGNORECASE,
)
# "Table 2. Continued" / "Table 2 (cont.)": a PDF page break, not a caption.
# The rows under it belong to the table whose caption was printed first.
_CONTINUED_CAPTION_RE = re.compile(
    r"^(?:#+\s*)?(?:\*\*)?\s*((?:e|supplementa\w+\s+|supp\.?\s+|online\s+)?table\s+"
    r"[a-z]?\d+[a-z]?)\s*[\.:]?\s*[\(\[]?\s*cont(?:inued|d|\.)?\s*[\)\]]?\s*\.?\s*(?:\*\*)?\s*$",
    re.IGNORECASE,
)
_BARE_TABLE_LABEL_RE = re.compile(
    r"^(?:#+\s*)?(?:\*\*)?\s*(?:e|supplementa\w+\s+|supp\.?\s+|online\s+)?table\s+"
    r"[a-z]?\d+[a-z]?\s*[\.:]?\s*(?:\*\*)?\s*$",
    re.IGNORECASE,
)

# The people the lab counts as affected in a case series.
_CASE_NOUN_RE = re.compile(
    r"\b(?:patients?|probands?|cases?|index[- ]cases?|index[- ]patients?)\b",
    re.IGNORECASE,
)
# Weaker people nouns: acceptable only when the caption names the disease.
_PEOPLE_NOUN_RE = re.compile(
    r"\b(?:individuals?|subjects?|persons?|people|participants?)\b", re.IGNORECASE
)
_CONTROL_RE = re.compile(
    r"\b(?:controls?|healthy|ostensibly\s+healthy|reference\s+alleles?|"
    r"reference\s+(?:individuals?|subjects?|population|samples?)|"
    r"normal\s+(?:individuals?|subjects?|volunteers?|controls?))\b",
    re.IGNORECASE,
)
_ALLELE_UNIT_RE = re.compile(r"\balleles?\b|\bchromosomes?\b", re.IGNORECASE)

# Cardiac and oncology disease names plus the generic "syndrome". Sudden-death
# vocabulary is deliberately absent: a molecular-autopsy victim is not a
# diagnosed case of the target channelopathy (RYR2 22677073).
_BUILTIN_DISEASE_RE = re.compile(
    r"\b(?:brugada|brs\d?|long[- ]qt|lqts?\d*|romano[- ]ward|jervell|"
    r"cpvt\d?|catecholaminergic|sqts?\d?|short[- ]qt|"
    r"hypertrophic\s+cardiomyopathy|hcm|dilated\s+cardiomyopathy|dcm|"
    r"arvc|arvd|arrhythmogenic|cardiomyopathy|conduction\s+(?:disease|defect|disorder)|"
    r"sick\s+sinus|atrial\s+fibrillation|"
    r"pulmonary\s+arterial\s+hypertension|pah|"
    r"breast\s+cancer|ovarian\s+cancer|cancer|carcinoma|tumou?rs?|syndrome)\b",
    re.IGNORECASE,
)
_DISEASE_STOPWORDS = frozenset(
    {
        "with",
        "and",
        "the",
        "type",
        "disease",
        "disorder",
        "familial",
        "hereditary",
        "syndrome",
        "cardiac",
        "heart",
        "gene",
        "variant",
        "mutation",
    }
)

# Caption vocabulary that says the table is NOT a plain case series.
_CAPTION_EXCLUDE_RE = re.compile(
    r"\b(?:famil(?:y|ies|ial)|kindreds?|pedigrees?|relatives?|segregat\w*|cascade|"
    r"carriers?|genotype[- ]positive|mutation[- ]positive|asymptomatic|symptomatic|"
    r"unaffected|penetran\w*|sub-?clinical|"
    r"clinical\s+(?:characteristics|features|data|findings|presentation|course)|"
    r"ecgs?|electrocardiograph\w*|qtc|"
    r"population|biobank|exomes?|genomes?|sequencing|wes|wgs|gnomad|exac|esp|"
    r"1000\s+genomes|topmed|dbsnp|clinvar|hgmd|"
    r"alleles?|allelic|allele\s+frequenc\w*|maf|minor\s+allele|"
    r"gwas|association\s+stud\w*|odds\s+ratio|"
    r"in[- ]vitro|functional|electrophysiolog\w*|expression|trafficking|"
    r"current\s+density|patch[- ]clamp|hek\s*293|xenopus|oocytes?|cho\s+cells|cells?|"
    r"transfect\w*|assays?|primers?|oligonucleotides?|"
    r"sudden\s+(?:cardiac|unexplained|arrhythmic|infant)?\s*deaths?|autops\w*|"
    r"post-?mortem|scd|suds?|sudy|sids|decedents?|victims?|"
    r"literature|published|previously\s+reported|reported\s+(?:in|by)|"
    r"meta-?analys\w*|review|incidental|secondary\s+findings?|"
    r"characteristics\s+of|"
    # Two patient groups, or a group defined by NOT having the phenotype
    # ("CTEPH patients and PE without PH patients"), are not one case series.
    r"without|versus|vs\.?|compared\s+(?:with|to))\b",
    re.IGNORECASE,
)
# Header vocabulary that marks a case-control layout: a bound "Cases" count
# next to a "Controls" count is a two-arm table, and the row's people count
# is one arm's denominator-side reading, not a phenotype partition.
_CASE_HEADER_RE = re.compile(r"\bcases?\b|\bpatients?\b|\bprobands?\b", re.IGNORECASE)
_CONTROL_HEADER_RE = re.compile(r"\bcontrols?\b|\bhealthy\b", re.IGNORECASE)
# A second column that carries clinical structure means the table is not a
# bare count list; leave it to the explicit-column and patient-row lanes.
_HEADER_EXCLUDE_RE = re.compile(
    r"affected|unaffected|symptom|phenotyp|diagnos|status|ecg|qtc|"
    r"relative|famil|kindred|pedigree|penetran|carrier|onset|event|"
    r"syncope|arrest|death|icd|therapy|treatment|followup|follow-up",
    re.IGNORECASE,
)
_COLUMN_EXCLUDE_TOKENS = (
    "famil",
    "kindred",
    "pedigree",
    "relative",
    "carrier",
    "screened",
    "tested",
    "cohort",
    "samplesize",
    "allele",
    "chromosome",
    "frequenc",
    "freq",
    "maf",
    "percent",
    "gnomad",
    "exac",
    "esp",
    "topmed",
    "clinvar",
    "dbsnp",
    "occurrence",
    "homozyg",
    "heterozyg",
    "genotype",
    "age",
    "year",
    "score",
    # The selected count header is skipped by the sibling-header check. It
    # must not turn a symptom/negative/assessed subset into the whole case
    # cohort merely because it also says "patients" or "controls". Use the
    # normalized label to cover converter-glued headers too.
    "affected",
    "symptom",
    "phenotyp",
    "diagnos",
    "ecg",
    "onset",
    "event",
    "syncope",
    "arrest",
    "death",
    "deceased",
    "therapy",
    "treatment",
    "followup",
    "without",
    "negative",
    "positive",
)
_COUNT_WORD_TOKENS = ("noof", "numberof", "number", "count", "nof", "total")
_BARE_COUNT_LABELS = frozenset({"n", "no", "number", "count", "counts", "total"})
_N_LABEL_RE = re.compile(r"^\s*n\s*(?:\(.*\))?\s*$|\(n\)|\bn\s*=", re.IGNORECASE)
# Column headers arrive glued by converters ("No. ofpatients",
# "No. of unrelatedindividuals"), so nouns in a header are matched as
# substrings of the normalized label, not with word boundaries.
_CASE_NOUN_TOKENS = ("patient", "proband", "indexcase", "case")
_PEOPLE_NOUN_TOKENS = ("individual", "subject", "person", "people", "participant")
_FIXED_WIDTH_COLUMN_GAP_RE = re.compile(r"\S {3,}\S")
_SENTENCE_END_RE = re.compile(r"[.!?](?=\s+[A-Z(])")
# A caption sentence carries ordinary lowercase words ("found in", "and",
# "variants"); a collapsed Title-Case header fragment does not.
_PROSE_WORD_RE = re.compile(r"(?:^|\s)[a-z][a-z-]{2,}\b|\.\s")
_CAPTION_HEADER_FRAGMENT_RE = re.compile(
    r"^(?:(?:patient|subject|participant)\s+(?:id|number|no\.?)(?:\s|$)|"
    r"(?:mutation|variant|mutation\s+or\s+(?:rare\s+)?variant)\s+"
    r"(?:site|location)(?:\s|$))",
    re.IGNORECASE,
)

_ASCERTAINMENT_RE = re.compile(
    r"\b(?:\d[\d,]*\s+)?(?:unrelated\s+|consecutive\s+|index\s+)?"
    r"(?:patients?|probands?|index\s+cases?|cases?)\s+"
    r"(?:with|referred\s+for|diagnosed\s+with|fulfilling|meeting|who\s+(?:were|had|met)|"
    r"carrying|affected\s+by)\b",
    re.IGNORECASE,
)
_TITLE_ASCERTAINMENT_RE = re.compile(
    r"\b(?:patients?|probands?|index\s+cases?|cases?)\s+"
    r"(?:with|referred\s+for|diagnosed\s+with|affected\s+by|fulfilling|meeting)\b",
    re.IGNORECASE,
)
_TITLE_COHORT_NOUN_RE = re.compile(
    r"\b(?:patients?|probands?|index\s+cases?|cases?)\b", re.IGNORECASE
)
# A caption that describes people rather than a mutation catalogue. One
# implicit carrier per row in such a table is a clinical roster (probands and
# relatives, screened individuals, follow-up subjects) whose phenotype must be
# read from the rows, never projected from the paper.
_ROSTER_CAPTION_RE = re.compile(
    r"\b(?:clinical|characteristics?|phenotyp\w*|demograph\w*|subjects?|"
    r"individuals?|persons?|participants?|carriers?|relatives?|famil\w*|kindreds?|"
    r"pedigrees?|cohorts?|registry|registries|ecgs?|symptoms?|diagnos\w*|"
    r"follow[- ]?up|outcomes?|genotype[- ]positive|screen\w*|tested)\b",
    re.IGNORECASE,
)
_CATALOGUE_NOUN_RE = re.compile(
    r"\b(?:mutations?|variants?|substitutions?|alterations?)\b", re.IGNORECASE
)
_PAPER_EXCLUDE_RE = re.compile(
    r"\b(?:population|biobank|general\s+population|community|"
    r"exome|genome|sequencing\s+of|family\s+members|relatives|cascade|"
    r"molecular\s+autopsy|sudden\s+(?:cardiac|unexplained)\s+death|decedents?)\b",
    re.IGNORECASE,
)
_SENTENCE_SPLIT_RE = re.compile(r"(?<=[.!?])\s+")
_PAPER_CONTEXT_CHARS = 8000


@dataclass(frozen=True)
class TableCohort:
    """Classification of one source table's cohort."""

    role: Optional[str]  # "case" | "control" | None
    tier: Optional[str]
    reason: str
    caption: str = ""
    count_label: str = ""
    quote: str = ""
    count_unit: str = "people"  # "people" | "alleles"

    def as_dict(self) -> dict[str, Any]:
        return {
            "role": self.role,
            "tier": self.tier,
            "reason": self.reason,
            "caption": self.caption[:600],
            "count_label": self.count_label[:200],
            "quote": self.quote[:600],
            "count_unit": self.count_unit,
        }


@dataclass
class PaperContext:
    """Paper-level ascertainment evidence for the optional second tier."""

    quote: str = ""
    excluded_by: str = ""
    disease_terms: list[str] = field(default_factory=list)
    from_title: bool = False


# --------------------------------------------------------------------------- #
# Settings
# --------------------------------------------------------------------------- #


def _setting(name: str, default: bool) -> bool:
    try:
        from config.settings import get_settings

        return bool(getattr(get_settings(), name, default))
    except Exception:  # pragma: no cover - settings unavailable in odd contexts
        return default


def table_cohort_phenotype_enabled() -> bool:
    return _setting("table_cohort_phenotype_enabled", True)


def paper_ascertainment_tier_enabled() -> bool:
    return _setting("table_cohort_paper_ascertainment_enabled", False)


def title_ascertainment_tier_enabled() -> bool:
    return _setting("table_cohort_title_ascertainment_enabled", True)


# --------------------------------------------------------------------------- #
# Text helpers
# --------------------------------------------------------------------------- #


def _normalize_label(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").lower())


def _squash(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "")).strip()


def _clean_caption_line(line: str) -> str:
    return _squash(line.lstrip("#").replace("**", "").replace("*", "").strip())


def _coerce_count(value: Any) -> Optional[int]:
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value if value >= 0 else None
    if isinstance(value, float) and value.is_integer() and value >= 0:
        return int(value)
    if isinstance(value, str):
        cleaned = value.strip().replace(",", "")
        if cleaned.isdigit():
            return int(cleaned)
    return None


def disease_pattern(disease: Optional[str]) -> Optional[re.Pattern[str]]:
    """Regex for the run's own disease phrase (whole phrase plus content words)."""
    text = _squash(disease)
    if not text:
        return None
    alternatives = [re.escape(text)]
    for token in re.findall(r"[A-Za-z][A-Za-z0-9-]{3,}", text):
        if token.lower() in _DISEASE_STOPWORDS:
            continue
        alternatives.append(re.escape(token))
    return re.compile(r"\b(?:" + "|".join(alternatives) + r")\b", re.IGNORECASE)


def _disease_hit(text: str, run_disease: Optional[re.Pattern[str]]) -> str:
    match = _BUILTIN_DISEASE_RE.search(text)
    if match:
        return match.group(0)
    if run_disease is not None:
        match = run_disease.search(text)
        if match:
            return match.group(0)
    return ""


def _is_deterministic_table_row(variant: dict[str, Any]) -> bool:
    """Only code-parsed table rows are eligible; model rows keep model rules."""
    if str(variant.get("source_layer") or "").strip().lower() == "regex_table":
        return True
    for holder in (variant, variant.get("patients")):
        if not isinstance(holder, dict):
            continue
        extra = holder.get("locator_extra")
        if isinstance(extra, dict) and str(extra.get("parser") or "").strip():
            return True
    notes = str(variant.get("additional_notes") or "").lower()
    return "deterministic" in notes and "parsed" in notes


# Extraction JSON written before 2026-09-08 by the fixed-width clinical
# mutation parser labelled an implicit one-carrier row with the "Coding Effect"
# text column. Archived runs keep that label; read it as the implicit-row label
# it always meant so replay and refresh classify those rows exactly as a fresh
# extraction would. A printed count ("Coding Effect count") is untouched.
_LEGACY_IMPLICIT_ROW_LABELS = frozenset(
    {("fixed_width_clinical_mutation", "coding effect")}
)


def _row_parser(variant: dict[str, Any]) -> str:
    for holder in (variant, variant.get("patients")):
        if not isinstance(holder, dict):
            continue
        extra = holder.get("locator_extra")
        if isinstance(extra, dict) and str(extra.get("parser") or "").strip():
            return str(extra["parser"]).strip().lower()
    return ""


def _count_label(variant: dict[str, Any]) -> str:
    label = _raw_count_label(variant)
    if (_row_parser(variant), label.casefold()) in _LEGACY_IMPLICIT_ROW_LABELS:
        penetrance = variant.get("penetrance_data")
        patients = variant.get("patients")
        carriers = _coerce_count(
            (penetrance or {}).get("total_carriers_observed")
            if isinstance(penetrance, dict)
            else None
        )
        if carriers is None and isinstance(patients, dict):
            carriers = _coerce_count(patients.get("count"))
        if carriers == 1:
            return IMPLICIT_ROW_CARRIER_LABEL
    return label


def _raw_count_label(variant: dict[str, Any]) -> str:
    provenance = variant.get("count_provenance")
    if isinstance(provenance, dict):
        label = _squash(provenance.get("carriers_column_label"))
        if label:
            return label
    patients = variant.get("patients")
    if isinstance(patients, dict):
        return _squash(patients.get("column_ref"))
    return ""


def _table_label(variant: dict[str, Any]) -> str:
    """The most descriptive table label the parser attached to this row."""
    candidates: list[str] = []
    for key in ("source_table", "source_table_caption", "source_ref"):
        candidates.append(_squash(variant.get(key)))
    patients = variant.get("patients")
    if isinstance(patients, dict):
        candidates.append(_squash(patients.get("source_ref")))
    location = _squash(variant.get("source_location"))
    location = _ROUTER_LOCATION_SUFFIX_RE.sub("", location)
    location = re.sub(r"\s*\((?:regex|router)[^)]*\)\s*$", "", location)
    candidates.append(location)
    candidates = [c for c in candidates if c and _TABLE_LABEL_RE.match(c)]
    if not candidates:
        return ""
    return max(candidates, key=len)


def _table_headers(variant: dict[str, Any]) -> list[str]:
    headers = variant.get("source_table_headers")
    if isinstance(headers, list):
        return [_squash(h) for h in headers if _squash(h)]
    return []


def _router_table_id(variant: dict[str, Any]) -> Optional[str]:
    for holder in (variant, variant.get("patients")):
        if not isinstance(holder, dict):
            continue
        extra = holder.get("locator_extra")
        if isinstance(extra, dict) and extra.get("table_id"):
            return str(extra["table_id"])
    return None


def _resolve_caption(
    label: str, source_text: str, *, max_lines: int = 3
) -> tuple[str, str]:
    """Resolve a caption only when every matching source anchor agrees.

    Parsers frequently keep only the anchor (``Table 2``) or the first wrapped
    line (``... Frequency in 406``). The descriptive sentence is what names the
    cohort, so look it up: find the label's occurrences, then append the
    following non-table lines until a blank line, a new table label, a heading,
    or ``max_lines`` lines.
    """
    label = _squash(label)
    if not label or not source_text:
        return label, "unresolved"
    lines = source_text.splitlines()
    wanted = _normalize_label(label)
    bare = bool(_BARE_TABLE_LABEL_RE.match(label))
    candidates: dict[str, str] = {}
    for index, raw in enumerate(lines):
        line = _clean_caption_line(raw)
        if not line or raw.lstrip().startswith("|"):
            continue
        normalized = _normalize_label(line)
        if not normalized.startswith(wanted):
            continue
        if bare and not _BARE_TABLE_LABEL_RE.match(line):
            # Accept an inline caption ("Table 2. ..."), but never a
            # different table number or a body sentence ("Table 2 shows...").
            anchor = _TABLE_LABEL_RE.match(line)
            if (
                anchor is None
                or _normalize_label(anchor.group(0)) != wanted
                or not anchor.group(0).rstrip().endswith((".", ":"))
            ):
                continue
        parts = [line]
        taken = 0
        cursor = index + 1
        while cursor < len(lines) and taken < max_lines:
            candidate = lines[cursor]
            stripped = candidate.strip()
            cursor += 1
            if not stripped:
                if taken == 0 and bare:
                    # "### Table 2" / blank / "Control variants ..."
                    continue
                break
            if stripped.startswith("|") or stripped.startswith("#"):
                break
            if _FIXED_WIDTH_COLUMN_GAP_RE.search(stripped):
                # Fixed-width text: the caption is followed by column header
                # fragments separated by runs of spaces, not by prose.
                break
            cleaned = _clean_caption_line(candidate)
            if not cleaned or _TABLE_LABEL_RE.match(cleaned):
                break
            if not _PROSE_WORD_RE.search(cleaned):
                # "History of Cardiac" / "Nucleotide Aborted Cardiac Events":
                # Title-Case header fragments whose spacing a converter has
                # collapsed. A caption sentence carries ordinary lowercase words.
                break
            if _CAPTION_HEADER_FRAGMENT_RE.search(cleaned):
                # Linearized headers can contain lowercase words, so the
                # prose-word check alone cannot distinguish a second copy's
                # "Mutation site ... Mean QTc" from a caption continuation.
                break
            parts.append(cleaned)
            taken += 1
            if bare or cleaned.endswith((".", ":")):
                break
        text = _squash(" ".join(parts))
        if not bare:
            # Keep the caption sentence only; a wrapped caption often runs on
            # into the table's descriptive paragraph.
            anchor_end = len(_squash(line))
            match = _SENTENCE_END_RE.search(text, max(anchor_end - 1, 0))
            if match:
                text = text[: match.end()]
        if bare and _normalize_label(text) == wanted:
            # An empty anchor/table-of-contents entry supplies no caption and
            # cannot contradict a descriptive occurrence later in the source.
            continue
        if _CONTINUED_CAPTION_RE.match(text):
            # A continuation heading repeats the label without describing the
            # table; it neither names a cohort nor contradicts the caption.
            continue
        candidates.setdefault(_normalize_label(text), text[:1200])
    if len(candidates) > 1:
        # Main text and supplements frequently reuse Table 1/2. A bare label
        # cannot choose between them, even if the first happens to name cases.
        return label, "ambiguous"
    if candidates:
        return next(iter(candidates.values())), "resolved"
    return label, "unresolved"


def expand_caption(label: str, source_text: str, *, max_lines: int = 3) -> str:
    """Expand a source-unique caption; preserve the label if it is ambiguous."""
    return _resolve_caption(label, source_text, max_lines=max_lines)[0]


def _is_bare_count_label(label: str) -> bool:
    normalized = _normalize_label(label)
    return normalized in _BARE_COUNT_LABELS or bool(_N_LABEL_RE.search(label or ""))


def _has_count_word(label: str) -> bool:
    normalized = _normalize_label(label)
    if any(token in normalized for token in _COUNT_WORD_TOKENS):
        return True
    return bool(_N_LABEL_RE.search(label or ""))


def _label_has_case_noun(label: str) -> bool:
    normalized = _normalize_label(label)
    return any(token in normalized for token in _CASE_NOUN_TOKENS)


def _label_has_people_noun(label: str) -> bool:
    normalized = _normalize_label(label)
    return any(token in normalized for token in _PEOPLE_NOUN_TOKENS)


def _column_excluded(label: str) -> str:
    normalized = _normalize_label(label)
    for token in _COLUMN_EXCLUDE_TOKENS:
        if token in normalized:
            return token
    if re.search(r"\bqtc\b", label, re.IGNORECASE):
        return "qtc"
    # Caption exclusions describe incompatible units/populations wherever
    # they are printed, including inside the selected count column. Keep
    # word boundaries here so ordinary words such as "individuals" survive.
    context_exclusion = _CAPTION_EXCLUDE_RE.search(label)
    if context_exclusion:
        return context_exclusion.group(0).lower()
    return ""


def _column_counts_people(label: str) -> bool:
    """A per-variant people count: a people/case noun with a count word, or N."""
    if not label:
        return False
    if _is_bare_count_label(label):
        return True
    if _label_has_case_noun(label) or _label_has_people_noun(label):
        return _has_count_word(label) or _is_bare_count_label(label)
    return False


def title_ascertains(
    title: Optional[str], run_disease: Optional[re.Pattern[str]] = None
) -> str:
    """Return the disease a paper title ascertains its cohort for, or ``""``.

    Accepts "probands with Brugada syndrome", "patients referred for long QT
    syndrome", or "long-QT syndrome patients": a people noun bound to a disease
    name within the title. A title that names an exome/genome/population/
    autopsy/relatives design is enrollment, not diagnosis, and never qualifies.
    """
    title_text = _squash(title)
    if not title_text or _PAPER_EXCLUDE_RE.search(title_text):
        return ""
    match = _TITLE_ASCERTAINMENT_RE.search(title_text)
    if match:
        disease = _disease_hit(title_text[match.end() : match.end() + 80], run_disease)
        if disease:
            return disease
    for noun in _TITLE_COHORT_NOUN_RE.finditer(title_text):
        disease = _disease_hit(
            title_text[max(0, noun.start() - 40) : noun.start()], run_disease
        )
        if disease:
            return disease
    return ""


def paper_context(
    source_text: str,
    *,
    title: Optional[str] = None,
    run_disease: Optional[re.Pattern[str]] = None,
    sentences: bool = True,
) -> PaperContext:
    """Find the title or one sentence that ascertains a disease case cohort."""
    context = PaperContext()
    head = _squash(source_text[:_PAPER_CONTEXT_CHARS])
    title_text = _squash(title)
    disease = title_ascertains(title_text, run_disease)
    if disease:
        context.quote = title_text[:400]
        context.disease_terms.append(disease)
        context.from_title = True
        return context
    if not sentences:
        return context
    for sentence in _SENTENCE_SPLIT_RE.split(head):
        if not _ASCERTAINMENT_RE.search(sentence):
            continue
        disease = _disease_hit(sentence, run_disease)
        if not disease:
            continue
        if _PAPER_EXCLUDE_RE.search(sentence):
            context.excluded_by = _PAPER_EXCLUDE_RE.search(sentence).group(0)  # type: ignore[union-attr]
            continue
        context.quote = _squash(sentence)[:400]
        context.disease_terms.append(disease)
        break
    return context


# --------------------------------------------------------------------------- #
# Classification
# --------------------------------------------------------------------------- #


def classify_table_cohort(
    caption: str,
    count_label: str,
    headers: Iterable[str] = (),
    *,
    run_disease: Optional[re.Pattern[str]] = None,
    paper: Optional[PaperContext] = None,
    allow_paper_tier: bool = False,
    allow_title_tier: bool = False,
) -> TableCohort:
    """Decide whether a count table enumerates cases, controls, or neither.

    ``caption`` is the printed caption (already expanded), ``count_label`` the
    header of the per-variant people-count column, ``headers`` every column
    header the parser saw (may be empty for fixed-width text).
    """
    caption = _squash(caption)
    count_label = _squash(count_label)
    header_list = [_squash(h) for h in headers if _squash(h)]
    base = {"caption": caption, "count_label": count_label}

    if count_label.casefold() == IMPLICIT_ROW_CARRIER_LABEL:
        return _classify_per_person_rows(
            caption,
            header_list,
            base,
            run_disease=run_disease,
            paper=paper,
            allow_title_tier=allow_title_tier,
        )
    excluded_token = _column_excluded(count_label)
    if excluded_token:
        return TableCohort(
            None, None, f"count_column_excluded:{excluded_token}", **base
        )

    control_caption = _CONTROL_RE.search(caption)
    control_column = _CONTROL_RE.search(count_label)
    case_in_caption = _CASE_NOUN_RE.search(caption)
    if control_column and (
        _label_has_case_noun(count_label)
        or _disease_hit(count_label, run_disease)
        or re.search(r"[+/&]|\band\b", count_label)
    ):
        # "BrS + LQT + Control" (SCN5A 25904541): a converter joined three
        # count columns into one label. A count that pools cases with controls
        # is neither phenotype class.
        return TableCohort(None, None, "count_column_mixes_cases_and_controls", **base)
    if control_column or control_caption:
        if case_in_caption or _disease_hit(caption, run_disease):
            return TableCohort(None, None, "caption_mixes_cases_and_controls", **base)
        if _ALLELE_UNIT_RE.search(caption + " " + count_label):
            # Preserve the explicit allele-unit abstention; it cannot stamp
            # either phenotype field. People-unit controls must pass the
            # same exclusions as cases below.
            return TableCohort(
                "control",
                TIER_COLUMN if control_column else TIER_CAPTION,
                "control_cohort",
                quote=(control_column or control_caption).group(0),
                count_unit="alleles",
                **base,
            )

    caption_exclusion = _CAPTION_EXCLUDE_RE.search(caption)
    if caption_exclusion:
        return TableCohort(
            None, None, f"caption_excluded:{caption_exclusion.group(0).lower()}", **base
        )
    count_key = _normalize_label(count_label)
    for header in header_list:
        if _normalize_label(header) == count_key:
            continue
        if _HEADER_EXCLUDE_RE.search(header):
            return TableCohort(
                None, None, f"clinical_column_present:{header[:60]}", **base
            )
    if any(_CASE_HEADER_RE.search(h) for h in header_list) and any(
        _CONTROL_HEADER_RE.search(h) for h in header_list
    ):
        return TableCohort(None, None, "case_and_control_columns_present", **base)

    if control_column or control_caption:
        # A control caption does not turn an arbitrary measurement column
        # into people. A bound "Controls" / "Control (8975)" header also
        # names the counted people directly.
        control_count_label = bool(control_column) and (
            _has_count_word(count_label)
            or bool(
                re.fullmatch(
                    r"controls?(?:\s*\(\s*(?:n\s*=?\s*)?\d+\s*\))?",
                    count_label,
                    re.IGNORECASE,
                )
            )
        )
        if not (_column_counts_people(count_label) or control_count_label):
            return TableCohort(None, None, "count_column_not_a_people_count", **base)
        return TableCohort(
            "control",
            TIER_COLUMN if control_column else TIER_CAPTION,
            "control_cohort",
            quote=(control_column or control_caption).group(0),
            **base,
        )

    if not _column_counts_people(count_label):
        return TableCohort(None, None, "count_column_not_a_people_count", **base)

    disease = _disease_hit(caption, run_disease)
    if _label_has_case_noun(count_label) and _has_count_word(count_label):
        if not (
            disease
            or _disease_hit(count_label, run_disease)
            or re.search(r"\bcases?\b", count_label, re.IGNORECASE)
        ):
            # "Number of Patients" in an uncaptioned mutation catalogue can
            # count a mixed-phenotype roster or a literature database. A
            # patient/proband noun alone says nothing about disease status.
            return TableCohort(
                None, None, "patient_count_without_disease_context", **base
            )
        return TableCohort(
            "case",
            TIER_COLUMN,
            "count_column_names_cases",
            quote=count_label,
            **base,
        )
    if disease and case_in_caption:
        return TableCohort(
            "case",
            TIER_CAPTION,
            "caption_names_disease_cases",
            quote=f"{disease}; {case_in_caption.group(0)}",
            **base,
        )
    if disease and (
        _label_has_case_noun(count_label) or _label_has_people_noun(count_label)
    ):
        return TableCohort(
            "case",
            TIER_CAPTION_DISEASE,
            "caption_names_disease_column_counts_people",
            quote=f"{disease}; {count_label}",
            **base,
        )
    if disease and _is_bare_count_label(count_label):
        return TableCohort(None, None, "disease_caption_without_people_noun", **base)

    title_evidence = bool(allow_title_tier and paper is not None and paper.from_title)
    if (allow_paper_tier or title_evidence) and paper is not None and paper.quote:
        if case_in_caption or _PEOPLE_NOUN_RE.search(caption):
            return TableCohort(
                None, None, "caption_names_cohort_without_disease", **base
            )
        return TableCohort(
            "case",
            TIER_TITLE if title_evidence else TIER_PAPER,
            "title_ascertainment" if title_evidence else "paper_ascertainment",
            quote=paper.quote,
            **base,
        )
    return TableCohort(None, None, "no_cohort_evidence", **base)


def _classify_per_person_rows(
    caption: str,
    header_list: list[str],
    base: dict[str, str],
    *,
    run_disease: Optional[re.Pattern[str]],
    paper: Optional[PaperContext],
    allow_title_tier: bool,
) -> TableCohort:
    """Classify a table whose parser inferred one carrier per row.

    A mutation catalogue (one proband's mutation per row, no clinical
    columns) may be projected when its own caption names the disease case
    series, or when the paper's title defines a disease proband cohort and the
    caption names nobody else. A clinical roster -- characteristics, relatives,
    carriers, screened or followed-up people, any phenotype/status column --
    always refuses: its rows carry their own phenotype and must be read.
    """
    refusal = "per_person_clinical_row"
    if _CAPTION_EXCLUDE_RE.search(caption) or _CONTROL_RE.search(caption):
        return TableCohort(None, None, refusal, **base)
    if _ROSTER_CAPTION_RE.search(caption) or not _CATALOGUE_NOUN_RE.search(caption):
        return TableCohort(None, None, refusal, **base)
    for header in header_list:
        if _HEADER_EXCLUDE_RE.search(header) or _CONTROL_HEADER_RE.search(header):
            return TableCohort(None, None, refusal, **base)
    disease = _disease_hit(caption, run_disease)
    if disease and _CASE_NOUN_RE.search(caption):
        return TableCohort(
            "case",
            TIER_CAPTION,
            "per_person_disease_case_catalogue",
            quote=f"{disease}; {caption}",
            **base,
        )
    if _CASE_NOUN_RE.search(caption) or _PEOPLE_NOUN_RE.search(caption):
        return TableCohort(
            None, None, "per_person_cohort_caption_without_disease", **base
        )
    if allow_title_tier and paper is not None and paper.from_title and paper.quote:
        return TableCohort(
            "case",
            TIER_TITLE,
            "per_person_proband_catalogue",
            quote=paper.quote,
            **base,
        )
    return TableCohort(None, None, refusal, **base)


# --------------------------------------------------------------------------- #
# Derivation
# --------------------------------------------------------------------------- #


def _identity(variant: dict[str, Any], index: int) -> str:
    for key in ("source_notation", "protein_notation", "cdna_notation"):
        value = variant.get(key)
        if value:
            return str(value)
    return f"variant_{index + 1}"


def _provenance_stamp(
    provenance: dict[str, Any], field_name: str, count_type: str, label: str
) -> None:
    provenance[f"{field_name}_column_label"] = label[:1000]
    provenance[f"{field_name}_count_type"] = count_type
    provenance[f"{field_name}_source"] = PROTOCOL_VERSION


def _append_fact(
    variant: dict[str, Any],
    *,
    fact_type: str,
    value: int,
    cohort: TableCohort,
    row: Any,
) -> None:
    facts = variant.get("fact_provenance")
    if not isinstance(facts, list):
        facts = []
        variant["fact_provenance"] = facts
    facts.append(
        {
            "fact_type": fact_type,
            "fact_value": str(value),
            "source_location": cohort.caption[:300] or None,
            "source_table": cohort.caption[:300] or None,
            "source_row": str(row) if row is not None else None,
            "source_column": cohort.count_label[:200] or None,
            "count_type": (
                "case"
                if fact_type == "affected_count" and cohort.role == "case"
                else "control"
                if fact_type == "affected_count"
                else "unaffected_control"
            ),
            "evidence_quote": (
                f"{cohort.caption} | column: {cohort.count_label} | "
                f"cohort: {cohort.role} ({cohort.tier}); {cohort.quote}"
            )[:2000],
        }
    )


def derive_table_cohort_phenotype_counts(
    extracted_data: dict[str, Any],
    source_text: str,
    *,
    gene_symbol: Optional[str] = None,
    disease: Optional[str] = None,
    title: Optional[str] = None,
    enabled: Optional[bool] = None,
    allow_paper_tier: Optional[bool] = None,
    allow_title_tier: Optional[bool] = None,
) -> dict[str, Any]:
    """Populate cohort-derived affected/unaffected on eligible table rows.

    Returns a deep copy; the input is not mutated. Rows that already carry a
    phenotype value are left alone, except that a parser-emitted control
    partition (``affected 0 / unaffected N``) or case copy (``affected N``)
    with no provenance is stamped when the table classifies the same way, so
    the always-on guard keeps a value the source actually supports.
    """
    if not isinstance(extracted_data, dict):
        return extracted_data
    variants = extracted_data.get("variants")
    if not isinstance(variants, list):
        return extracted_data
    if enabled is None:
        enabled = table_cohort_phenotype_enabled()
    if allow_paper_tier is None:
        allow_paper_tier = paper_ascertainment_tier_enabled()
    if allow_title_tier is None:
        allow_title_tier = title_ascertainment_tier_enabled()

    result = copy.deepcopy(extracted_data)
    metadata = result.setdefault("extraction_metadata", {})
    if not enabled:
        metadata["table_cohort_phenotype_derivation"] = {
            "protocol_version": PROTOCOL_VERSION,
            "attempted": False,
            "applied": False,
            "reason": "disabled_by_setting",
        }
        return result

    source_text = source_text or ""
    run_disease = disease_pattern(disease)
    paper = (
        paper_context(
            source_text,
            title=title,
            run_disease=run_disease,
            sentences=bool(allow_paper_tier),
        )
        if (allow_paper_tier or allow_title_tier)
        else None
    )
    router_tables: Optional[dict[str, Any]] = None
    table_cache: dict[tuple[str, str, tuple[str, ...]], TableCohort] = {}
    caption_cache: dict[str, tuple[str, str]] = {}
    tables_summary: dict[str, dict[str, Any]] = {}
    outcomes: list[dict[str, Any]] = []
    applied = 0
    stamped = 0

    def router_caption(table_id: str) -> tuple[str, list[str]]:
        nonlocal router_tables
        if router_tables is None:
            router_tables = {}
            try:
                from pipeline.table_router import enumerate_markdown_tables

                for table in enumerate_markdown_tables(
                    source_text, only_variant_like=False
                ):
                    router_tables[table.table_id] = table
            except Exception:  # pragma: no cover - defensive
                router_tables = {}
        table = router_tables.get(table_id)
        if table is None:
            return "", []
        return _squash(table.caption), [_squash(h) for h in table.header_cells]

    for index, variant in enumerate(result["variants"]):
        if not isinstance(variant, dict):
            continue
        identity = _identity(variant, index)
        if not _is_deterministic_table_row(variant):
            outcomes.append({"variant": identity, "status": "model_authored_row"})
            continue
        penetrance = variant.get("penetrance_data")
        if not isinstance(penetrance, dict):
            penetrance = {}
            variant["penetrance_data"] = penetrance
        patients = (
            variant.get("patients") if isinstance(variant.get("patients"), dict) else {}
        )
        carriers = _coerce_count(penetrance.get("total_carriers_observed"))
        if carriers is None:
            carriers = _coerce_count(patients.get("count"))
        if carriers is None or carriers <= 0:
            outcomes.append({"variant": identity, "status": "no_carrier_count"})
            continue
        provenance = variant.get("count_provenance")
        if not isinstance(provenance, dict):
            provenance = {}
            variant["count_provenance"] = provenance
        carrier_type = str(provenance.get("carriers_count_type") or "").strip().lower()
        if carrier_type not in {"", "per_variant_carrier"}:
            outcomes.append(
                {"variant": identity, "status": f"carrier_type_{carrier_type}"}
            )
            continue
        existing_affected = _coerce_count(penetrance.get("affected_count"))
        existing_unaffected = _coerce_count(penetrance.get("unaffected_count"))
        stampable_case_copy = (
            existing_affected == carriers and existing_unaffected is None
        )
        stampable_control = (
            existing_affected in (None, 0) and existing_unaffected == carriers
        )
        if (
            (existing_affected is not None or existing_unaffected is not None)
            and not stampable_case_copy
            and not stampable_control
        ):
            outcomes.append(
                {"variant": identity, "status": "phenotype_already_populated"}
            )
            continue
        count_label = _count_label(variant)
        label = _table_label(variant)
        headers = _table_headers(variant)
        caption = ""
        table_id = _router_table_id(variant)
        if table_id and (not label or _ROUTER_TABLE_ID_RE.match(label)):
            caption, router_headers = router_caption(table_id)
            headers = headers or router_headers
            label = label or f"Table {table_id}"
        if not caption:
            if label not in caption_cache:
                caption_cache[label] = _resolve_caption(label, source_text)
            caption, caption_status = caption_cache[label]
            if caption_status == "ambiguous":
                outcomes.append(
                    {"variant": identity, "status": "ambiguous_source_caption"}
                )
                continue
        continued = _CONTINUED_CAPTION_RE.match(caption) if caption else None
        if continued:
            # The parser kept the page-break heading; the table's cohort is
            # named by the caption printed at its first page.
            base_label = _squash(continued.group(1))
            if base_label not in caption_cache:
                caption_cache[base_label] = _resolve_caption(base_label, source_text)
            resolved, caption_status = caption_cache[base_label]
            if caption_status == "ambiguous":
                outcomes.append(
                    {"variant": identity, "status": "ambiguous_source_caption"}
                )
                continue
            if caption_status == "resolved":
                caption = resolved
        if not caption and not count_label:
            outcomes.append({"variant": identity, "status": "no_table_context"})
            continue
        key = (caption, count_label, tuple(headers))
        cohort = table_cache.get(key)
        if cohort is None:
            cohort = classify_table_cohort(
                caption,
                count_label,
                headers,
                run_disease=run_disease,
                paper=paper,
                allow_paper_tier=bool(allow_paper_tier),
                allow_title_tier=bool(allow_title_tier),
            )
            table_cache[key] = cohort
        summary = tables_summary.setdefault(
            caption or label or count_label,
            {**cohort.as_dict(), "rows_seen": 0, "rows_applied": 0, "rows_stamped": 0},
        )
        summary["rows_seen"] += 1

        affected = _coerce_count(penetrance.get("affected_count"))
        unaffected = _coerce_count(penetrance.get("unaffected_count"))
        aff_type = str(provenance.get("affected_count_type") or "").strip().lower()
        una_type = str(provenance.get("unaffected_count_type") or "").strip().lower()
        row = variant.get("source_row") or patients.get("row_ordinal")

        if cohort.role is None:
            outcomes.append(
                {"variant": identity, "status": f"not_derived:{cohort.reason}"}
            )
            continue
        if count_label.casefold() == IMPLICIT_ROW_CARRIER_LABEL and carriers != 1:
            # The parser asserted one carrier per row; any other count means
            # the row is not one proband and the catalogue rule does not apply.
            outcomes.append(
                {
                    "variant": identity,
                    "status": "not_derived:per_person_row_count_not_one",
                }
            )
            continue

        # The replay must undo only this lane, including when it certifies
        # a parser value the ordinary guard would already have kept.
        guarded_counts_without_projection = {
            "affected": affected,
            "unaffected": unaffected,
        }
        for cleared in phenotype_fields_to_clear(variant):
            if cleared.field in guarded_counts_without_projection:
                guarded_counts_without_projection[cleared.field] = None

        if cohort.role == "case":
            if unaffected is not None or (
                affected is not None and affected != carriers
            ):
                outcomes.append(
                    {"variant": identity, "status": "phenotype_already_populated"}
                )
                continue
            if (
                affected == carriers
                and aff_type
                and aff_type
                not in {
                    "per_variant_carrier",
                    "unknown",
                }
            ):
                outcomes.append(
                    {"variant": identity, "status": "phenotype_already_populated"}
                )
                continue
            copied = affected == carriers
            penetrance["affected_count"] = carriers
            _provenance_stamp(
                provenance,
                "affected",
                "case",
                f"{cohort.caption} | {cohort.count_label}",
            )
            _append_fact(
                variant,
                fact_type="affected_count",
                value=carriers,
                cohort=cohort,
                row=row,
            )
            if isinstance(variant.get("patients"), dict) and not variant[
                "patients"
            ].get("phenotype"):
                variant["patients"]["phenotype"] = (
                    f"{gene_symbol}-associated disease"
                    if gene_symbol
                    else "disease case"
                )
            status = "case_copy_stamped" if copied else "applied"
        else:  # control
            if affected not in (None, 0) or (
                unaffected is not None and unaffected != carriers
            ):
                outcomes.append(
                    {"variant": identity, "status": "phenotype_already_populated"}
                )
                continue
            if cohort.count_unit != "people":
                # "2,600 reference alleles": an allele count is not a person
                # count, so certify nothing. The parser's own row-level control
                # reading (if any) is left exactly as it was.
                outcomes.append(
                    {
                        "variant": identity,
                        "status": "not_derived:control_count_unit_is_alleles",
                        "source_table": cohort.caption[:200],
                    }
                )
                continue
            if (
                unaffected == carriers
                and una_type
                and una_type not in {"per_variant_carrier", "unknown"}
                and str(provenance.get("unaffected_source") or "")
            ):
                outcomes.append(
                    {"variant": identity, "status": "phenotype_already_populated"}
                )
                continue
            copied = unaffected == carriers
            penetrance["unaffected_count"] = carriers
            _provenance_stamp(
                provenance,
                "unaffected",
                "unaffected_control",
                f"{cohort.caption} | {cohort.count_label}",
            )
            _append_fact(
                variant,
                fact_type="unaffected_count",
                value=carriers,
                cohort=cohort,
                row=row,
            )
            penetrance["affected_count"] = 0
            _provenance_stamp(
                provenance,
                "affected",
                "control",
                f"{cohort.caption} | {cohort.count_label}",
            )
            _append_fact(
                variant,
                fact_type="affected_count",
                value=0,
                cohort=cohort,
                row=row,
            )
            if isinstance(variant.get("patients"), dict):
                variant["patients"]["phenotype"] = "unaffected control"
            status = "control_partition_stamped" if copied else "applied"

        variant["phenotype_derivation"] = {
            "protocol_version": PROTOCOL_VERSION,
            "method": METHOD,
            "cohort_role": cohort.role,
            "tier": cohort.tier,
            "source_table": cohort.caption[:600],
            "count_column": cohort.count_label[:200],
            "count_unit": cohort.count_unit,
            "evidence_quote": cohort.quote[:600],
            "guarded_counts_without_projection": guarded_counts_without_projection,
            "operational_rule": (
                "A case-series count column counts disease-ascertained people, "
                "so the per-variant count is the affected count; unaffected is "
                "left unassessed. A control table's count is unaffected; its "
                "affected is zero only when the unit is people."
            ),
        }
        if status == "applied":
            applied += 1
            summary["rows_applied"] += 1
        else:
            stamped += 1
            summary["rows_stamped"] += 1
        outcomes.append(
            {
                "variant": identity,
                "status": status,
                "cohort_role": cohort.role,
                "tier": cohort.tier,
                "source_table": cohort.caption[:200],
            }
        )

    metadata["table_cohort_phenotype_derivation"] = {
        "protocol_version": PROTOCOL_VERSION,
        "attempted": True,
        "applied": bool(applied or stamped),
        "paper_tier_enabled": bool(allow_paper_tier),
        "title_tier_enabled": bool(allow_title_tier),
        "paper_ascertainment_quote": (paper.quote if paper else ""),
        "paper_ascertainment_from_title": bool(paper.from_title) if paper else False,
        "paper_variant_count": len(variants),
        "applied_variant_count": applied,
        "stamped_variant_count": stamped,
        "tables": tables_summary,
        "outcomes": outcomes,
    }
    return result
