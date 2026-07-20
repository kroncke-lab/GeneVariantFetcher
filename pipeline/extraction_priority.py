"""Deterministic pre-LLM prioritization for extraction candidates.

The scorer is intentionally transparent and cheap. It does not decide whether a
paper is scientifically valid; it only orders papers so original variant/count
evidence is extracted before reviews, protocols, and weakly matched records.
"""

from __future__ import annotations

import csv
import json
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional

from utils.gene_metadata import gene_alias_regex
from utils.pmid_utils import extract_pmid_from_filename


MAX_SOURCE_CHARS = 300_000


VARIANT_PATTERNS = (
    re.compile(r"\bc\.\d+[A-Za-z0-9_>*+\-?]+", re.IGNORECASE),
    re.compile(r"\bp\.[A-Za-z]{1,3}\d+[A-Za-z*]{1,3}\b", re.IGNORECASE),
    re.compile(r"\brs\d{4,}\b", re.IGNORECASE),
    re.compile(r"\bIVS\d+[+\-]\d+[A-Z]>[A-Z]\b", re.IGNORECASE),
    re.compile(r"\b\d{2,6}(?:del|ins|dup)[A-Za-z0-9]*\b", re.IGNORECASE),
    re.compile(r"\b(?:epsilon|eps|e|E|[εe])\s*[234]\b"),
)

ORIGINAL_DATA_RE = re.compile(
    r"\b("
    r"cohort|case[- ]control|patients?|probands?|famil(?:y|ies)|"
    r"index cases?|unrelated|consecutive|screen(?:ed|ing)|sequenc(?:ed|ing)|"
    r"genotyp(?:ed|ing)|mutation analysis|variant analysis|study population|"
    r"participants?|subjects?|enrolled|recruited|ascertained|clinical characteristics"
    r")\b",
    re.IGNORECASE,
)

CARRIER_COUNT_RE = re.compile(
    r"\b("
    r"carriers?|non[- ]?carriers?|mutation carriers?|variant carriers?|"
    r"allele frequency|carrier frequency|genotype frequency|minor allele frequency|"
    r"heterozyg(?:ote|ous|osity)|homozyg(?:ote|ous|osity)|affected|unaffected|"
    r"cases?|controls?|odds ratio|OR"
    r")\b|"
    r"\bn\s*[=:]\s*\d+\b|\b\d+\s*/\s*\d+\b",
    re.IGNORECASE,
)

# Carrier-first / penetrance study signal: papers that ASCERTAIN carriers and
# track who develops disease (segregation, cascade, prospective follow-up) are
# the primary source of UNAFFECTED-carrier data, which patient-first case series
# under-report. Boost them so they are not out-prioritized by high-volume case
# series (criticisms 3 and 6).
PENETRANCE_RE = re.compile(
    r"\b("
    r"penetrance|co[- ]?segregation|segregat(?:ion|ed|es)|"
    r"cascade\s+(?:screening|testing)|unaffected\s+carriers?|"
    r"asymptomatic\s+carriers?|prospective\s+(?:cohort|follow[- ]?up)|"
    r"family\s+segregation|genotype[- ]first"
    r")\b",
    re.IGNORECASE,
)

TABLE_RE = re.compile(
    r"\b(?:table|supplementary table|supplemental table|e[- ]table)\s+[A-Za-z0-9]+",
    re.IGNORECASE,
)

REVIEW_TITLE_RE = re.compile(
    r"\b("
    r"review|systematic review|meta[- ]analysis|editorial|commentary|"
    r"guideline|consensus|position statement|protocol|study design"
    r")\b",
    re.IGNORECASE,
)

METHODS_RESULTS_RE = re.compile(
    r"\b(methods?|results?|materials and methods|statistical analysis|"
    r"patient characteristics|clinical data)\b",
    re.IGNORECASE,
)


def _gene_regex(gene_symbol: str) -> re.Pattern[str]:
    return gene_alias_regex(gene_symbol, include_query_aliases=True)


@dataclass
class ExtractionCandidate:
    """A scored candidate paper/source for LLM extraction."""

    pmid: str
    source_kind: str
    source_file: str
    score: int
    title: str = ""
    journal: str = ""
    year: str = ""
    rank: int = 0
    selected: bool = False
    text_chars: int = 0
    reasons: list[str] = field(default_factory=list)
    signals: dict[str, Any] = field(default_factory=dict)


@dataclass
class PriorityResult:
    """Selected sources plus an auditable ranking table."""

    selected_markdown_files: list[Path]
    selected_abstract_papers: list[tuple[str, str]]
    candidates: list[ExtractionCandidate]
    report_dir: Optional[Path] = None

    @property
    def selected_candidates(self) -> list[ExtractionCandidate]:
        return [candidate for candidate in self.candidates if candidate.selected]


def _load_abstract_record(
    record_path: Optional[str | Path],
) -> tuple[dict[str, Any], str]:
    if not record_path:
        return {}, ""
    path = Path(record_path)
    if not path.exists():
        return {}, ""
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}, ""
    metadata = data.get("metadata") if isinstance(data, dict) else {}
    abstract = data.get("abstract") if isinstance(data, dict) else ""
    return (metadata or {}, abstract or "")


def _read_source_sample(path: Path, max_chars: int = MAX_SOURCE_CHARS) -> str:
    size = path.stat().st_size
    if size <= max_chars:
        return path.read_text(encoding="utf-8", errors="replace")

    head = max_chars * 2 // 3
    tail = max_chars - head
    with path.open("rb") as handle:
        head_bytes = handle.read(head)
        handle.seek(max(0, size - tail))
        tail_bytes = handle.read(tail)
    return (
        head_bytes.decode("utf-8", errors="replace")
        + "\n\n[...middle omitted for prioritization...]\n\n"
        + tail_bytes.decode("utf-8", errors="replace")
    )


def _count_pattern(pattern: re.Pattern[str], text: str, cap: int) -> int:
    count = 0
    for _ in pattern.finditer(text):
        count += 1
        if count >= cap:
            return cap
    return count


def _variant_mentions(text: str, cap: int = 20) -> int:
    total = 0
    for pattern in VARIANT_PATTERNS:
        total += _count_pattern(pattern, text, cap - total)
        if total >= cap:
            return cap
    return total


def _line_cooccurrence_count(
    text: str,
    gene_pattern: re.Pattern[str],
    signal_patterns: Iterable[re.Pattern[str]],
    cap: int,
) -> int:
    count = 0
    for line in text.splitlines():
        if not gene_pattern.search(line):
            continue
        if any(pattern.search(line) for pattern in signal_patterns):
            count += 1
            if count >= cap:
                return cap
    return count


def _normalized_terms(terms: Optional[Iterable[str]]) -> list[str]:
    normalized: list[str] = []
    seen: set[str] = set()
    for term in terms or []:
        clean = " ".join(str(term).lower().split())
        if len(clean) < 3 or clean in seen:
            continue
        seen.add(clean)
        normalized.append(clean)
        ascii_apostrophe = clean.replace("'", "")
        if ascii_apostrophe != clean and ascii_apostrophe not in seen:
            seen.add(ascii_apostrophe)
            normalized.append(ascii_apostrophe)
    return normalized


def _disease_hits(text: str, disease_terms: list[str]) -> int:
    lower = text.lower()
    return sum(1 for term in disease_terms if term in lower)


def _artifact_table_count(source_file: Path, pmid: str) -> int:
    artifact_file = source_file.parent / f"{pmid}_artifacts.json"
    if not artifact_file.exists():
        return 0
    try:
        data = json.loads(artifact_file.read_text(encoding="utf-8"))
    except Exception:
        return 0
    main_text = data.get("main_text") if isinstance(data, dict) else {}
    if not isinstance(main_text, dict):
        return 0
    return int(main_text.get("table_captions_count") or 0) + int(
        main_text.get("supplement_descriptions_count") or 0
    )


def _score_candidate(
    *,
    pmid: str,
    source_kind: str,
    source_file: str,
    title: str,
    journal: str,
    year: str,
    body_text: str,
    abstract: str,
    disease_terms: list[str],
    gene_pattern: re.Pattern[str],
    source_size_bytes: int = 0,
) -> ExtractionCandidate:
    title_text = title or ""
    combined = "\n".join(part for part in (title_text, abstract, body_text) if part)
    source_path = Path(source_file)
    text_chars = len(body_text)

    title_disease_hits = _disease_hits(title_text, disease_terms)
    abstract_disease_hits = _disease_hits(abstract, disease_terms)
    body_disease_hits = _disease_hits(body_text, disease_terms)
    disease_hit_score = (
        title_disease_hits * 3 + abstract_disease_hits * 2 + min(body_disease_hits, 3)
    )

    variant_mentions = _variant_mentions(combined)
    original_mentions = _count_pattern(ORIGINAL_DATA_RE, combined, 20)
    carrier_mentions = _count_pattern(CARRIER_COUNT_RE, combined, 20)
    penetrance_mentions = _count_pattern(PENETRANCE_RE, combined, 20)
    table_mentions = _count_pattern(TABLE_RE, combined, 20)
    if source_kind == "fulltext":
        table_mentions += _artifact_table_count(source_path, pmid)
    methods_results = bool(METHODS_RESULTS_RE.search(body_text))
    gene_title_mentions = _count_pattern(gene_pattern, title_text, 5)
    gene_abstract_mentions = _count_pattern(gene_pattern, abstract, 10)
    gene_body_mentions = _count_pattern(gene_pattern, body_text, 20)
    gene_variant_lines = _line_cooccurrence_count(
        combined, gene_pattern, VARIANT_PATTERNS, 20
    )
    gene_carrier_lines = _line_cooccurrence_count(
        combined, gene_pattern, (CARRIER_COUNT_RE,), 20
    )

    likely_full_text = source_kind == "fulltext" and text_chars >= 5_000
    rich_full_text = source_kind == "fulltext" and text_chars >= 40_000
    oversized_source = source_kind == "fulltext" and source_size_bytes > 20_000_000
    review_title = bool(REVIEW_TITLE_RE.search(title_text))
    review_abstract = bool(REVIEW_TITLE_RE.search(abstract[:2_000]))
    protocol_like = bool(
        re.search(r"\b(protocol|study design|trial design)\b", combined, re.IGNORECASE)
    )

    score = 0
    reasons: list[str] = []

    if likely_full_text:
        score += 25
        reasons.append("usable full text")
    elif source_kind == "abstract":
        score -= 30
        reasons.append("abstract only")
    if rich_full_text:
        score += 8
        reasons.append("rich full text")
    if disease_hit_score:
        score += min(25, disease_hit_score * 5)
        reasons.append("disease match")
    if gene_title_mentions:
        score += 35
        reasons.append("gene in title")
    if gene_abstract_mentions:
        score += min(25, gene_abstract_mentions * 5)
        reasons.append("gene in abstract")
    if gene_body_mentions:
        score += min(20, gene_body_mentions)
        reasons.append("gene in full text")
    if gene_variant_lines:
        score += min(60, gene_variant_lines * 6)
        reasons.append("gene-specific variant lines")
    if gene_carrier_lines:
        score += min(45, gene_carrier_lines * 5)
        reasons.append("gene-specific carrier/count lines")
    if variant_mentions:
        score += min(20, variant_mentions * 2)
        reasons.append("explicit variant terms")
    if carrier_mentions:
        score += min(20, carrier_mentions * 2)
        reasons.append("carrier/count terms")
    if penetrance_mentions:
        score += min(24, penetrance_mentions * 4)
        reasons.append("penetrance/segregation study signal")
    if original_mentions:
        score += min(20, original_mentions * 2)
        reasons.append("original cohort/case terms")
    if table_mentions:
        score += min(20, table_mentions * 3)
        reasons.append("table/supplement signal")
    if methods_results:
        score += 6
        reasons.append("methods/results sections")
    if oversized_source:
        score -= 100
        reasons.append("oversized source over 20MB")
    if not (gene_title_mentions or gene_abstract_mentions or gene_body_mentions):
        score -= 45
        reasons.append("no gene-specific mention")

    if review_title:
        score -= 55
        reasons.append("review/editorial/protocol title")
    elif review_abstract:
        score -= 25
        reasons.append("review/editorial abstract")
    if protocol_like and variant_mentions == 0:
        score -= 15
        reasons.append("protocol/design without variants")

    signals = {
        "likely_full_text": likely_full_text,
        "rich_full_text": rich_full_text,
        "source_size_bytes": source_size_bytes,
        "oversized_source": oversized_source,
        "disease_term_hits": title_disease_hits
        + abstract_disease_hits
        + body_disease_hits,
        "title_disease_hits": title_disease_hits,
        "abstract_disease_hits": abstract_disease_hits,
        "body_disease_hits": body_disease_hits,
        "variant_mentions": variant_mentions,
        "gene_title_mentions": gene_title_mentions,
        "gene_abstract_mentions": gene_abstract_mentions,
        "gene_body_mentions": gene_body_mentions,
        "gene_variant_lines": gene_variant_lines,
        "gene_carrier_lines": gene_carrier_lines,
        "carrier_count_mentions": carrier_mentions,
        "penetrance_study_mentions": penetrance_mentions,
        "original_data_mentions": original_mentions,
        "table_mentions": table_mentions,
        "methods_results": methods_results,
        "review_title_penalty": review_title,
        "review_abstract_penalty": review_abstract,
        "protocol_like": protocol_like,
    }

    return ExtractionCandidate(
        pmid=pmid,
        source_kind=source_kind,
        source_file=source_file,
        score=score,
        title=title,
        journal=journal,
        year=str(year or ""),
        text_chars=text_chars,
        reasons=reasons,
        signals=signals,
    )


def _year_sort_value(year: str) -> int:
    match = re.search(r"\d{4}", str(year or ""))
    return int(match.group(0)) if match else 0


def _sort_key(candidate: ExtractionCandidate) -> tuple[Any, ...]:
    signals = candidate.signals
    return (
        -candidate.score,
        candidate.source_kind != "fulltext",
        -int(signals.get("gene_title_mentions") or 0),
        -int(signals.get("gene_variant_lines") or 0),
        -int(signals.get("gene_carrier_lines") or 0),
        -int(signals.get("gene_abstract_mentions") or 0),
        -int(signals.get("variant_mentions") or 0),
        -int(signals.get("carrier_count_mentions") or 0),
        -int(signals.get("table_mentions") or 0),
        -int(signals.get("disease_term_hits") or 0),
        -_year_sort_value(candidate.year),
        candidate.pmid,
    )


def prioritize_extraction_sources(
    *,
    markdown_files: list[Path],
    abstract_papers: list[tuple[str, str]],
    gene_symbol: str,
    disease_terms: Optional[Iterable[str]] = None,
    abstract_records: Optional[dict[str, str]] = None,
    top_n: Optional[int] = None,
    offset: int = 0,
    report_dir: Optional[Path] = None,
) -> PriorityResult:
    """Rank and optionally limit extraction sources before LLM submission."""

    terms = _normalized_terms(disease_terms)
    gene_pattern = _gene_regex(gene_symbol)
    abstract_records = abstract_records or {}
    candidates: list[ExtractionCandidate] = []

    for source_path in markdown_files:
        pmid = extract_pmid_from_filename(source_path)
        if not pmid:
            continue
        metadata, abstract = _load_abstract_record(abstract_records.get(pmid))
        body_text = _read_source_sample(source_path)
        candidates.append(
            _score_candidate(
                pmid=pmid,
                source_kind="fulltext",
                source_file=str(source_path),
                title=str(metadata.get("title") or ""),
                journal=str(metadata.get("journal") or ""),
                year=str(metadata.get("year") or ""),
                body_text=body_text,
                abstract=abstract,
                disease_terms=terms,
                gene_pattern=gene_pattern,
                source_size_bytes=source_path.stat().st_size,
            )
        )

    fulltext_pmids = {candidate.pmid for candidate in candidates}
    for pmid, record_path in abstract_papers:
        if pmid in fulltext_pmids:
            continue
        metadata, abstract = _load_abstract_record(record_path)
        candidates.append(
            _score_candidate(
                pmid=pmid,
                source_kind="abstract",
                source_file=str(record_path),
                title=str(metadata.get("title") or ""),
                journal=str(metadata.get("journal") or ""),
                year=str(metadata.get("year") or ""),
                body_text="",
                abstract=abstract,
                disease_terms=terms,
                gene_pattern=gene_pattern,
                source_size_bytes=Path(record_path).stat().st_size
                if Path(record_path).exists()
                else 0,
            )
        )

    candidates.sort(key=_sort_key)
    start = max(0, offset)
    end = start + top_n if top_n and top_n > 0 else len(candidates)
    selected_pmids = {candidate.pmid for candidate in candidates[start:end]}
    for rank, candidate in enumerate(candidates, start=1):
        candidate.rank = rank
        candidate.selected = candidate.pmid in selected_pmids

    selected_markdown_files = [
        path
        for path in markdown_files
        if (pmid := extract_pmid_from_filename(path)) and pmid in selected_pmids
    ]
    selected_abstract_papers = [
        (pmid, record_path)
        for pmid, record_path in abstract_papers
        if pmid in selected_pmids
    ]

    # Preserve priority order for worker submission.
    rank_by_pmid = {candidate.pmid: candidate.rank for candidate in candidates}
    selected_markdown_files.sort(
        key=lambda path: rank_by_pmid.get(extract_pmid_from_filename(path) or "", 10**9)
    )
    selected_abstract_papers.sort(key=lambda item: rank_by_pmid.get(item[0], 10**9))

    result = PriorityResult(
        selected_markdown_files=selected_markdown_files,
        selected_abstract_papers=selected_abstract_papers,
        candidates=candidates,
        report_dir=report_dir,
    )
    if report_dir:
        write_priority_reports(result, report_dir)
    return result


def write_priority_reports(result: PriorityResult, report_dir: Path) -> None:
    """Write ranking artifacts used to audit top-N extraction runs."""

    report_dir.mkdir(parents=True, exist_ok=True)

    json_path = report_dir / "priority_candidates.json"
    json_path.write_text(
        json.dumps([asdict(candidate) for candidate in result.candidates], indent=2),
        encoding="utf-8",
    )

    pmids_path = report_dir / "priority_pmids.txt"
    pmids_path.write_text(
        "\n".join(candidate.pmid for candidate in result.selected_candidates) + "\n",
        encoding="utf-8",
    )

    tsv_path = report_dir / "priority_candidates.tsv"
    fieldnames = [
        "rank",
        "selected",
        "pmid",
        "source_kind",
        "score",
        "year",
        "title",
        "journal",
        "text_chars",
        "source_size_bytes",
        "likely_full_text",
        "disease_term_hits",
        "variant_mentions",
        "gene_title_mentions",
        "gene_abstract_mentions",
        "gene_body_mentions",
        "gene_variant_lines",
        "gene_carrier_lines",
        "carrier_count_mentions",
        "original_data_mentions",
        "table_mentions",
        "review_title_penalty",
        "reasons",
        "source_file",
    ]
    with tsv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for candidate in result.candidates:
            signals = candidate.signals
            writer.writerow(
                {
                    "rank": candidate.rank,
                    "selected": int(candidate.selected),
                    "pmid": candidate.pmid,
                    "source_kind": candidate.source_kind,
                    "score": candidate.score,
                    "year": candidate.year,
                    "title": candidate.title,
                    "journal": candidate.journal,
                    "text_chars": candidate.text_chars,
                    "source_size_bytes": signals.get("source_size_bytes", 0),
                    "likely_full_text": int(bool(signals.get("likely_full_text"))),
                    "disease_term_hits": signals.get("disease_term_hits", 0),
                    "variant_mentions": signals.get("variant_mentions", 0),
                    "gene_title_mentions": signals.get("gene_title_mentions", 0),
                    "gene_abstract_mentions": signals.get("gene_abstract_mentions", 0),
                    "gene_body_mentions": signals.get("gene_body_mentions", 0),
                    "gene_variant_lines": signals.get("gene_variant_lines", 0),
                    "gene_carrier_lines": signals.get("gene_carrier_lines", 0),
                    "carrier_count_mentions": signals.get("carrier_count_mentions", 0),
                    "original_data_mentions": signals.get("original_data_mentions", 0),
                    "table_mentions": signals.get("table_mentions", 0),
                    "review_title_penalty": int(
                        bool(signals.get("review_title_penalty"))
                    ),
                    "reasons": "; ".join(candidate.reasons),
                    "source_file": candidate.source_file,
                }
            )
