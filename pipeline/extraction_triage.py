"""Cheap pre-extraction triage for priority-ranked papers."""

from __future__ import annotations

import csv
import json
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional

from pipeline.extraction_priority import (
    CARRIER_COUNT_RE,
    ORIGINAL_DATA_RE,
    TABLE_RE,
    VARIANT_PATTERNS,
    ExtractionCandidate,
    PriorityResult,
    _gene_regex,
    _read_source_sample,
)
from utils.llm_utils import BaseLLMCaller
from utils.pmid_utils import extract_pmid_from_filename

logger = logging.getLogger(__name__)

TRIAGE_DECISIONS = {"extract_now", "defer", "skip"}


@dataclass
class TriageDecision:
    pmid: str
    decision: str
    confidence: float
    priority_rank: int
    source_kind: str
    source_file: str
    model: str = "deterministic"
    deterministic_decision: str = ""
    deterministic_confidence: float = 0.0
    llm_decision: str = ""
    llm_confidence: float = 0.0
    useful_info_types: list[str] = field(default_factory=list)
    false_positive_risks: list[str] = field(default_factory=list)
    reasons: list[str] = field(default_factory=list)


@dataclass
class TriageResult:
    decisions: list[TriageDecision]
    report_dir: Optional[Path] = None

    @property
    def extract_pmids(self) -> set[str]:
        return {d.pmid for d in self.decisions if d.decision == "extract_now"}

    @property
    def defer_pmids(self) -> set[str]:
        return {d.pmid for d in self.decisions if d.decision == "defer"}


def deterministic_triage(candidate: ExtractionCandidate) -> TriageDecision:
    """Cheap rule-based triage calibrated to avoid APOE-like association noise."""

    signals = candidate.signals
    reasons: list[str] = []
    risks: list[str] = []
    useful: list[str] = []

    gene_mentions = (
        int(signals.get("gene_title_mentions") or 0)
        + int(signals.get("gene_abstract_mentions") or 0)
        + int(signals.get("gene_body_mentions") or 0)
    )
    gene_variant_lines = int(signals.get("gene_variant_lines") or 0)
    gene_carrier_lines = int(signals.get("gene_carrier_lines") or 0)
    table_mentions = int(signals.get("table_mentions") or 0)
    original_mentions = int(signals.get("original_data_mentions") or 0)
    variant_mentions = int(signals.get("variant_mentions") or 0)
    carrier_mentions = int(signals.get("carrier_count_mentions") or 0)

    if signals.get("oversized_source"):
        return TriageDecision(
            pmid=candidate.pmid,
            decision="skip",
            confidence=0.98,
            priority_rank=candidate.rank,
            source_kind=candidate.source_kind,
            source_file=candidate.source_file,
            deterministic_decision="skip",
            deterministic_confidence=0.98,
            false_positive_risks=["oversized_source"],
            reasons=["source is too large for reliable cheap triage/extraction"],
        )

    if signals.get("review_title_penalty") and gene_variant_lines == 0:
        return TriageDecision(
            pmid=candidate.pmid,
            decision="skip",
            confidence=0.9,
            priority_rank=candidate.rank,
            source_kind=candidate.source_kind,
            source_file=candidate.source_file,
            deterministic_decision="skip",
            deterministic_confidence=0.9,
            false_positive_risks=["review_or_protocol"],
            reasons=[
                "review/editorial/protocol title without target-gene variant evidence"
            ],
        )

    if not gene_mentions:
        return TriageDecision(
            pmid=candidate.pmid,
            decision="skip",
            confidence=0.9,
            priority_rank=candidate.rank,
            source_kind=candidate.source_kind,
            source_file=candidate.source_file,
            deterministic_decision="skip",
            deterministic_confidence=0.9,
            false_positive_risks=["no_target_gene_mention"],
            reasons=["no target-gene mention in title/abstract/sample"],
        )

    if (
        gene_variant_lines
        and gene_carrier_lines
        and (table_mentions or original_mentions)
    ):
        useful.extend(["target_gene_variants", "carrier_or_genotype_counts"])
        if table_mentions:
            useful.append("tables_or_supplements")
        return TriageDecision(
            pmid=candidate.pmid,
            decision="extract_now",
            confidence=0.88,
            priority_rank=candidate.rank,
            source_kind=candidate.source_kind,
            source_file=candidate.source_file,
            deterministic_decision="extract_now",
            deterministic_confidence=0.88,
            useful_info_types=useful,
            reasons=[
                "target gene co-occurs with variant and carrier/count evidence",
                "original/table signal present",
            ],
        )

    if gene_variant_lines and (table_mentions or carrier_mentions):
        return TriageDecision(
            pmid=candidate.pmid,
            decision="extract_now",
            confidence=0.74,
            priority_rank=candidate.rank,
            source_kind=candidate.source_kind,
            source_file=candidate.source_file,
            deterministic_decision="extract_now",
            deterministic_confidence=0.74,
            useful_info_types=["target_gene_variants"],
            reasons=[
                "target gene co-occurs with variant evidence and table/count signal"
            ],
        )

    if variant_mentions and not gene_variant_lines:
        risks.append("variant_terms_not_tied_to_target_gene")
    if gene_mentions and not (gene_variant_lines or gene_carrier_lines):
        risks.append("gene_mentioned_without_variant_count_evidence")

    decision = "defer" if gene_mentions else "skip"
    confidence = 0.62 if decision == "defer" else 0.8
    return TriageDecision(
        pmid=candidate.pmid,
        decision=decision,
        confidence=confidence,
        priority_rank=candidate.rank,
        source_kind=candidate.source_kind,
        source_file=candidate.source_file,
        deterministic_decision=decision,
        deterministic_confidence=confidence,
        false_positive_risks=risks,
        reasons=["ambiguous target-gene extraction value"],
    )


def _evidence_snippets(
    candidate: ExtractionCandidate, gene_symbol: str, max_chars: int
) -> str:
    path = Path(candidate.source_file)
    if not path.exists():
        return ""
    text = _read_source_sample(path, max_chars=max_chars)
    gene_pattern = _gene_regex(gene_symbol)
    lines = []
    for line in text.splitlines():
        lower_line = line.lower()
        has_gene = bool(gene_pattern.search(line))
        has_signal = (
            any(pattern.search(line) for pattern in VARIANT_PATTERNS)
            or CARRIER_COUNT_RE.search(line)
            or TABLE_RE.search(line)
            or ORIGINAL_DATA_RE.search(line)
        )
        if has_gene and has_signal:
            lines.append(line.strip())
        elif has_gene and len(lines) < 8:
            lines.append(line.strip())
        elif (
            has_signal
            and ("table" in lower_line or "carrier" in lower_line)
            and len(lines) < 12
        ):
            lines.append(line.strip())
        if len(lines) >= 18:
            break
    return "\n".join(line[:500] for line in lines if line)[:max_chars]


class ExtractionTriageClassifier(BaseLLMCaller):
    SYSTEM_MESSAGE = (
        "You are a biomedical literature triage classifier. Decide whether a paper "
        "should receive expensive variant/carrier extraction for the target gene. "
        "Be skeptical of review articles, protocols, functional-only papers, and "
        "association papers where variants belong to another gene. A paper is useful "
        "when it likely contains original target-gene variant/genotype/allele data, "
        "carrier counts, case/control counts, family/proband data, or extractable "
        "tables for the requested disease context. Return JSON only."
    )

    PROMPT = """Target gene: {gene}
Disease context: {disease}

PMID: {pmid}
Title: {title}
Journal/year: {journal} / {year}
Priority score: {score}
Deterministic signal summary:
{signals}

Evidence snippets from abstract/full text/table-like lines:
{snippets}

Classify this paper for expensive extraction.

Return a JSON object:
{{
  "decision": "extract_now" | "defer" | "skip",
  "confidence": 0.0-1.0,
  "useful_info_types": ["target_gene_variants", "carrier_or_genotype_counts", "case_control_counts", "family_or_proband_data", "tables_or_supplements"],
  "false_positive_risks": ["review_or_protocol", "variant_terms_not_tied_to_target_gene", "functional_only", "wrong_disease", "no_counts"],
  "reason": "short explanation"
}}

Use "extract_now" only when target-gene data is likely enough to justify the expensive extractor.
Use "defer" when it may be useful but evidence is weak.
Use "skip" when it is unlikely to contain extractable target-gene carrier/variant evidence.
"""

    def classify(
        self,
        candidate: ExtractionCandidate,
        *,
        gene_symbol: str,
        disease: str = "",
        max_snippet_chars: int = 6000,
    ) -> TriageDecision:
        signals = {
            key: candidate.signals.get(key)
            for key in (
                "gene_title_mentions",
                "gene_abstract_mentions",
                "gene_body_mentions",
                "gene_variant_lines",
                "gene_carrier_lines",
                "variant_mentions",
                "carrier_count_mentions",
                "original_data_mentions",
                "table_mentions",
                "review_title_penalty",
                "oversized_source",
            )
        }
        prompt = self.PROMPT.format(
            gene=gene_symbol,
            disease=disease or "(none)",
            pmid=candidate.pmid,
            title=candidate.title,
            journal=candidate.journal,
            year=candidate.year,
            score=candidate.score,
            signals=json.dumps(signals, indent=2),
            snippets=_evidence_snippets(candidate, gene_symbol, max_snippet_chars),
        )
        data = self.call_llm_json(prompt, system_message=self.SYSTEM_MESSAGE)
        decision = str(data.get("decision") or "").strip().lower()
        if decision not in TRIAGE_DECISIONS:
            decision = "defer"
        try:
            confidence = float(data.get("confidence", 0.5))
        except (TypeError, ValueError):
            confidence = 0.5
        confidence = max(0.0, min(1.0, confidence))
        useful = data.get("useful_info_types") or []
        risks = data.get("false_positive_risks") or []
        return TriageDecision(
            pmid=candidate.pmid,
            decision=decision,
            confidence=confidence,
            priority_rank=candidate.rank,
            source_kind=candidate.source_kind,
            source_file=candidate.source_file,
            model=self.model,
            llm_decision=decision,
            llm_confidence=confidence,
            useful_info_types=[str(item) for item in useful if str(item).strip()],
            false_positive_risks=[str(item) for item in risks if str(item).strip()],
            reasons=[str(data.get("reason") or "LLM triage decision")],
        )


def _merge_decisions(
    deterministic: TriageDecision,
    llm_decision: Optional[TriageDecision],
    *,
    mode: str,
) -> TriageDecision:
    if llm_decision is None or mode == "deterministic":
        return deterministic
    if mode == "llm":
        llm_decision.deterministic_decision = deterministic.decision
        llm_decision.deterministic_confidence = deterministic.confidence
        return llm_decision

    # Hybrid: trust strong deterministic skips for obvious non-target/review
    # cases, otherwise let the cheap model make the final call.
    if deterministic.decision == "skip" and deterministic.confidence >= 0.9:
        deterministic.llm_decision = llm_decision.decision
        deterministic.llm_confidence = llm_decision.confidence
        deterministic.model = llm_decision.model
        deterministic.reasons.extend(
            [
                f"LLM would choose {llm_decision.decision}: {', '.join(llm_decision.reasons)}"
            ]
        )
        return deterministic
    llm_decision.deterministic_decision = deterministic.decision
    llm_decision.deterministic_confidence = deterministic.confidence
    return llm_decision


def triage_priority_result(
    priority_result: PriorityResult,
    *,
    gene_symbol: str,
    disease: str = "",
    mode: str = "hybrid",
    model: Optional[str] = None,
    max_llm_candidates: Optional[int] = None,
    report_dir: Optional[Path] = None,
    max_snippet_chars: int = 6000,
) -> TriageResult:
    """Run cheap triage over the currently selected priority candidates."""

    mode = (mode or "hybrid").lower()
    if mode not in {"deterministic", "llm", "hybrid"}:
        raise ValueError(f"Unsupported triage mode: {mode}")

    classifier: Optional[ExtractionTriageClassifier] = None
    if mode in {"llm", "hybrid"}:
        if model is None:
            from config.settings import get_settings

            model = get_settings().get_tier2_model()
        classifier = ExtractionTriageClassifier(
            model=model,
            temperature=0.0,
            max_tokens=1200,
        )

    decisions: list[TriageDecision] = []
    selected = priority_result.selected_candidates
    llm_calls = 0
    for candidate in selected:
        deterministic = deterministic_triage(candidate)
        llm_decision = None
        should_call_llm = classifier is not None
        if (
            mode == "hybrid"
            and deterministic.decision == "skip"
            and deterministic.confidence >= 0.9
        ):
            should_call_llm = False
        if max_llm_candidates is not None and llm_calls >= max_llm_candidates:
            should_call_llm = False
        if should_call_llm and classifier is not None:
            try:
                llm_decision = classifier.classify(
                    candidate,
                    gene_symbol=gene_symbol,
                    disease=disease,
                    max_snippet_chars=max_snippet_chars,
                )
                llm_calls += 1
            except Exception as exc:
                logger.warning(
                    "Triage LLM failed for PMID %s: %s; using deterministic decision",
                    candidate.pmid,
                    exc,
                )
        decisions.append(_merge_decisions(deterministic, llm_decision, mode=mode))

    result = TriageResult(decisions=decisions, report_dir=report_dir)
    if report_dir:
        write_triage_reports(result, report_dir)
    return result


def apply_triage_filter(
    priority_result: PriorityResult,
    triage_result: TriageResult,
    *,
    include_defer: bool = False,
) -> PriorityResult:
    """Filter a priority result to papers triage selected for extraction."""

    allowed = set(triage_result.extract_pmids)
    if include_defer:
        allowed |= triage_result.defer_pmids
    for candidate in priority_result.candidates:
        candidate.selected = candidate.pmid in allowed

    selected_markdown = [
        path
        for path in priority_result.selected_markdown_files
        if (pmid := extract_pmid_from_filename(path)) and pmid in allowed
    ]
    selected_abstracts = [
        (pmid, path)
        for pmid, path in priority_result.selected_abstract_papers
        if pmid in allowed
    ]
    rank_by_pmid = {
        candidate.pmid: candidate.rank for candidate in priority_result.candidates
    }
    selected_markdown.sort(
        key=lambda path: rank_by_pmid.get(extract_pmid_from_filename(path) or "", 10**9)
    )
    selected_abstracts.sort(key=lambda item: rank_by_pmid.get(item[0], 10**9))
    return PriorityResult(
        selected_markdown_files=selected_markdown,
        selected_abstract_papers=selected_abstracts,
        candidates=priority_result.candidates,
        report_dir=priority_result.report_dir,
    )


def write_triage_reports(result: TriageResult, report_dir: Path) -> None:
    report_dir.mkdir(parents=True, exist_ok=True)
    decisions = sorted(result.decisions, key=lambda d: d.priority_rank)

    (report_dir / "triage_candidates.json").write_text(
        json.dumps([asdict(decision) for decision in decisions], indent=2),
        encoding="utf-8",
    )
    (report_dir / "triage_pmids.txt").write_text(
        "\n".join(
            decision.pmid
            for decision in decisions
            if decision.decision == "extract_now"
        )
        + "\n",
        encoding="utf-8",
    )

    summary = {
        "total": len(decisions),
        "extract_now": sum(1 for d in decisions if d.decision == "extract_now"),
        "defer": sum(1 for d in decisions if d.decision == "defer"),
        "skip": sum(1 for d in decisions if d.decision == "skip"),
    }
    (report_dir / "triage_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    fieldnames = [
        "priority_rank",
        "pmid",
        "decision",
        "confidence",
        "deterministic_decision",
        "deterministic_confidence",
        "llm_decision",
        "llm_confidence",
        "model",
        "useful_info_types",
        "false_positive_risks",
        "reasons",
        "source_kind",
        "source_file",
    ]
    with (report_dir / "triage_candidates.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for decision in decisions:
            writer.writerow(
                {
                    "priority_rank": decision.priority_rank,
                    "pmid": decision.pmid,
                    "decision": decision.decision,
                    "confidence": decision.confidence,
                    "deterministic_decision": decision.deterministic_decision,
                    "deterministic_confidence": decision.deterministic_confidence,
                    "llm_decision": decision.llm_decision,
                    "llm_confidence": decision.llm_confidence,
                    "model": decision.model,
                    "useful_info_types": "; ".join(decision.useful_info_types),
                    "false_positive_risks": "; ".join(decision.false_positive_risks),
                    "reasons": "; ".join(decision.reasons),
                    "source_kind": decision.source_kind,
                    "source_file": decision.source_file,
                }
            )
