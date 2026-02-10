"""Structured failure taxonomy logging for extraction pipeline.

Tracks why extractions fail with categorized failure codes,
enabling systematic analysis and improvement of extraction recall.
"""

import json
import logging
from dataclasses import asdict, dataclass, field
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


class FailureCode(Enum):
    """Categorized failure codes for extraction failures."""

    PAYWALL = "PAYWALL"
    NO_FULL_TEXT = "NO_FULL_TEXT"
    XML_PARSE_ERROR = "XML_PARSE_ERROR"
    SUPPLEMENT_MISSING = "SUPPLEMENT_MISSING"
    LLM_TIMEOUT = "LLM_TIMEOUT"
    NO_VARIANTS_FOUND = "NO_VARIANTS_FOUND"
    LLM_INPUT_ERROR = "LLM_INPUT_ERROR"
    INCORRECT_JSON = "INCORRECT_JSON"
    TEXT_TOO_SHORT = "TEXT_TOO_SHORT"
    METADATA_ERROR = "METADATA_ERROR"
    HTML_GARBAGE = "HTML_GARBAGE"


@dataclass
class ExtractionFailure:
    """A single extraction failure record."""

    pmid: str
    failure_code: FailureCode
    description: str
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())
    model_used: Optional[str] = None
    additional_info: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["failure_code"] = self.failure_code.value
        return d

    @classmethod
    def from_dict(cls, data: dict) -> "ExtractionFailure":
        data = dict(data)
        data["failure_code"] = FailureCode(data["failure_code"])
        return cls(**data)


class FailureLogger:
    """Logs and queries structured extraction failures.

    Failures are persisted to a JSON file for cross-run analysis.
    """

    def __init__(self, output_path: Optional[Path] = None):
        if output_path is None:
            output_path = Path("gvf_output") / "extraction_failures.json"
        self.output_path = Path(output_path)
        self._failures: List[ExtractionFailure] = []
        self._load()

    def _load(self) -> None:
        """Load existing failures from disk."""
        if self.output_path.exists():
            try:
                data = json.loads(self.output_path.read_text())
                self._failures = [ExtractionFailure.from_dict(f) for f in data]
                logger.debug("Loaded %d existing failures", len(self._failures))
            except (json.JSONDecodeError, KeyError) as e:
                logger.warning("Could not load failure log: %s", e)
                self._failures = []

    def _save(self) -> None:
        """Persist failures to disk."""
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        self.output_path.write_text(
            json.dumps([f.to_dict() for f in self._failures], indent=2)
        )

    def log_failure(
        self,
        pmid: str,
        code: FailureCode,
        description: str,
        **kwargs: Any,
    ) -> ExtractionFailure:
        """Log a new extraction failure.

        Args:
            pmid: PubMed ID of the failed paper.
            code: Categorized failure code.
            description: Human-readable description of what went wrong.
            **kwargs: Optional fields - model_used, additional_info, or
                      any extra keys folded into additional_info.
        """
        model_used = kwargs.pop("model_used", None)
        additional_info = kwargs.pop("additional_info", {})
        # Fold remaining kwargs into additional_info
        additional_info.update(kwargs)

        failure = ExtractionFailure(
            pmid=pmid,
            failure_code=code,
            description=description,
            model_used=model_used,
            additional_info=additional_info,
        )
        self._failures.append(failure)
        self._save()
        logger.info(
            "Logged failure for PMID %s: %s - %s", pmid, code.value, description
        )
        return failure

    def get_failures_by_code(self, code: FailureCode) -> List[ExtractionFailure]:
        """Return all failures matching a specific code."""
        return [f for f in self._failures if f.failure_code == code]

    def get_failure_stats(self) -> Dict[str, int]:
        """Return count of failures grouped by failure code."""
        stats: Dict[str, int] = {}
        for f in self._failures:
            key = f.failure_code.value
            stats[key] = stats.get(key, 0) + 1
        return stats

    def get_failures_for_pmid(self, pmid: str) -> List[ExtractionFailure]:
        """Return all failures for a specific PMID."""
        return [f for f in self._failures if f.pmid == pmid]

    @property
    def failures(self) -> List[ExtractionFailure]:
        """All recorded failures."""
        return list(self._failures)
