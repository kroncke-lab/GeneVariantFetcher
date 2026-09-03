"""The authoritative finalized per-paper source ledger for a run.

Why this exists
---------------
Run artifacts used to be written from whatever a pipeline *step* returned, not
from the finalized list of papers that extraction actually consumed. On a
resume path (``--skip extract`` / ``--resume-dir``) the download step
contributes nothing, so ``download_fulltext`` classified every paper as
abstract-only via ``abstract_only = [p for p in attempted if p not in
successfully_downloaded]``. One real run directory therefore shipped four
mutually contradictory statements about the same 50 papers:

* per-paper extraction records: ``source_type='fulltext'`` for 50/50
* ``source_completeness.json``: ``papers_with_fulltext=0``, 50 abstract-only
* ``RUN_STATUS.json``: ``full_text=50``, ``abstract_only=0``
* ``workflow_summary``: ``papers_downloaded=0`` and ``papers_from_fulltext=50``

A curator reading that cannot tell whether they are reviewing full-paper
extraction, abstract extraction, or a mix.

The rule this module enforces
-----------------------------
The denominator is the **finalized list**: the papers that produced an
extraction record. Each row states what extraction actually parsed (file,
digest, byte size) rather than what a step intended to fetch.

The declared class is not trusted on its own. ``source_type='fulltext'`` only
means extraction *ran* in full-text mode; it is not proof that the bytes were a
whole article. Every row therefore carries both the declared class and a class
re-derived from the recorded bytes, and any disagreement is surfaced as an
explicit ``discrepancies`` entry instead of being silently resolved. When the
recorded source file is no longer readable the row says ``unverified``; it
never upgrades to ``fulltext`` on the strength of the declaration alone.

Nothing here mutates a run. It reads finalized records and reports.
"""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional

from pipeline.source_quality import is_abstract_only_fallback_text

logger = logging.getLogger(__name__)

__all__ = [
    "SourceLedgerRow",
    "SourceLedger",
    "build_source_ledger",
    "ledger_from_extraction_dir",
]

#: Byte floor below which recorded source cannot be a full article body. Well
#: under any real paper; this only catches empty/truncated writes.
_EMPTY_SOURCE_FLOOR_BYTES = 512

CLASS_FULLTEXT = "fulltext"
CLASS_ABSTRACT_ONLY = "abstract_only"
CLASS_UNVERIFIED = "unverified"


@dataclass
class SourceLedgerRow:
    """One finalized paper: what extraction actually read."""

    pmid: str
    declared_class: str = CLASS_UNVERIFIED
    observed_class: str = CLASS_UNVERIFIED
    source_file: Optional[str] = None
    source_sha256: Optional[str] = None
    source_size_bytes: Optional[int] = None
    source_present: bool = False
    variants: int = 0
    counted_variants: int = 0
    max_carriers: Optional[int] = None
    discrepancy: Optional[str] = None

    @property
    def effective_class(self) -> str:
        """The class to report. Never better than what the bytes support."""
        if self.observed_class == CLASS_UNVERIFIED:
            return CLASS_UNVERIFIED
        return self.observed_class

    def as_dict(self) -> dict[str, Any]:
        return {
            "pmid": self.pmid,
            "declared_class": self.declared_class,
            "observed_class": self.observed_class,
            "effective_class": self.effective_class,
            "source_file": self.source_file,
            "source_sha256": self.source_sha256,
            "source_size_bytes": self.source_size_bytes,
            "source_present": self.source_present,
            "variants": self.variants,
            "counted_variants": self.counted_variants,
            "max_carriers": self.max_carriers,
            "discrepancy": self.discrepancy,
        }


@dataclass
class SourceLedger:
    """Finalized source ledger for one run. Lists are never truncated."""

    rows: list[SourceLedgerRow] = field(default_factory=list)
    requested_pmids: list[str] = field(default_factory=list)

    @property
    def finalized_pmids(self) -> list[str]:
        return [row.pmid for row in self.rows]

    @property
    def missing_pmids(self) -> list[str]:
        """Requested but with no extraction record. Real coverage loss."""
        finalized = set(self.finalized_pmids)
        return [p for p in self.requested_pmids if p not in finalized]

    def _by_class(self, wanted: str) -> list[str]:
        return [row.pmid for row in self.rows if row.effective_class == wanted]

    @property
    def fulltext_pmids(self) -> list[str]:
        return self._by_class(CLASS_FULLTEXT)

    @property
    def abstract_only_pmids(self) -> list[str]:
        return self._by_class(CLASS_ABSTRACT_ONLY)

    @property
    def unverified_pmids(self) -> list[str]:
        return self._by_class(CLASS_UNVERIFIED)

    @property
    def zero_variant_pmids(self) -> list[str]:
        return [row.pmid for row in self.rows if row.variants == 0]

    @property
    def single_carrier_pmids(self) -> list[str]:
        """Papers whose every counted variant reports at most one carrier.

        The predicate this replaces read ``variant["carriers_total"]`` and
        ``variant["total_carriers"]``. Neither key exists in the extraction
        schema, so the expression collapsed to its ``or 1`` default and every
        paper in every run was flagged. A paper with no counted variant at all
        is not evidence of a missed cohort table and is excluded.
        """
        return [
            row.pmid
            for row in self.rows
            if row.counted_variants > 0
            and row.max_carriers is not None
            and row.max_carriers <= 1
        ]

    @property
    def discrepancies(self) -> list[dict[str, str]]:
        return [
            {"pmid": row.pmid, "discrepancy": row.discrepancy}
            for row in self.rows
            if row.discrepancy
        ]

    def as_completeness_dict(self) -> dict[str, Any]:
        """Render the ledger. Every count is ``len()`` of the list beside it."""
        total = len(self.rows)
        fulltext = self.fulltext_pmids
        abstract_only = self.abstract_only_pmids
        unverified = self.unverified_pmids
        return {
            "denominator_source": "finalized_extraction_records",
            "papers_finalized": total,
            "papers_requested": len(self.requested_pmids),
            "papers_missing_extraction": len(self.missing_pmids),
            "missing_extraction_pmids": self.missing_pmids,
            "papers_with_fulltext": len(fulltext),
            "fulltext_pmids": fulltext,
            "papers_abstract_only": len(abstract_only),
            "abstract_only_pmids": abstract_only,
            "papers_source_unverified": len(unverified),
            "source_unverified_pmids": unverified,
            "zero_variant_papers": len(self.zero_variant_pmids),
            "zero_variant_pmids": self.zero_variant_pmids,
            "single_carrier_papers": len(self.single_carrier_pmids),
            "single_carrier_pmids": self.single_carrier_pmids,
            "fulltext_coverage_pct": (
                round(len(fulltext) / total * 100, 1) if total else 0.0
            ),
            "source_class_discrepancies": self.discrepancies,
            "rows": [row.as_dict() for row in self.rows],
        }


def _coerce_int(value: Any) -> Optional[int]:
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float) and value.is_integer():
        return int(value)
    if isinstance(value, str) and value.strip().lstrip("-").isdigit():
        return int(value.strip())
    return None


def _classify_recorded_bytes(
    path: Optional[Path], recorded_size: Optional[int]
) -> tuple[str, bool]:
    """Classify from the bytes extraction recorded, not from its declaration.

    Returns ``(observed_class, source_present)``. An unreadable or absent file
    yields ``unverified``: the declaration alone must never promote a row.
    """
    if path is None:
        return CLASS_UNVERIFIED, False
    try:
        if not path.is_file():
            return CLASS_UNVERIFIED, False
        head = path.read_text(encoding="utf-8", errors="replace")[:8192]
        size = path.stat().st_size
    except OSError:
        return CLASS_UNVERIFIED, False
    if size < _EMPTY_SOURCE_FLOOR_BYTES:
        return CLASS_ABSTRACT_ONLY, True
    if is_abstract_only_fallback_text(head):
        return CLASS_ABSTRACT_ONLY, True
    if recorded_size is not None and recorded_size < _EMPTY_SOURCE_FLOOR_BYTES:
        return CLASS_ABSTRACT_ONLY, True
    return CLASS_FULLTEXT, True


def _resolve_source_path(raw: Any, search_dirs: Iterable[Path]) -> Optional[Path]:
    if not raw:
        return None
    candidate = Path(str(raw))
    if candidate.is_file():
        return candidate
    for base in search_dirs:
        for probe in (base / candidate.name, base / candidate):
            try:
                if probe.is_file():
                    return probe
            except OSError:
                continue
    return None


def build_source_ledger(
    records: Iterable[dict[str, Any]],
    *,
    requested_pmids: Optional[Iterable[str]] = None,
    search_dirs: Optional[Iterable[Path]] = None,
) -> SourceLedger:
    """Build the finalized ledger from per-paper extraction payloads.

    ``records`` are extraction JSON payloads (``paper_metadata`` +
    ``extraction_metadata`` + ``variants``). One row per finalized paper.
    """
    from pipeline.phenotype_count_guard import read_phenotype_counts

    dirs = list(search_dirs or [])
    ledger = SourceLedger(requested_pmids=[str(p) for p in (requested_pmids or [])])
    seen: set[str] = set()

    for record in records:
        if not isinstance(record, dict):
            continue
        meta = record.get("extraction_metadata")
        meta = meta if isinstance(meta, dict) else {}
        paper = record.get("paper_metadata")
        paper = paper if isinstance(paper, dict) else {}
        pmid = str(paper.get("pmid") or meta.get("pmid") or "").strip()
        if not pmid or pmid in seen:
            continue
        seen.add(pmid)

        declared = str(meta.get("source_type") or "").strip().lower()
        if meta.get("abstract_only"):
            declared = CLASS_ABSTRACT_ONLY
        if declared not in {CLASS_FULLTEXT, CLASS_ABSTRACT_ONLY}:
            declared = CLASS_UNVERIFIED

        recorded_size = _coerce_int(meta.get("source_size_bytes"))
        source_path = _resolve_source_path(meta.get("source_file"), dirs)
        observed, present = _classify_recorded_bytes(source_path, recorded_size)

        variants = record.get("variants")
        variants = variants if isinstance(variants, list) else []
        counted = 0
        max_carriers: Optional[int] = None
        for variant in variants:
            if not isinstance(variant, dict):
                continue
            carriers = read_phenotype_counts(variant).get("carriers")
            if carriers is None:
                continue
            counted += 1
            max_carriers = (
                carriers if max_carriers is None else max(max_carriers, carriers)
            )

        discrepancy = None
        if (
            declared != CLASS_UNVERIFIED
            and observed != CLASS_UNVERIFIED
            and declared != observed
        ):
            discrepancy = (
                f"extraction declared {declared!r} but recorded source bytes "
                f"read as {observed!r}"
            )
        elif observed == CLASS_UNVERIFIED and declared != CLASS_UNVERIFIED:
            # Applies to a declared abstract-only paper too: whichever class was
            # declared, an unreadable recorded source means the run can no
            # longer show what it parsed, and that must be visible rather than
            # quietly landing in the unverified bucket.
            discrepancy = (
                f"extraction declared {declared!r} but its recorded source file "
                "is not readable now, so the class cannot be confirmed"
            )

        ledger.rows.append(
            SourceLedgerRow(
                pmid=pmid,
                declared_class=declared,
                observed_class=observed,
                source_file=str(meta.get("source_file") or "") or None,
                source_sha256=str(meta.get("source_sha256") or "") or None,
                source_size_bytes=recorded_size,
                source_present=present,
                variants=len(variants),
                counted_variants=counted,
                max_carriers=max_carriers,
                discrepancy=discrepancy,
            )
        )

    ledger.rows.sort(key=lambda row: row.pmid)
    return ledger


def ledger_from_extraction_dir(
    extraction_dir: Path,
    *,
    requested_pmids: Optional[Iterable[str]] = None,
    search_dirs: Optional[Iterable[Path]] = None,
) -> SourceLedger:
    """Build the ledger by reading a run's finalized extraction JSON files."""
    records: list[dict[str, Any]] = []
    if extraction_dir and Path(extraction_dir).is_dir():
        for path in sorted(Path(extraction_dir).glob("*.json")):
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError) as exc:
                logger.warning("source ledger: unreadable extraction %s: %s", path, exc)
                continue
            if isinstance(payload, dict):
                records.append(payload)
    return build_source_ledger(
        records, requested_pmids=requested_pmids, search_dirs=search_dirs
    )


def file_digest(path: Path) -> Optional[str]:
    """SHA-256 of a file, or None when unreadable. Never raises."""
    try:
        hasher = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                hasher.update(chunk)
        return hasher.hexdigest()
    except OSError:
        return None
