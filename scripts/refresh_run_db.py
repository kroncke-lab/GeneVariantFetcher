#!/usr/bin/env python3
"""Refresh an existing GVF run after extractor/source-handling fixes.

This is the safe alternative to patching SQLite rows in place:

1. Select PMIDs whose extraction JSON is stale or under-counted relative to
   the current source files and deterministic parsers.
2. Re-extract only those PMIDs from source markdown into canonical
   ``<GENE>_PMID_<PMID>.json`` files, backing up replaced JSONs.
3. Rebuild a fresh SQLite DB from the complete extraction directory.
4. Optionally rerun recovery layers against the fresh DB. Gold is optional; if
   no gold exists, recovery still runs in DB-PMID mode and scoring is skipped.

The script is gene-agnostic and does not use gold standards to choose PMIDs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from harvesting.migrate_to_sqlite import (  # noqa: E402
    create_database_schema,
    migrate_extraction_directory,
)
from pipeline.extraction import ExpertExtractor  # noqa: E402
from pipeline.source_quality import is_usable_fulltext_source  # noqa: E402
from utils.models import Paper  # noqa: E402

logger = logging.getLogger("refresh_run_db")

DETERMINISTIC_ABSOLUTE_LIFT_OVERRIDE = 50


@dataclass
class ReplayCandidate:
    pmid: str
    source_file: Path
    output_file: Path
    current_variants: int
    deterministic_variants: int
    reasons: list[str] = field(default_factory=list)
    # Distinct codon positions the gene-scoped deterministic TABLE parse found in
    # this source — the trusted structural baseline the gate scores a re-extraction
    # against (covering them = faithful; under-covering = a lossy under-extraction).
    # Empty when the selector found no deterministic table (e.g. prose-only papers).
    deterministic_positions: frozenset[int] = frozenset()


def _json_load(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _read_pmid_file(path: Path) -> set[str]:
    pmids: set[str] = set()
    with path.open(encoding="utf-8") as f:
        for line in f:
            value = line.strip()
            if not value or value.startswith("#"):
                continue
            pmids.add(value)
    return pmids


def _split_model_args(values: list[str]) -> list[str]:
    models: list[str] = []
    for value in values:
        for model in str(value or "").split(","):
            model = model.strip()
            if model:
                models.append(model)
    return models


def load_report_pmids(
    *,
    report: Path,
    gene: str,
    failure_classes: set[str],
    min_missing_rows: int,
    max_row_recall: float | None = None,
) -> set[str]:
    """Load PMIDs from a paper disagreement report for targeted replay."""
    selected: set[str] = set()
    with report.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if not _report_row_selected(
                row,
                gene=gene,
                failure_classes=failure_classes,
                min_missing_rows=min_missing_rows,
                max_row_recall=max_row_recall,
            ):
                continue
            pmid = str(row.get("pmid") or "").strip()
            if pmid:
                selected.add(pmid)
    return selected


def _report_row_selected(
    row: dict[str, str],
    *,
    gene: str,
    failure_classes: set[str],
    min_missing_rows: int,
    max_row_recall: float | None,
) -> bool:
    if str(row.get("gene") or "").upper() != gene.upper():
        return False
    failure_class = str(row.get("failure_class") or "")
    if failure_classes and failure_class not in failure_classes:
        return False
    try:
        missing_rows = int(row.get("missing_rows") or 0)
    except ValueError:
        missing_rows = 0
    if missing_rows < min_missing_rows:
        return False
    if max_row_recall is not None:
        try:
            row_recall = float(row.get("row_recall") or 0)
        except ValueError:
            row_recall = 0
        if row_recall > max_row_recall:
            return False
    return True


def load_report_available_contexts(
    *,
    report: Path,
    gene: str,
    failure_classes: set[str],
    min_missing_rows: int,
    max_row_recall: float | None = None,
    context_search_roots: list[Path] | None = None,
) -> dict[str, Path]:
    """Load PMID -> available_context_path from a disagreement report."""
    selected: dict[str, Path] = {}
    with report.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if not _report_row_selected(
                row,
                gene=gene,
                failure_classes=failure_classes,
                min_missing_rows=min_missing_rows,
                max_row_recall=max_row_recall,
            ):
                continue
            pmid = str(row.get("pmid") or "").strip()
            raw_path = str(
                row.get("available_context_path") or row.get("context_path") or ""
            ).strip()
            if not pmid:
                continue
            path = _largest_context_path(
                pmid=pmid,
                raw_path=raw_path,
                search_roots=context_search_roots or [],
            )
            if path is None:
                continue
            current = selected.get(pmid)
            current_size = (
                current.stat().st_size if current and current.exists() else -1
            )
            if path.stat().st_size > current_size:
                selected[pmid] = path
    return selected


def load_source_override_csv(path: Path) -> dict[str, Path]:
    """Load PMID -> source path from a route-filtered acquisition worklist.

    The CSV can be the output of
    ``scripts/recall_audit/build_acquisition_worklist.py``. If it has an
    ``action`` column, only ``refresh_replay`` rows are used; otherwise every
    row with a PMID and usable source path is considered.
    """
    selected: dict[str, Path] = {}
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            return selected
        lower_to_field = {name.strip().lower(): name for name in reader.fieldnames}
        pmid_col = lower_to_field.get("pmid")
        if pmid_col is None:
            raise ValueError(f"No PMID column found in {path}")
        action_col = lower_to_field.get("action")
        path_col = None
        for candidate in (
            "available_context_path",
            "source_file",
            "context_path",
            "full_context_path",
        ):
            if candidate in lower_to_field:
                path_col = lower_to_field[candidate]
                break
        if path_col is None:
            raise ValueError(f"No source path column found in {path}")
        for row in reader:
            action = str(row.get(action_col) or "").strip() if action_col else ""
            if action_col and action != "refresh_replay":
                continue
            pmid = str(row.get(pmid_col) or "").strip()
            raw_path = str(row.get(path_col) or "").strip()
            if not pmid or not raw_path:
                continue
            source_path = Path(raw_path).expanduser()
            if not source_path.exists() or not is_usable_fulltext_source(source_path):
                continue
            current = selected.get(pmid)
            current_size = (
                current.stat().st_size if current and current.exists() else -1
            )
            if source_path.stat().st_size > current_size:
                selected[pmid] = source_path
    return selected


def _largest_context_path(
    *,
    pmid: str,
    raw_path: str,
    search_roots: list[Path],
) -> Path | None:
    candidates: list[Path] = []
    if raw_path:
        candidates.append(Path(raw_path).expanduser())
    for root in search_roots:
        root = root.expanduser()
        if root.exists():
            candidates.extend(root.rglob(f"{pmid}_FULL_CONTEXT.md"))
    usable = [
        path for path in candidates if path.exists() and is_usable_fulltext_source(path)
    ]
    if not usable:
        return None
    return max(usable, key=lambda path: path.stat().st_size)


def _variant_count(data: dict[str, Any]) -> int:
    variants = data.get("variants")
    if isinstance(variants, list):
        return len(variants)
    meta_count = (data.get("extraction_metadata") or {}).get("total_variants_found")
    try:
        return int(meta_count)
    except (TypeError, ValueError):
        return 0


def _metadata_mentions_abstract_only(metadata: dict[str, Any]) -> bool:
    fields: list[str] = []
    for key in ("source_type", "source_kind", "notes"):
        value = metadata.get(key)
        if value:
            fields.append(str(value))
    challenges = metadata.get("challenges")
    if isinstance(challenges, list):
        fields.extend(str(item) for item in challenges)
    elif challenges:
        fields.append(str(challenges))
    joined = " ".join(fields).lower()
    return (
        bool(metadata.get("abstract_only"))
        or "abstract-only" in joined
        or "abstract only" in joined
        or "full text not available" in joined
        or "full text could not be retrieved" in joined
    )


def _metadata_source_is_unbound(
    metadata: dict[str, Any],
    *,
    output_file: Path,
    extraction_dir: Path,
) -> bool:
    """True when metadata records the JSON artifact instead of its source text."""
    source_file = metadata.get("source_file")
    if not source_file:
        return False

    source_path = Path(str(source_file)).expanduser()
    if source_path.suffix.lower() == ".json":
        return True

    try:
        source_path.resolve().relative_to(extraction_dir.resolve())
        return True
    except (OSError, ValueError):
        pass

    try:
        return source_path.resolve() == output_file.resolve()
    except OSError:
        return False


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _source_metadata(source_file: Path) -> dict[str, Any]:
    stat = source_file.stat()
    return {
        "source_file": str(source_file),
        "source_size_bytes": stat.st_size,
        "source_sha256": _sha256(source_file),
        "source_type": "fulltext",
        "abstract_only": False,
    }


def fold_on_disk_supplements(harvest_dir: Path) -> set[str]:
    """Fold any on-disk ``{pmid}_supplements/`` into ``{pmid}_FULL_CONTEXT.md``.

    The re-extraction path otherwise never re-reads supplement files, so a paper
    whose supplement tables were downloaded but not folded (stale binding, thin
    harvest-time fold, side-directory recovery) silently loses those variants on
    replay. The fold is non-destructive (``.pre_fold_bak``), sentinel-delimited,
    and idempotent. Returns the set of PMIDs whose FULL_CONTEXT carries folded
    supplements (so discovery can avoid a stale condensed source for them).
    """
    from harvesting.supplement_fold import fold_supplements_into_full_context

    folded: set[str] = set()
    for supp_dir in sorted(harvest_dir.glob("*_supplements")):
        if not supp_dir.is_dir():
            continue
        pmid = supp_dir.name.replace("_supplements", "")
        try:
            out = fold_supplements_into_full_context(pmid, harvest_dir)
        except Exception as exc:  # noqa: BLE001
            logger.warning("supplement fold failed for %s: %s", pmid, exc)
            continue
        if out is not None:
            folded.add(pmid)
    if folded:
        logger.info(
            "folded on-disk supplements into FULL_CONTEXT for %d PMID(s)", len(folded)
        )
    return folded


def _has_folded_supplements(path: Path) -> bool:
    """True if *path*'s text contains the supplement-fold sentinel."""
    from harvesting.supplement_fold import FOLD_BEGIN

    try:
        return FOLD_BEGIN in path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return False


def discover_source_files(harvest_dir: Path) -> dict[str, Path]:
    """Return PMID -> preferred extraction source, matching pipeline priority."""
    data_zones = {
        f.name.replace("_DATA_ZONES.md", ""): f
        for f in harvest_dir.glob("*_DATA_ZONES.md")
        if is_usable_fulltext_source(f)
    }
    cleaned = {
        f.name.replace("_CLEANED.md", ""): f
        for f in harvest_dir.glob("*_CLEANED.md")
        if is_usable_fulltext_source(f)
    }
    full_context = {
        f.name.replace("_FULL_CONTEXT.md", ""): f
        for f in harvest_dir.glob("*_FULL_CONTEXT.md")
        if is_usable_fulltext_source(f)
    }

    sources: dict[str, Path] = {}
    for pmid in sorted(set(data_zones) | set(cleaned) | set(full_context)):
        if pmid in data_zones:
            chosen = data_zones[pmid]
        elif pmid in cleaned:
            chosen = cleaned[pmid]
        else:
            chosen = full_context[pmid]
        # If supplements were folded into FULL_CONTEXT but the preferred condensed
        # source predates the fold (lacks the sentinel), that condensed form is
        # stale and would drop the supplement tables — fall back to the grown
        # FULL_CONTEXT for just this PMID.
        fc = full_context.get(pmid)
        if (
            fc is not None
            and chosen is not fc
            and _has_folded_supplements(fc)
            and not _has_folded_supplements(chosen)
        ):
            chosen = fc
        sources[pmid] = chosen
    return sources


# Fraction of the prior extraction's distinct codon positions a re-extraction
# must still cover to be accepted when it has FEWER total rows. This is the
# no-gold safety net: it rejects a re-extraction that pairs a few rows while
# silently dropping real variants, on genes that have no manually curated answer
# (the cardiac-gene --land path catches this with a gold rescore; new genes do
# not). Intrinsic signal only — no gene/PMID/gold reference.
_POSITION_COVERAGE_TOL = float(os.getenv("GVF_REPLAY_POSITION_COVERAGE_TOL", "0.85"))


def _dedup_variants(variants: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Variants deduplicated by (cdna, protein, genomic) identity, first kept."""
    seen: set[tuple[str, str, str]] = set()
    out: list[dict[str, Any]] = []
    for v in variants:
        if not (
            v.get("cdna_notation")
            or v.get("protein_notation")
            or v.get("genomic_position")
        ):
            continue
        key = (
            str(v.get("cdna_notation") or "").lower(),
            str(v.get("protein_notation") or "").lower(),
            str(v.get("genomic_position") or "").lower(),
        )
        if key in seen:
            continue
        seen.add(key)
        out.append(v)
    return out


def _dedup_count(variants: list[dict[str, Any]]) -> int:
    return len(_dedup_variants(variants))


def _variant_positions(variants: list) -> set[int]:
    """Distinct codon positions across a variant list, unifying cDNA + protein.

    A protein residue number and a cDNA implied codon ((nt + 2) // 3) for the
    same variant collapse to one position, so a cDNA-only row and a protein-only
    row for one variant are counted once. Gene- and gold-agnostic; used to detect
    recall loss (dropped positions) and quality lift (re-pairing) without a gold
    standard, so it generalizes to genes with no curated answer.
    """
    positions: set[int] = set()
    for variant in variants or []:
        if not isinstance(variant, dict):
            continue
        protein = str(variant.get("protein_notation") or "")
        pm = re.search(r"\d{1,5}", protein)
        if pm:
            positions.add(int(pm.group(0)))
        cdna = str(variant.get("cdna_notation") or "")
        cm = re.search(r"\d+", cdna)
        if cm:
            positions.add((int(cm.group(0)) + 2) // 3)
    return positions


def _position_coverage(new_variants: list, prior_variants: list) -> float:
    """Fraction of the prior extraction's distinct positions still present.

    1.0 when the prior had no positional variants (nothing to lose).
    """
    return _coverage_of_positions(new_variants, _variant_positions(prior_variants))


def _coverage_of_positions(variants: list, target_positions: set[int]) -> float:
    """Fraction of ``target_positions`` present in ``variants`` (1.0 if empty).

    Used to score a re-extraction against the deterministic TABLE parse of its
    own source — the trusted structural baseline. Covering the table parse means
    the re-extraction is faithful (an over-counted prior can be dropped); falling
    below it means the re-extraction under-extracted real table rows. No gold
    standard required, so it generalizes to genes with no curated answer.
    """
    if not target_positions:
        return 1.0
    return len(_variant_positions(variants) & target_positions) / len(target_positions)


def _paired_count(variants: list) -> int:
    """Distinct variants carrying BOTH a cDNA and a protein notation."""
    n = 0
    for v in _dedup_variants([x for x in (variants or []) if isinstance(x, dict)]):
        if (v.get("cdna_notation") or "").strip() and (
            v.get("protein_notation") or ""
        ).strip():
            n += 1
    return n


def deterministic_variant_list(
    extractor: ExpertExtractor,
    source_file: Path,
    gene: str,
) -> list[dict[str, Any]]:
    """Deduplicated deterministic-parser variants for a source (cheap, no LLM).

    Proxy for what a re-extraction would yield; used by candidate selection to
    detect both a count lift and a notation-pairing (quality) lift.
    """
    if not is_usable_fulltext_source(source_file):
        return []
    text = source_file.read_text(encoding="utf-8", errors="replace")
    text = extractor._augment_pdf_linearized_tables(text)
    variants: list[dict[str, Any]] = []
    variants.extend(extractor._parse_markdown_table_variants(text, gene))
    fixed_width_variants = extractor._parse_fixed_width_table_variants(text, gene)
    if fixed_width_variants:
        pmid_match = re.search(r"\d{6,}", source_file.name)
        logger.info(
            "deterministic_fixed_width_parser_candidate pmid=%s gene=%s source=%s variants=%d",
            pmid_match.group(0) if pmid_match else "unknown",
            gene,
            source_file,
            len(fixed_width_variants),
        )
    variants.extend(fixed_width_variants)
    variants.extend(extractor._parse_vertical_gene_table_variants(text, gene))
    return _dedup_variants(variants)


def deterministic_variant_count(
    extractor: ExpertExtractor,
    source_file: Path,
    gene: str,
) -> int:
    return len(deterministic_variant_list(extractor, source_file, gene))


def select_replay_candidates(
    *,
    gene: str,
    harvest_dir: Path,
    extraction_dir: Path,
    min_deterministic_variants: int,
    min_deterministic_lift: int,
    deterministic_lift_ratio: float,
    include_source_newer: bool,
    replay_missing_fingerprint: bool,
    replay_unbound_source: bool,
    force_pmids: set[str] | None = None,
    source_overrides: dict[str, Path] | None = None,
) -> list[ReplayCandidate]:
    sources = discover_source_files(harvest_dir)
    for pmid, source_file in (source_overrides or {}).items():
        if is_usable_fulltext_source(source_file):
            sources[str(pmid)] = source_file
    extractor = ExpertExtractor(models=["noop"], tier_threshold=0)
    candidates: list[ReplayCandidate] = []

    for pmid, source_file in sources.items():
        output_file = extraction_dir / f"{gene}_PMID_{pmid}.json"
        data = _json_load(output_file) if output_file.exists() else {}
        metadata = data.get("extraction_metadata") or {}
        current_count = _variant_count(data)
        det_variants = deterministic_variant_list(extractor, source_file, gene)
        deterministic_count = len(det_variants)

        reasons: list[str] = []
        if not output_file.exists():
            reasons.append("missing_extraction")
        if force_pmids and pmid in force_pmids:
            reasons.append("forced_pmid")
        if metadata and _metadata_mentions_abstract_only(metadata):
            reasons.append("stale_abstract_only")
        if metadata and replay_unbound_source:
            if _metadata_source_is_unbound(
                metadata,
                output_file=output_file,
                extraction_dir=extraction_dir,
            ):
                reasons.append("unbound_source_metadata")

        existing_sha = metadata.get("source_sha256")
        if existing_sha and existing_sha != _sha256(source_file):
            reasons.append("source_fingerprint_mismatch")
        elif replay_missing_fingerprint and output_file.exists() and not existing_sha:
            reasons.append("missing_source_fingerprint")

        if include_source_newer and output_file.exists():
            if source_file.stat().st_mtime > output_file.stat().st_mtime + 1:
                reasons.append("source_newer_than_extraction")

        deterministic_lift = deterministic_count - current_count
        ratio_ok = (
            deterministic_count >= current_count * deterministic_lift_ratio
            if current_count
            else deterministic_count >= min_deterministic_variants
        )
        absolute_lift_ok = deterministic_lift >= max(
            min_deterministic_lift, DETERMINISTIC_ABSOLUTE_LIFT_OVERRIDE
        )
        if (
            deterministic_count >= min_deterministic_variants
            and deterministic_lift >= min_deterministic_lift
            and (ratio_ok or absolute_lift_ok)
        ):
            reason = (
                "deterministic_parser_lift"
                if ratio_ok
                else "deterministic_parser_absolute_lift"
            )
            reasons.append(f"{reason}:{deterministic_count}>{current_count}")

        # Quality lift: the gene-scoped deterministic table parse finds MORE paired
        # (cDNA+protein) variants than the current extraction has -> re-extraction
        # would improve notation completeness (a clean paired table replacing a
        # stale over-counted cDNA-only extraction the count-only selector skips).
        # We deliberately do NOT require the deterministic parse to cover the
        # current extraction's positions: an over-counted current is exactly the
        # case to re-extract, and the gate (not the selector) decides keep/revert.
        # Pairing, not raw count or any gene/gold reference, so it generalizes.
        current_variants_list = data.get("variants") or []
        det_paired = _paired_count(det_variants)
        current_paired = _paired_count(current_variants_list)
        if (
            deterministic_count >= min_deterministic_variants
            and det_paired > current_paired
        ):
            reasons.append(
                f"deterministic_quality_lift:{det_paired}>{current_paired}_paired"
            )

        if reasons:
            candidates.append(
                ReplayCandidate(
                    pmid=pmid,
                    source_file=source_file,
                    output_file=output_file,
                    current_variants=current_count,
                    deterministic_variants=deterministic_count,
                    deterministic_positions=frozenset(_variant_positions(det_variants)),
                    reasons=reasons,
                )
            )
    return candidates


def write_candidates_csv(candidates: list[ReplayCandidate], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "pmid",
                "source_file",
                "output_file",
                "current_variants",
                "deterministic_variants",
                "reasons",
            ],
        )
        writer.writeheader()
        for c in candidates:
            writer.writerow(
                {
                    "pmid": c.pmid,
                    "source_file": str(c.source_file),
                    "output_file": str(c.output_file),
                    "current_variants": c.current_variants,
                    "deterministic_variants": c.deterministic_variants,
                    "reasons": ";".join(c.reasons),
                }
            )


def _is_variant_explosion(
    prior_variant_count: Optional[int],
    new_variant_count: int,
    *,
    ratio: float,
    min_new: int,
    min_delta: int,
) -> bool:
    """True when a replay's variant count blows up suspiciously vs the prior JSON.

    The regression gate only catches ``new < prior``. This catches the opposite
    failure mode: re-binding to a larger but garbage / wrong-paper / multi-gene
    source (or OCR noise) that explodes the row count. We require ALL of a large
    multiple of the prior count, a large absolute count, and a large absolute
    delta — the signature of garbage/leakage rather than a legitimate supplement
    recovery (e.g. 5 -> 28 stays well under the absolute floor). Pure
    internal-consistency check; no gold standard. A recovery from a 0/None prior
    is never gated (nothing to compare against).
    """
    if prior_variant_count is None or prior_variant_count <= 0:
        return False
    return (
        new_variant_count > prior_variant_count * ratio
        and new_variant_count >= min_new
        and (new_variant_count - prior_variant_count) >= min_delta
    )


def replay_candidates(
    *,
    candidates: list[ReplayCandidate],
    gene: str,
    harvest_dir: Path,
    backup_dir: Path,
    tier_threshold: int,
    dry_run: bool,
    gate_regressions: bool = True,
    gate_explosions: bool = True,
    explosion_ratio: float = 10.0,
    explosion_min_new: int = 400,
    explosion_min_delta: int = 300,
    replay_models: list[str] | None = None,
) -> dict[str, Any]:
    if dry_run:
        return {
            "attempted": 0,
            "successful": 0,
            "failed": 0,
            "gated": 0,
            "attempted_pmids": [],
            "successful_pmids": [],
            "failed_pmids": [],
            "gated_pmids": [],
            "gated_regressions": [],
            "gated_explosions": [],
            "errors": [],
            "replay_models": replay_models or [],
        }

    backup_dir.mkdir(parents=True, exist_ok=True)
    extractor = ExpertExtractor(
        models=replay_models,
        tier_threshold=tier_threshold,
        fulltext_dir=str(harvest_dir),
    )
    attempted = 0
    successful = 0
    gated = 0
    successful_pmids: list[str] = []
    failed_pmids: list[str] = []
    gated_regressions: list[dict[str, Any]] = []
    gated_explosions: list[dict[str, Any]] = []
    errors: list[dict[str, str]] = []

    for candidate in candidates:
        attempted += 1
        backup_file = backup_dir / candidate.output_file.name
        if candidate.output_file.exists() and not backup_file.exists():
            shutil.copy2(candidate.output_file, backup_file)

        try:
            if not is_usable_fulltext_source(candidate.source_file):
                raise RuntimeError("source is abstract-only fallback or empty")

            text = candidate.source_file.read_text(encoding="utf-8", errors="replace")
            result = extractor.extract(
                Paper(pmid=candidate.pmid, full_text=text, gene_symbol=gene)
            )
            if not result.success or not result.extracted_data:
                raise RuntimeError(result.error or "empty extraction result")
            new_variants = result.extracted_data.get("variants", []) or []
            new_variant_count = len(new_variants)
            prior_variant_count = _backup_variant_count(backup_file)
            # Quality-aware regression gate: a re-extraction with fewer total
            # rows is NOT a regression when it carries more gold-matchable
            # (paired cDNA+protein) content. This stops a stale over-counted
            # cDNA-only extraction from blocking a cleaner paired one.
            # A re-extraction with FEWER rows than the prior is accepted only if
            # (1) it does not lose paired (cDNA+protein) variants, AND (2) it stays
            # faithful to the gene-scoped deterministic TABLE parse of this source
            # (covers >=85% of its positions). (2) is the no-gold discriminator
            # between replacing an OVER-counted prior (prior >> table; the table is
            # the trusted structure, so the excess is droppable) and a lossy
            # UNDER-extraction (new misses table rows the deterministic parser
            # found). A leaky multi-gene deterministic parse only makes this
            # conservatively keep the prior — it never corrupts.
            new_paired = _paired_count(new_variants)
            prior_paired = _paired_count(_backup_variants(backup_file))
            if candidate.deterministic_positions:
                # Trusted table baseline exists: the re-extraction must be faithful
                # to it (>=85% of its positions) and not lose paired variants.
                quality_held = (
                    new_paired >= prior_paired
                    and _coverage_of_positions(
                        new_variants, candidate.deterministic_positions
                    )
                    >= _POSITION_COVERAGE_TOL
                )
            else:
                # No deterministic table (prose-only source): fall back to the prior
                # extraction — require a STRICT pairing improvement that still keeps
                # the prior's positions, so a plain count drop is not waved through.
                quality_held = (
                    new_paired > prior_paired
                    and _position_coverage(new_variants, _backup_variants(backup_file))
                    >= _POSITION_COVERAGE_TOL
                )
            if (
                gate_regressions
                and prior_variant_count is not None
                and new_variant_count < prior_variant_count
                and not quality_held
            ):
                gated += 1
                failed_pmids.append(candidate.pmid)
                gated_regressions.append(
                    {
                        "pmid": candidate.pmid,
                        "prior_variant_count": prior_variant_count,
                        "new_variant_count": new_variant_count,
                        "delta": new_variant_count - prior_variant_count,
                    }
                )
                logger.warning(
                    "PMID %s replay gated: prior=%d new=%d Δ=%d (kept backup)",
                    candidate.pmid,
                    prior_variant_count,
                    new_variant_count,
                    new_variant_count - prior_variant_count,
                )
                if backup_file.exists():
                    shutil.copy2(backup_file, candidate.output_file)
                continue
            if gate_explosions and _is_variant_explosion(
                prior_variant_count,
                new_variant_count,
                ratio=explosion_ratio,
                min_new=explosion_min_new,
                min_delta=explosion_min_delta,
            ):
                gated += 1
                failed_pmids.append(candidate.pmid)
                gated_explosions.append(
                    {
                        "pmid": candidate.pmid,
                        "prior_variant_count": prior_variant_count,
                        "new_variant_count": new_variant_count,
                        "delta": new_variant_count - prior_variant_count,
                    }
                )
                logger.warning(
                    "PMID %s replay gated (variant explosion): prior=%s new=%d "
                    "Δ=+%d (kept backup; re-run with --no-gate-explosions or a "
                    "higher --explosion-min-new to accept)",
                    candidate.pmid,
                    prior_variant_count,
                    new_variant_count,
                    new_variant_count - prior_variant_count,
                )
                if backup_file.exists():
                    shutil.copy2(backup_file, candidate.output_file)
                continue
            metadata = result.extracted_data.setdefault("extraction_metadata", {})
            metadata.update(_source_metadata(candidate.source_file))
            metadata["model_used"] = result.model_used
            candidate.output_file.write_text(
                json.dumps(result.extracted_data, indent=2),
                encoding="utf-8",
            )
            successful += 1
            successful_pmids.append(candidate.pmid)
            logger.info(
                "replayed PMID %s: %s variants via %s",
                candidate.pmid,
                new_variant_count,
                result.model_used,
            )
        except Exception as exc:  # noqa: BLE001
            errors.append({"pmid": candidate.pmid, "error": str(exc)})
            failed_pmids.append(candidate.pmid)
            logger.warning("PMID %s replay failed: %s", candidate.pmid, exc)
            if backup_file.exists():
                shutil.copy2(backup_file, candidate.output_file)

    return {
        "attempted": attempted,
        "successful": successful,
        "failed": len(errors),
        "gated": gated,
        "attempted_pmids": [candidate.pmid for candidate in candidates],
        "successful_pmids": successful_pmids,
        "failed_pmids": sorted(set(failed_pmids)),
        "gated_pmids": sorted(
            {
                row["pmid"]
                for row in [*gated_regressions, *gated_explosions]
                if row.get("pmid")
            }
        ),
        "gated_regressions": gated_regressions,
        "gated_explosions": gated_explosions,
        "errors": errors,
        "backup_dir": str(backup_dir),
        "replay_models": replay_models or [],
    }


def _backup_variant_count(backup_file: Path) -> Optional[int]:
    """Return the number of variants in the backup extraction JSON, or None.

    Used by `replay_candidates` to gate against per-PMID regressions: if the
    re-extraction produces fewer variants than the prior JSON (preserved as
    backup at the start of the replay), the backup is restored. Pure
    internal-consistency check; no gold standard required.
    """
    if not backup_file.exists():
        return None
    try:
        data = json.loads(backup_file.read_text(encoding="utf-8"))
    except Exception:
        return None
    if not isinstance(data, dict):
        return None
    return _variant_count(data)


def _backup_variants(backup_file: Path) -> list[dict[str, Any]]:
    """Variant list from the backup extraction JSON (empty if absent/unreadable).

    Used for the gate's position-coverage check against the prior extraction.
    """
    if not backup_file.exists():
        return []
    try:
        data = json.loads(backup_file.read_text(encoding="utf-8"))
    except Exception:
        return []
    if not isinstance(data, dict):
        return []
    variants = data.get("variants")
    return variants if isinstance(variants, list) else []


def rebuild_db(
    extraction_dir: Path, output_db: Path, *, dry_run: bool
) -> dict[str, Any]:
    if dry_run:
        return {"output_db": str(output_db), "skipped": True}
    if output_db.exists():
        raise FileExistsError(f"Output DB already exists: {output_db}")
    conn = create_database_schema(str(output_db))
    try:
        stats = migrate_extraction_directory(conn, extraction_dir)
    finally:
        conn.close()
    stats["output_db"] = str(output_db)
    return stats


def run_recovery_layers(
    *,
    gene: str,
    db: Path,
    run_dir: Path,
    gold: Optional[Path],
    outdir: Path,
    skip_layers: list[str],
    dry_run: bool,
) -> Optional[dict[str, Any]]:
    if dry_run:
        return {"skipped": True, "outdir": str(outdir)}
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "recall_recovery" / "run_all_layers.py"),
        "--gene",
        gene,
        "--db",
        str(db),
        "--pmc-dir",
        str(run_dir / "pmc_fulltext"),
        "--outdir",
        str(outdir),
    ]
    if gold:
        cmd.extend(["--gold", str(gold)])
    for layer in skip_layers:
        cmd.extend(["--skip", layer])
    logger.info("recovery → %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(result.stderr[-1000:] or result.stdout[-1000:])
    summary_path = outdir / "progression.json"
    if summary_path.exists():
        return json.loads(summary_path.read_text(encoding="utf-8"))
    return {"outdir": str(outdir), "stdout_tail": result.stdout[-1000:]}


def _default_gold(gene: str) -> Optional[Path]:
    path = (
        REPO_ROOT
        / "gene_variant_fetcher_gold_standard"
        / "normalized"
        / f"{gene}_recall_input.csv"
    )
    return path if path.exists() else None


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gene", required=True, help="Gene symbol")
    p.add_argument("--run-dir", required=True, type=Path, help="Existing GVF run dir")
    p.add_argument("--harvest-dir", type=Path, default=None)
    p.add_argument("--extraction-dir", type=Path, default=None)
    p.add_argument(
        "--no-supplement-fold",
        action="store_true",
        help=(
            "Skip folding on-disk {pmid}_supplements/ into FULL_CONTEXT before "
            "discovery. By default the fold runs (non-destructive, idempotent) so "
            "downloaded supplement tables are visible to re-extraction."
        ),
    )
    p.add_argument(
        "--stage-extractions",
        action="store_true",
        help=(
            "Copy the extraction directory into the refresh directory and replay "
            "against that copy. This is safer for experiments because the active "
            "run's extraction JSONs are left untouched."
        ),
    )
    p.add_argument("--output-db", type=Path, default=None)
    p.add_argument(
        "--replace-db",
        action="store_true",
        help="After refresh/recovery, back up <run-dir>/<GENE>.db and replace it.",
    )
    p.add_argument("--gold", type=Path, default=None)
    p.add_argument(
        "--layers-outdir",
        type=Path,
        default=None,
        help=(
            "Directory for recovery-layer outputs. Default: "
            "<refresh_dir>/layers, so each refresh keeps independent metrics."
        ),
    )
    p.add_argument("--skip-recovery", action="store_true")
    p.add_argument(
        "--skip-layer",
        action="append",
        default=[],
        choices=("clinvar", "pubtator", "figures"),
        help="Recovery layer to skip. Repeatable.",
    )
    p.add_argument("--tier-threshold", type=int, default=1)
    p.add_argument(
        "--replay-model",
        action="append",
        default=[],
        help=(
            "Override Tier 3 model(s) for replay extraction only. Repeatable "
            "or comma-separated. Candidate selection still uses deterministic "
            "parsers only."
        ),
    )
    p.add_argument("--min-deterministic-variants", type=int, default=20)
    p.add_argument("--min-deterministic-lift", type=int, default=5)
    p.add_argument("--deterministic-lift-ratio", type=float, default=1.2)
    p.add_argument("--include-source-newer", action="store_true")
    p.add_argument("--replay-missing-fingerprint", action="store_true")
    p.add_argument(
        "--replay-unbound-source",
        action="store_true",
        help=(
            "Replay extractions whose metadata source_file points at the JSON "
            "artifact/extractions directory instead of a usable source document."
        ),
    )
    p.add_argument(
        "--pmids-file",
        action="append",
        type=Path,
        default=[],
        help=(
            "Force replay for PMIDs listed one per line when a usable source exists. "
            "Repeatable."
        ),
    )
    p.add_argument(
        "--candidate-report",
        type=Path,
        default=None,
        help=(
            "Force replay for PMIDs selected from a paper_disagreement_report.csv. "
            "This is evaluation-aided and should not be used for cold-start claims."
        ),
    )
    p.add_argument(
        "--report-class",
        action="append",
        default=[],
        help=(
            "Failure class to include from --candidate-report. Repeatable. "
            "Default: include all classes."
        ),
    )
    p.add_argument(
        "--report-min-missing-rows",
        type=int,
        default=1,
        help="Minimum missing_rows required for --candidate-report PMIDs.",
    )
    p.add_argument(
        "--report-max-row-recall",
        type=float,
        default=None,
        help=(
            "Optional maximum existing row_recall for --candidate-report PMIDs. "
            "Use this to avoid replaying already high-recall papers."
        ),
    )
    p.add_argument(
        "--use-report-available-context",
        action="store_true",
        help=(
            "When --candidate-report is used, replay each selected PMID from "
            "available_context_path instead of only using the run directory's "
            "normal source priority."
        ),
    )
    p.add_argument(
        "--report-context-search-root",
        action="append",
        type=Path,
        default=[],
        help=(
            "When --use-report-available-context is set, also search this root "
            "for larger <PMID>_FULL_CONTEXT.md files and use the largest usable "
            "context. Repeatable."
        ),
    )
    p.add_argument(
        "--source-override-csv",
        action="append",
        type=Path,
        default=[],
        help=(
            "CSV with PMID and source path columns to force replay from explicit "
            "source files. If an action column is present, only refresh_replay "
            "rows are used. Repeatable."
        ),
    )
    p.add_argument(
        "--only-forced-pmids",
        action="store_true",
        help=(
            "After candidate discovery, keep only PMIDs explicitly supplied by "
            "--pmids-file, --candidate-report, or --source-override-csv."
        ),
    )
    p.add_argument(
        "--no-gate-regressions",
        action="store_true",
        help=(
            "Disable per-PMID acceptance gating. By default, when a replay "
            "produces fewer variants than the prior extraction JSON for the "
            "same PMID, the backup is restored and the new extraction is "
            "discarded. This flag overrides that protection and overwrites "
            "the prior JSON unconditionally; use only when you know the new "
            "extraction is authoritative even if it has fewer variants."
        ),
    )
    p.add_argument(
        "--no-gate-explosions",
        action="store_true",
        help=(
            "Disable the variant-explosion gate. By default, when a replay "
            "produces a suspiciously larger variant count than the prior JSON "
            "(a large multiple AND a large absolute count AND a large delta), "
            "the backup is restored and the explosion is recorded for audit. "
            "This guards against re-binding to a garbage / wrong-paper / "
            "multi-gene source. Use this flag to accept the larger extraction."
        ),
    )
    p.add_argument(
        "--explosion-ratio",
        type=float,
        default=10.0,
        help="Variant-explosion gate: trip when new > prior * this. Default 10.0.",
    )
    p.add_argument(
        "--explosion-min-new",
        type=int,
        default=400,
        help=(
            "Variant-explosion gate: new count must reach this absolute floor "
            "to trip. Default 400 (keeps legitimate large recoveries)."
        ),
    )
    p.add_argument(
        "--explosion-min-delta",
        type=int,
        default=300,
        help="Variant-explosion gate: new - prior must reach this to trip. Default 300.",
    )
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--verbose", "-v", action="store_true")
    return p


def main() -> int:
    args = build_parser().parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    gene = args.gene.upper()
    run_dir = args.run_dir.expanduser().resolve()
    harvest_dir = (args.harvest_dir or run_dir / "pmc_fulltext").expanduser()
    extraction_dir = (args.extraction_dir or run_dir / "extractions").expanduser()
    gold = args.gold.expanduser() if args.gold else _default_gold(gene)

    if not run_dir.is_dir():
        sys.exit(f"Run dir not found: {run_dir}")
    if not harvest_dir.is_dir():
        sys.exit(f"Harvest dir not found: {harvest_dir}")
    if not extraction_dir.is_dir():
        sys.exit(f"Extraction dir not found: {extraction_dir}")

    # Fold any on-disk supplements into FULL_CONTEXT so re-extraction sees them
    # (the discovery glob never re-reads {pmid}_supplements/ on its own).
    if not args.no_supplement_fold:
        fold_on_disk_supplements(harvest_dir)
    if args.stage_extractions and args.replace_db:
        sys.exit("--stage-extractions cannot be combined with --replace-db")

    force_pmids: set[str] = set()
    source_overrides: dict[str, Path] = {}
    for pmids_file in args.pmids_file:
        force_pmids.update(_read_pmid_file(pmids_file.expanduser()))
    if args.candidate_report:
        report = args.candidate_report.expanduser()
        report_classes = {str(item) for item in args.report_class}
        force_pmids.update(
            load_report_pmids(
                report=report,
                gene=gene,
                failure_classes=report_classes,
                min_missing_rows=args.report_min_missing_rows,
                max_row_recall=args.report_max_row_recall,
            )
        )
        if args.use_report_available_context:
            source_overrides.update(
                load_report_available_contexts(
                    report=report,
                    gene=gene,
                    failure_classes=report_classes,
                    min_missing_rows=args.report_min_missing_rows,
                    max_row_recall=args.report_max_row_recall,
                    context_search_roots=[
                        path.expanduser() for path in args.report_context_search_root
                    ],
                )
            )
    for override_csv in args.source_override_csv:
        overrides_from_csv = load_source_override_csv(override_csv.expanduser())
        source_overrides.update(overrides_from_csv)
        force_pmids.update(overrides_from_csv)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    refresh_dir = run_dir / f"refresh_{timestamp}"
    refresh_dir.mkdir(parents=True, exist_ok=True)
    original_extraction_dir = extraction_dir
    if args.stage_extractions:
        staged_extraction_dir = refresh_dir / "staged_extractions"
        shutil.copytree(extraction_dir, staged_extraction_dir)
        extraction_dir = staged_extraction_dir
        logger.info(
            "Staged extraction directory %s → %s",
            original_extraction_dir,
            extraction_dir,
        )

    candidates = select_replay_candidates(
        gene=gene,
        harvest_dir=harvest_dir,
        extraction_dir=extraction_dir,
        min_deterministic_variants=args.min_deterministic_variants,
        min_deterministic_lift=args.min_deterministic_lift,
        deterministic_lift_ratio=args.deterministic_lift_ratio,
        include_source_newer=args.include_source_newer,
        replay_missing_fingerprint=args.replay_missing_fingerprint,
        replay_unbound_source=args.replay_unbound_source,
        force_pmids=force_pmids or None,
        source_overrides=source_overrides or None,
    )
    if args.only_forced_pmids:
        if not force_pmids:
            sys.exit(
                "--only-forced-pmids requires --pmids-file, --candidate-report, "
                "or --source-override-csv"
            )
        candidates = [c for c in candidates if c.pmid in force_pmids]
    candidates_csv = refresh_dir / "replay_candidates.csv"
    write_candidates_csv(candidates, candidates_csv)
    logger.info("Selected %d replay candidates → %s", len(candidates), candidates_csv)
    replay_models = _split_model_args(args.replay_model)

    replay_stats = replay_candidates(
        candidates=candidates,
        gene=gene,
        harvest_dir=harvest_dir,
        backup_dir=refresh_dir / "extraction_json_backup",
        tier_threshold=args.tier_threshold,
        dry_run=args.dry_run,
        gate_regressions=not args.no_gate_regressions,
        gate_explosions=not args.no_gate_explosions,
        explosion_ratio=args.explosion_ratio,
        explosion_min_new=args.explosion_min_new,
        explosion_min_delta=args.explosion_min_delta,
        replay_models=replay_models or None,
    )

    output_db = args.output_db
    if output_db is None:
        output_db = run_dir / f"{gene}.refresh_{timestamp}.db"
    output_db = output_db.expanduser()
    db_stats = rebuild_db(extraction_dir, output_db, dry_run=args.dry_run)

    recovery_summary = None
    if not args.skip_recovery:
        layers_outdir = (args.layers_outdir or refresh_dir / "layers").expanduser()
        recovery_summary = run_recovery_layers(
            gene=gene,
            db=output_db,
            run_dir=run_dir,
            gold=gold,
            outdir=layers_outdir,
            skip_layers=args.skip_layer,
            dry_run=args.dry_run,
        )

    replace_info = None
    active_db = output_db
    current_db = run_dir / f"{gene}.db"
    if args.replace_db and not args.dry_run:
        backup_db = run_dir / f"{gene}.db.before_refresh_{timestamp}.db"
        if current_db.exists():
            shutil.copy2(current_db, backup_db)
        shutil.copy2(output_db, current_db)
        active_db = current_db
        replace_info = {"backup_db": str(backup_db), "active_db": str(current_db)}

    summary = {
        "gene": gene,
        "run_dir": str(run_dir),
        "harvest_dir": str(harvest_dir),
        "extraction_dir": str(extraction_dir),
        "original_extraction_dir": str(original_extraction_dir),
        "staged_extractions": bool(args.stage_extractions),
        "gold": str(gold) if gold else None,
        "candidate_count": len(candidates),
        "forced_pmid_count": len(force_pmids),
        "replay_models": replay_models,
        "candidates_csv": str(candidates_csv),
        "replay": replay_stats,
        "db_rebuild": db_stats,
        "recovery": recovery_summary,
        "replace": replace_info,
        "active_db": str(active_db),
    }
    summary_path = refresh_dir / "refresh_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
