#!/usr/bin/env python3
"""Audit fixed collaborator extraction runs without pretending QC is precision.

The audit verifies pinned cohort membership, workflow/trace/VF gates, source
coverage, identity hygiene, count invariants, trust quarantine, and reference
enrichment.  It deliberately leaves publication on HOLD: structural checks do
not replace a source-grounded precision sample or a trusted-field importer.
"""

from __future__ import annotations

import argparse
import json
import re
import sqlite3
from collections import Counter
from pathlib import Path
from typing import Any

from pipeline.filters import names_nonhuman_ortholog
from utils.gene_metadata import gene_alias_regex

# A VariantFeatures miss can still be a defensible paper identity: truly novel
# protein alleles and transcript-bound cDNA variants are expected not to exist
# in the reference feature table.  Every other live miss is held for manual
# adjudication; high-confidence wrong-residue/out-of-range calls must live only
# in ``quarantined_variants``.
_ALLOWED_UNMATCHED_VF_CLASSES = frozenset(
    {
        "novel_in_range",
        "cdna_only_unmatched",
        "known_isoform_offset",
        "legacy_source_notation",
    }
)


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _read_pmids(path: Path) -> list[str]:
    return [
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def _parse_mapping(values: list[str], label: str) -> dict[str, Path]:
    parsed: dict[str, Path] = {}
    for raw in values:
        if "=" not in raw:
            raise SystemExit(f"{label} must use GENE=PATH: {raw}")
        gene, path = raw.split("=", 1)
        parsed[gene.strip().upper()] = Path(path).expanduser().resolve()
    return parsed


def _scalar(con: sqlite3.Connection, sql: str, params: tuple[Any, ...] = ()) -> int:
    row = con.execute(sql, params).fetchone()
    return int(row[0] or 0)


def _group_counts(con: sqlite3.Connection, sql: str) -> dict[str, int]:
    return {str(key or "unknown"): int(value) for key, value in con.execute(sql)}


def _extraction_json_issues(gene: str, run_dir: Path, expected: list[str]) -> list[str]:
    """Validate that each claimed extraction is a real, paper-bound payload."""
    issues: list[str] = []
    for pmid in expected:
        path = run_dir / "extractions" / f"{gene}_PMID_{pmid}.json"
        if not path.is_file():
            continue
        try:
            payload = _read_json(path)
        except (OSError, json.JSONDecodeError) as exc:
            issues.append(f"PMID {pmid}: unreadable extraction JSON ({exc})")
            continue
        if not isinstance(payload, dict) or not payload:
            issues.append(f"PMID {pmid}: empty extraction payload")
            continue
        paper_metadata = payload.get("paper_metadata")
        if not isinstance(paper_metadata, dict) and payload.get("pmid") is not None:
            # Deterministic recovery artifacts use the older flat paper schema.
            paper_metadata = {
                "pmid": payload.get("pmid"),
                "title": payload.get("title"),
            }
        extraction_metadata = payload.get("extraction_metadata")
        variants = payload.get("variants")
        if not isinstance(paper_metadata, dict):
            issues.append(f"PMID {pmid}: paper_metadata is missing")
        elif str(paper_metadata.get("pmid") or "") != pmid:
            issues.append(f"PMID {pmid}: paper_metadata PMID does not match filename")
        if not isinstance(extraction_metadata, dict):
            issues.append(f"PMID {pmid}: extraction_metadata is missing")
        if not isinstance(variants, list):
            issues.append(f"PMID {pmid}: variants must be a list")
            continue
        for index, variant in enumerate(variants):
            if not isinstance(variant, dict):
                continue
            headers = variant.get("source_table_headers")
            quote = str(variant.get("evidence_quote") or "").strip()
            if not isinstance(headers, list) or not quote.startswith("|"):
                continue
            gene_columns = [
                column
                for column, header in enumerate(headers)
                if re.sub(r"[^a-z0-9]+", " ", str(header).lower()).strip()
                in {"gene", "gene symbol", "gene name"}
            ]
            cells = [cell.strip() for cell in quote.strip("|").split("|")]
            for column in gene_columns:
                if column >= len(cells):
                    continue
                source_gene = cells[column]
                if source_gene and not gene_alias_regex(
                    gene, include_query_aliases=False
                ).search(source_gene):
                    issues.append(
                        f"PMID {pmid} variant {index}: source Gene column "
                        f"names {source_gene!r}, not {gene}"
                    )
                    break
    return issues


def _effective_trace(status: dict[str, Any], run_dir: Path) -> dict[str, Any]:
    """Return the strongest trace evidence retained in a resumed run directory.

    A validation-only ``gvf-run --resume-dir`` correctly creates a fresh empty
    trace session and replaces ``RUN_STATUS.json``.  That must not hide the
    original extraction trace in ``llm_traces/trace_manifest.json``.  Prefer
    whichever record has more captured calls/events, while keeping the current
    status when it contains the substantive trace.
    """

    status_trace = dict(status.get("llm_trace") or {})
    status_trace.setdefault("evidence_source", "RUN_STATUS.json")
    primary_manifest = run_dir / "llm_traces" / "trace_manifest.json"
    if not primary_manifest.exists():
        return status_trace
    try:
        manifest = _read_json(primary_manifest)
    except (OSError, json.JSONDecodeError):
        return status_trace
    verification = manifest.get("verification") or {}
    manifest_trace = {
        "run_id": manifest.get("run_id"),
        "trace_root": manifest.get("trace_root"),
        "llm_call_count": int(manifest.get("llm_call_count") or 0),
        "decision_event_count": int(manifest.get("decision_event_count") or 0),
        "integrity_level": verification.get("level"),
        "integrity_errors": verification.get("errors") or [],
        "missing_decision_links": [],
        "evidence_source": str(primary_manifest),
    }
    status_records = int(status_trace.get("llm_call_count") or 0) + int(
        status_trace.get("decision_event_count") or 0
    )
    manifest_records = (
        manifest_trace["llm_call_count"] + manifest_trace["decision_event_count"]
    )
    return manifest_trace if manifest_records > status_records else status_trace


def _nonhuman_gene_title_links(
    con: sqlite3.Connection, gene: str, run_dir: Path
) -> list[dict[str, Any]]:
    """Find conservatively species-scoped target-gene papers with live links.

    The shared extractor predicate recognizes strong ortholog constructions in
    either order.  It deliberately preserves papers that explicitly assay a
    human gene in a model organism.
    """

    findings: list[dict[str, Any]] = []
    for pmid, db_title, link_count in con.execute(
        """SELECT p.pmid, COALESCE(p.title, ''), COUNT(vp.variant_id)
           FROM papers p
           JOIN variant_papers vp ON vp.pmid = p.pmid
           GROUP BY p.pmid, p.title
           ORDER BY p.pmid"""
    ):
        title = str(db_title or "")
        if not title.strip():
            abstract_path = run_dir / "abstract_json" / f"{pmid}.json"
            if abstract_path.exists():
                try:
                    abstract_data = _read_json(abstract_path)
                    title = str(
                        (abstract_data.get("metadata") or {}).get("title") or ""
                    )
                except (OSError, json.JSONDecodeError):
                    title = ""
        source_head = ""
        for suffix in ("FULL_CONTEXT.md", "CLEANED.md"):
            source_path = run_dir / "pmc_fulltext" / f"{pmid}_{suffix}"
            if source_path.is_file():
                try:
                    source_head = "\n".join(source_path.read_text().splitlines()[:8])
                except OSError:
                    source_head = ""
                break
        if names_nonhuman_ortholog(title, gene) or names_nonhuman_ortholog(
            source_head, gene
        ):
            findings.append(
                {"pmid": str(pmid), "title": str(title), "variant_links": link_count}
            )
    return findings


def audit_gene(gene: str, run_dir: Path, manifest_path: Path) -> dict[str, Any]:
    expected = _read_pmids(manifest_path)
    staged_manifest_path = run_dir / f"{gene}_pmids.txt"
    staged = _read_pmids(staged_manifest_path) if staged_manifest_path.exists() else []
    extracted = sorted(
        match.group(1)
        for path in (run_dir / "extractions").glob(f"{gene}_PMID_*.json")
        if (match := re.search(r"_PMID_(\d+)\.json$", path.name))
    )
    extraction_json_issues = _extraction_json_issues(gene, run_dir, expected)

    status_path = run_dir / "RUN_STATUS.json"
    status = _read_json(status_path) if status_path.exists() else {}
    source_path = run_dir / "source_completeness.json"
    source = _read_json(source_path) if source_path.exists() else {}
    trace = _effective_trace(status, run_dir)

    db_path = run_dir / f"{gene}.db"
    con = sqlite3.connect(db_path)
    try:
        variant_columns = {row[1] for row in con.execute("PRAGMA table_info(variants)")}
        legacy_identity_sql = (
            "AND NULLIF(TRIM(COALESCE(legacy_notation, '')), '') IS NULL"
            if "legacy_notation" in variant_columns
            else ""
        )
        variants = {
            "variants": _scalar(con, "SELECT COUNT(*) FROM variants"),
            "paper_variant_links": _scalar(con, "SELECT COUNT(*) FROM variant_papers"),
            "papers_in_db": _scalar(con, "SELECT COUNT(*) FROM papers"),
            "placeholder_paper_titles": _scalar(
                con,
                """SELECT COUNT(*) FROM papers
                   WHERE NULLIF(TRIM(COALESCE(title, '')), '') IS NULL
                      OR LOWER(TRIM(title)) IN ('unknown', 'unknown title')
                      OR LOWER(TRIM(title)) = 'paper ' || LOWER(TRIM(pmid))""",
            ),
            "nameless": _scalar(
                con,
                f"""SELECT COUNT(*) FROM variants
                   WHERE NULLIF(TRIM(COALESCE(cdna_notation, '')), '') IS NULL
                     AND NULLIF(TRIM(COALESCE(protein_notation, '')), '') IS NULL
                     AND NULLIF(TRIM(COALESCE(genomic_position, '')), '') IS NULL
                     AND NULLIF(TRIM(COALESCE(structural_description, '')), '') IS NULL
                     {legacy_identity_sql}""",
            ),
            "wrong_gene_symbol": _scalar(
                con,
                f"SELECT COUNT(*) FROM variants WHERE UPPER(gene_symbol) != '{gene}'",
            ),
        }
        db_pmids = {
            str(row[0]).strip()
            for row in con.execute("SELECT pmid FROM papers")
            if str(row[0] or "").strip()
        }
        counts = {
            "rows": _scalar(con, "SELECT COUNT(*) FROM penetrance_data"),
            "carriers_supplied": _scalar(
                con,
                "SELECT COUNT(*) FROM penetrance_data WHERE total_carriers_observed IS NOT NULL",
            ),
            "affected_supplied": _scalar(
                con,
                "SELECT COUNT(*) FROM penetrance_data WHERE affected_count IS NOT NULL",
            ),
            "unaffected_supplied": _scalar(
                con,
                "SELECT COUNT(*) FROM penetrance_data WHERE unaffected_count IS NOT NULL",
            ),
            "penetrance_percentage_supplied": _scalar(
                con,
                "SELECT COUNT(*) FROM penetrance_data WHERE penetrance_percentage IS NOT NULL",
            ),
            "age_penetrance_percentage_supplied": _scalar(
                con,
                "SELECT COUNT(*) FROM age_dependent_penetrance WHERE penetrance_percentage IS NOT NULL",
            ),
            "negative_rows": _scalar(
                con,
                """SELECT COUNT(*) FROM penetrance_data
                   WHERE total_carriers_observed < 0 OR affected_count < 0
                      OR unaffected_count < 0 OR uncertain_count < 0""",
            ),
            "impossible_partitions_raw": _scalar(
                con,
                """SELECT COUNT(*) FROM penetrance_data
                   WHERE total_carriers_observed IS NOT NULL
                     AND affected_count IS NOT NULL
                     AND unaffected_count IS NOT NULL
                     AND affected_count + unaffected_count > total_carriers_observed""",
            ),
            "impossible_partitions_unquarantined": _scalar(
                con,
                """SELECT COUNT(*) FROM penetrance_data
                   WHERE total_carriers_observed IS NOT NULL
                     AND affected_count IS NOT NULL
                     AND unaffected_count IS NOT NULL
                     AND affected_count + unaffected_count > total_carriers_observed
                     AND COALESCE(trust_tier, 'trusted') != 'quarantine'""",
            ),
        }
        trust = {
            "tiers": _group_counts(
                con,
                "SELECT COALESCE(trust_tier, 'unknown'), COUNT(*) "
                "FROM penetrance_data GROUP BY COALESCE(trust_tier, 'unknown')",
            ),
            "quarantined_fields": Counter(),
        }
        for (raw,) in con.execute(
            "SELECT field_trust FROM penetrance_data WHERE field_trust IS NOT NULL"
        ):
            try:
                fields = json.loads(raw)
            except (TypeError, json.JSONDecodeError):
                continue
            for field, decision in fields.items():
                if decision != "trusted":
                    trust["quarantined_fields"][field] += 1
        trust["quarantined_fields"] = dict(trust["quarantined_fields"])

        vf = {
            "rows": _scalar(
                con,
                "SELECT COUNT(*) FROM vf_enrichment e JOIN variants v USING(variant_id)",
            ),
            "matched": _scalar(
                con,
                "SELECT COUNT(*) FROM vf_enrichment e JOIN variants v USING(variant_id) "
                "WHERE e.matched = 1",
            ),
            "residual_classes": _group_counts(
                con,
                "SELECT COALESCE(fp_class, 'none'), COUNT(*) "
                "FROM vf_enrichment e JOIN variants v USING(variant_id) "
                "GROUP BY COALESCE(fp_class, 'none')",
            ),
            "quarantined_variants": _group_counts(
                con,
                "SELECT COALESCE(fp_class, 'unknown'), COUNT(*) "
                "FROM quarantined_variants GROUP BY COALESCE(fp_class, 'unknown')",
            ),
        }
        vf_uncovered = _scalar(
            con,
            """SELECT COUNT(*) FROM variants v
               LEFT JOIN vf_enrichment e USING(variant_id)
               WHERE e.variant_id IS NULL""",
        )
        vf_unmatched = _scalar(
            con,
            """SELECT COUNT(*) FROM variants v
               LEFT JOIN vf_enrichment e USING(variant_id)
               WHERE COALESCE(e.matched, 0) != 1""",
        )
        allowed_placeholders = ",".join(
            "?" for _ in sorted(_ALLOWED_UNMATCHED_VF_CLASSES)
        )
        vf_blocked = _scalar(
            con,
            f"""SELECT COUNT(*) FROM variants v
                LEFT JOIN vf_enrichment e USING(variant_id)
                WHERE COALESCE(e.matched, 0) != 1
                  AND COALESCE(e.fp_class, '') NOT IN ({allowed_placeholders})""",
            tuple(sorted(_ALLOWED_UNMATCHED_VF_CLASSES)),
        )
        vf_gate_passed = vf_uncovered == 0 and vf_blocked == 0
        vf["uncovered_live_variants"] = vf_uncovered
        vf["unmatched_live_variants"] = vf_unmatched
        vf["blocked_live_variants"] = vf_blocked
        vf["allowed_unmatched_classes"] = sorted(_ALLOWED_UNMATCHED_VF_CLASSES)
        vf["gate_evidence"] = (
            "every live variant is VariantFeatures-matched or belongs to an "
            "explicitly allowed novel/cDNA-only class; ambiguous and rejected "
            "identities remain in quarantined_variants only"
        )
        paper_scope = {
            "nonhuman_target_title_links": _nonhuman_gene_title_links(
                con, gene, run_dir
            )
        }
    finally:
        con.close()

    failures: list[str] = []
    if not expected:
        failures.append("pinned PMID manifest is empty")
    if status.get("status") != "completed" or status.get("exit_code") != 0:
        failures.append("run did not complete with exit 0")
    if status.get("stage_failures"):
        failures.append(f"stage failures: {status['stage_failures']}")
    if staged != expected:
        failures.append("staged PMID manifest differs in membership or order")
    if set(extracted) != set(expected):
        failures.append("extraction JSON PMID set differs from pinned manifest")
    if extraction_json_issues:
        failures.append(f"invalid extraction JSON payloads: {extraction_json_issues}")
    if db_pmids != set(expected):
        failures.append("papers table PMID set differs from pinned manifest")
    if variants["variants"] == 0 or variants["paper_variant_links"] == 0:
        failures.append("run contains no publishable variant evidence")
    if vf["uncovered_live_variants"]:
        failures.append(
            "mandatory VariantFeatures enrichment/quarantine evidence missing"
        )
    if vf["blocked_live_variants"]:
        failures.append(
            "mandatory VariantFeatures enrichment/quarantine gate leaves "
            f"{vf['blocked_live_variants']} ambiguous live variant identities"
        )
    if trace.get("integrity_errors"):
        failures.append("LLM trace integrity errors present")
    if trace.get("integrity_level") != "write_time_verified":
        failures.append("LLM trace is not write-time verified")
    if int(trace.get("llm_call_count") or 0) <= 0:
        failures.append("LLM trace contains no model calls")
    if trace.get("missing_decision_links"):
        failures.append("LLM trace has missing decision links")
    if variants["nameless"]:
        failures.append("nameless variant rows remain")
    if variants["wrong_gene_symbol"]:
        failures.append("variant rows carry a conflicting gene symbol")
    if variants["placeholder_paper_titles"]:
        failures.append("papers retain missing or placeholder bibliography titles")
    if counts["negative_rows"]:
        failures.append("negative patient counts remain")
    if counts["impossible_partitions_unquarantined"]:
        failures.append("an unquarantined affected + unaffected exceeds carriers")
    if paper_scope["nonhuman_target_title_links"]:
        failures.append("non-human target-gene papers retain downstream variant links")

    fulltext = int(source.get("papers_with_fulltext") or 0)
    discovered = int(source.get("total_pmids_discovered") or 0)
    abstract_only = int(source.get("papers_abstract_only") or 0)
    if discovered != len(expected):
        failures.append("source-completeness total differs from pinned manifest")
    if fulltext + abstract_only != discovered:
        failures.append("source-completeness classes do not account for every paper")
    active_db = str(status.get("active_db") or "")
    if Path(active_db).name != f"{gene}.db":
        failures.append("RUN_STATUS active_db does not name the audited gene DB")
    return {
        "gene": gene,
        "run_dir": str(run_dir),
        "manifest": str(manifest_path),
        "run": {
            "status": status.get("status"),
            "exit_code": status.get("exit_code"),
            "stage_failures": status.get("stage_failures") or [],
            "stage_warnings": status.get("stage_warnings") or [],
            "duration_seconds": status.get("duration_seconds"),
        },
        "cohort": {
            "expected": len(expected),
            "staged": len(staged),
            "extracted": len(extracted),
            "exact_order": staged == expected,
            "missing_extractions": sorted(set(expected) - set(extracted)),
            "extra_extractions": sorted(set(extracted) - set(expected)),
            "invalid_extractions": extraction_json_issues,
            "missing_db_papers": sorted(set(expected) - db_pmids),
            "extra_db_papers": sorted(db_pmids - set(expected)),
        },
        "source": {
            "fulltext": fulltext,
            "abstract_only": abstract_only,
            "total": discovered,
            "fulltext_fraction": fulltext / discovered if discovered else None,
            "abstract_only_pmids": source.get("abstract_only_pmids") or [],
        },
        "trace": {
            "run_id": trace.get("run_id"),
            "evidence_source": trace.get("evidence_source"),
            "integrity_level": trace.get("integrity_level"),
            "integrity_errors": trace.get("integrity_errors") or [],
            "llm_calls": trace.get("llm_call_count"),
            "decision_events": trace.get("decision_event_count"),
            "missing_decision_links": trace.get("missing_decision_links") or [],
        },
        "vf_gate_passed": vf_gate_passed,
        "variants": variants,
        "counts": counts,
        "trust": trust,
        "variant_features": vf,
        "paper_scope": paper_scope,
        "structural_gate_passed": not failures,
        "structural_gate_failures": failures,
        "publication_recommendation": (
            "HOLD_PENDING_SOURCE_ADJUDICATION_AND_TRUSTED_IMPORT"
        ),
    }


def _markdown(report: dict[str, Any]) -> str:
    lines = [
        "# Collaborator extraction readiness audit",
        "",
        "Structural QC is not an empirical precision estimate. Publication remains "
        "on hold until a source-grounded sample is adjudicated and the Variant "
        "Browser imports trusted fields/quarantines.",
        "",
        "| Gene | Gate | Papers | Full text | Variants | Nameless | Species-scope links | VF matched | VF held | VF quarantined |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for item in report["genes"]:
        vf = item["variant_features"]
        vf_pct = vf["matched"] / vf["rows"] if vf["rows"] else 0
        quarantined = sum(vf["quarantined_variants"].values())
        lines.append(
            f"| {item['gene']} | {'PASS' if item['structural_gate_passed'] else 'FAIL'} "
            f"| {item['cohort']['extracted']}/{item['cohort']['expected']} "
            f"| {item['source']['fulltext']}/{item['source']['total']} "
            f"| {item['variants']['variants']} | {item['variants']['nameless']} "
            f"| {sum(row['variant_links'] for row in item['paper_scope']['nonhuman_target_title_links'])} "
            f"| {vf['matched']}/{vf['rows']} ({vf_pct:.1%}) "
            f"| {vf['blocked_live_variants']} | {quarantined} |"
        )
    lines.extend(["", "## Gate failures", ""])
    for item in report["genes"]:
        failures = item["structural_gate_failures"] or ["none"]
        lines.append(f"- **{item['gene']}:** {'; '.join(failures)}")
    lines.extend(
        [
            "",
            "## Publication decision",
            "",
            "**HOLD.** Passing these invariants is necessary but cannot establish "
            "paper-level precision. Keep public publication disabled pending a "
            "source-grounded sample and a trusted-field-only Variant Browser import.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="append", required=True, help="GENE=RUN_DIR")
    parser.add_argument(
        "--manifest", action="append", required=True, help="GENE=PMID_FILE"
    )
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--out-md", type=Path, required=True)
    args = parser.parse_args()

    runs = _parse_mapping(args.run, "--run")
    manifests = _parse_mapping(args.manifest, "--manifest")
    if set(runs) != set(manifests):
        raise SystemExit("--run and --manifest genes must match exactly")
    report = {
        "schema_version": 1,
        "scope": "collaborator_fixed_manifests",
        "genes": [audit_gene(gene, runs[gene], manifests[gene]) for gene in runs],
        "publication_recommendation": "HOLD_PENDING_SOURCE_ADJUDICATION_AND_TRUSTED_IMPORT",
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_md.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    args.out_md.write_text(_markdown(report))


if __name__ == "__main__":
    main()
