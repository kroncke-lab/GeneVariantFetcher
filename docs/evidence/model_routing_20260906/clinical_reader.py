"""Independent Astra reading plus a source-validated, missing-field-only overlay.

This campaign tool never updates production databases. Its two diagnostic lanes
are locked and scored externally; the raw lane is explicitly unvalidated.
"""

import argparse
import copy
import concurrent.futures
import json
import os
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from campaign import ROOT, GOLD, HARNESS, PYTHON, FROZEN, environment, run, run_dir, sha

sys.path.insert(0, str(ROOT))
from budget_guard import install
from pipeline.count_recovery import (
    PaperGap,
    VariantGap,
    validate_paper_response,
    result_to_dict,
    _notation_matches_target,
)
from utils.llm_trace import (
    configure_llm_tracing,
    llm_trace_scope,
    record_trace_event,
    last_llm_trace,
)
from utils.llm_utils import litellm_completion, parse_llm_json_response

FIELDS = ("carriers", "affected", "unaffected")
SIGNAL = re.compile(
    r"\b(?:carriers?|patients?|probands?|symptomatic|asymptomatic|affected|unaffected)\b",
    re.I,
)


def event(stage, data):
    return record_trace_event(stage, data, stage=stage)


def identity_parts(notation):
    # Match existing explicit HGVS components only; never infer a protein from a codon.
    parts = re.findall(r"(?:p\.|c\.)[^\s]+", notation)
    return tuple(parts or [notation])


def reconcile(paper, payload, source):
    rows = copy.deepcopy(paper["variants"])
    raw_rows = copy.deepcopy(rows)
    targets = [
        VariantGap(
            i,
            row["variant"],
            [f for f in FIELDS if row.get(f) is None],
            parts=identity_parts(row["variant"]),
            paper_derived=True,
        )
        for i, row in enumerate(rows)
    ]
    audit = {"accepted": [], "rejected": [], "contradictions": [], "raw_fills": []}
    candidates = {}
    for proposal in payload.get("variants", []):
        if not isinstance(proposal, dict):
            continue
        notation = str(proposal.get("variant") or "")
        proposal_parts = identity_parts(notation)
        matched = [
            v
            for v in targets
            if any(
                _notation_matches_target(part, v, paper["gene"])
                for part in proposal_parts
            )
        ]
        if len(matched) != 1:
            audit["rejected"].append(
                {
                    "variant": notation,
                    "reason": "no unique retained baseline identity",
                    "matching_rows": len(matched),
                    "proposal": proposal,
                }
            )
            continue
        target = matched[0]
        for field in FIELDS:
            value = proposal.get(field)
            if value is None:
                continue
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                audit["rejected"].append(
                    {
                        "variant": notation,
                        "field": field,
                        "value": value,
                        "reason": "invalid count",
                    }
                )
                continue
            evidence = proposal.get(field + "_evidence") or {}
            if not isinstance(evidence, dict):
                evidence = {}
            candidate = {
                "variant": target.notation,
                field: value,
                "quote": evidence.get("quote"),
                "count_role": evidence.get("count_role"),
                "evidence_locator": evidence.get("evidence_locator"),
            }
            item = {
                "target": target.variant_id,
                "variant": notation,
                "field": field,
                "value": value,
                "entry": candidate,
                "evidence": evidence,
            }
            if field not in target.missing:
                if rows[target.variant_id].get(field) != value:
                    audit["contradictions"].append(item)
                continue
            candidates.setdefault((target.variant_id, field), []).append(item)
    for (idx, field), items in candidates.items():
        values = {x["value"] for x in items}
        if len(values) != 1:
            audit["rejected"].append(
                {
                    "variant": targets[idx].notation,
                    "field": field,
                    "reason": "conflicting independent proposals",
                    "proposals": items,
                }
            )
            continue
        item = items[0]
        raw_rows[idx][field] = item["value"]
        raw_rows[idx]["evidence"] = (
            str(raw_rows[idx].get("evidence") or "")
            + "\nUNVALIDATED independent proposal: "
            + str(item["evidence"].get("quote") or "no quote supplied")
        )
        audit["raw_fills"].append(item)
        evidence = item["evidence"]
        if (
            evidence.get("basis") != "explicit"
            or evidence.get("scope") != "current_study"
        ):
            audit["rejected"].append(
                dict(
                    item,
                    reason="acceptance requires explicit current-study count; derived/prior/unclear retained only in raw diagnostic",
                )
            )
            continue
        result = validate_paper_response(
            PaperGap(paper["gene"], str(paper["pmid"]), [targets[idx]]),
            [item["entry"]],
            source,
        )
        detail = result_to_dict(result)
        audit["accepted"].extend(detail["accepted"])
        audit["rejected"].extend(detail["rejected"])
        for accepted in result.accepted:
            rows[accepted.variant_id][accepted.field] = accepted.value
    # The source validator checks individual fields. This overlay additionally
    # refuses new fills that contradict the retained per-variant person total.
    inconsistent = set()
    for index, row in enumerate(rows):
        carrier = row.get("carriers")
        supplied = [row.get(f) for f in ("affected", "unaffected")]
        if carrier is not None and (
            any(v is not None and v > carrier for v in supplied)
            or (all(v is not None for v in supplied) and sum(supplied) > carrier)
        ):
            inconsistent.add(index)
    retained = []
    for fill in audit["accepted"]:
        index = fill["variant_id"]
        if index in inconsistent:
            rows[index][fill["field"]] = paper["variants"][index].get(fill["field"])
            audit["rejected"].append(
                dict(
                    fill,
                    reason="row person totals conflict; cannot choose which field is wrong, so hold all new fills",
                )
            )
        else:
            retained.append(fill)
            rows[index]["evidence"] = (
                str(rows[index].get("evidence") or "")
                + "\nIndependent Astra evidence: "
                + fill["quote"]
            )
    audit["accepted"] = retained
    for row in rows + raw_rows:
        row.setdefault(
            "inclusion_rationale",
            "Identity retained from the locked baseline trusted paper projection.",
        )
        row["count_rationale"] = str(row.get("count_rationale") or "") + (
            "\n"
            "Original non-null fields retained; null fields changed only by the declared accepted or raw diagnostic overlay. See locked per-paper audit."
        )
    return rows, raw_rows, audit


def main(base_arm):
    base = run_dir(base_arm)
    assert (base / "LOCK.json").is_file(), "Lock baseline before independent reading"
    folder = HARNESS / "runs" / f"20260906_model12_{base_arm}_astra_clinical_verified"
    assert not folder.exists(), folder
    run(
        [
            PYTHON,
            HARNESS / "run_eval.py",
            "prepare",
            "--run-id",
            folder.name,
            "--seed",
            "2026090601",
            "--paper-manifest",
            HERE / "paper_manifest.tsv",
            "--corpus-root",
            FROZEN,
            "--gold-root",
            GOLD,
            "--eligibility-mode",
            "variant",
            "--minimum-chars",
            "1000",
        ],
        HERE / "clinical_prepare.log",
    )
    selection = json.loads((folder / "selection.json").read_text())
    baseline = json.loads((base / "predictions.json").read_text())
    baseline_lock = json.loads((base / "LOCK.json").read_text())
    assert sha(base / "predictions.json") == baseline_lock["predictions_sha256"]
    assert sha(base / "selection.json") == baseline_lock["selection_sha256"]
    by_key = {(p["gene"], str(p["pmid"])): p for p in baseline["papers"]}
    assert {(p["gene"], str(p["pmid"])) for p in selection["papers"]} == set(by_key)
    template = (HERE / "clinical_reader_prompt.txt").read_text()
    triggers = []
    for source_paper in selection["papers"]:
        paper = by_key[(source_paper["gene"], str(source_paper["pmid"]))]
        source = Path(source_paper["source"]).read_text()
        assert sha(Path(source_paper["source"])) == sha(
            base
            / "frozen_corpus"
            / paper["gene"]
            / str(paper["pmid"])
            / f"{paper['pmid']}_FULL_CONTEXT.md"
        )
        missing = sum(row.get(f) is None for row in paper["variants"] for f in FIELDS)
        triggers.append(
            {
                "gene": paper["gene"],
                "pmid": paper["pmid"],
                "missing_fields": missing,
                "clinical_signal": bool(SIGNAL.search(source)),
                "triggered": bool(missing and SIGNAL.search(source)),
                "source_sha256": sha(Path(source_paper["source"])),
            }
        )
    (folder / "trigger_manifest.json").write_text(json.dumps(triggers, indent=2) + "\n")
    setup = {
        "classification": "Opened-calibration independent clinical reader, accepted additive overlay and unvalidated raw diagnostic",
        "base_run": str(base),
        "base_lock_sha256": sha(base / "LOCK.json"),
        "prompt_sha256": sha(HERE / "clinical_reader_prompt.txt"),
        "script_sha256": sha(Path(__file__)),
        "trigger_manifest_sha256": sha(folder / "trigger_manifest.json"),
        "no_baseline_answers_sent_to_reader": True,
        "gold_use": "prepare checks PMID eligibility only; no gold identities or values supplied to reader or overlay",
    }
    (folder / "analysis_setup.json").write_text(json.dumps(setup, indent=2) + "\n")
    os.environ.update(environment("astra_medium_verified"))
    os.environ["GVF_EXPERIMENT_ARM"] = base_arm + "_astra_clinical_verified"
    install()
    configure_llm_tracing(folder / "llm_traces", run_id=folder.name)
    output = copy.deepcopy(baseline)
    output.update(
        schema_version=2,
        run_id=folder.name,
        strategy="independent_astra_clinical_additive_overlay",
        primary_score_lane="source_validated_additive",
        comparison_score_lanes=["reader_raw_additive"],
        papers=[],
    )
    output["clinical_reader_setup"] = setup
    output["clinical_reader_audits"] = []
    start = time.monotonic()

    def read_paper(pair):
        source_paper, trigger = pair
        gene, pmid = source_paper["gene"], str(source_paper["pmid"])
        paper = copy.deepcopy(by_key[(gene, pmid)])
        source = Path(source_paper["source"]).read_text()
        assert sha(Path(source_paper["source"])) == trigger["source_sha256"]
        t = time.monotonic()
        with llm_trace_scope(gene=gene, pmid=pmid):
            refs = [
                event(
                    "representation_route",
                    {
                        "method": "fixed source-only text; no model routing decision",
                        "trigger": trigger,
                    },
                ),
                event(
                    "representation_route_decision",
                    {
                        "tool": "text",
                        "source_sha256": trigger["source_sha256"],
                        "source_chars": len(source),
                    },
                ),
            ]
            payload = {"variants": []}
            usage = (
                None
                if trigger["triggered"]
                else {"prompt_tokens": 0, "completion_tokens": 0, "total_tokens": 0}
            )
            error = None
            if trigger["triggered"]:
                prompt = (
                    template.replace("{gene}", gene)
                    .replace("{pmid}", pmid)
                    .replace("{source}", source)
                )
                try:
                    with llm_trace_scope(
                        stage="paper_curation", component="independent_clinical_reader"
                    ):
                        response = litellm_completion(
                            model="azure_ai/gpt-6-astra",
                            messages=[{"role": "user", "content": prompt}],
                            max_tokens=32000,
                            reasoning_effort="medium",
                            temperature=0,
                            response_format={"type": "json_object"},
                            timeout=1200,
                            num_retries=0,
                        )
                        refs.append(last_llm_trace())
                        usage = response.usage.model_dump()
                        if response.choices[0].finish_reason != "stop":
                            raise RuntimeError(
                                "Independent response incomplete: "
                                + str(response.choices[0].finish_reason)
                            )
                        payload = parse_llm_json_response(
                            response.choices[0].message.content
                        )
                        if not isinstance(payload, dict):
                            raise RuntimeError("Independent JSON was not an object")
                except Exception as exc:
                    error = type(exc).__name__ + ": " + str(exc)[:600]
                    payload = {"variants": []}
                    refs.append(
                        event("paper_curation", {"status": "failed", "error": error})
                    )
            else:
                refs.append(
                    event(
                        "paper_curation",
                        {"status": "not_triggered", "retained_baseline": True},
                    )
                )
            accepted, raw, audit = reconcile(paper, payload, source)
            audit.update(gene=gene, pmid=pmid, trigger=trigger, error=error)
            refs.append(event("paper_curation_decision", audit))
        paper.update(
            variants=accepted,
            comparison_variants={"reader_raw_additive": raw},
            tool="text",
            tool_rationale="Frozen source-only independent reading when clinical gaps trigger; otherwise retain baseline.",
            source_completeness=source_paper.get("source_completeness")
            or paper.get("source_completeness"),
            elapsed_seconds=time.monotonic() - t,
            curation_rationale="Missing-field-only overlay; accepted values pass existing literal source validator. Original counts and identities stay fixed.",
            llm_trace_refs=[r for r in refs if r],
            token_usage={
                "telemetry_available": usage is not None,
                "input_tokens": usage.get("prompt_tokens") if usage else None,
                "output_tokens": usage.get("completion_tokens") if usage else None,
                "total_tokens": usage.get("total_tokens") if usage else None,
                "note": "Incremental reader usage only; unknown failed-call usage is null, with reservation retained in campaign ledger.",
            },
        )
        print(
            gene,
            pmid,
            "trigger",
            trigger["triggered"],
            "accepted",
            len(audit["accepted"]),
            "raw",
            len(audit["raw_fills"]),
            "error",
            error,
            flush=True,
        )
        for row_set in ("external_linkage_variants", "unattributed_variants"):
            for row in paper.get(row_set, []):
                row["inclusion_rationale"] = (
                    "Retained baseline secondary provenance row; excluded from primary overlay."
                )
                row["count_rationale"] = (
                    "No independent clinical changes applied to this secondary row."
                )
        return paper, audit

    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        futures = [
            pool.submit(read_paper, pair) for pair in zip(selection["papers"], triggers)
        ]
        for future in concurrent.futures.as_completed(futures):
            paper, audit = future.result()
            output["papers"].append(paper)
            output["clinical_reader_audits"].append(audit)
            (folder / "clinical_progress.json").write_text(
                json.dumps(output, indent=2) + "\n"
            )
    output["papers"].sort(key=lambda p: (p["gene"], str(p["pmid"])))
    output["completed_at"] = datetime.now(timezone.utc).isoformat()
    output["extraction_elapsed_seconds"] = time.monotonic() - start
    output["token_usage"]["note"] = (
        "Baseline production usage only; incremental clinical usage is exact per paper and in campaign budget ledger."
    )
    (folder / "predictions.json").write_text(json.dumps(output, indent=2) + "\n")
    run(
        [PYTHON, HARNESS / "run_eval.py", "lock", "--run-dir", folder],
        HERE / "clinical_lock.log",
    )
    print("LOCKED CLINICAL", folder, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base-arm", default="grok46_verified", choices=["grok43", "grok46_verified"]
    )
    args = parser.parse_args()
    main(args.base_arm)
