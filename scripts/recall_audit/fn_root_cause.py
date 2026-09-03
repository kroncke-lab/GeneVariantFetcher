#!/usr/bin/env python3
"""Root-cause every paper-derived false negative of a scored production run.

For each gold row the locked report lists under ``missed_gold``, walk the
pipeline stages the run left on disk and find the first one where the variant
disappears. Every check is a string search for the gold row's notation forms,
so the tool never rewrites a score; it only says *where* a miss happened so the
owning fix is unambiguous.

Stages, in pipeline order
-------------------------
``corpus``      any on-disk representation under ``corpus/<GENE>/<PMID>/``
                (from the source-presence sweep's per-row class)
``run_text``    the run-local ``pmc_fulltext/<PMID>_FULL_CONTEXT.md`` /
                ``_CLEANED.md`` / ``_DATA_ZONES.md`` the run staged
``llm_request`` the user message of any ``paper_variant_extraction`` LLM call
                (what the model actually saw)
``llm_response`` the raw ``output_text`` of any such call (what the model said)
``extraction``  the per-paper extraction JSON (after parsing/merging)
``db``          a live ``variants`` row joined to this PMID, with its
                ``source_layer`` origins
``predictions`` the locked ``predictions.json`` paper-derived lane; the
                linkage/unattributed lanes are reported separately

Leaf (owner) is the first missing stage:

``acquisition``       absent from corpus (sweep class in the hard ceiling)
``unknown_notation``  absent from all text but notation is not string-searchable,
                      or figures are on disk (sweep unknown classes)
``source_selection``  in corpus, not in the text the run staged
``condensing``        in run text, not in any LLM request (scout/data-zones cut,
                      or a deterministic table route that dropped the row)
``model_missed``      in the request, not in any response
``parser_dropped``    in a response, not in the extraction JSON
``postprocess_dropped`` in extraction JSON, not a live DB row (verifier/trust/migrate)
``paper_row_lost_to_linkage_origin`` live in the DB and in the linkage lane only:
                      the paper layer found it but the surviving row is attributed
                      to ClinVar/PubTator, so the primary lane never saw it
``projection_dropped`` live in DB, in no lane
``matcher``           in the paper-derived lane but the scorer did not pair it
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from benchmarks.codex_paper_eval.run_eval import matches  # noqa: E402
from scripts.recall_audit.gold_source_presence_sweep import (  # noqa: E402
    HARD_CEILING,
    WIDE_CEILING,
    expand_forms,
)
from scripts.recall_audit.tier_source_reachability import FormIndex, _squash  # noqa: E402

NOTATION_KEYS = (
    "cdna_notation",
    "protein_notation",
    "legacy_notation",
    "source_notation",
)


def _production_run_dir(run_dir: Path, gene: str) -> Path | None:
    candidates = sorted((run_dir / "production_runs" / gene).glob("*/"))
    return candidates[-1] if candidates else None


def _read(path: Path) -> str:
    try:
        return path.read_text(errors="replace")
    except OSError:
        return ""


class RunPaper:
    """Lazy, cached view of everything the run left on disk for one paper."""

    def __init__(self, run_dir: Path, gene: str, pmid: str):
        self.gene, self.pmid = gene, pmid
        self.prod = _production_run_dir(run_dir, gene)
        self._text: dict[str, str] = {}
        self._traces: list[dict] | None = None
        self._extraction: list[str] | None = None
        self._db: tuple[list[str], list[str], list[tuple[str, str]]] | None = None

    def run_text(self) -> str:
        if "run_text" not in self._text:
            parts = []
            if self.prod:
                for suffix in ("_FULL_CONTEXT.md", "_CLEANED.md", "_DATA_ZONES.md"):
                    parts.append(
                        _read(self.prod / "pmc_fulltext" / f"{self.pmid}{suffix}")
                    )
                    parts.append(
                        _read(self.prod / "scout_output" / f"{self.pmid}{suffix}")
                    )
            self._text["run_text"] = _squash("\n".join(parts))
        return self._text["run_text"]

    def traces(self) -> list[dict]:
        if self._traces is None:
            out = []
            if self.prod:
                pattern = str(
                    self.prod / "llm_traces" / self.gene / self.pmid / "*.json"
                )
                for f in glob.glob(pattern):
                    try:
                        out.append(json.load(open(f)))
                    except Exception:  # noqa: BLE001
                        continue
            self._traces = out
        return self._traces

    def llm_request(self) -> str:
        if "llm_request" not in self._text:
            parts = []
            for t in self.traces():
                if t.get("record_type") != "llm_call":
                    continue
                payload = (t.get("request") or {}).get("payload") or {}
                for m in payload.get("messages") or []:
                    content = m.get("content")
                    if isinstance(content, list):
                        content = " ".join(
                            c.get("text", "") for c in content if isinstance(c, dict)
                        )
                    parts.append(str(content or ""))
            self._text["llm_request"] = _squash("\n".join(parts))
        return self._text["llm_request"]

    def llm_response(self) -> str:
        if "llm_response" not in self._text:
            parts = [
                str((t.get("response") or {}).get("output_text") or "")
                for t in self.traces()
                if t.get("record_type") == "llm_call"
            ]
            self._text["llm_response"] = _squash("\n".join(parts))
        return self._text["llm_response"]

    def extraction_variants(self) -> list[str]:
        if self._extraction is None:
            names: list[str] = []
            if self.prod:
                pattern = str(self.prod / "extractions" / f"*PMID_{self.pmid}.json")
                for f in glob.glob(pattern):
                    try:
                        d = json.load(open(f))
                    except Exception:  # noqa: BLE001
                        continue
                    for v in d.get("variants") or []:
                        names.extend(str(v[k]) for k in NOTATION_KEYS if v.get(k))
            self._extraction = names
        return self._extraction

    def db_variants(self) -> tuple[list[str], list[str], list[tuple[str, str]]]:
        """(live notations, quarantined notations, (notation, source_layer))."""
        if self._db is None:
            live: list[str] = []
            quarantined: list[str] = []
            layers: list[tuple[str, str]] = []
            db_path = None
            if self.prod:
                try:
                    db_name = json.load(open(self.prod / "RUN_STATUS.json")).get(
                        "active_db"
                    )
                    db_path = self.prod / db_name if db_name else None
                except Exception:  # noqa: BLE001
                    db_path = None
            if db_path and db_path.is_file():
                con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
                try:
                    rows = con.execute(
                        "select v.cdna_notation, v.protein_notation, v.legacy_notation, "
                        "vp.source_notation, vp.source_layer "
                        "from variants v join variant_papers vp on vp.variant_id=v.variant_id "
                        "where vp.pmid=?",
                        (self.pmid,),
                    ).fetchall()
                    for row in rows:
                        for x in row[:4]:
                            if x:
                                live.append(str(x))
                                layers.append((str(x), str(row[4] or "")))
                    try:
                        qrows = con.execute(
                            "select * from quarantined_variants where pmid=?",
                            (self.pmid,),
                        ).fetchall()
                        quarantined = [
                            str(x) for row in qrows for x in row if isinstance(x, str)
                        ]
                    except sqlite3.Error:
                        quarantined = []
                finally:
                    con.close()
            self._db = (live, quarantined, layers)
        return self._db


def _any_form_in(forms: frozenset[str], squashed_text: str) -> bool:
    return any(f in squashed_text for f in forms)


def _any_match(gold: str, names: list[str], gene: str) -> bool:
    for n in names:
        try:
            if matches(gold, n, gene):
                return True
        except Exception:  # noqa: BLE001
            continue
    return False


def classify_fn(
    gene: str,
    pmid: str,
    gold_variant: str,
    paper: RunPaper,
    paper_lane_names: list[str],
    linkage_names: list[str],
    sweep_class: str,
    forms: FormIndex,
) -> dict:
    fset = expand_forms(forms.get(gene, gold_variant))
    stages: dict[str, bool | None] = {
        "corpus": (sweep_class not in WIDE_CEILING) if sweep_class else None,
        "run_text": _any_form_in(fset, paper.run_text()),
        "llm_request": _any_form_in(fset, paper.llm_request()),
        "llm_response": _any_form_in(fset, paper.llm_response()),
        "extraction": _any_match(gold_variant, paper.extraction_variants(), gene),
    }
    live, quarantined, layers = paper.db_variants()
    stages["db"] = _any_match(gold_variant, live, gene)
    stages["db_quarantined"] = _any_match(gold_variant, quarantined, gene)
    db_layers = sorted(
        {layer for name, layer in layers if _any_match(gold_variant, [name], gene)}
    )
    stages["paper_lane"] = _any_match(gold_variant, paper_lane_names, gene)
    stages["linkage_lane"] = _any_match(gold_variant, linkage_names, gene)

    if sweep_class in HARD_CEILING:
        leaf = "acquisition"
    elif sweep_class in WIDE_CEILING and not stages["run_text"]:
        leaf = "unknown_notation"
    elif not stages["run_text"]:
        leaf = "source_selection"
    elif not stages["llm_request"]:
        leaf = "condensing"
    elif not stages["llm_response"]:
        leaf = "model_missed"
    elif not stages["extraction"]:
        leaf = "parser_dropped"
    elif not stages["db"]:
        leaf = "postprocess_dropped" + (
            "(quarantined)" if stages["db_quarantined"] else ""
        )
    elif not stages["paper_lane"]:
        leaf = (
            "paper_row_lost_to_linkage_origin"
            if stages["linkage_lane"]
            else "projection_dropped"
        )
    else:
        leaf = "matcher"
    return {
        "gene": gene,
        "pmid": pmid,
        "variant": gold_variant,
        "sweep_class": sweep_class,
        "leaf": leaf,
        "db_layers": ",".join(db_layers),
        **{f"in_{k}": v for k, v in stages.items()},
    }


def _lane_names(predictions: dict) -> tuple[dict, dict]:
    paper_lane: dict[tuple[str, str], list[str]] = defaultdict(list)
    linkage: dict[tuple[str, str], list[str]] = defaultdict(list)
    for p in predictions.get("papers", []):
        key = (p["gene"], str(p["pmid"]))
        for v in p.get("variants") or []:
            paper_lane[key].append(str(v.get("variant") or ""))
        for lane in ("external_linkage_variants", "unattributed_variants"):
            for v in p.get(lane) or []:
                linkage[key].append(str(v.get("variant") or ""))
    return paper_lane, linkage


def render(report: dict, rows: list[dict]) -> str:
    leaves = Counter(r["leaf"] for r in rows)
    by_paper: dict = defaultdict(Counter)
    for r in rows:
        by_paper[(r["gene"], r["pmid"])][r["leaf"]] += 1
    lines = [
        f"# FN root cause — `{report.get('run_id')}`",
        "",
        f"{len(rows)} paper-derived false negatives.",
        "",
        "| leaf | rows |",
        "|---|---:|",
    ]
    for leaf, n in leaves.most_common():
        lines.append(f"| `{leaf}` | {n} |")
    lines += ["", "| gene | PMID | FN | leaves |", "|---|---:|---:|---|"]
    for (g, p), c in sorted(by_paper.items(), key=lambda kv: -sum(kv[1].values())):
        lines.append(
            f"| {g} | {p} | {sum(c.values())} | "
            + ", ".join(f"{k}={v}" for k, v in c.most_common())
            + " |"
        )
    lines += [
        "",
        "## Rows the reading protocol could have found",
        "",
        "| gene | PMID | variant | leaf | db layers | run_text | request | response | extraction | db | paper lane | linkage lane |",
        "|---|---:|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        if r["leaf"] in ("acquisition", "unknown_notation"):
            continue
        lines.append(
            f"| {r['gene']} | {r['pmid']} | `{r['variant']}` | {r['leaf']} | {r['db_layers']} | "
            f"{r['in_run_text']} | {r['in_llm_request']} | {r['in_llm_response']} | "
            f"{r['in_extraction']} | {r['in_db']} | {r['in_paper_lane']} | {r['in_linkage_lane']} |"
        )
    lines += [
        "",
        "## Notation-unknown rows (probe could not search the notation)",
        "",
        "| gene | PMID | variant | sweep class |",
        "|---|---:|---|---|",
    ]
    for r in rows:
        if r["leaf"] == "unknown_notation":
            lines.append(
                f"| {r['gene']} | {r['pmid']} | `{r['variant']}` | {r['sweep_class']} |"
            )
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--run-dir", type=Path, required=True)
    ap.add_argument(
        "--sweep-rows",
        type=Path,
        default=REPO_ROOT
        / "docs/evidence/gold_source_presence_sweep_20260903/gold_rows.tsv",
        help="per-row output of gold_source_presence_sweep.py",
    )
    ap.add_argument("--out-dir", type=Path, default=None)
    args = ap.parse_args()

    report = json.load(open(args.run_dir / "report.json"))
    predictions = json.load(open(args.run_dir / "predictions.json"))
    paper_lane, linkage = _lane_names(predictions)

    sweep: dict[tuple[str, str, str], str] = {}
    if args.sweep_rows.is_file():
        for r in csv.DictReader(open(args.sweep_rows), delimiter="\t"):
            sweep[(r["gene"], r["pmid"], r["variant"])] = r["class"]

    forms = FormIndex()
    rows: list[dict] = []
    for paper in report["papers"]:
        gene, pmid = paper["gene"], str(paper["pmid"])
        missed = paper.get("missed_gold") or []
        if not missed:
            continue
        rp = RunPaper(args.run_dir, gene, pmid)
        key = (gene, pmid)
        for gold_variant in missed:
            rows.append(
                classify_fn(
                    gene,
                    pmid,
                    gold_variant,
                    rp,
                    paper_lane.get(key, []),
                    linkage.get(key, []),
                    sweep.get((gene, pmid, gold_variant), ""),
                    forms,
                )
            )

    out_dir = args.out_dir or (args.run_dir / "diagnostics" / "fn_root_cause")
    out_dir.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys()) if rows else ["gene", "pmid", "variant", "leaf"]
    with (out_dir / "fn_root_cause.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    leaves = Counter(r["leaf"] for r in rows)
    by_paper: dict = defaultdict(Counter)
    for r in rows:
        by_paper[f"{r['gene']}:{r['pmid']}"][r["leaf"]] += 1
    summary = {
        "run_id": report.get("run_id"),
        "fn_total": len(rows),
        "by_leaf": dict(leaves),
        "by_paper": {
            k: dict(c)
            for k, c in sorted(by_paper.items(), key=lambda kv: -sum(kv[1].values()))
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=1))
    (out_dir / "summary.md").write_text(render(report, rows))
    print(json.dumps(summary["by_leaf"], indent=1))
    print(f"wrote {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
