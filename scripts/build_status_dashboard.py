#!/usr/bin/env python3
"""Build the GVF strategy/status dashboard (docs/dashboard/index.html).

This renders a SELF-CONTAINED, offline HTML dashboard from the metrics the
recall suite already produces. It has two clearly separated halves:

1. DETERMINISTIC DATA (this script): recall / precision / MAE across every gene,
   the gap-to-90% scorecard, per-source-layer precision, the failure-mode
   ("barriers") split, and a recall trajectory sparkline. Pulled straight from
   ``recall_metrics/<run>/summary.json`` + its ``paper_disagreement_report.csv``.
   This always renders, even with no LLM involvement.

2. STRATEGY NARRATIVE (optional sidecar ``docs/dashboard/strategy.json``): the
   overnight reasoning pass writes barriers-to-success, things-to-try, and
   candidate-solutions here, optionally with Grok / Claude second opinions.
   If the sidecar is missing, the section shows a placeholder.

Usage:
    python scripts/build_status_dashboard.py \
        [--metrics-dir recall_metrics/<run>] \
        [--repo-root .] \
        [--out docs/dashboard/index.html]

With no --metrics-dir it auto-selects the newest valid summary.json (one that
carries both aggregate_recall and gene_results), so a fresh overnight recompute
is picked up automatically.

No network. No third-party deps. Charts are inline SVG.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import html
import json
import sys
from pathlib import Path
from typing import Any, Optional

TARGET = 0.90
METRIC_ORDER = [
    ("pmids", "PMIDs"),
    ("variant_rows", "Variant rows"),
    ("unique_variants", "Unique variants"),
    ("patients", "Patients/carriers"),
    ("affected", "Affected"),
    ("unaffected", "Unaffected"),
]


def esc(x: Any) -> str:
    return html.escape("" if x is None else str(x))


def pct(x: Optional[float]) -> str:
    return "-" if x is None else f"{x * 100:.1f}%"


def band(recall: Optional[float]) -> str:
    """CSS class by distance to the 90% target."""
    if recall is None:
        return "na"
    if recall >= TARGET:
        return "good"
    if recall >= 0.80:
        return "warn"
    return "bad"


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def _valid_summary(p: Path) -> Optional[dict]:
    try:
        d = json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return None
    if isinstance(d, dict) and "aggregate_recall" in d and "gene_results" in d:
        return d
    return None


def find_metrics(repo_root: Path, override: Optional[str]) -> tuple[dict, Path]:
    if override:
        d = repo_root / override
        p = d / "summary.json" if d.is_dir() else d
        s = _valid_summary(p)
        if not s:
            sys.exit(f"error: {p} is not a valid recall summary.json")
        return s, p.parent
    candidates = list((repo_root / "recall_metrics").rglob("summary.json"))
    scored: list[tuple[str, Path, dict]] = []
    for p in candidates:
        s = _valid_summary(p)
        if s:
            ts = s.get("generated_at_utc", "")
            # fall back to mtime when timestamp missing
            key = ts or _dt.datetime.utcfromtimestamp(p.stat().st_mtime).isoformat()
            scored.append((key, p, s))
    if not scored:
        sys.exit("error: no valid recall summary.json under recall_metrics/")
    scored.sort(key=lambda t: t[0])
    _, p, s = scored[-1]
    return s, p.parent


def load_failure_split(metrics_dir: Path) -> list[tuple[str, int]]:
    csv_path = metrics_dir / "paper_disagreement_report.csv"
    if not csv_path.exists():
        return []
    counts: dict[str, int] = {}
    try:
        with csv_path.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                fc = (row.get("failure_class") or "").strip()
                if fc:
                    counts[fc] = counts.get(fc, 0) + 1
    except Exception:
        return []
    return sorted(counts.items(), key=lambda kv: kv[1], reverse=True)


# Human-readable gloss for the raw failure_class labels.
FAILURE_GLOSS = {
    "source_missing_or_stub": "Paper/source never landed or only a stub landed.",
    "source_abstract_only": "Abstract available, but mutation tables/body missing.",
    "available_source_underextraction": "Usable source exists but extraction missed rows.",
    "source_missing_table_bodies": "Full text landed without the relevant tables.",
    "partial_underextraction": "Some rows extracted, table not exhausted.",
    "count_semantics": "Variant present but carrier/affected/unaffected semantics wrong.",
    "overinclusive_extraction": "DB has many extra rows for the PMID.",
    "pipeline_only_extras": "Rows the pipeline added that are not in the counted gold packet.",
    "matched_or_low_priority": "Matched or low-priority residual; not a recall lever.",
    "phenotype_or_count_direction": "Phenotype/count-direction disagreement.",
}
# Failure classes that represent a genuine recall loss worth acting on.
RECALL_BARRIERS = {
    "source_missing_or_stub",
    "source_abstract_only",
    "available_source_underextraction",
    "source_missing_table_bodies",
    "partial_underextraction",
    "count_semantics",
}


def detect_nogold_genes(repo_root: Path, scored: set[str]) -> list[str]:
    out: set[str] = set()
    results = repo_root / "results"
    if results.is_dir():
        for child in results.iterdir():
            name = child.name
            if child.is_dir() and name.isupper() and name.isalnum() and len(name) >= 3:
                if name not in scored and name not in {"PURPOSE"}:
                    out.add(name)
    return sorted(out)


def load_strategy(dashboard_dir: Path) -> Optional[dict]:
    p = dashboard_dir / "strategy.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return None


def gather_trend(repo_root: Path, limit: int = 20) -> list[dict]:
    """Aggregate unique_variants + pmids recall over time from all summaries."""
    pts: list[dict] = []
    for p in (repo_root / "recall_metrics").rglob("summary.json"):
        s = _valid_summary(p)
        if not s:
            continue
        ts = s.get("generated_at_utc")
        if not ts:
            continue
        agg = s.get("aggregate_recall", {})
        uv = (agg.get("unique_variants") or {}).get("recall")
        pm = (agg.get("pmids") or {}).get("recall")
        if uv is None:
            continue
        pts.append({"ts": ts, "unique_variants": uv, "pmids": pm})
    # de-dup by timestamp, sort, keep the last `limit`
    seen: dict[str, dict] = {}
    for pt in pts:
        seen[pt["ts"]] = pt
    ordered = sorted(seen.values(), key=lambda d: d["ts"])
    return ordered[-limit:]


# ---------------------------------------------------------------------------
# SVG helpers (inline, no deps)
# ---------------------------------------------------------------------------
def svg_bar(recall: Optional[float], width: int = 160) -> str:
    if recall is None:
        return "<span class='muted'>n/a</span>"
    filled = int(round(recall * width))
    target_x = int(round(TARGET * width))
    cls = band(recall)
    return (
        f"<svg class='bar' width='{width}' height='14' viewBox='0 0 {width} 14'>"
        f"<rect x='0' y='2' width='{width}' height='10' rx='3' fill='var(--track)'/>"
        f"<rect x='0' y='2' width='{filled}' height='10' rx='3' class='fill-{cls}'/>"
        f"<line x1='{target_x}' y1='0' x2='{target_x}' y2='14' stroke='var(--target)' "
        f"stroke-width='2' stroke-dasharray='2,1'/></svg>"
    )


def svg_trend(points: list[dict], width: int = 640, height: int = 150) -> str:
    if len(points) < 2:
        return "<p class='muted'>Not enough history yet for a trajectory line.</p>"
    pad = 28
    xs = [pad + i * (width - 2 * pad) / (len(points) - 1) for i in range(len(points))]

    def y(v: float) -> float:
        lo, hi = 0.4, 1.0  # fixed recall window for readability
        v = max(lo, min(hi, v))
        return height - pad - (v - lo) / (hi - lo) * (height - 2 * pad)

    def path(key: str, cls: str) -> str:
        pts = [
            (xs[i], y(points[i][key]))
            for i in range(len(points))
            if points[i].get(key) is not None
        ]
        if len(pts) < 2:
            return ""
        d = "M " + " L ".join(f"{x:.1f},{yy:.1f}" for x, yy in pts)
        dots = "".join(
            f"<circle cx='{x:.1f}' cy='{yy:.1f}' r='2.5' class='{cls}'/>"
            for x, yy in pts
        )
        return f"<path d='{d}' fill='none' class='{cls}' stroke-width='2'/>{dots}"

    target_y = y(TARGET)
    grid = "".join(
        f"<line x1='{pad}' y1='{y(v):.1f}' x2='{width - pad}' y2='{y(v):.1f}' class='grid'/>"
        f"<text x='4' y='{y(v) + 3:.1f}' class='axis'>{int(v * 100)}%</text>"
        for v in (0.5, 0.7, 0.9)
    )
    return (
        f"<svg width='100%' viewBox='0 0 {width} {height}' class='trend'>"
        f"{grid}"
        f"<line x1='{pad}' y1='{target_y:.1f}' x2='{width - pad}' y2='{target_y:.1f}' class='targetline'/>"
        f"{path('unique_variants', 'line-uv')}"
        f"{path('pmids', 'line-pm')}"
        f"</svg>"
        f"<div class='legend'><span class='k uv'></span>Unique-variant recall "
        f"<span class='k pm'></span>PMID recall "
        f"<span class='k tg'></span>90% target</div>"
    )


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------
def render_strategy(strategy: Optional[dict]) -> str:
    if not strategy:
        return (
            "<div class='card placeholder'><h3>Strategy narrative</h3>"
            "<p class='muted'>No <code>strategy.json</code> yet. The scheduled overnight "
            "run writes barriers, things to try, and candidate solutions here "
            "(optionally with Grok / Claude second opinions), then rebuilds this page.</p></div>"
        )
    gen = strategy.get("generated_at", "")
    models = strategy.get("models") or []
    parts = ["<div class='card strategy'><h3>Strategy narrative</h3>"]
    if gen or models:
        tag = f"generated {esc(gen)}" if gen else ""
        if models:
            tag += (" · " if tag else "") + "models: " + esc(", ".join(models))
        parts.append(f"<p class='muted small'>{tag}</p>")

    def block(title: str, key: str) -> str:
        items = strategy.get(key) or []
        if not items:
            return ""
        lis = "".join(
            f"<li>{esc(it)}</li>"
            if isinstance(it, str)
            else f"<li><b>{esc(it.get('title', ''))}</b> — {esc(it.get('detail', ''))}</li>"
            for it in items
        )
        return f"<h4>{esc(title)}</h4><ul>{lis}</ul>"

    parts.append(block("Barriers to success", "barriers"))
    parts.append(block("Things to try", "ideas"))
    parts.append(block("Candidate solutions", "solutions"))

    notes = strategy.get("model_notes") or {}
    if notes:
        parts.append("<h4>Model second opinions</h4><div class='notes'>")
        for model, text in notes.items():
            parts.append(
                f"<div class='note'><div class='note-h'>{esc(model)}</div>"
                f"<div class='note-b'>{esc(text)}</div></div>"
            )
        parts.append("</div>")
    parts.append("</div>")
    return "".join(parts)


def render(
    summary: dict, metrics_dir: Path, repo_root: Path, strategy: Optional[dict]
) -> str:
    agg = summary["aggregate_recall"]
    gaps = summary.get("aggregate_target_gaps", {})
    mae = summary.get("aggregate_mae", {})
    prec = summary.get("aggregate_precision", {})
    genes = summary["gene_results"]
    scored_names = {g["gene"] for g in genes}
    nogold = detect_nogold_genes(repo_root, scored_names)
    failures = load_failure_split(metrics_dir)
    trend = gather_trend(repo_root)
    gen_utc = summary.get("generated_at_utc", "")
    now = _dt.datetime.now().strftime("%Y-%m-%d %H:%M")

    # --- aggregate scorecard ---
    score_cards = []
    for key, label in METRIC_ORDER:
        m = agg.get(key, {})
        r = m.get("recall")
        gap = gaps.get(key)
        gap_txt = (
            f"+{gap} to 90%"
            if isinstance(gap, (int, float)) and gap > 0
            else "at target"
        )
        score_cards.append(
            f"<div class='sc {band(r)}'><div class='sc-l'>{esc(label)}</div>"
            f"<div class='sc-v'>{pct(r)}</div>"
            f"<div class='sc-s'>{esc(m.get('matched', '?'))}/{esc(m.get('gold', '?'))} · {esc(gap_txt)}</div>"
            f"{svg_bar(r, 150)}</div>"
        )

    # --- MAE cards ---
    mae_cards = []
    for key in ("carriers", "affected", "unaffected"):
        m = mae.get(key, {})
        v = m.get("mae")
        mae_cards.append(
            f"<div class='mc'><div class='sc-l'>MAE {esc(key)}</div>"
            f"<div class='sc-v'>{'-' if v is None else f'{v:.3f}'}</div>"
            f"<div class='sc-s'>n={esc(m.get('n_matched', '?'))}</div></div>"
        )

    # --- per-gene table ---
    rows = []
    for g in sorted(genes, key=lambda g: g["gene"]):
        s = g.get("summary", {})
        rec = s.get("recall", {})
        gmae = s.get("mae", {}).get("carriers", {}).get("mae")
        cells = [f"<td class='gene'>{esc(g['gene'])}</td>"]
        for key, _ in METRIC_ORDER:
            r = (rec.get(key) or {}).get("recall")
            cells.append(f"<td class='{band(r)}'>{pct(r)}<br>{svg_bar(r, 90)}</td>")
        cells.append(f"<td>{'-' if gmae is None else f'{gmae:.3f}'}</td>")
        rows.append(f"<tr>{''.join(cells)}</tr>")
    for gn in nogold:
        cells = [f"<td class='gene'>{esc(gn)}</td>"]
        cells += ["<td class='na'>no gold</td>"] * len(METRIC_ORDER)
        cells.append("<td class='na'>—</td>")
        rows.append(f"<tr class='nogold'>{''.join(cells)}</tr>")
    head = "".join(f"<th>{esc(lbl)}</th>" for _, lbl in METRIC_ORDER)

    # --- precision by source layer ---
    layer_rows = []
    for layer, d in sorted((prec.get("by_source_layer") or {}).items()):
        if not isinstance(d, dict):
            continue
        layer_rows.append(
            f"<tr><td>{esc(layer)}</td>"
            f"<td>{esc(d.get('matched_db', '?'))}</td>"
            f"<td>{esc(d.get('extra_rows', d.get('extra', '?')))}</td>"
            f"<td>{esc(d.get('counted_extra_rows', d.get('counted_extra', '?')))}</td>"
            f"<td>{pct(d.get('precision_vs_counted_gold_pmids'))}</td></tr>"
        )
    headline_prec = prec.get("precision_vs_counted_gold_pmids")

    # --- barriers (failure split) ---
    bar_rows = []
    total_missing = sum(c for _, c in failures) or 1
    for fc, c in failures:
        is_lever = fc in RECALL_BARRIERS
        w = int(round(c / total_missing * 100))
        bar_rows.append(
            f"<tr class='{'lever' if is_lever else 'nonlever'}'>"
            f"<td>{esc(fc)}</td><td class='num'>{c}</td>"
            f"<td><div class='hbar'><div style='width:{w}%'></div></div></td>"
            f"<td class='gloss'>{esc(FAILURE_GLOSS.get(fc, ''))}</td></tr>"
        )

    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GVF Recall Status Dashboard</title>
<style>
:root{{
  --bg:#0f1216; --panel:#171b21; --panel2:#1d232b; --ink:#e7edf3; --muted:#95a1af;
  --line:#2a323c; --track:#232b34; --target:#f4c430;
  --good:#3fb26b; --warn:#e0a13a; --bad:#e0573f; --uv:#5aa9e6; --pm:#a98be6;
}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--ink);
  font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}}
.wrap{{max-width:1100px;margin:0 auto;padding:28px 20px 60px}}
h1{{font-size:22px;margin:0 0 2px}} h2{{font-size:15px;margin:30px 0 12px;color:var(--muted);
  text-transform:uppercase;letter-spacing:.06em;font-weight:600}}
h3{{margin:0 0 10px;font-size:16px}} h4{{margin:16px 0 6px;font-size:13px;color:var(--muted);
  text-transform:uppercase;letter-spacing:.04em}}
.sub{{color:var(--muted);font-size:12.5px;margin-bottom:4px}}
.grid-sc{{display:grid;grid-template-columns:repeat(auto-fit,minmax(165px,1fr));gap:10px}}
.sc,.mc{{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:12px 14px}}
.sc-l{{color:var(--muted);font-size:12px}} .sc-v{{font-size:24px;font-weight:700;margin:2px 0}}
.sc-s{{color:var(--muted);font-size:11.5px;margin-bottom:8px}}
.sc.good{{border-left:3px solid var(--good)}} .sc.warn{{border-left:3px solid var(--warn)}}
.sc.bad{{border-left:3px solid var(--bad)}} .sc.na{{opacity:.6}}
.grid-mae{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin-top:10px}}
.bar .fill-good{{fill:var(--good)}} .bar .fill-warn{{fill:var(--warn)}} .bar .fill-bad{{fill:var(--bad)}}
table{{width:100%;border-collapse:collapse;background:var(--panel);border:1px solid var(--line);
  border-radius:10px;overflow:hidden;font-size:13px}}
th,td{{padding:8px 10px;text-align:left;border-bottom:1px solid var(--line)}}
th{{background:var(--panel2);color:var(--muted);font-weight:600;font-size:12px}}
td.good{{color:#8fe0ad}} td.warn{{color:#f0cd8a}} td.bad{{color:#f0a08f}} td.na,.muted{{color:var(--muted)}}
td.gene{{font-weight:700}} tr.nogold td{{opacity:.55}}
.num{{text-align:right;font-variant-numeric:tabular-nums}}
.hbar{{background:var(--track);border-radius:4px;height:10px;width:100%}}
.hbar>div{{background:var(--uv);height:10px;border-radius:4px}}
tr.lever td:first-child{{border-left:3px solid var(--bad)}}
tr.nonlever{{opacity:.7}} .gloss{{color:var(--muted);font-size:12px}}
.card{{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:16px 18px;margin-top:12px}}
.card.placeholder{{border-style:dashed}}
.strategy ul{{margin:4px 0 10px;padding-left:18px}} .strategy li{{margin:3px 0}}
.notes{{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:10px}}
.note{{background:var(--panel2);border:1px solid var(--line);border-radius:8px;padding:10px}}
.note-h{{font-weight:700;font-size:12px;color:var(--target);margin-bottom:4px}}
.note-b{{font-size:12.5px;white-space:pre-wrap}}
.trend .grid{{stroke:var(--line)}} .trend .axis{{fill:var(--muted);font-size:9px}}
.trend .targetline{{stroke:var(--target);stroke-dasharray:4,3;stroke-width:1.5}}
.trend .line-uv{{stroke:var(--uv)}} .trend .line-uv circle{{fill:var(--uv)}}
.trend .line-pm{{stroke:var(--pm)}} .line-uv{{stroke:var(--uv)}} .line-pm{{stroke:var(--pm)}}
circle.line-uv{{fill:var(--uv)}} circle.line-pm{{fill:var(--pm)}}
.legend{{color:var(--muted);font-size:12px;margin-top:6px}}
.legend .k{{display:inline-block;width:10px;height:10px;border-radius:2px;margin:0 4px 0 12px;vertical-align:middle}}
.legend .k.uv{{background:var(--uv)}} .legend .k.pm{{background:var(--pm)}} .legend .k.tg{{background:var(--target)}}
.small{{font-size:12px}} code{{background:var(--panel2);padding:1px 5px;border-radius:4px}}
.foot{{color:var(--muted);font-size:11.5px;margin-top:34px;border-top:1px solid var(--line);padding-top:12px}}
</style></head>
<body><div class="wrap">
<h1>GeneVariantFetcher — Recall Status</h1>
<div class="sub">Target {int(TARGET * 100)}% across every metric · scored run generated {esc(gen_utc)} · page built {esc(now)}</div>
<div class="sub">Source: <code>{esc(str(metrics_dir.relative_to(repo_root)) if metrics_dir.is_relative_to(repo_root) else metrics_dir)}</code></div>

<h2>Aggregate recall — gap to 90%</h2>
<div class="grid-sc">{"".join(score_cards)}</div>
<div class="grid-mae">{"".join(mae_cards)}</div>

<h2>Trajectory</h2>
<div class="card">{svg_trend(trend)}</div>

<h2>Per-gene recall</h2>
<div class="sub">Bars show recall; the dashed gold line marks the 90% target. Genes without a gold standard are extraction-only.</div>
<table><thead><tr><th>Gene</th>{head}<th>carriers MAE</th></tr></thead><tbody>{"".join(rows)}</tbody></table>

<h2>Precision</h2>
<div class="sub">Headline precision (counted extras on gold PMIDs): <b>{pct(headline_prec)}</b>. Per source layer:</div>
<table><thead><tr><th>Source layer</th><th>Matched DB rows</th><th>Extra rows</th><th>Counted extra</th><th>Precision (counted)</th></tr></thead>
<tbody>{"".join(layer_rows) or "<tr><td colspan=5 class=muted>No per-layer precision in this run.</td></tr>"}</tbody></table>

<h2>Barriers — failure-mode split</h2>
<div class="sub">Red-marked classes are genuine recall losses (the levers). Others are precision/labeling residuals.</div>
<table><thead><tr><th>Failure class</th><th class="num">PMIDs</th><th>Share</th><th>What it means</th></tr></thead>
<tbody>{"".join(bar_rows) or "<tr><td colspan=4 class=muted>No paper_disagreement_report.csv in this run.</td></tr>"}</tbody></table>

<h2>Strategy</h2>
{render_strategy(strategy)}

<div class="foot">Generated by <code>scripts/build_status_dashboard.py</code>. Metrics are read-only from the recall suite;
the strategy narrative is written by the scheduled overnight run into <code>docs/dashboard/strategy.json</code>.
Authoritative numbers live in <code>docs/RECALL_STATUS.md</code>.</div>
</div></body></html>
"""


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--metrics-dir",
        default=None,
        help="recall_metrics/<run> dir or a summary.json path; default = newest valid",
    )
    ap.add_argument("--repo-root", default=".", help="repo root (default .)")
    ap.add_argument(
        "--out",
        default="docs/dashboard/index.html",
        help="output HTML path (default docs/dashboard/index.html)",
    )
    args = ap.parse_args(argv)

    repo_root = Path(args.repo_root).resolve()
    summary, metrics_dir = find_metrics(repo_root, args.metrics_dir)
    dashboard_dir = (repo_root / args.out).parent
    dashboard_dir.mkdir(parents=True, exist_ok=True)
    strategy = load_strategy(dashboard_dir)

    html_out = render(summary, metrics_dir, repo_root, strategy)
    out_path = repo_root / args.out
    out_path.write_text(html_out, encoding="utf-8")
    print(
        f"wrote {out_path}  (metrics: {metrics_dir}, strategy: {'yes' if strategy else 'none'})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
