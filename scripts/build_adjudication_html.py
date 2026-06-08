#!/usr/bin/env python3
"""Build a self-contained HTML review packet for adjudicating a cold-start run.

Produces a folder (zipped) a reviewer opens locally — no server, no accounts,
no spreadsheet. ``index.html`` lists every sampled paper with the AI's extracted
variants and the provenance behind each (counts, count_provenance, source
location, verbatim quote, model notes, germline/somatic flag), a link to the
paper on PubMed AND to its full text, and per-variant inputs to record a verdict
(correct / wrong_count / not_in_paper / unsure) + corrected counts, plus an
"add a variant the AI missed" control. Progress autosaves to the browser; an
Export button downloads the reviewer's answers as a CSV that
``score_curation_packet.py`` consumes.

Pairs with ``build_adjudication_sheet.py`` (same data) but is the reviewer-
facing UI instead of a CSV.

Usage:
    python scripts/build_adjudication_html.py \
        --run-dir results/BRCA2/<ts>/BRCA2/<ts> \
        --manifest curation_packets/BRCA2_gold_50/manifest.csv \
        --gene BRCA2 --out curation_packets/BRCA2_review
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.build_adjudication_sheet import _load, build_rows  # noqa: E402

PAPER_TEMPLATE = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>{pmid} — full text</title>
<style>body{{font:15px/1.55 -apple-system,Segoe UI,Roboto,sans-serif;max-width:860px;
margin:24px auto;padding:0 18px;color:#1a1a1a}}h1{{font-size:18px}}
a{{color:#1558d6}}pre{{white-space:pre-wrap;word-wrap:break-word;font:13px/1.5
ui-monospace,Menlo,monospace;background:#f6f8fa;border:1px solid #e1e4e8;
border-radius:8px;padding:16px}}</style></head><body>
<h1>PMID {pmid} — full text (as the pipeline ingested it)</h1>
<p><a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" target="_blank">Open on PubMed ↗</a>
&nbsp;|&nbsp; <a href="../index.html">← back to review</a></p>
<pre>{body}</pre></body></html>"""


def _full_text(run_dir: Path, pmid: str) -> str:
    p = run_dir / "pmc_fulltext" / f"{pmid}_FULL_CONTEXT.md"
    if p.exists():
        try:
            return p.read_text(encoding="utf-8")
        except OSError:
            return ""
    return ""


def build_data(records: list[dict]) -> list[dict]:
    """Group the flat adjudication rows into per-paper structures for the UI."""
    rows = build_rows(records)
    by_pmid: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_pmid[r["pmid"]].append(r)
    papers = []
    for rec in records:
        pmid = rec["pmid"]
        rrows = by_pmid.get(pmid, [])
        variants = [
            {
                "variant": r.get("pipeline_variant", ""),
                "clin": r.get("clinical_significance", ""),
                "carriers": r.get("carriers", ""),
                "affected": r.get("affected", ""),
                "unaffected": r.get("unaffected", ""),
                "gs": r.get("germline_or_somatic", ""),
                "loc": r.get("source_location", ""),
                "quote": r.get("key_quote", ""),
                "prov": r.get("count_provenance", ""),
                "note": r.get("pipeline_notes", ""),
            }
            for r in rrows
            if r.get("pipeline_variant") != "(none extracted)"
        ]
        papers.append(
            {
                "pmid": pmid,
                "title": rec.get("title", ""),
                "pubmed": rec.get("url", f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"),
                "fulltext": f"papers/{pmid}.html",
                "variants": variants,
            }
        )
    return papers


def render_index(papers: list[dict], gene: str) -> str:
    data_json = json.dumps({"gene": gene, "papers": papers}).replace("</", "<\\/")
    return _INDEX_TEMPLATE.replace("__GENE__", html.escape(gene)).replace(
        "__DATA__", data_json
    )


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--no-zip", action="store_true")
    args = ap.parse_args()

    run_dir = args.run_dir.expanduser().resolve()
    gene = args.gene.upper()
    records = _load(run_dir, args.manifest.expanduser())
    papers = build_data(records)

    out = args.out.expanduser().resolve()
    if out.exists():
        shutil.rmtree(out)
    (out / "papers").mkdir(parents=True)
    (out / "index.html").write_text(render_index(papers, gene), encoding="utf-8")
    for rec in records:
        pmid = rec["pmid"]
        body = html.escape(_full_text(run_dir, pmid) or "(full text not on disk)")
        (out / "papers" / f"{pmid}.html").write_text(
            PAPER_TEMPLATE.format(pmid=html.escape(pmid), body=body), encoding="utf-8"
        )

    n_var = sum(len(p["variants"]) for p in papers)
    archive = ""
    if not args.no_zip:
        archive = shutil.make_archive(str(out), "zip", out.parent, out.name)
    print(f"HTML review packet: {out}")
    print(f"  {len(papers)} papers, {n_var} AI-extracted variants, full text per paper")
    if archive:
        print(f"  zip -> {archive}")
    print("  Reviewer: open index.html, record verdicts, click Export.")
    return 0


# --- the single-file review UI (static HTML + JS; data injected at __DATA__) ---
_INDEX_TEMPLATE = r"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__GENE__ — verify the AI's variant extractions</title>
<style>
:root{--bd:#e1e4e8;--mut:#586069;--ok:#e6f4ea;--bad:#fce8e6;--warn:#fff4e5}
*{box-sizing:border-box}body{font:15px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;
margin:0;color:#1a1a1a;background:#fafbfc}
header{position:sticky;top:0;background:#fff;border-bottom:1px solid var(--bd);
padding:12px 20px;z-index:5;display:flex;gap:16px;align-items:center;flex-wrap:wrap}
header h1{font-size:17px;margin:0}.grow{flex:1}
button{font:14px inherit;padding:8px 14px;border:1px solid var(--bd);border-radius:8px;
background:#1558d6;color:#fff;cursor:pointer}button.sec{background:#fff;color:#1558d6}
main{max-width:1100px;margin:0 auto;padding:18px}
.paper{background:#fff;border:1px solid var(--bd);border-radius:10px;margin:0 0 18px;
padding:14px 16px}.paper h2{font-size:15px;margin:0 0 4px}
.links a{color:#1558d6;text-decoration:none;margin-right:14px;font-size:13px}
.muted{color:var(--mut);font-size:12px}
table{width:100%;border-collapse:collapse;margin-top:10px;font-size:13px}
th,td{border:1px solid var(--bd);padding:6px 8px;vertical-align:top;text-align:left}
th{background:#f6f8fa;position:sticky}.prov{color:var(--mut);font-size:12px;max-width:340px}
.q{font-style:italic;color:#444}
select,input[type=text],input[type=number]{font:13px inherit;padding:4px;border:1px solid var(--bd);
border-radius:6px;width:100%}
input[type=number]{width:64px}
tr.correct{background:var(--ok)}tr.wrong{background:var(--warn)}tr.absent{background:var(--bad)}
.none{color:var(--mut);font-style:italic;padding:8px 0}
.addbtn{margin-top:8px;background:#fff;color:#1558d6}
.prog{font-size:13px;color:var(--mut)}
</style></head><body>
<header>
  <h1>__GENE__ — verify the AI's extractions</h1>
  <span class="prog" id="prog"></span>
  <span class="grow"></span>
  <button class="sec" onclick="toggleHelp()">How to use</button>
  <button onclick="exportCSV()">⬇ Export my answers (CSV)</button>
</header>
<main>
<div id="help" class="paper" style="display:none">
  <b>What to do.</b> For each paper, open its <b>full text</b> (or PubMed) and check the
  AI's extracted variants. For each row pick a <b>verdict</b>; if a count is wrong, type the
  correct number. If the AI <b>missed</b> a variant, click "Add a variant the AI missed".
  Counts: <b>carriers</b>=individuals with the variant, <b>affected</b>=of those, how many
  had cancer, <b>unaffected</b>=cancer-free. Leave blank if not reported. Work saves in this
  browser automatically; click <b>Export</b> when done and send the file back.
</div>
<div id="papers"></div>
</main>
<script>
const DATA = __DATA__;
const KEY = "adjud_" + DATA.gene;
const VERDICTS = ["", "correct", "wrong_count", "wrong_variant", "not_in_paper", "unsure"];

function load(){ try { return JSON.parse(localStorage.getItem(KEY)) || {}; } catch(e){ return {}; } }
function save(){
  const st = {};
  document.querySelectorAll("[data-k]").forEach(el => { if(el.value) st[el.dataset.k] = el.value; });
  // also persist dynamically-added "missed" rows' identity
  st.__missed = missed;
  localStorage.setItem(KEY, JSON.stringify(st));
  updateProg();
}
let store = load();
let missed = store.__missed || {};  // pmid -> count of extra rows

function cell(k, kind){
  const v = store[k] ? ` value="${store[k].replace(/"/g,'&quot;')}"` : "";
  if(kind==="verdict"){
    return `<select data-k="${k}" onchange="onVerdict(this)">` +
      VERDICTS.map(o=>`<option${store[k]===o?" selected":""}>${o}</option>`).join("") + `</select>`;
  }
  if(kind==="num") return `<input type="number" data-k="${k}"${v} oninput="save()">`;
  return `<input type="text" data-k="${k}"${v} oninput="save()">`;
}

function variantRow(pmid, i, vt){
  const base = `${pmid}|${i}`;
  return `<tr data-row="${base}">
   <td><b>${esc(vt.variant)}</b><div class="muted">${esc(vt.clin||"")} · ${esc(vt.gs||"")}</div></td>
   <td>${esc(vt.carriers)}</td><td>${esc(vt.affected)}</td><td>${esc(vt.unaffected)}</td>
   <td class="prov"><div>${esc(vt.loc||"")}</div>${vt.quote?`<div class="q">“${esc(vt.quote)}”</div>`:""}
     ${vt.prov?`<div class="muted">counts: ${esc(vt.prov)}</div>`:""}
     ${vt.note?`<div class="muted">${esc(vt.note)}</div>`:""}</td>
   <td>${cell(base+"|verdict","verdict")}</td>
   <td>${cell(base+"|c","num")}</td><td>${cell(base+"|a","num")}</td><td>${cell(base+"|u","num")}</td>
   <td>${cell(base+"|note","text")}</td></tr>`;
}

function missedRow(pmid, j){
  const base = `${pmid}|M${j}`;
  return `<tr data-row="${base}" class="wrong">
   <td>${cell(base+"|variant","text")}<div class="muted">missed by AI</div></td>
   <td colspan="3" class="muted">— AI did not extract this —</td>
   <td class="muted">enter the truth →</td>
   <td><input type="text" data-k="${base}|verdict" value="missed" readonly></td>
   <td>${cell(base+"|c","num")}</td><td>${cell(base+"|a","num")}</td><td>${cell(base+"|u","num")}</td>
   <td>${cell(base+"|note","text")}</td></tr>`;
}

function paperHTML(p){
  const head = `<h2>PMID <a href="${p.pubmed}" target="_blank">${esc(p.pmid)}</a> — ${esc(p.title||"")}</h2>
    <div class="links"><a href="${p.fulltext}" target="_blank">📄 Open full text ↗</a>
    <a href="${p.pubmed}" target="_blank">PubMed ↗</a>
    <span class="muted">${p.variants.length} AI-extracted variant(s)</span></div>`;
  let body;
  if(!p.variants.length){
    body = `<div class="none">The AI extracted no variants from this paper. If that's correct, leave it.
      If it missed some, add them below.</div>`;
  } else {
    body = `<table><thead><tr>
      <th>AI variant</th><th>carriers</th><th>affected</th><th>unaffected</th>
      <th>provenance (where / why)</th><th>verdict</th>
      <th>true carriers</th><th>true affected</th><th>true unaffected</th><th>notes</th>
      </tr></thead><tbody>` +
      p.variants.map((vt,i)=>variantRow(p.pmid,i,vt)).join("") + `</tbody></table>`;
  }
  const m = missed[p.pmid]||0;
  let extra = "";
  for(let j=0;j<m;j++) extra += `<table style="margin-top:6px"><tbody>${missedRow(p.pmid,j)}</tbody></table>`;
  return `<section class="paper" id="p_${p.pmid}">${head}${body}${extra}
    <button class="addbtn" onclick="addMissed('${p.pmid}')">+ Add a variant the AI missed</button></section>`;
}

function render(){ document.getElementById("papers").innerHTML = DATA.papers.map(paperHTML).join(""); updateProg(); }
function addMissed(pmid){ missed[pmid]=(missed[pmid]||0)+1; save(); render(); }
function onVerdict(sel){
  const tr = sel.closest("tr"); tr.className="";
  if(sel.value==="correct") tr.className="correct";
  else if(sel.value==="not_in_paper") tr.className="absent";
  else if(sel.value) tr.className="wrong";
  save();
}
function updateProg(){
  const tot = document.querySelectorAll('select[data-k$="|verdict"]').length;
  let done=0; document.querySelectorAll('select[data-k$="|verdict"]').forEach(s=>{if(s.value)done++;});
  document.getElementById("prog").textContent = `${done}/${tot} AI variants reviewed`;
}
function esc(s){ return String(s==null?"":s).replace(/[&<>]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;"}[c])); }
function toggleHelp(){ const h=document.getElementById("help"); h.style.display = h.style.display==="none"?"":"none"; }

function exportCSV(){
  const rows = [["pmid","variant","germline_or_somatic","carriers","affected","unaffected","verdict","reviewer_notes"]];
  const g = (k)=> (document.querySelector(`[data-k="${CSS.escape(k)}"]`)||{}).value || "";
  DATA.papers.forEach(p=>{
    p.variants.forEach((vt,i)=>{
      const b=`${p.pmid}|${i}`, verdict=g(b+"|verdict");
      if(!verdict) return; // only export reviewed rows
      rows.push([p.pmid, vt.variant, vt.gs, g(b+"|c")||vt.carriers, g(b+"|a")||vt.affected,
                 g(b+"|u")||vt.unaffected, verdict, g(b+"|note")]);
    });
    const m=missed[p.pmid]||0;
    for(let j=0;j<m;j++){ const b=`${p.pmid}|M${j}`, v=g(b+"|variant"); if(!v) continue;
      rows.push([p.pmid, v, "", g(b+"|c"), g(b+"|a"), g(b+"|u"), "missed", g(b+"|note")]); }
  });
  const csv = rows.map(r=>r.map(x=>{x=String(x==null?"":x); return /[",\n]/.test(x)?'"'+x.replace(/"/g,'""')+'"':x;}).join(",")).join("\n");
  const a=document.createElement("a");
  a.href=URL.createObjectURL(new Blob([csv],{type:"text/csv"}));
  a.download=DATA.gene+"_adjudication_FILLED.csv"; a.click();
}
render();
</script></body></html>"""


if __name__ == "__main__":
    raise SystemExit(main())
