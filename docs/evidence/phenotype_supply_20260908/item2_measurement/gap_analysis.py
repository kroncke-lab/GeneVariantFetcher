"""Where does the remaining positive-gold affected gap sit after the projection?"""

import copy
import json
import sys
import collections
import importlib.util
from pathlib import Path

REPO = Path("/Users/kronckbm/GitRepos/GeneVariantFetcher")
sys.path.insert(0, str(REPO))
spec = importlib.util.spec_from_file_location(
    "replay", REPO / "scripts/replay_table_cohort_phenotype.py"
)
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)
from benchmarks.codex_paper_eval.production_run import resolve_active_gene_run

run_eval = replay.load_run_eval()
cat = collections.Counter()
per_paper = collections.Counter()
captions = collections.Counter()
wrong = collections.Counter()
total_pos = 0
supplied_exact = 0
for run_id in (
    "20260905_protocol_cont120_02_candidate",
    "20260905_protocol_cont120_03_candidate",
):
    run_dir = REPO / "benchmarks/codex_paper_eval/runs" / run_id
    selection = json.loads((run_dir / "selection.json").read_text())
    predictions = json.loads((run_dir / "predictions.json").read_text())
    setup = json.loads((run_dir / "setup.json").read_text())
    gold_root = Path(setup["cohort"]["gold_root"])
    prod = run_dir / "production_runs"
    gene_runs = {}
    for gene in sorted({p["gene"] for p in selection["papers"]}):
        try:
            gene_runs[gene] = resolve_active_gene_run(prod, gene)[0]
        except Exception as e:
            print(gene, e, file=sys.stderr)
    papers = {
        (p["gene"], str(p["pmid"])): copy.deepcopy(p) for p in predictions["papers"]
    }
    for paper in selection["papers"]:
        gene, pmid = paper["gene"], str(paper["pmid"])
        key = (gene, pmid)
        if gene not in ("KCNH2", "KCNQ1", "SCN5A", "RYR2") or key not in papers:
            continue
        gr = gene_runs.get(gene)
        if gr is None:
            continue
        path = gr / "extractions" / f"{gene}_PMID_{pmid}.json"
        extraction = (
            json.loads(path.read_text()) if path.is_file() else {"variants": []}
        )
        src = replay.source_text_for(extraction, gr, pmid) if path.is_file() else ""
        records, meta = (
            replay.derived_records(
                extraction,
                src,
                gene=gene,
                disease=replay.disease_for(gr),
                mode="derive",
                paper_tier=False,
            )
            if path.is_file()
            else ([], {})
        )
        if records:
            replay.patch_paper(papers[key], records, "derive")
        gold = run_eval.load_gold(gold_root, gene, pmid)
        score = run_eval.score_one(gene, pmid, papers[key], gold)
        assigned = replay.scored_gold_rows(score, gold)
        # index extraction variants by tokens
        ext = [
            (replay.record_tokens(v), v)
            for v in extraction.get("variants", [])
            if isinstance(v, dict)
        ]
        outcomes = {
            str(o.get("variant")): o.get("status") for o in (meta.get("outcomes") or [])
        }
        pred_rows = {r.get("variant"): r for r in papers[key].get("variants", [])}
        for pred_name, grow in assigned.items():
            ga = grow.get("affected")
            if ga is None or int(ga) <= 0:
                continue
            total_pos += 1
            prow = pred_rows.get(pred_name) or {}
            if prow.get("affected") is not None:
                if int(prow["affected"]) == int(ga):
                    supplied_exact += 1
                else:
                    wrong[(gene, pmid)] += 1
                continue
            toks = replay.prediction_tokens(prow)
            match = next((v for t, v in ext if t & toks), None)
            if match is None:
                c = "no_extraction_row_match"
                cap = ""
            else:
                pat = match.get("patients") or {}
                lay = str(match.get("source_layer") or "")
                parser = str(
                    ((match.get("locator_extra") or {}).get("parser"))
                    or (pat.get("locator_extra") or {}).get("parser")
                    or ""
                )
                ident = next(
                    (
                        str(match.get(k))
                        for k in (
                            "source_notation",
                            "protein_notation",
                            "cdna_notation",
                        )
                        if match.get(k)
                    ),
                    "",
                )
                status = outcomes.get(ident, "?")
                if pat.get("column_ref") == "implicit one carrier per clinical row":
                    c = f"per_person_row[{parser}]"
                elif lay == "regex_table" or parser:
                    c = f"deterministic_table[{status}]"
                elif (match.get("penetrance_data") or {}).get(
                    "total_carriers_observed"
                ) is None:
                    c = f"model_row_no_carriers"
                else:
                    c = f"model_row_with_carriers[{status}]"
                cap = str(pat.get("source_ref") or match.get("source_table") or "")[:70]
            cat[c] += 1
            per_paper[(gene, pmid, c)] += 1
            if cap:
                captions[(gene, pmid, c, cap)] += 1
print(
    "positive-gold affected rows matched:",
    total_pos,
    "| supplied exact:",
    supplied_exact,
    "| wrong:",
    sum(wrong.values()),
    dict(wrong),
)
print("\nSTILL NULL by category:")
[print(f"{n:5} {c}") for c, n in cat.most_common()]
print("\nTop (paper,category):")
[print(f"{n:5} {k}") for k, n in per_paper.most_common(25)]
print("\nTop captions:")
[print(f"{n:5} {k}") for k, n in captions.most_common(25)]
