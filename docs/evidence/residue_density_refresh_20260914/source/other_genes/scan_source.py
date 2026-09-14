#!/usr/bin/env python3
"""Bounded source-context sanity screen; no automatic count corrections."""

from __future__ import annotations

import hashlib
import json
import re
import sqlite3
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
PRIOR = REPO / "docs/evidence/brca2_distance_prior_audit_20260913/source"
ORIGINAL = (
    REPO.parent / "BayesianPenetranceEstimator/iterations/grant_e2e_20260909/results"
)
CLASS = (
    REPO
    / "docs/evidence/class_matched_penetrance_20260912/analysis/empirical_posteriors"
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(frame, name):
    kwargs = dict(index=False, lineterminator="\n")
    if name.endswith(".gz"):
        kwargs["compression"] = dict(method="gzip", mtime=0)
    frame.to_csv(HERE / name, **kwargs)
    assert (HERE / name).stat().st_size < 1_200_000


def main():
    rows, coverages, hashes = [], [], {}
    manifest = pd.read_csv(PRIOR / "cross_gene_db_provenance.csv")
    for gene in ["HNF1A", "LDLR", "KCNQ1"]:
        obs_path = ORIGINAL / f"{gene}_protocol/observations.csv"
        obs = pd.read_csv(obs_path, dtype={"pmid": str})
        hashes[str(obs_path)] = sha(obs_path)
        top_pmids = obs.groupby("pmid").affected.sum().nlargest(3).index
        class_paths = sorted(CLASS.glob(f"{gene}.part*.csv.gz"))
        model = pd.concat([pd.read_csv(p) for p in class_paths])
        keys = set(
            model.loc[model.variant_type.eq("missense"), "literature_key"].dropna()
        )
        gene_manifest = manifest.loc[manifest.gene.eq(gene)]
        for source in gene_manifest.itertuples():
            path = Path(source.path)
            assert sha(path) == source.sha256
            hashes[str(path)] = source.sha256
            with sqlite3.connect(f"file:{path}?mode=ro", uri=True) as db:
                sub = obs.loc[
                    (obs.variant_id.astype(int) // 10_000_000).eq(source.db_index)
                ]
                for r in sub.itertuples():
                    local_id = int(r.variant_id) % 10_000_000
                    paper = db.execute(
                        "SELECT title FROM papers WHERE CAST(pmid AS TEXT)=?", (r.pmid,)
                    ).fetchone()
                    title = paper[0] if paper else ""
                    result = db.execute(
                        "SELECT source_location,key_quotes,count_provenance FROM variant_papers WHERE variant_id=? AND CAST(pmid AS TEXT)=?",
                        (local_id, r.pmid),
                    ).fetchall()
                    assert result, (gene, r.key, r.pmid)
                    quotes = "\n".join(str(x[1] or "") for x in result)
                    location = "\n".join(str(x[0] or "") for x in result)
                    provenance = "\n".join(str(x[2] or "") for x in result)
                    source_text = title + "\n" + location + "\n" + quotes
                    patterns = {
                        "HNF1A": r"adenoma|adenomatosis|hepatocellular|somatic|renal cell carcinoma",
                        "LDLR": r"somatic|tumor|tumour|carcinoma|cancer",
                        "KCNQ1": r"atrial fibrillation|short[ -]?qt|Jervell|deafness",
                    }
                    hits = sorted(
                        set(re.findall(patterns[gene], source_text, flags=re.I))
                    )
                    catalog = bool(
                        re.search(r"VariationID|Canonical SPDI", source_text)
                    )
                    rows.append(
                        dict(
                            gene=gene,
                            key=r.key,
                            pmid=r.pmid,
                            variant_id=int(r.variant_id),
                            vclass=r.vclass,
                            affected=r.affected,
                            unaffected=r.unaffected,
                            canonical_missense_key=r.key in keys,
                            top_three_affected_source=r.pmid in top_pmids,
                            keyword_flag=bool(hits),
                            keyword_hits=";".join(hits),
                            catalog_header=catalog,
                            title=title,
                            source_location=location,
                            key_quotes=quotes,
                            count_provenance=provenance,
                            db_path=str(path),
                            source_db_variant_id=local_id,
                            primary_url=f"https://pubmed.ncbi.nlm.nih.gov/{r.pmid}/"
                            if r.pmid.isdigit()
                            else "",
                        )
                    )
        frame = pd.DataFrame([r for r in rows if r["gene"] == gene])
        for scope, subset in [
            ("all_frozen_observations", frame),
            ("canonical_missense_keys", frame.loc[frame.canonical_missense_key]),
        ]:
            for coverage, chosen in [
                ("all", subset),
                ("top_three_A_sources", subset.loc[subset.top_three_affected_source]),
                (
                    "keyword_flagged",
                    subset.loc[subset.keyword_flag | subset.catalog_header],
                ),
                (
                    "top_three_or_keyword_flagged",
                    subset.loc[
                        subset.top_three_affected_source
                        | subset.keyword_flag
                        | subset.catalog_header
                    ],
                ),
            ]:
                coverages.append(
                    dict(
                        gene=gene,
                        scope=scope,
                        coverage=coverage,
                        rows=len(chosen),
                        A=chosen.affected.sum(),
                        U=chosen.unaffected.sum(),
                        papers=chosen.pmid.nunique(),
                    )
                )
    all_rows = pd.DataFrame(rows)
    save(all_rows, "source_metadata_screen.csv.gz")
    save(
        all_rows.loc[
            all_rows.top_three_affected_source
            | all_rows.keyword_flag
            | all_rows.catalog_header
        ],
        "reviewed_source_scope.csv.gz",
    )
    save(
        all_rows.loc[all_rows.keyword_flag | all_rows.catalog_header],
        "keyword_flagged_rows.csv.gz",
    )
    save(pd.DataFrame(coverages), "coverage.csv")
    by_paper = all_rows.groupby(["gene", "pmid"], as_index=False).agg(
        rows=("key", "size"),
        A=("affected", "sum"),
        U=("unaffected", "sum"),
        top_three=("top_three_affected_source", "max"),
        keyword_rows=("keyword_flag", "sum"),
        title=("title", "first"),
    )
    save(
        by_paper.loc[by_paper.top_three | by_paper.keyword_rows.gt(0)],
        "papers_in_scope.csv",
    )
    (HERE / "screen_checks.json").write_text(
        json.dumps(
            dict(
                passed=True,
                original_DB_hashes_match=True,
                hashes=hashes,
                method="Top three PMIDs by frozen affected count per gene plus explicit source-metadata warning keywords; all variant types and canonical-missense-key coverage shown separately.",
                limitation="Screening metadata and top-source context is not source re-extraction or validation of every count; sparse quotes may miss problems. Keyword matches are warnings, not automatic exclusions. Counts and plot inputs unchanged.",
            ),
            indent=2,
        )
        + "\n"
    )
    print(
        by_paper.loc[by_paper.top_three | by_paper.keyword_rows.gt(0)].to_string(
            index=False
        )
    )


if __name__ == "__main__":
    main()
