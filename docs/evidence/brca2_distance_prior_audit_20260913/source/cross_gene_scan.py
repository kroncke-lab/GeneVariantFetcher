#!/usr/bin/env python3
"""Read-only, bounded search for the proven catalog/implicit-carrier signature.

This is a signature screen of selected frozen observations, not clinical validation.
The known BRCA2 block is a positive control. Empty quotes do not rule out a catalog.
Requires the original frozen local DBs; does not query live services or mutate them.
"""

from __future__ import annotations

import hashlib
import json
import sqlite3
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
SOURCE = (
    REPO.parent / "BayesianPenetranceEstimator/iterations/grant_e2e_20260909/results"
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(frame, filename):
    kwargs = dict(index=False, lineterminator="\n")
    if filename.endswith(".gz"):
        kwargs["compression"] = dict(method="gzip", mtime=0)
    frame.to_csv(HERE / filename, **kwargs)
    assert (HERE / filename).stat().st_size < 1_200_000


def main():
    all_rows, db_manifest = [], []
    for gene in ["HNF1A", "GCK", "LDLR", "KCNQ1", "BRCA2"]:
        folder = SOURCE / f"{gene}_protocol"
        obs = pd.read_csv(folder / "observations.csv", dtype={"pmid": str})
        summary = json.loads((folder / "observation_summary.json").read_text())
        paths = summary["gvf_db"]
        expected = summary["gvf_db_sha256"]
        if isinstance(paths, str):
            paths, expected = [paths], [expected]
        obs["db_index"] = obs.variant_id.astype(int) // 10_000_000
        for db_index, path in enumerate(paths):
            actual_hash = sha(path)
            assert actual_hash == expected[db_index], (gene, path)
            db_manifest.append(
                dict(gene=gene, db_index=db_index, path=path, sha256=actual_hash)
            )
            with sqlite3.connect(f"file:{path}?mode=ro", uri=True) as db:
                for row in obs.loc[obs.db_index.eq(db_index)].itertuples():
                    local_id = int(row.variant_id) % 10_000_000
                    records = db.execute(
                        "SELECT source_location,key_quotes,count_provenance,source_layer "
                        "FROM variant_papers WHERE variant_id=? AND CAST(pmid AS TEXT)=?",
                        (local_id, row.pmid),
                    ).fetchall()
                    quotes = "\n".join(str(r[1] or "") for r in records)
                    provenance = "\n".join(str(r[2] or "") for r in records)
                    catalog_header = (
                        "VariationID" in quotes or "Canonical SPDI" in quotes
                    )
                    implicit = "implicit one carrier per clinical row" in provenance
                    all_rows.append(
                        dict(
                            gene=gene,
                            key=row.key,
                            pmid=row.pmid,
                            variant_id=int(row.variant_id),
                            db_index=db_index,
                            db_variant_id=local_id,
                            affected=row.affected,
                            unaffected=row.unaffected,
                            source_records_found=len(records),
                            nonempty_source_quotes=bool(quotes.strip("\n []")),
                            catalog_header=catalog_header,
                            implicit_carrier=implicit,
                            exact_signature=catalog_header and implicit,
                            source_location="; ".join(str(r[0] or "") for r in records),
                            source_layer="; ".join(str(r[3] or "") for r in records),
                        )
                    )
    rows = pd.DataFrame(all_rows)
    rows["target_pmid"] = rows.pmid.eq("40664060")
    summary = rows.groupby("gene", sort=False).agg(
        frozen_observations=("key", "size"),
        target_pmid_observations=("target_pmid", "sum"),
        implicit_carrier_rows=("implicit_carrier", "sum"),
        exact_catalog_signature_rows=("exact_signature", "sum"),
        rows_with_nonempty_source_quotes=("nonempty_source_quotes", "sum"),
        rows_with_source_record=("source_records_found", lambda s: int(s.gt(0).sum())),
    )
    assert summary.loc["BRCA2", "exact_catalog_signature_rows"] == 4070
    assert summary.drop(index="BRCA2").target_pmid_observations.eq(0).all()
    save(rows, "cross_gene_selected_source_signature.csv.gz")
    save(summary.reset_index(), "cross_gene_signature_summary.csv")
    save(pd.DataFrame(db_manifest), "cross_gene_db_provenance.csv")
    checks = dict(
        passed=True,
        scope="Selected frozen observation source metadata; exact known catalog header plus implicit carrier inference only.",
        limitation="A negative signature is not validation. Empty/truncated quotes can hide catalog context; implicit clinical rows can also be legitimate patients. No source counts changed by this screen.",
        original_observation_hashes={
            gene: sha(SOURCE / f"{gene}_protocol/observations.csv")
            for gene in summary.index
        },
        all_db_hashes_match_frozen=True,
    )
    (HERE / "cross_gene_checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(summary.to_string())


if __name__ == "__main__":
    main()
