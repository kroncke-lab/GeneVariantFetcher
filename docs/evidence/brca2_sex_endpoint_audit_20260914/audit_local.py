"""Audit BRCA2 sex metadata and endpoint scope without rewriting clinical counts."""

from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import re
import sqlite3

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
EVIDENCE = HERE.parent
OLD = EVIDENCE / "brca2_distance_prior_audit_20260913/source"
REFRESH = EVIDENCE / "residue_density_refresh_20260914"
HASHES = {}


def read(path):
    HASHES[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pd.read_csv(path, low_memory=False)


def save(df, name):
    blob = df.to_csv(index=False, lineterminator="\n").encode()
    if name.endswith(".gz"):
        blob = gzip.compress(blob, mtime=0)
    (HERE / name).write_bytes(blob)


def main():
    obs = read(OLD / "original_observations.csv.gz")
    original_count = len(obs)
    obs["original_row_index"] = np.arange(len(obs))
    duplicated_source_keys = int(obs.duplicated(["key", "pmid", "variant_id"]).sum())
    for path in [
        OLD / "BRCA2_catalog_quarantine_observations.csv",
        REFRESH / "source/BRCA2/BRCA2_observation_decisions.csv.gz",
    ]:
        decisions = read(path)
        for row in decisions.itertuples():
            mask = (
                obs.key.eq(row.original_key)
                & obs.pmid.eq(row.pmid)
                & obs.variant_id.eq(row.original_observation_variant_id)
            )
            assert mask.sum() == 1, (row.original_key, row.pmid)
            obs.loc[mask, "affected"] -= row.A_to_remove
            obs.loc[mask, "unaffected"] -= row.U_to_remove
    assert obs[["affected", "unaffected"]].ge(0).all().all()
    obs = obs.loc[(obs.affected + obs.unaffected).gt(0)].copy()
    restored = read(OLD / "BRCA2_clinical_restoration_observations.csv")
    obs = pd.concat(
        [
            obs,
            pd.DataFrame(
                {
                    "key": restored.protein_key,
                    "pmid": restored.pmid,
                    "variant_id": -np.arange(1, len(restored) + 1),
                    "affected": restored.A_to_restore,
                    "unaffected": restored.U_to_restore,
                    "vclass": "missense",
                }
            ),
        ],
        ignore_index=True,
    )
    clinical = read(REFRESH / "source/BRCA2/BRCA2_literature_identity_corrected.csv.gz")
    by_key = obs.groupby("key")[["affected", "unaffected"]].sum()
    np.testing.assert_allclose(
        by_key.reindex(clinical.key).fillna(0).to_numpy(),
        clinical[["affected_literature", "unaffected_literature"]].to_numpy(),
    )
    units = pd.concat(
        [
            read(p)
            for p in sorted(
                (REFRESH / "analysis/BRCA2").glob("missense_posteriors.part*.csv.gz")
            )
        ],
        ignore_index=True,
    )
    selected = obs.loc[obs.key.isin(units.literature_key.dropna())].copy()
    assert selected.affected.sum() == units.affected.sum()
    assert selected.unaffected.sum() == units.unaffected_literature.sum()
    dbs = read(OLD / "cross_gene_db_provenance.csv")
    dbs = dbs.loc[dbs.gene.eq("BRCA2")].set_index("db_index")
    conns = {}
    for index, row in dbs.iterrows():
        p = Path(row.path)
        assert hashlib.sha256(p.read_bytes()).hexdigest() == row.sha256
        HASHES[str(p)] = row.sha256
        conns[index] = sqlite3.connect(f"file:{p}?mode=ro", uri=True)
    records = []
    for row in selected.itertuples():
        data = dict(
            key=row.key,
            pmid=int(row.pmid),
            observation_variant_id=int(row.variant_id),
            affected=row.affected,
            unaffected=row.unaffected,
            original_row_index=row.original_row_index,
        )
        if row.variant_id < 0:
            data.update(
                title="Previously restored primary breast-cancer table",
                source_sex_status="needs_source_review",
                individual_record_count=0,
                individual_sex_counts="{}",
                title_male_or_prostate_signal=False,
            )
            records.append(data)
            continue
        index, local_id = divmod(int(row.variant_id), 10_000_000)
        con = conns[index]
        paper = con.execute(
            "SELECT title FROM papers WHERE pmid=?", (str(row.pmid),)
        ).fetchone()
        fields = con.execute(
            "SELECT sex,affected_status FROM individual_records WHERE variant_id=? AND pmid=?",
            (local_id, str(row.pmid)),
        ).fetchall()
        sexes = Counter(str(r[0] or "unknown") for r in fields)
        title = paper[0] if paper else ""
        if not title or title.startswith("Paper "):
            abstract = (
                Path(dbs.loc[index, "path"]).parent
                / "abstract_json"
                / f"{row.pmid}.json"
            )
            if abstract.exists():
                HASHES[str(abstract)] = hashlib.sha256(
                    abstract.read_bytes()
                ).hexdigest()
                title = (
                    (
                        json.loads(abstract.read_text())
                        .get("metadata", {})
                        .get("title", title)
                    )
                    or title
                    or ""
                )
        data.update(
            title=title,
            source_sex_status="structured_rows_present_not_reconciled_to_aggregate"
            if fields
            else "no_linked_individual_sex_records",
            individual_record_count=len(fields),
            individual_sex_counts=json.dumps(sexes, sort_keys=True),
            title_male_or_prostate_signal=bool(
                re.search(r"\b(male|men|prostat\w*)\b", title, flags=re.I)
            ),
            db_index=index,
            db_variant_id=local_id,
        )
        records.append(data)
    for con in conns.values():
        con.close()
    records = pd.DataFrame(records)
    save(records, "retained_missense_observation_scope.csv.gz")
    papers = (
        records.groupby("pmid")
        .agg(
            title=("title", "first"),
            observations=("key", "size"),
            affected=("affected", "sum"),
            unaffected=("unaffected", "sum"),
            records_with_individual_metadata=(
                "individual_record_count",
                lambda s: int(s.gt(0).sum()),
            ),
            title_male_or_prostate_signal=("title_male_or_prostate_signal", "max"),
        )
        .sort_values("affected", ascending=False)
        .reset_index()
    )
    save(papers, "retained_paper_inventory.csv")
    save(
        records.loc[records.title_male_or_prostate_signal.eq(True)],
        "male_endpoint_review_queue.csv",
    )
    columns = [
        "key",
        "affected",
        "unaffected_literature",
        "gnomad_carriers",
        "posterior_mean",
        "member_alleles",
    ]
    selected_keys = set(units.nlargest(10, "posterior_mean").key) | {
        "D2723H",
        "R3052W",
        "W2626C",
        "K2729N",
        "N372H",
        "I2675V",
    }
    save(
        units.loc[units.key.isin(selected_keys), columns], "selected_variant_counts.csv"
    )
    save(records.loc[records.key.isin(selected_keys)], "selected_variant_sources.csv")
    checks = dict(
        original_observations=original_count,
        duplicate_source_key_rows_preserved=duplicated_source_keys,
        corrected_key_totals_reproduced=True,
        final_missense_units=len(units),
        missense_literature_observations=len(records),
        missense_papers=len(papers),
        affected=float(records.affected.sum()),
        literature_unaffected=float(records.unaffected.sum()),
        gnomad_all_sex_unaffected=int(units.gnomad_carriers.sum()),
        observation_rows_with_linked_individual_records=int(
            records.individual_record_count.gt(0).sum()
        ),
        raw_all_affected_units=int(
            ((units.affected > 0) & (units.unaffected == 0)).sum()
        ),
        raw_all_affected_singletons=int(
            ((units.affected == 1) & (units.unaffected == 0)).sum()
        ),
        posterior_at_least_90_percent=int(units.posterior_mean.ge(0.9).sum()),
        sex_filter_present_in_frozen_count_model=False,
        metadata_screen_is_not_a_sex_specific_count_partition=True,
        verified_source_database_snapshots=len(dbs),
        title_male_or_prostate_signal_observations=int(
            records.title_male_or_prostate_signal.sum()
        ),
        title_male_or_prostate_signal_affected=float(
            records.loc[
                records.title_male_or_prostate_signal.eq(True), "affected"
            ].sum()
        ),
    )
    (HERE / "local_checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    HASHES[str(Path(__file__))] = hashlib.sha256(
        Path(__file__).read_bytes()
    ).hexdigest()
    (HERE / "input_hashes.json").write_text(json.dumps(HASHES, indent=2) + "\n")
    print(json.dumps(checks, indent=2))
    print(papers.head(8).to_string(index=False))


if __name__ == "__main__":
    main()
