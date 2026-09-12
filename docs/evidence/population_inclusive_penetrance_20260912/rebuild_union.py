"""Build a population-inclusive, auditable variant union and empirical priors.

Population carriers are assumed unaffected. Actual observed, QC-passing alleles
enter independently of literature inclusion. Hypothetical unobserved variants
never become negative observations. All original clinical rows remain in a
separate join ledger, including identity quarantines.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

import numpy as np
import pandas as pd
from scipy.stats import beta as beta_dist

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / "structural_density_plan_20260912"
GENES = ["HNF1A", "GCK", "LDLR", "BRCA2", "KCNQ1"]
POINT_CLASSES = {"missense", "synonymous", "nonsense", "start_lost", "stop_lost"}


def save(frame, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(
        path,
        index=False,
        compression={"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None,
    )


def text(value):
    return "" if pd.isna(value) else str(value).strip()


def normalized_cdna(value):
    return unquote(text(value)).split(":")[-1].replace(" ", "")


def cdna_matches_transcript(value, transcript):
    """Bare cDNA inherits a verified row transcript; explicit prefixes must agree."""
    value = unquote(text(value)).replace(" ", "")
    transcript = text(transcript)
    if not transcript or not value:
        return False
    prefix, separator, change = value.rpartition(":")
    if separator and prefix != transcript:
        return False
    return (change if separator else value).startswith("c.")


def canonical_annotation_compatible(row):
    canonical = text(row.get("canonical_transcript_id"))
    return bool(canonical) and text(row.get("transcript_id")) == canonical


def validate_counts(frame, fields, label):
    for field in fields:
        values = pd.to_numeric(frame[field], errors="coerce").to_numpy(dtype=float)
        if not (
            np.isfinite(values).all()
            and (values >= 0).all()
            and np.equal(values, np.floor(values)).all()
        ):
            raise ValueError(f"{label} {field} must contain nonnegative integer counts")


def normalized_protein(value):
    """Match equivalent point annotations: R186*/R186X and G258G/G258=."""
    value = text(value).removeprefix("p.")
    match = re.fullmatch(r"([A-Z])(\d+)([A-Z*=])", value)
    if not match:
        return value
    ref, pos, alt = match.groups()
    if alt == ref:
        alt = "="
    if alt == "*":
        alt = "X"
    return f"{ref}{pos}{alt}"


def build_union(population, literature):
    """Assign each eligible population allele and clinical aggregate at most once.

    A protein-only clinical aggregate absorbs all population alleles producing
    that exact canonical substitution. Unmatched population alleles remain
    genomic units. Shared population matches across clinical keys are an
    unresolved alias component: clinical counts are quarantined, never summed.
    """
    if literature.duplicated(["gene", "key"]).any():
        raise ValueError("Clinical gene/key identities must be unique")
    validate_counts(
        literature, ["affected_literature", "unaffected_literature"], "clinical"
    )
    eligible_population = population.loc[population.population_eligible]
    validate_counts(
        eligible_population, ["gnomad_carriers", "joint_ac", "joint_hom"], "population"
    )
    if not (
        (eligible_population.joint_ac >= 2 * eligible_population.joint_hom).all()
        and (
            eligible_population.gnomad_carriers
            == eligible_population.joint_ac - eligible_population.joint_hom
        ).all()
    ):
        raise ValueError(
            "Population carrier counts must equal coherent joint AC minus homozygotes"
        )
    units, ledger, membership = [], [], []
    for gene in GENES:
        pop = population.loc[
            population.gene.eq(gene) & population.population_eligible
        ].copy()
        pop = pop.sort_values("variant_id").reset_index(drop=True)
        if pop.variant_id.duplicated().any():
            raise ValueError(f"Duplicate population allele in {gene}")
        canonical_transcripts = {
            text(value) for value in pop.canonical_transcript_id if text(value)
        }
        if len(canonical_transcripts) > 1:
            raise ValueError(f"Multiple canonical transcripts supplied for {gene}")
        by_protein, by_cdna = defaultdict(set), defaultdict(set)
        for i, row in pop.iterrows():
            if not canonical_annotation_compatible(row):
                # Preserve the observed genomic allele but do not promote an
                # alternate-transcript annotation into canonical protein evidence.
                if text(row.protein_key) or text(row.hgvsc) or text(row.aa_ref):
                    pop.loc[i, ["protein_key", "aa_ref", "aa_alt"]] = ""
                    pop.loc[i, "aa_pos"] = np.nan
                    pop.loc[i, "canonical_wt_status"] = "unavailable"
                    pop.loc[i, "coding_or_splice"] = False
                    pop.loc[i, "vclass"] = "unresolved_canonical_annotation"
                continue
            if text(row.protein_key) and row.canonical_wt_status == "match":
                by_protein[normalized_protein(row.protein_key)].add(i)
            if cdna_matches_transcript(row.hgvsc, row.canonical_transcript_id):
                by_cdna[normalized_cdna(row.hgvsc)].add(i)
        clinical = (
            literature.loc[literature.gene.eq(gene)].copy().reset_index(drop=True)
        )
        candidates, reasons, join_methods = {}, {}, {}
        for i, row in clinical.iterrows():
            reasons[i] = (
                text(row.identity_quarantine_reason)
                if row.identity_quarantine_recommended
                else ""
            )
            if reasons[i]:
                candidates[i], join_methods[i] = set(), "identity_quarantine"
                continue
            clinical_transcript = text(row.feature_transcript)
            if (
                canonical_transcripts
                and clinical_transcript not in canonical_transcripts
            ):
                candidates[i], join_methods[i] = set(), "identity_quarantine"
                reasons[i] = "clinical_feature_transcript_not_canonical"
                continue
            cdna_prefix, cdna_separator, _ = (
                unquote(text(row.cdna)).replace(" ", "").rpartition(":")
            )
            if cdna_separator and cdna_prefix != clinical_transcript:
                candidates[i], join_methods[i] = set(), "identity_quarantine"
                reasons[i] = "clinical_cdna_transcript_conflict"
                continue
            pkey = (
                normalized_protein(row.canonical_point_key)
                if row.protein_point_join_allowed
                else ""
            )
            pcandidates = by_protein[pkey] if pkey else set()
            ccandidates = (
                by_cdna[normalized_cdna(row.cdna)] if text(row.cdna) else set()
            )
            # A cDNA label cannot override contradictory canonical protein evidence.
            if pkey and ccandidates:
                mismatch = {
                    j
                    for j in ccandidates
                    if normalized_protein(pop.iloc[j].protein_key) != pkey
                }
                if mismatch:
                    candidates[i], join_methods[i] = set(), "identity_quarantine"
                    reasons[i] = "canonical_cdna_and_protein_conflict"
                    continue
            if pkey and pcandidates:
                candidates[i], join_methods[i] = (
                    set(pcandidates),
                    "canonical_protein_aggregate",
                )
            elif ccandidates:
                if len(ccandidates) > 1:
                    candidates[i], join_methods[i] = (
                        set(ccandidates),
                        "identity_quarantine",
                    )
                    reasons[i] = "canonical_cdna_matches_multiple_genomic_alleles"
                    continue
                candidates[i], join_methods[i] = (
                    set(ccandidates),
                    "exact_canonical_cdna",
                )
            else:
                candidates[i], join_methods[i] = set(), "literature_only_unmatched"
        owners = defaultdict(list)
        for i, indices in candidates.items():
            if reasons[i]:
                continue
            for j in indices:
                owners[j].append(i)
        ambiguous = {i for values in owners.values() if len(values) > 1 for i in values}
        for i in ambiguous:
            reasons[i] = "multiple_literature_keys_match_same_population_allele"
            join_methods[i] = "identity_quarantine"
        consumed = set()
        for i, row in clinical.iterrows():
            match = sorted(candidates[i])
            record = row.to_dict()
            record.update(
                {
                    "join_method": join_methods[i],
                    "join_exclusion_reason": reasons[i],
                    "prior_literature_eligible": not bool(reasons[i]),
                    "candidate_population_members": ";".join(
                        pop.iloc[match].variant_id
                    ),
                }
            )
            if reasons[i]:
                record["unit_id"] = ""
                ledger.append(record)
                continue
            members = pop.iloc[match]
            unit_id = f"{gene}|lit:{row.key}"
            record["unit_id"] = unit_id
            ledger.append(record)
            if consumed.intersection(match):
                raise ValueError(
                    "A population allele was attached to multiple clinical units"
                )
            consumed.update(match)
            molecular = normalized_protein(row.canonical_point_key)
            resolved = members.iloc[0] if len(members) else None
            vclass = row.vclass
            aa_pos, aa_ref, aa_alt = row.aa_pos, row.aa_ref, row.aa_alt
            if not molecular and resolved is not None and text(resolved.protein_key):
                molecular = normalized_protein(resolved.protein_key)
                vclass = resolved.vclass
                aa_pos, aa_ref, aa_alt = (
                    resolved.aa_pos,
                    resolved.aa_ref,
                    resolved.aa_alt,
                )
            coding = (
                bool(members.coding_or_splice.any())
                if len(members)
                else vclass not in {"intron", "utr", "noncoding"}
            )
            units.append(
                {
                    "gene": gene,
                    "unit_id": unit_id,
                    "key": row.key,
                    "literature_key": row.key,
                    "origin": "literature_and_population"
                    if len(members)
                    else "literature_only",
                    "unit_grain": "clinical_protein_aggregate"
                    if len(members) > 1
                    else "clinical_variant_aggregate",
                    "member_alleles": ";".join(members.variant_id),
                    "n_member_alleles": len(members),
                    "protein_key": molecular,
                    "vclass": vclass,
                    "aa_pos": aa_pos,
                    "aa_ref": aa_ref,
                    "aa_alt": aa_alt,
                    "canonical_wt_status": row.canonical_wt_status
                    if text(row.aa_ref)
                    else (
                        resolved.canonical_wt_status
                        if resolved is not None
                        else "unavailable"
                    ),
                    "coding_or_splice": coding,
                    "affected": row.affected_literature,
                    "unaffected_literature": row.unaffected_literature,
                    "gnomad_carriers": members.gnomad_carriers.sum(),
                    "gnomad_ac": members.joint_ac.sum(),
                    "gnomad_homozygotes": members.joint_hom.sum(),
                    "join_method": join_methods[i],
                }
            )
            for allele in members.variant_id:
                membership.append(
                    {"gene": gene, "variant_id": allele, "unit_id": unit_id}
                )
        for j, row in pop.iterrows():
            if j in consumed:
                continue
            unit_id = f"{gene}|g:{row.variant_id}"
            units.append(
                {
                    "gene": gene,
                    "unit_id": unit_id,
                    "key": row.variant_id,
                    "literature_key": "",
                    "origin": "population_only",
                    "unit_grain": "genomic_allele",
                    "member_alleles": row.variant_id,
                    "n_member_alleles": 1,
                    "protein_key": normalized_protein(row.protein_key),
                    "vclass": row.vclass,
                    "aa_pos": row.aa_pos,
                    "aa_ref": row.aa_ref,
                    "aa_alt": row.aa_alt,
                    "canonical_wt_status": row.canonical_wt_status,
                    "coding_or_splice": bool(row.coding_or_splice),
                    "affected": 0,
                    "unaffected_literature": 0,
                    "gnomad_carriers": row.gnomad_carriers,
                    "gnomad_ac": row.joint_ac,
                    "gnomad_homozygotes": row.joint_hom,
                    "join_method": "no_unambiguous_clinical_match",
                }
            )
            membership.append(
                {"gene": gene, "variant_id": row.variant_id, "unit_id": unit_id}
            )
    units = pd.DataFrame(units)
    ledger = pd.DataFrame(
        ledger,
        columns=[
            *literature.columns,
            "join_method",
            "join_exclusion_reason",
            "prior_literature_eligible",
            "candidate_population_members",
            "unit_id",
        ],
    )
    membership = pd.DataFrame(membership, columns=["gene", "variant_id", "unit_id"])
    if units.empty:
        raise ValueError("No eligible clinical or observed population units")
    assert units.unit_id.is_unique
    assert not membership.duplicated(["gene", "variant_id"]).any()
    assert len(membership) == int(population.population_eligible.sum())
    units["unaffected"] = units.unaffected_literature + units.gnomad_carriers
    units["n"] = units.affected + units.unaffected
    assert (units.n > 0).all()
    for field in ["affected", "unaffected", "n", "gnomad_carriers"]:
        assert np.isfinite(units[field]).all() and (units[field] >= 0).all()
        assert np.equal(units[field], np.floor(units[field])).all()
    for gene in GENES:
        pop = population.loc[population.gene.eq(gene) & population.population_eligible]
        u = units.loc[units.gene.eq(gene)]
        clinical = ledger.loc[ledger.gene.eq(gene) & ledger.prior_literature_eligible]
        assert u.gnomad_carriers.sum() == pop.gnomad_carriers.sum()
        assert u.affected.sum() == clinical.affected_literature.sum()
        assert u.unaffected_literature.sum() == clinical.unaffected_literature.sum()
    return units, ledger, membership


def fit_empirical(units, variance_mode="historical"):
    y = units.affected.to_numpy() / units.n.to_numpy()
    weights = 1 - 1 / (units.n.to_numpy() + 0.01)
    mean = float(np.average(y, weights=weights))
    mse = float(np.mean(weights * (y - mean) ** 2))
    normalized_mse = float(np.average((y - mean) ** 2, weights=weights))
    variance = mse if variance_mode == "historical" else normalized_mse
    if not 0 < mean < 1 or not 0 < variance < mean * (1 - mean):
        raise ValueError(f"Invalid empirical Beta moments: {mean=}, {variance=}")
    strength = mean * (1 - mean) / variance - 1
    alpha, beta = mean * strength, (1 - mean) * strength
    return {
        "mean": mean,
        "variance": variance,
        "historical_mse": mse,
        "normalized_weighted_mse": normalized_mse,
        "weighted_mae": float(np.average(abs(y - mean), weights=weights)),
        "alpha_empirical": alpha,
        "beta_empirical": beta,
        "strength": strength,
        "pooled_carrier_fraction": float(units.affected.sum() / units.n.sum()),
        "equal_variant_fraction_mean": float(y.mean()),
    }


def posterior_table(units, parameters):
    out = units.copy()
    out["alpha_empirical"] = parameters["alpha_empirical"]
    out["beta_empirical"] = parameters["beta_empirical"]
    out["posterior_alpha"] = out.alpha_empirical + out.affected
    out["posterior_beta"] = out.beta_empirical + out.unaffected
    out["posterior_mean"] = out.posterior_alpha / (
        out.posterior_alpha + out.posterior_beta
    )
    out["posterior_variance"] = (
        out.posterior_mean
        * (1 - out.posterior_mean)
        / (out.posterior_alpha + out.posterior_beta + 1)
    )
    out["posterior_lower_95"] = beta_dist.ppf(
        0.025, out.posterior_alpha, out.posterior_beta
    )
    out["posterior_upper_95"] = beta_dist.ppf(
        0.975, out.posterior_alpha, out.posterior_beta
    )
    assert np.allclose(out.posterior_alpha - out.alpha_empirical, out.affected)
    assert np.allclose(
        out.posterior_beta - out.beta_empirical,
        out.unaffected_literature + out.gnomad_carriers,
    )
    return out


def load_population():
    """Exact full-span allele inventory plus original boundary-padding alleles.

    Region-only alleles keep missing canonical annotations. Overlapping counts
    must agree before transferring annotations from the canonical snapshot.
    """
    footprint_path = HERE / "population/population_variants.csv.gz"
    footprint = pd.read_csv(footprint_path)
    inputs = [footprint_path]
    tables = []
    for gene in GENES:
        manifest_path = HERE / f"population/full_locus/{gene}_provenance.json"
        manifest = json.loads(manifest_path.read_text())
        inputs.append(manifest_path)
        parts = []
        for source in manifest["output_files"]:
            path = manifest_path.parent / source["file"]
            if hashlib.sha256(path.read_bytes()).hexdigest() != source["sha256"]:
                raise ValueError(f"Population source hash mismatch: {path}")
            inputs.append(path)
            parts.append(pd.read_csv(path, low_memory=False))
        region = pd.concat(parts, ignore_index=True)
        old = footprint.loc[footprint.gene.eq(gene)].set_index("variant_id")
        region = region.set_index("variant_id")
        if not region.index.is_unique:
            raise ValueError(f"Duplicate full-locus allele: {gene}")
        interval = manifest["interval"]
        absent = old.loc[~old.index.isin(region.index)]
        if absent.pos.between(interval["start"], interval["stop"]).any():
            raise ValueError(
                f"Full-locus source missing an inside-span footprint allele: {gene}"
            )
        overlap = region.index.intersection(old.index)
        for column in [
            "joint_ac",
            "joint_an",
            "joint_hom",
            "gnomad_carriers",
            "qc_pass",
            "population_eligible",
            "exome_ac",
            "exome_an",
            "exome_hom",
            "exome_filters",
            "genome_ac",
            "genome_an",
            "genome_hom",
            "genome_filters",
            "joint_filters",
            "reconstructed_joint_filter",
        ]:
            left, right = region.loc[overlap, column], old.loc[overlap, column]
            if not (left.eq(right) | (left.isna() & right.isna())).all():
                raise ValueError(f"Population overlap changed: {gene} {column}")
        # Only annotate exact observed alleles, never a nearby or protein proxy.
        annotation_cols = [c for c in old if c not in region]
        region = region.join(old[annotation_cols], how="left")
        region["in_coding_footprint"] = region.index.isin(old.index)
        region["in_gene_span"] = True
        outside = old.loc[~old.index.isin(region.index)].copy()
        outside["in_coding_footprint"] = True
        outside["in_gene_span"] = False
        combined = pd.concat([region, outside]).reset_index()
        for column, default in {
            "vclass": "region_only_unannotated",
            "canonical_wt_status": "unavailable",
            "canonical_annotation_status": "outside_canonical_query_footprint",
            "coding_or_splice": False,
        }.items():
            combined[column] = combined[column].fillna(default)
        combined["coding_or_splice"] = combined.coding_or_splice.astype(bool)
        tables.append(combined)
    return pd.concat(tables, ignore_index=True), inputs


def save_shards(frame, directory, stem):
    """Bounded, deterministic files; keep the complete analysis reviewable in Git."""
    directory.mkdir(parents=True, exist_ok=True)
    # This script owns these exact generated shard names.
    for path in directory.glob(f"{stem}*.part*.csv.gz"):
        path.unlink()
    for gene, group in frame.groupby("gene", sort=False):
        for part, start in enumerate(range(0, len(group), 5000), 1):
            path = directory / f"{stem}{gene}.part{part:03d}.csv.gz"
            save(group.iloc[start : start + 5000], path)
            if path.stat().st_size > 1100000:
                raise ValueError(f"Reduce shard size: {path}")


def main():
    population, inputs = load_population()
    clinical_path = HERE / "audit/literature_identity_flags.csv.gz"
    clinical = pd.read_csv(clinical_path)
    if len(clinical) != 8317 or clinical.duplicated(["gene", "key"]).any():
        raise ValueError("Frozen clinical universe changed")
    units, ledger, membership = build_union(population, clinical)
    footprint_ids = set(
        population.loc[population.in_coding_footprint, "gene"]
        + "|g:"
        + population.loc[population.in_coding_footprint, "variant_id"]
    )
    units["in_footprint_comparison"] = units.origin.ne(
        "population_only"
    ) | units.unit_id.isin(footprint_ids)
    out = HERE / "analysis"
    for name, frame in [("union_counts", units), ("population_membership", membership)]:
        save_shards(frame, out / name, "")
        obsolete = out / f"{name}.csv.gz"
        if obsolete.exists():
            obsolete.unlink()
    save(ledger, out / "literature_join_ledger.csv.gz")
    moments, posteriors = [], []
    for gene in GENES:
        frame = units.loc[units.gene.eq(gene)].copy()
        subsets = {
            "all_observed_full_locus": frame,
            "all_observed_footprint": frame.loc[frame.in_footprint_comparison],
            "coding_or_splice": frame.loc[frame.coding_or_splice],
            "literature_only_refreshed": frame.loc[frame.origin.ne("population_only")],
        }
        for scope, subset in subsets.items():
            parameters = fit_empirical(subset)
            moments.append(
                {
                    "gene": gene,
                    "scope": scope,
                    "variants": len(subset),
                    "population_only": int(subset.origin.eq("population_only").sum()),
                    "affected": int(subset.affected.sum()),
                    "unaffected_literature": int(subset.unaffected_literature.sum()),
                    "unaffected_gnomad": int(subset.gnomad_carriers.sum()),
                    **parameters,
                }
            )
            if scope == "all_observed_full_locus":
                posteriors.append(posterior_table(frame, parameters))
        ac_frame = frame.copy()
        ac_frame["unaffected"] = ac_frame.unaffected_literature + ac_frame.gnomad_ac
        ac_frame["n"] = ac_frame.affected + ac_frame.unaffected
        # The primary keeps distinct genomic alleles unless a clinical aggregate
        # cannot resolve its constituent DNA alleles. Quantify that mixed grain.
        grain = frame.copy()
        grain["aggregate_key"] = np.where(
            grain.protein_key.fillna("").ne("") & grain.canonical_wt_status.eq("match"),
            "p:" + grain.protein_key.fillna(""),
            grain.unit_id,
        )
        grain = grain.groupby("aggregate_key", sort=False).agg(
            affected=("affected", "sum"),
            unaffected_literature=("unaffected_literature", "sum"),
            gnomad_carriers=("gnomad_carriers", "sum"),
            n=("n", "sum"),
            origin=(
                "origin",
                lambda x: (
                    "population_only"
                    if x.eq("population_only").all()
                    else "contains_literature"
                ),
            ),
        )
        grain["unaffected"] = grain.unaffected_literature + grain.gnomad_carriers
        for label, data, mode in [
            ("all_allele_count_proxy", ac_frame, "historical"),
            ("all_normalized_mse", frame, "normalized"),
            ("all_protein_aggregate_grain", grain, "historical"),
        ]:
            moments.append(
                {
                    "gene": gene,
                    "scope": label,
                    "variants": len(data),
                    "population_only": int(data.origin.eq("population_only").sum()),
                    "affected": int(data.affected.sum()),
                    "unaffected_literature": int(data.unaffected_literature.sum()),
                    "unaffected_gnomad": int(
                        (data.unaffected - data.unaffected_literature).sum()
                    ),
                    **fit_empirical(data, mode),
                }
            )
    save(pd.DataFrame(moments), out / "empirical_prior_comparison.csv")
    save_shards(
        pd.concat(posteriors, ignore_index=True), out / "empirical_posteriors", ""
    )
    obsolete = out / "empirical_posteriors.csv.gz"
    if obsolete.exists():
        obsolete.unlink()
    summary = []
    for gene in GENES:
        frame = units.loc[units.gene.eq(gene)]
        lit = ledger.loc[ledger.gene.eq(gene)]
        raw = population.loc[population.gene.eq(gene)]
        summary.append(
            {
                "gene": gene,
                "source_population_alleles": len(raw),
                "source_gene_span_alleles": int(raw.in_gene_span.sum()),
                "boundary_padding_alleles": int((~raw.in_gene_span).sum()),
                "eligible_population_alleles": int(raw.population_eligible.sum()),
                "clinical_input_keys": len(lit),
                "quarantined_clinical_keys": int(
                    (~lit.prior_literature_eligible).sum()
                ),
                "quarantined_affected_observations": int(
                    lit.loc[~lit.prior_literature_eligible].affected_literature.sum()
                ),
                "clinical_units": int(frame.origin.ne("population_only").sum()),
                "population_only_units": int(frame.origin.eq("population_only").sum()),
                "union_units": len(frame),
                "multi_allele_clinical_units": int((frame.n_member_alleles > 1).sum()),
                "affected": int(frame.affected.sum()),
                "unaffected_literature": int(frame.unaffected_literature.sum()),
                "gnomad_unaffected_carrier_observations": int(
                    frame.gnomad_carriers.sum()
                ),
                "alpha_gets_affected_beta_gets_all_unaffected": True,
                "each_population_allele_once": True,
                "each_eligible_clinical_count_once": True,
            }
        )
    save(pd.DataFrame(summary), out / "union_coverage.csv")
    inputs.append(clinical_path)
    (out / "union_checks.json").write_text(
        json.dumps(
            {
                "input_sha256": {
                    str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in inputs
                },
                "population_assumption": "all included gnomAD carriers unaffected",
                "population_count": "joint AC minus joint homozygotes; never summed with exome/genome counts",
                "primary_scope": "all observed QC-passing small variants in full genomic gene span plus original CDS+75bp boundary padding and eligible clinical rows",
                "missing_population_members_are_not_imputed_from_hypothetical_variants": True,
                "unresolved_clinical_identity_rows_retained_in_ledger": True,
            },
            indent=2,
        )
        + "\n"
    )
    print(
        pd.DataFrame(moments)
        .loc[lambda d: d.scope.eq("all_observed_full_locus")]
        .to_string(index=False),
        flush=True,
    )


if __name__ == "__main__":
    main()
