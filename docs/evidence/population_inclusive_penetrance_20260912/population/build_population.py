"""Freeze complete gnomAD gene-query inventories with canonical annotations.

Uses only previously downloaded public files; never opens the feature warehouse.
Raw inventory includes AC0 and QC failures, marked ineligible for prior updates.
"""

from collections import Counter
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re

from fetch_annotations import API, HERE, RAW, REPO, TRANSCRIPTS

ACCESSIONS = {
    "HNF1A": "P20823",
    "GCK": "P35557",
    "LDLR": "P01130",
    "BRCA2": "P51587",
    "KCNQ1": "P51787",
}
AA = dict(
    zip(
        (
            "Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val Ter"
        ).split(),
        "ARNDCQEGHILKMFPSTWYV*",
        strict=True,
    )
)
CODING = {
    "missense_variant",
    "synonymous_variant",
    "frameshift_variant",
    "stop_gained",
    "start_lost",
    "stop_lost",
    "stop_retained_variant",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
    "coding_sequence_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "splice_region_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "start_retained_variant",
}
CLASS = {
    "missense_variant": "missense",
    "synonymous_variant": "synonymous",
    "frameshift_variant": "frameshift",
    "stop_gained": "stop_gained",
    "start_lost": "start_loss",
    "stop_lost": "stop_loss",
    "stop_retained_variant": "synonymous",
    "inframe_insertion": "inframe",
    "inframe_deletion": "inframe",
    "protein_altering_variant": "other_coding",
    "coding_sequence_variant": "other_coding",
}


def load_gzip(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def check_counts(variant):
    for assay in ("exome", "genome", "joint"):
        data = variant.get(assay)
        if data is None:
            continue
        ac, an, hom = (data.get(field) for field in ("ac", "an", "homozygote_count"))
        if any(not isinstance(value, int) for value in (ac, an, hom)):
            raise ValueError(
                f"Missing/noninteger population counts: {variant['variant_id']}"
            )
        if not 0 <= 2 * hom <= ac <= an:
            raise ValueError(f"Inconsistent autosomal counts: {variant['variant_id']}")
        if not isinstance(data.get("filters"), list):
            raise ValueError("Unknown QC state")


def joint_filter(variant):
    failed = [
        assay
        for assay in ("exome", "genome")
        if variant.get(assay) is not None and variant[assay]["filters"]
    ]
    if len(failed) == 2:
        return "BOTH_FILTERED"
    if failed:
        return f"{failed[0].upper()}S_FILTERED"
    return "PASS"


def protein_fields(hgvsp, sequence):
    out = {
        "aa_ref": "",
        "aa_pos": "",
        "aa_alt": "",
        "protein_key": "",
        "canonical_wt_status": "not_applicable",
    }
    if not hgvsp:
        return out
    change = hgvsp.removeprefix("p.").strip("()")
    anchors = re.findall(r"(?:^|_)([A-Z][a-z]{2})(\d+)", change)
    if not anchors:
        out["canonical_wt_status"] = "unparsed_reference"
        return out
    first, pos = anchors[0]
    out.update(aa_ref=AA.get(first, ""), aa_pos=int(pos))
    for ref, position in anchors:
        expected = AA.get(ref)
        position = int(position)
        if expected == "*" and position == len(sequence) + 1:
            continue
        if (
            expected is None
            or not 1 <= position <= len(sequence)
            or sequence[position - 1] != expected
        ):
            out["canonical_wt_status"] = "mismatch"
            return out
    out["canonical_wt_status"] = "match"
    simple = re.fullmatch(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=)", change)
    if simple:
        ref, pos, alt = simple.groups()
        out["aa_alt"] = AA[ref] if alt == "=" else AA.get(alt, "")
        out["protein_key"] = f"{AA[ref]}{pos}{out['aa_alt']}"
    return out


def transcript_label(variant):
    if not variant or not variant.get("transcript_id"):
        return ""
    return str(variant["transcript_id"]) + (
        "." + str(variant["transcript_version"])
        if variant.get("transcript_version")
        else ""
    )


def main():
    rows, summaries, sources, contexts = [], [], [], {}
    for gene, expected_transcript in TRANSCRIPTS.items():
        population_path = RAW / f"{gene}_population_20260911.json.gz"
        annotation_path = RAW / f"{gene}_annotations.json.gz"
        uniprot_path = RAW / f"{gene}_UniProt.json"
        source = load_gzip(population_path)
        annotation = load_gzip(annotation_path)
        response = annotation["response"]
        if response.get("errors"):
            raise ValueError("Partial/error annotation response")
        data = response["data"]
        seqdata = json.loads(uniprot_path.read_text())
        if source["gene"] != gene or source["release"] != "4.1.1":
            raise ValueError("Population source scope/release changed")
        if (
            seqdata["primaryAccession"] != ACCESSIONS[gene]
            or seqdata["organism"]["taxonId"] != 9606
            or not any(
                g.get("geneName", {}).get("value") == gene for g in seqdata["genes"]
            )
        ):
            raise ValueError("UniProt protein/gene identity mismatch")
        sequence = seqdata["sequence"]["value"]
        (HERE / f"{gene}_canonical.fasta").write_text(
            f">{gene} {ACCESSIONS[gene]} canonical displayed sequence\n{sequence}\n"
        )
        gene_data, transcript_data = data["gene"], data["transcript"]
        if (
            gene_data["canonical_transcript_id"] != expected_transcript
            or transcript_data["transcript_id"] != expected_transcript
            or gene_data["gene_id"] != transcript_data["gene_id"]
            or gene_data["gene_id"] != source["population"]["gene_id"]
        ):
            raise ValueError("Canonical transcript identity mismatch")
        canonical_label = (
            expected_transcript + "." + transcript_data["transcript_version"]
        )
        original = {v["variant_id"]: v for v in source["population"]["variants"]}
        gene_variants = {v["variant_id"]: v for v in gene_data["variants"]}
        canonical_variants = {v["variant_id"]: v for v in transcript_data["variants"]}
        if (
            len(original) != len(source["population"]["variants"])
            or len(gene_variants) != len(gene_data["variants"])
            or len(canonical_variants) != len(transcript_data["variants"])
        ):
            raise ValueError("Duplicate genomic alleles in endpoint response")
        if (
            original.keys() != gene_variants.keys()
            or not canonical_variants.keys() <= original.keys()
        ):
            raise ValueError("Annotation and population inventories differ")
        gene_rows = []
        for key in sorted(
            original, key=lambda value: (int(value.split("-")[1]), value)
        ):
            var = original[key]
            gene_annotation = gene_variants[key]
            canonical_annotation = canonical_variants.get(key)
            for assay in ("exome", "genome", "joint"):
                for field in ("ac", "an", "homozygote_count", "filters"):
                    if (var.get(assay) or {}).get(field) != (
                        gene_annotation.get(assay) or {}
                    ).get(field):
                        raise ValueError(
                            f"Population counts or QC drifted: {gene} {key}"
                        )
            check_counts(var)
            chrom, pos, ref, alt = key.split("-")
            if chrom != gene_data["chrom"] or chrom not in {
                str(x) for x in range(1, 23)
            }:
                raise ValueError("Carrier formula requires verified autosomal gene")
            canonical = canonical_annotation or {}
            consequence = canonical.get("consequence") or ""
            csq = canonical.get("transcript_consequence") or {}
            if canonical and (
                canonical["transcript_id"] != expected_transcript
                or csq.get("gene_id") != gene_data["gene_id"]
            ):
                raise ValueError("Wrong canonical transcript consequence")
            row = {
                "gene": gene,
                "variant_id": key,
                "chrom": chrom,
                "pos": int(pos),
                "ref": ref,
                "alt": alt,
                "variant_type": "SNV"
                if len(ref) == len(alt) == 1
                else "MNV"
                if len(ref) == len(alt)
                else "indel",
                "gnomad_release": "4.1.1",
                "canonical_transcript_id": canonical_label,
                "transcript_id": transcript_label(canonical),
                "hgvsc": canonical.get("hgvsc") or "",
                "hgvsp": canonical.get("hgvsp") or "",
                "consequence": consequence,
                "consequence_terms": ";".join(csq.get("consequence_terms") or []),
                "canonical_annotation_status": "canonical_transcript_endpoint"
                if canonical
                else "outside_canonical_transcript_endpoint",
                **protein_fields(canonical.get("hgvsp"), sequence),
                "vclass": CLASS.get(
                    consequence,
                    "splice"
                    if consequence.startswith("splice_")
                    else "noncoding"
                    if consequence
                    else "other_transcript",
                ),
                "coding_or_splice": bool(
                    set(csq.get("consequence_terms") or []) & CODING
                ),
                "gene_annotation_transcript_id": transcript_label(gene_annotation),
                "gene_annotation_hgvsc": gene_annotation.get("hgvsc") or "",
                "gene_annotation_hgvsp": gene_annotation.get("hgvsp") or "",
                "gene_annotation_consequence": gene_annotation.get("consequence") or "",
                "lof": csq.get("lof") or "",
                "lof_filter": csq.get("lof_filter") or "",
                "lof_flags": csq.get("lof_flags") or "",
            }
            for assay in ("exome", "genome", "joint"):
                value = var.get(assay)
                row[f"{assay}_present"] = value is not None
                for output, field in (
                    ("ac", "ac"),
                    ("an", "an"),
                    ("hom", "homozygote_count"),
                ):
                    row[f"{assay}_{output}"] = value[field] if value is not None else ""
                row[f"{assay}_filters"] = (
                    ";".join(value["filters"]) if value is not None else "ABSENT"
                )
            joint = var.get("joint")
            row.update(
                joint_api_pass=joint is not None and not joint["filters"],
                reconstructed_joint_filter=joint_filter(var),
                qc_pass=joint is not None and joint_filter(var) == "PASS",
                gnomad_carriers=joint["ac"] - joint["homozygote_count"]
                if joint is not None
                else "",
                observed=joint is not None and joint["ac"] > 0,
            )
            row["population_eligible"] = row["observed"] and row["qc_pass"]
            gene_rows.append(row)
        rows.extend(gene_rows)
        eligible = [r for r in gene_rows if r["population_eligible"]]
        coding = [r for r in eligible if r["coding_or_splice"]]
        summary = {
            "gene": gene,
            "all_source_rows": len(gene_rows),
            "canonical_annotation_rows": len(canonical_variants),
            "observed_joint_ac_positive": sum(r["observed"] for r in gene_rows),
            "observed_reconstructed_qc_pass": len(eligible),
            "coding_splice_qc_pass": len(coding),
            "canonical_missense_qc_pass": sum(
                r["vclass"] == "missense" for r in eligible
            ),
            "canonical_missense_wt_valid_qc_pass": sum(
                r["vclass"] == "missense" and r["canonical_wt_status"] == "match"
                for r in eligible
            ),
            "observed_qc_fail": sum(
                r["observed"] and not r["qc_pass"] for r in gene_rows
            ),
            "joint_ac_total_qc_pass": sum(r["joint_ac"] for r in eligible),
            "joint_hom_total_qc_pass": sum(r["joint_hom"] for r in eligible),
            "joint_carriers_total_qc_pass": sum(r["gnomad_carriers"] for r in eligible),
            "source_indels": sum(r["variant_type"] == "indel" for r in gene_rows),
            "qc_pass_indels": sum(r["variant_type"] == "indel" for r in eligible),
            "canonical_wt_mismatch_rows": sum(
                r["canonical_wt_status"] == "mismatch" for r in gene_rows
            ),
            "allele_inventory_and_all_cohort_counts_match_prior_snapshot": True,
        }
        summaries.append(summary)
        contexts[gene] = {
            "gene_metadata": {k: v for k, v in gene_data.items() if k != "variants"},
            "canonical_transcript_metadata": {
                k: v for k, v in transcript_data.items() if k != "variants"
            },
            "canonical_protein_accession": ACCESSIONS[gene],
            "canonical_protein_length": len(sequence),
            "population_fetched_at": source["fetched_at"],
            "annotations_fetched_at": annotation["fetched_at"],
            "major_consequence_counts_all_canonical_rows": dict(
                Counter(r["consequence"] for r in gene_rows)
            ),
        }
        for path, url in (
            (population_path, API),
            (annotation_path, API),
            (
                uniprot_path,
                f"https://rest.uniprot.org/uniprotkb/{ACCESSIONS[gene]}.json",
            ),
        ):
            sources.append(
                {
                    "path": str(path.relative_to(REPO)),
                    "url": url,
                    "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                    "bytes": path.stat().st_size,
                }
            )
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    (HERE / "population_variants.csv.gz").write_bytes(
        gzip.compress(output.getvalue().encode(), mtime=0)
    )
    with (HERE / "population_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(summaries[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(summaries)
    report = {
        "gnomad_release": "4.1.1",
        "dataset_selector": "gnomad_r4",
        "reference_genome": "GRCh38",
        "gene_query_boundary": "union of all gene CDS intervals plus 75bp padding on each side; merged overlaps; gene_id filter; paginated complete endpoint result; not entire gene locus",
        "canonical_query_boundary": "canonical transcript CDS intervals plus 75bp padding; exact transcript_id filter",
        "carrier_policy": "official joint AC minus joint homozygote_count, once per genomic allele; all genes autosomal; gnomAD carriers assumed unaffected",
        "qc_policy": "Reconstruct documented joint VCF categories from both assay filter lists. PASS requires every present assay PASS; absent assay permitted. Joint API filters retained but not trusted as joint VCF FILTER.",
        "qc_source": "https://discuss.gnomad.broadinstitute.org/t/is-joint-combined-genome-exome-faf-unreliable-if-either-genome-or-exome-fails-filters/88/3",
        "browser_source_commit": "5c14bfea3298f46854df14781adc1ccead58da3b",
        "browser_source_url": "https://github.com/broadinstitute/gnomad-browser/blob/5c14bfea3298f46854df14781adc1ccead58da3b/graphql-api/src/queries/variant-datasets/gnomad-v4-variant-queries.ts",
        "browser_source_findings": "getFilteredRegions/fetchVariantsByGene/fetchVariantsByTranscript define CDS+75bp boundary; shapeVariantSummary reads joint.filter singular despite selected joint.filters and may append AC0 from exome AC",
        "coding_splice_subset": sorted(CODING),
        "population_only_outcome_assumption": "No literature match in downstream union means A=0 under the user's gnomAD-unaffected assumption; collector itself does not assign literature overlap.",
        "noncoding_scope": "All returned observed QC-pass rows retained; coding/splice subset is separately identifiable. No claim that deep intronic, whole-locus regulatory, CNV or SV inventories are complete.",
        "sources": sources,
        "gene_contexts": contexts,
        "summaries": summaries,
    }
    (HERE / "population_provenance.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
