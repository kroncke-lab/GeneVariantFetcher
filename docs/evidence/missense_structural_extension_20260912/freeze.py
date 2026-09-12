"""Freeze provenance and verify artifacts locally or against committed Git blobs."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DESTINATION = HERE / "manifest.json"
PPA = REPO.parent / "ProteinProximityAnalysis"


def digest(blob):
    return hashlib.sha256(blob).hexdigest()


def files():
    return [
        p
        for p in sorted(HERE.rglob("*"))
        if p.is_file() and p != DESTINATION and "__pycache__" not in p.parts
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--git", action="store_true")
    args = parser.parse_args()
    if args.verify or args.git:
        manifest = json.loads(DESTINATION.read_text())
        if args.git:
            committed_manifest = subprocess.check_output(
                ["git", "show", f"HEAD:{DESTINATION.relative_to(REPO)}"], cwd=REPO
            )
            assert committed_manifest == DESTINATION.read_bytes(), (
                "Manifest differs from committed blob"
            )
        for relative, expected in manifest["artifacts"].items():
            path = HERE / relative
            if args.git:
                blob = subprocess.check_output(
                    ["git", "show", f"HEAD:{path.relative_to(REPO)}"], cwd=REPO
                )
            else:
                blob = path.read_bytes()
            assert digest(blob) == expected, relative
        for relative, expected in manifest["source_snapshots"].items():
            path = HERE.parent / relative
            blob = (
                subprocess.check_output(
                    ["git", "show", f"HEAD:{path.relative_to(REPO)}"], cwd=REPO
                )
                if args.git
                else path.read_bytes()
            )
            assert digest(blob) == expected, relative
        for relative, expected in manifest["ppa_modules"].items():
            blob = (
                subprocess.check_output(
                    ["git", "show", f"{manifest['ppa_commit']}:{relative}"], cwd=PPA
                )
                if args.git
                else (PPA / relative).read_bytes()
            )
            assert digest(blob) == expected, relative
        print(
            f"Verified {len(manifest['artifacts'])} artifacts, {len(manifest['source_snapshots'])} source snapshots and 2 engine modules"
        )
        return
    source_paths = set()
    for gene in ["HNF1A", "LDLR", "KCNQ1", "BRCA2"]:
        check = json.loads((HERE / f"analysis/{gene}/checks.json").read_text())
        source_paths.update(HERE.parent / relative for relative in check["inputs"])
    source_paths.update(
        [
            HERE.parent / "gck_structural_pilot_20260912/pilot_statistics.py",
            HERE.parent
            / "class_matched_penetrance_20260912/structural/GCK_primary_variant_density.csv",
            HERE.parent
            / "class_matched_penetrance_20260912/structural/GCK_variant_loo_predictions.csv.gz",
        ]
    )
    artifact_paths = files()
    assert all(p.stat().st_size < 1_200_000 for p in artifact_paths)
    manifest = {
        "schema_version": 1,
        "analysis_date": "2026-09-12",
        "contract": "Gene-specific missense priors; alpha affected; beta literature unaffected plus all gnomAD; fixed priors and variant-only outer LOO; positive-tail sigmoid; 3D/polymer sources separate.",
        "ppa_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=PPA, text=True
        ).strip(),
        "ppa_modules": {
            relative: digest((PPA / relative).read_bytes())
            for relative in [
                "src/alphafold_rin/empirical_density.py",
                "src/alphafold_rin/empirical_context.py",
            ]
        },
        "source_snapshots": {
            str(p.relative_to(HERE.parent)): digest(p.read_bytes())
            for p in sorted(source_paths)
        },
        "artifacts": {
            str(p.relative_to(HERE)): digest(p.read_bytes()) for p in artifact_paths
        },
    }
    DESTINATION.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Froze {len(artifact_paths)} artifacts")


if __name__ == "__main__":
    main()
