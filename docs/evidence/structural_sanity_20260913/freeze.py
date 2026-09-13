"""Freeze this correction and verify bytes, source provenance and Git preservation."""

import argparse
import ast
import gzip
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
EVIDENCE = HERE.parent
REPOSITORIES = REPO.parent
DESTINATION = HERE / "manifest.json"
REPO_NAMES = [
    "GeneVariantFetcher",
    "ProteinProximityAnalysis",
    "BayesianPenetranceEstimator",
]


def digest(blob):
    return hashlib.sha256(blob).hexdigest()


def file_digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def read(relative):
    return json.loads((HERE / relative).read_text())


def artifact_paths():
    return [
        path
        for path in sorted(HERE.rglob("*"))
        if path.is_file() and path != DESTINATION and "__pycache__" not in path.parts
    ]


def verify_artifact_format(path, blob):
    assert len(blob) < 1_200_000, f"Oversized artifact: {path}"
    if path.name.endswith(".csv.gz"):
        blob = gzip.decompress(blob)
    if path.name.endswith((".csv", ".csv.gz")):
        assert b"\r" not in blob, f"CSV payload must use LF: {path}"
        assert blob.endswith(b"\n"), path


def verify_formatting():
    transitions = {}
    for row in read("formatting_receipt.json")["files"]:
        path = HERE / row["file"]
        assert file_digest(path) == row["formatted_source_sha256"], path
        tree = ast.dump(ast.parse(path.read_bytes()), include_attributes=False)
        assert digest(tree.encode()) == row["identical_ast_sha256"], path
        transitions[(path, row["executed_source_sha256"])] = row[
            "formatted_source_sha256"
        ]
    executed = read("analysis/BRCA2/checks.json")["outer_loo_sha256"]
    current = file_digest(HERE / "outer_loo.py")
    assert (
        executed == current or transitions[(HERE / "outer_loo.py", executed)] == current
    )


def input_paths():
    inputs = {}

    def add(path, expected=None):
        path = Path(path)
        actual = file_digest(path)
        assert expected is None or actual == expected, f"Source changed: {path}"
        if not path.is_relative_to(HERE):
            inputs[path] = actual

    for relative, expected in read("analysis/BRCA2/checks.json")[
        "source_hashes"
    ].items():
        add(EVIDENCE / relative, expected)
    for relative, expected in read("gck/checks.json")["source_hashes"].items():
        add(REPO / relative, expected)
    for absolute, expected in read("analysis/allgene_audit_receipt.json")[
        "sources"
    ].items():
        add(absolute, expected)
    for relative, expected in read("reviews/local_source_receipt.json")[
        "sources"
    ].items():
        add(REPO / relative, expected)
    for source in read("geometry/BRCA2/geometry_manifest.json")["sources"]:
        relative = source["file"]
        base = (
            REPOSITORIES if relative.startswith("ProteinProximityAnalysis/") else REPO
        )
        add(base / relative, source["sha256"])
    provenance = read("gck/frozen_endpoint_provenance.json")
    add(provenance["source_path"], provenance["source_sha256"])
    for source in read("gck/primary_sources.json"):
        if source["raw_cached_path"]:
            add(REPO / source["raw_cached_path"], source["raw_sha256"])
    cache = read("analysis/BRCA2/cache_manifest.json")
    raw = REPO / "results/structural_sanity_20260913/BRCA2"
    add(raw / "all_global_exclusions.float64", cache["exclusion_cache_sha256"])
    for relative, expected in cache["raw_weight_shards"].items():
        add(raw / relative, expected)
    for relative in [
        "gck_structural_pilot_20260912/pilot_statistics.py",
        "missense_structural_extension_20260912/run_structure.py",
    ]:
        add(EVIDENCE / relative)
    for path in (
        EVIDENCE / "missense_structural_extension_20260912/analysis/BRCA2"
    ).glob("primary_density.part*.csv.gz"):
        add(path)
    for module in ["empirical_density.py", "empirical_context.py"]:
        add(REPOSITORIES / "ProteinProximityAnalysis/src/alphafold_rin" / module)
    return inputs


def git_blob(repo, revision, relative):
    return subprocess.check_output(["git", "show", f"{revision}:{relative}"], cwd=repo)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--git", action="store_true")
    args = parser.parse_args()
    verify_formatting()
    files = artifact_paths()
    if args.verify or args.git:
        manifest = read("manifest.json")
        assert set(manifest["artifacts"]) == {
            str(path.relative_to(HERE)) for path in files
        }
        if args.git:
            assert (
                git_blob(REPO, "HEAD", DESTINATION.relative_to(REPO))
                == DESTINATION.read_bytes()
            )
        for relative, expected in manifest["artifacts"].items():
            path = HERE / relative
            blob = (
                git_blob(REPO, "HEAD", path.relative_to(REPO))
                if args.git
                else path.read_bytes()
            )
            assert digest(blob) == expected, relative
            verify_artifact_format(path, blob)
        verified = skipped = 0
        for relative, record in manifest["inputs"].items():
            name, within = relative.split("/", 1)
            repo = REPOSITORIES / name
            if args.git and record["storage"] == "local_cache":
                skipped += 1
                continue
            if args.git:
                revision = (
                    "HEAD"
                    if name == REPO.name
                    else manifest["repository_commits"][name]
                )
                actual = digest(git_blob(repo, revision, within))
            else:
                actual = file_digest(repo / within)
            assert actual == record["sha256"], relative
            verified += 1
        print(
            f"Verified {len(files)} artifacts and {verified} inputs; {skipped} local caches outside Git"
        )
        return
    inputs = input_paths()
    tracked, commits = {}, {}
    for name in REPO_NAMES:
        repo = REPOSITORIES / name
        tracked[name] = set(
            subprocess.check_output(["git", "ls-files", "-z"], cwd=repo)
            .decode()
            .split("\0")
        )
        commits[name] = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=repo, text=True
        ).strip()
    for path in files:
        verify_artifact_format(path, path.read_bytes())
    records = {}
    for path, sha in sorted(inputs.items()):
        relative = path.relative_to(REPOSITORIES)
        name, within = str(relative).split("/", 1)
        assert name in REPO_NAMES, path
        records[str(relative)] = {
            "sha256": sha,
            "storage": "git" if within in tracked[name] else "local_cache",
        }
    manifest = {
        "schema_version": 1,
        "analysis_date": "2026-09-13",
        "contract": "Fixed gene-by-missense empirical priors; alpha affected; beta literature unaffected plus all gnomAD assumed unaffected. Variant-only global outer LOO retaining other same-residue units. Positive-tail sigmoid; separate 3D frames and same-segment canonical polymer. Neighborhood scores are features, not target disease risks.",
        "repository_commits": commits,
        "inputs": records,
        "artifacts": {str(path.relative_to(HERE)): file_digest(path) for path in files},
        "cache_note": "Large numeric caches and original external structure/paper downloads stay outside Git. Their hashes are verified locally; compact geometry, all prediction rows and independent audits are committed.",
        "formatting_note": "Original execution byte hashes retained. formatting_receipt.json proves AST-identical Ruff formatting; independent output audits repeated using formatted sources.",
    }
    DESTINATION.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Froze {len(files)} artifacts and {len(records)} inputs")


if __name__ == "__main__":
    main()
