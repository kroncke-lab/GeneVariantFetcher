"""Pin audit artifacts and implementation; verify working, staged or HEAD bytes."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
MANIFEST = HERE / "manifest.json"
IMPLEMENTATION = ["TASKS.md"]


def digest(blob):
    return hashlib.sha256(blob).hexdigest()


def artifact_paths():
    return [
        p
        for p in sorted(HERE.rglob("*"))
        if p.is_file() and p != MANIFEST and "__pycache__" not in p.parts
    ]


def validate_format(path, blob):
    assert len(blob) < 1_200_000, (path, len(blob))
    if path.name.endswith(".csv.gz"):
        blob = gzip.decompress(blob)
    if path.name.endswith((".csv", ".csv.gz")):
        assert b"\r" not in blob and blob.endswith(b"\n"), path


def read_blob(path, tree):
    if tree == "working":
        return path.read_bytes()
    relative = str(path.relative_to(REPO))
    ref = f":{relative}" if tree == "staged" else f"HEAD:{relative}"
    return subprocess.check_output(["git", "show", ref], cwd=REPO)


def external_sources():
    sources = {}
    manifests = sorted((HERE / "analysis").glob("*/input_hashes.json"))
    manifests.append(HERE / "plot_input_hashes.json")
    for manifest in manifests:
        for relative, expected in json.loads(manifest.read_text()).items():
            path = (HERE.parent / relative).resolve()
            assert digest(path.read_bytes()) == expected, path
            if not path.is_relative_to(HERE):
                sources[str(path)] = expected
    return sources


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true")
    parser.add_argument(
        "--tree", choices=["working", "staged", "HEAD"], default="working"
    )
    args = parser.parse_args()
    if args.build:
        assert args.tree == "working"
        artifacts = {}
        for path in artifact_paths():
            blob = path.read_bytes()
            validate_format(path, blob)
            artifacts[str(path.relative_to(HERE))] = dict(
                bytes=len(blob), sha256=digest(blob)
            )
        manifest = dict(
            schema=1,
            artifacts=artifacts,
            implementation={
                relative: digest((REPO / relative).read_bytes())
                for relative in IMPLEMENTATION
            },
            external_sources=external_sources(),
        )
        MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    manifest = json.loads(MANIFEST.read_text())
    if args.tree != "working":
        assert read_blob(MANIFEST, args.tree) == MANIFEST.read_bytes(), MANIFEST
    for relative, record in manifest["artifacts"].items():
        path = HERE / relative
        blob = read_blob(path, args.tree)
        validate_format(path, blob)
        assert len(blob) == record["bytes"] and digest(blob) == record["sha256"], path
    for relative, expected in manifest["implementation"].items():
        assert digest(read_blob(REPO / relative, args.tree)) == expected, relative
    for absolute, expected in manifest["external_sources"].items():
        assert digest(Path(absolute).read_bytes()) == expected, absolute
    if args.tree == "working":
        assert {str(p.relative_to(HERE)) for p in artifact_paths()} == set(
            manifest["artifacts"]
        )
    print(
        json.dumps(
            dict(
                status="pass",
                tree=args.tree,
                artifacts=len(manifest["artifacts"]),
                implementation_files=len(manifest["implementation"]),
                external_sources=len(manifest["external_sources"]),
            )
        )
    )


if __name__ == "__main__":
    main()
