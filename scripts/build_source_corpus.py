#!/usr/bin/env python3
"""Build / incrementally update the consolidated source corpus.

The corpus is the single deduplicated home for fetched source:

    corpus/<GENE>/<PMID>/{PMID}_FULL_CONTEXT.md   # best usable full text
                         {PMID}_CLEANED.md          # (if present)
                         {PMID}_artifacts.json      # (if present)
                         {PMID}_figures/            # richest figure set
                         {PMID}_supplements/        # richest supplement set
    corpus/INDEX.json / INDEX.csv                   # gene -> pmid -> status
    corpus/MANIFEST.sha256                          # integrity for zip/transfer

Idempotent + incremental + never-downgrade. It scans EVERY ``*_FULL_CONTEXT.md``
under the source roots (not just ``pmc_fulltext/`` dirs, so recovery/acquisition
dirs are picked up too), treats the current corpus copy as one more candidate,
and for each (gene, PMID) keeps:

* full text: the **usable** copy (``is_usable_fulltext_source``) with the most
  bytes — so a real body always beats a paywall/abstract stub;
* figures / supplements: the copy with the most files.

A category is only rewritten when a strictly better candidate appears (e.g. a
new publisher key unlocks a body, or a supplement finally downloads), so:

* re-running with no new source is a no-op;
* changing a key + re-running only ADDS new papers and UPGRADES compromised
  categories, never touching already-complete ones.

Default is a dry-run; pass ``--apply`` to write. ``--verify`` checks corpus
integrity against MANIFEST.sha256. ``--rebuild`` wipes the corpus first.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from pipeline.source_quality import is_usable_fulltext_source  # noqa: E402

DEFAULT_ROOTS = ["results", "validation_runs"]
KNOWN_GENES = {
    "KCNH2",
    "KCNQ1",
    "SCN5A",
    "RYR2",
    "KCNE1",
    "KCNE2",
    "KCNJ2",
    "CACNA1C",
    "GJA5",
    "KCNA5",
    "NPPA",
    "PITX2",
}
SUFFIX = "_FULL_CONTEXT.md"
IMG_EXT = {".png", ".jpg", ".jpeg", ".gif", ".tif", ".tiff", ".webp", ".bmp"}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def looks_like_gene_dir(name: str) -> bool:
    """Return True for path components that plausibly name a gene."""
    cand = name.upper()
    return cand in KNOWN_GENES or (
        cand == name and any(ch.isalpha() for ch in cand) and len(cand) <= 12
    )


def infer_gene(
    p: Path,
    pmid: str,
    corpus_pmid_gene: dict[str, str],
    source_bases: list[Path],
) -> str | None:
    for base in source_bases:
        try:
            rel = p.relative_to(base)
        except ValueError:
            continue
        for part in rel.parts[:-1]:
            if looks_like_gene_dir(part):
                return part.upper()
    up = p.as_posix().upper()
    for g in KNOWN_GENES:
        if f"/{g}/" in up or f"_{g}_" in up or f"_{g}." in up:
            return g
    # Already-known PMID -> reuse its gene (handles gene-less manual dirs).
    if pmid in corpus_pmid_gene:
        return corpus_pmid_gene[pmid]
    # Last resort: scan the body for a single known gene token.
    try:
        head = p.read_text(encoding="utf-8", errors="replace")[:8192].upper()
        hits = {g for g in KNOWN_GENES if g in head}
        if len(hits) == 1:
            return next(iter(hits))
    except OSError:
        pass
    return None


def count_files(d: Path | None, exts: set[str] | None = None) -> int:
    if not d or not d.is_dir():
        return 0
    return sum(
        1
        for f in d.rglob("*")
        if f.is_file() and (exts is None or f.suffix.lower() in exts)
    )


def candidate(ft: Path):
    """Build a candidate record for one {PMID}_FULL_CONTEXT.md file."""
    d, pmid = ft.parent, ft.name[: -len(SUFFIX)]
    figd, supd = d / f"{pmid}_figures", d / f"{pmid}_supplements"
    return {
        "ft": ft,
        "ft_size": ft.stat().st_size,
        "usable": is_usable_fulltext_source(ft),
        "cleaned": (d / f"{pmid}_CLEANED.md")
        if (d / f"{pmid}_CLEANED.md").exists()
        else None,
        "art": (d / f"{pmid}_artifacts.json")
        if (d / f"{pmid}_artifacts.json").exists()
        else None,
        "figs": figd if figd.is_dir() else None,
        "nfigs": count_files(figd, IMG_EXT),
        "suppl": supd if supd.is_dir() else None,
        "nsuppl": count_files(supd),
    }


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--apply", action="store_true", help="Write changes (default: dry-run)."
    )
    ap.add_argument(
        "--rebuild", action="store_true", help="Wipe corpus/ before merging."
    )
    ap.add_argument(
        "--verify",
        action="store_true",
        help="Check corpus against MANIFEST.sha256 and exit.",
    )
    ap.add_argument(
        "--out",
        default=str(REPO / "corpus"),
        help="Corpus dir (default: <repo>/corpus).",
    )
    ap.add_argument(
        "--roots",
        nargs="*",
        default=DEFAULT_ROOTS,
        help=(
            "Source roots to scan (default: results validation_runs). "
            "Pass legacy/ad-hoc roots explicitly when needed."
        ),
    )
    args = ap.parse_args()
    out = Path(args.out).expanduser().resolve()
    dry = not args.apply
    source_bases = [(REPO / r).resolve() for r in args.roots]

    if args.verify:
        return verify(out)

    if args.rebuild and args.apply and out.exists():
        shutil.rmtree(out)

    # Map existing corpus PMID -> gene (for gene-less source dirs) and current entries.
    corpus_pmid_gene: dict[str, str] = {}
    if out.exists():
        for ft in out.rglob("*" + SUFFIX):
            corpus_pmid_gene.setdefault(ft.parts[-2], ft.parts[-3])

    # Gather candidates per (gene, pmid) from all source roots.
    cand: dict[tuple[str, str], list[dict]] = defaultdict(list)
    skipped_no_gene: list[str] = []
    for base in source_bases:
        if not base.exists():
            continue
        for ft in base.rglob("*" + SUFFIX):
            if not ft.is_file() or "/corpus/" in ft.as_posix():
                continue  # skip broken symlinks / dirs
            pmid = ft.name[: -len(SUFFIX)]
            if not pmid.isdigit():
                continue
            gene = infer_gene(ft.resolve(), pmid, corpus_pmid_gene, source_bases)
            if not gene:
                skipped_no_gene.append(ft.as_posix())
                continue
            cand[(gene, pmid)].append(candidate(ft))

    # Add the current corpus copy as a candidate so we never downgrade.
    for ft in out.rglob("*" + SUFFIX) if out.exists() else []:
        if ft.is_file():
            cand[(ft.parts[-3], ft.parts[-2])].append(
                {**candidate(ft), "_corpus": True}
            )

    actions = Counter()
    index: list[dict] = []
    for (gene, pmid), recs in sorted(cand.items()):
        best_ft = max(recs, key=lambda r: (r["usable"], r["ft_size"]))
        best_fig = max(recs, key=lambda r: r["nfigs"])
        best_sup = max(recs, key=lambda r: r["nsuppl"])
        dest = out / gene / pmid
        cur = dest / best_ft["ft"].name
        cur_sha = sha256_file(cur) if cur.exists() else None
        new_sha = sha256_file(best_ft["ft"])

        ft_change = new_sha != cur_sha and not best_ft.get("_corpus")
        cur_nfig = (
            count_files(dest / f"{pmid}_figures", IMG_EXT) if dest.exists() else 0
        )
        cur_nsup = count_files(dest / f"{pmid}_supplements") if dest.exists() else 0
        fig_change = best_fig["nfigs"] > cur_nfig and not best_fig.get("_corpus")
        sup_change = best_sup["nsuppl"] > cur_nsup and not best_sup.get("_corpus")

        if not dest.exists():
            action = "add"
        elif ft_change or fig_change or sup_change:
            action = "upgrade"
        else:
            action = "noop"
        actions[action] += 1

        if not dry and action != "noop":
            dest.mkdir(parents=True, exist_ok=True)
            if ft_change or action == "add":
                shutil.copy2(best_ft["ft"], dest / best_ft["ft"].name)
                for k in ("cleaned", "art"):
                    if best_ft[k]:
                        shutil.copy2(best_ft[k], dest / best_ft[k].name)
            for key, changed in (
                ("figs", fig_change or action == "add"),
                ("suppl", sup_change or action == "add"),
            ):
                src = best_fig[key] if key == "figs" else best_sup[key]
                if src and changed:
                    dst = dest / src.name
                    if dst.exists():
                        shutil.rmtree(dst)
                    shutil.copytree(src, dst)

        index.append(
            {
                "gene": gene,
                "pmid": pmid,
                "fulltext": str((dest / best_ft["ft"].name).relative_to(REPO)),
                "fulltext_bytes": best_ft["ft_size"],
                "source_sha256": new_sha,
                "full_text_status": "ok" if best_ft["usable"] else "stub",
                "n_figures": max(best_fig["nfigs"], cur_nfig),
                "n_supplement_files": max(best_sup["nsuppl"], cur_nsup),
                "n_source_copies": sum(1 for r in recs if not r.get("_corpus")),
                "chosen_from": str(best_ft["ft"].parent.relative_to(REPO))
                if not best_ft.get("_corpus")
                else "corpus",
            }
        )

    if not dry and index:
        write_index(out, index)
        write_manifest(out)

    # Summary
    per_gene = Counter(r["gene"] for r in index)
    stubs = [r for r in index if r["full_text_status"] == "stub"]
    print(
        "MODE:",
        "DRY-RUN (nothing written; pass --apply)" if dry else f"APPLIED -> {out}",
    )
    print(
        f"actions: add={actions['add']} upgrade={actions['upgrade']} noop={actions['noop']}"
    )
    print(
        f"corpus (gene,PMID) entries: {len(index)}  "
        + " ".join(f"{g}={c}" for g, c in sorted(per_gene.items()))
    )
    print(f"full_text ok: {len(index) - len(stubs)}  stub/compromised: {len(stubs)}")
    print(
        f"with figures: {sum(1 for r in index if r['n_figures'])}  "
        f"with supplements: {sum(1 for r in index if r['n_supplement_files'])}"
    )
    if skipped_no_gene:
        print(
            f"WARN: {len(skipped_no_gene)} FULL_CONTEXT files skipped (gene not inferable), e.g. {skipped_no_gene[:3]}"
        )
    if stubs and not dry:
        print(
            f"NOTE: {len(stubs)} entries are still stub-only (no usable full text found anywhere)."
        )
    return 0


def write_index(out: Path, index: list[dict]) -> None:
    out.mkdir(parents=True, exist_ok=True)
    by_gene: dict[str, dict] = defaultdict(dict)
    for r in index:
        by_gene[r["gene"]][r["pmid"]] = {
            k: r[k]
            for k in (
                "fulltext",
                "fulltext_bytes",
                "source_sha256",
                "full_text_status",
                "n_figures",
                "n_supplement_files",
                "n_source_copies",
                "chosen_from",
            )
        }
    (out / "INDEX.json").write_text(
        json.dumps(
            {
                "description": "Consolidated, deduplicated, quality-gated GVF source corpus; best usable copy per (gene, PMID).",
                "genes": by_gene,
            },
            indent=2,
        )
    )
    with (out / "INDEX.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(index[0].keys()))
        w.writeheader()
        w.writerows(index)


def write_manifest(out: Path) -> None:
    """sha256 of every FULL_CONTEXT.md + a corpus-wide file/byte summary."""
    lines, total_files, total_bytes = [], 0, 0
    for f in sorted(out.rglob("*")):
        if f.is_file() and f.name not in {"MANIFEST.sha256"}:
            total_files += 1
            total_bytes += f.stat().st_size
            if f.name.endswith(SUFFIX):
                lines.append(f"{sha256_file(f)}  {f.relative_to(out)}")
    header = [f"# total_files {total_files}", f"# total_bytes {total_bytes}"]
    (out / "MANIFEST.sha256").write_text("\n".join(header + lines) + "\n")


def verify(out: Path) -> int:
    man = out / "MANIFEST.sha256"
    if not man.exists():
        print("ERROR: no MANIFEST.sha256; run with --apply first.")
        return 1
    exp_files = exp_bytes = None
    ft_sha: dict[str, str] = {}
    for line in man.read_text().splitlines():
        if line.startswith("# total_files"):
            exp_files = int(line.split()[-1])
        elif line.startswith("# total_bytes"):
            exp_bytes = int(line.split()[-1])
        elif line and not line.startswith("#"):
            sha, rel = line.split("  ", 1)
            ft_sha[rel] = sha
    bad = [
        rel
        for rel, sha in ft_sha.items()
        if not (out / rel).exists() or sha256_file(out / rel) != sha
    ]
    act_files = sum(
        1 for f in out.rglob("*") if f.is_file() and f.name != "MANIFEST.sha256"
    )
    act_bytes = sum(
        f.stat().st_size
        for f in out.rglob("*")
        if f.is_file() and f.name != "MANIFEST.sha256"
    )
    ok = not bad and act_files == exp_files and act_bytes == exp_bytes
    print(f"full-text files checked: {len(ft_sha)}  mismatched/missing: {len(bad)}")
    print(f"total files: {act_files}/{exp_files}  total bytes: {act_bytes}/{exp_bytes}")
    print("VERIFY OK" if ok else "VERIFY FAILED")
    if bad:
        print("  e.g.", bad[:5])
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
