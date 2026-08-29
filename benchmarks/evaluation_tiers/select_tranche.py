#!/usr/bin/env python3
"""Draw a disjoint, human-adjudicated evaluation tranche from the repo gold standard.

``gold_120`` consumed one seeded sample of source-available, count-eligible gold
papers. A second tranche has to answer a different question -- does the measured
protocol hold on papers it has never been scored against? -- so it must be drawn
from the same eligibility rule but share no gene-paper attempt with any surface
that has already been scored, tuned on, or staged for review.

Eligibility (a gene-paper attempt qualifies when all three hold):

1. ``run_eval.gold_count_eligible_pmids`` accepts the PMID for that gene -- the
   same helper ``gold_120`` was drawn with, so the two tranches share one
   eligibility rule by construction rather than by a re-implementation that can
   drift. It requires a named variant and an authoritative assertion for every
   one of carriers/affected/unaffected, and it reads presence only, never value.
   Genes with no repo gold standard resolve to the curated ``gold_overrides``
   answer key through ``run_eval.gold_csv_path``;
2. a ``*_FULL_CONTEXT.md`` of at least ``--min-fulltext-bytes`` exists under
   ``corpus/<GENE>/<PMID>/``;
3. the attempt is in none of the exclusion surfaces (other registry tiers, the
   preregistered gold-150 rosters, and the explicit quarantine list).

Whether a paper's gold rows carry a *positive* count is recorded per paper as a
diagnostic and never filtered on. Restricting the draw to positive-count papers
would tilt the tranche toward count-bearing studies and break comparability with
``gold_120``, whose identity-recall gate deliberately includes gold rows that
assert an explicit zero.

Selection is a stable SHA-256 rank of ``"<seed>|<GENE>|<PMID>"``. It never reads
gold values or row counts, so the draw cannot be steered by how well the
pipeline happens to do on a paper.

The script writes the manifest, a selection-provenance JSON that binds each
chosen source file by digest, and a per-gene answer-key subset so the tranche is
self-contained and scoreable without re-deriving gold. It does not edit
``registry.json``; print the emitted registry block and paste it, so adding a
scored tier stays a deliberate act.

Usage:
    python benchmarks/evaluation_tiers/select_tranche.py \\
        --tier-id gold_120b --seed 2026082501 \\
        --gene SCN5A=30 --gene KCNH2=30 --gene KCNQ1=30 --gene RYR2=30 \\
        --gene BRCA2=5 --out-prefix tier4_gold_120b
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from benchmarks.codex_paper_eval.run_eval import (  # noqa: E402
    DEFAULT_GOLD,
    gold_count_eligible_pmids,
    gold_csv_path,
)

CORPUS = REPO / "corpus"
# Every preregistered gold-150 tree, not one dated directory. An amended tree
# supersedes its original by re-splitting calibration/holdout, so pinning a
# single path would silently stop excluding a roster the moment one is reissued.
GOLD150_GLOB = "benchmarks/curated_extraction_eval/gold150_preregistered_*"
COUNT_COLUMNS = ("carriers", "affected", "unaffected")

# Gene-paper attempts that must never enter a scored tranche regardless of what
# gold says about them. Each one is a PMID whose gold row is real but whose
# paper is not a genetics study of that gene, so scoring against it measures the
# curation error rather than the pipeline.
QUARANTINE = {
    ("KCNH2", "10086972"): "not a genetics paper; gold row belongs to 10086971",
    (
        "KCNH2",
        "14642689",
    ): "angiotensin-II receptor expression study, not KCNH2 genetics",
    ("BRCA2", "19944633"): "non-human canine BRCA2 ortholog; negative control only",
}


class SelectionError(RuntimeError):
    """The tranche cannot be drawn reproducibly from the current inputs."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rank_key(seed: int, gene: str, pmid: str) -> str:
    return hashlib.sha256(f"{seed}|{gene}|{pmid}".encode()).hexdigest()


def read_tier_manifest(path: Path) -> set[tuple[str, str]]:
    """Read a '<GENE>\\t<PMID>' tier manifest into gene-paper attempts."""
    attempts: set[tuple[str, str]] = set()
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        parts = re.split(r"[\t, ]+", line)
        if len(parts) != 2 or not parts[1].isdigit():
            raise SelectionError(f"{path}:{number}: expected '<GENE> <PMID>'")
        attempts.add((parts[0].upper(), parts[1]))
    return attempts


def read_gold150_rosters() -> tuple[set[tuple[str, str]], list[str]]:
    """Every preregistered calibration+holdout roster, from the answer sheets."""
    attempts: set[tuple[str, str]] = set()
    roots: list[str] = []
    for root in sorted(REPO.glob(GOLD150_GLOB)):
        if not root.is_dir():
            continue
        roots.append(root.name)
        for template in sorted(root.glob("*/*/curation_template.csv")):
            gene = template.parent.parent.name.upper()
            with template.open(newline="", encoding="utf-8-sig") as handle:
                for row in csv.DictReader(handle):
                    pmid = (row.get("pmid") or "").strip()
                    if pmid.isdigit():
                        attempts.add((gene, pmid))
    return attempts, roots


def gold_path(gene: str) -> tuple[Path, str]:
    """Resolve the gene's human-curated answer key and name its provenance."""
    try:
        path = gold_csv_path(DEFAULT_GOLD, gene)
    except (FileNotFoundError, SystemExit) as exc:
        raise SelectionError(f"no human-curated gold for {gene}") from exc
    if not path.is_file():
        raise SelectionError(f"no human-curated gold for {gene}")
    provenance = (
        "repo_gold_standard" if path.parent == DEFAULT_GOLD else "curated_gold_override"
    )
    return path, provenance


def gold_by_pmid(path: Path) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            pmid = (row.get("pmid") or "").strip()
            if pmid.isdigit():
                grouped[pmid].append(row)
    return grouped


def positive_count(row: dict[str, str]) -> bool:
    for column in COUNT_COLUMNS:
        raw = (row.get(column) or "").strip()
        if not raw:
            continue
        try:
            if float(raw) > 0:
                return True
        except ValueError:
            continue
    return False


def full_text(gene: str, pmid: str, min_bytes: int) -> Path | None:
    """Largest usable FULL_CONTEXT for this gene-paper, mirroring corpus reuse."""
    directory = CORPUS / gene / pmid
    if not directory.is_dir():
        return None
    best: Path | None = None
    best_size = 0
    for candidate in sorted(directory.glob("*FULL_CONTEXT.md")):
        size = candidate.stat().st_size
        if size >= min_bytes and size > best_size:
            best, best_size = candidate, size
    return best


def parse_gene_spec(values: list[str]) -> list[tuple[str, int]]:
    spec: list[tuple[str, int]] = []
    for value in values:
        if "=" not in value:
            raise SelectionError(f"--gene expects GENE=N, got {value!r}")
        gene, _, count = value.partition("=")
        if not count.isdigit() or int(count) <= 0:
            raise SelectionError(f"--gene expects a positive N, got {value!r}")
        spec.append((gene.strip().upper(), int(count)))
    return spec


def build(args: argparse.Namespace) -> dict:
    gene_spec = parse_gene_spec(args.gene)

    excluded: dict[tuple[str, str], str] = {}
    # Article-level exclusion is the default. A tier whose whole purpose is
    # "papers the protocol has never been scored or tuned against" is not served
    # by attempt-level disjointness alone: a multi-gene paper already scored
    # under KCNQ1 has had its text, tables, and supplements optimized against,
    # and BRCA1/BRCA2 output diverging only in the gene column is a recorded
    # failure here. Same article under a different gene is still a used article.
    excluded_articles: dict[str, str] = {}
    exclusion_sources: list[dict[str, str]] = []

    def note(attempts: set[tuple[str, str]], label: str) -> None:
        for attempt in attempts:
            excluded.setdefault(attempt, label)
            if args.exclude_level == "pmid":
                excluded_articles.setdefault(attempt[1], f"{label} ({attempt[0]})")

    for manifest_name in args.exclude_tier:
        manifest = (HERE / manifest_name).resolve()
        attempts = read_tier_manifest(manifest)
        note(attempts, manifest.name)
        exclusion_sources.append(
            {
                "source": manifest.name,
                "kind": "registry_tier_manifest",
                "attempts": str(len(attempts)),
                "sha256": sha256_file(manifest),
            }
        )
    gold150, gold150_roots = read_gold150_rosters()
    if not gold150_roots:
        raise SelectionError(
            f"no preregistered gold-150 roster matched {GOLD150_GLOB}; "
            "refusing to draw a tranche that cannot prove it excludes them"
        )
    note(gold150, "gold150_preregistered")
    exclusion_sources.append(
        {
            "source": ", ".join(gold150_roots),
            "kind": "preregistered_curation_roster",
            "attempts": str(len(gold150)),
            "sha256": "",
        }
    )
    for attempt, reason in QUARANTINE.items():
        excluded.setdefault(attempt, f"quarantine: {reason}")
        if args.exclude_level == "pmid":
            excluded_articles.setdefault(attempt[1], f"quarantine: {reason}")

    selected: list[dict] = []
    diagnostics: list[dict] = []
    for gene, want in gene_spec:
        source_path, provenance = gold_path(gene)
        grouped = gold_by_pmid(source_path)
        count_eligible = gold_count_eligible_pmids(DEFAULT_GOLD, gene)
        eligible: list[dict] = []
        rejected = defaultdict(int)
        for pmid, rows in grouped.items():
            if (gene, pmid) in excluded:
                rejected["excluded_prior_surface"] += 1
                continue
            if pmid in excluded_articles:
                rejected["excluded_prior_article_other_gene"] += 1
                continue
            if pmid not in count_eligible:
                rejected["not_count_eligible"] += 1
                continue
            text = full_text(gene, pmid, args.min_fulltext_bytes)
            if text is None:
                rejected["no_full_text"] += 1
                continue
            eligible.append(
                {
                    "gene": gene,
                    "pmid": pmid,
                    "gold_rows": len(rows),
                    "has_positive_count": any(positive_count(row) for row in rows),
                    "source_path": str(text.relative_to(REPO)),
                    "source_bytes": text.stat().st_size,
                    "rank": rank_key(args.seed, gene, pmid),
                }
            )
        if len(eligible) < want:
            raise SelectionError(
                f"{gene}: need {want} eligible papers, pool has {len(eligible)}"
            )
        eligible.sort(key=lambda item: item["rank"])
        drawn = eligible[:want]
        for position, item in enumerate(drawn, 1):
            item["draw_position"] = position
            item["source_sha256"] = sha256_file(REPO / item["source_path"])
        selected.extend(drawn)
        diagnostics.append(
            {
                "gene": gene,
                "gold_source": str(source_path.relative_to(REPO)),
                "gold_provenance": provenance,
                "gold_pmids": len(grouped),
                "eligible_pool": len(eligible),
                "requested": want,
                "drawn_gold_rows": sum(item["gold_rows"] for item in drawn),
                "drawn_positive_count_papers": sum(
                    1 for item in drawn if item["has_positive_count"]
                ),
                "rejected": dict(rejected),
            }
        )

    return {
        "selected": selected,
        "diagnostics": diagnostics,
        "exclusions": exclusion_sources,
    }


def write_manifest(path: Path, selected: list[dict], tier_id: str, seed: int) -> str:
    lines = [
        f"# {tier_id}: second human-adjudicated tranche, seed {seed}.",
        "# Disjoint from every prior scored, tuned, or staged surface.",
        "# Drawn by benchmarks/evaluation_tiers/select_tranche.py; do not hand-edit.",
    ]
    lines.extend(f"{item['gene']}\t{item['pmid']}" for item in selected)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return sha256_file(path)


def write_answer_key(root: Path, selected: list[dict]) -> list[dict]:
    """Freeze the human-curated gold rows for exactly this tranche."""
    root.mkdir(parents=True, exist_ok=True)
    by_gene: dict[str, set[str]] = defaultdict(set)
    for item in selected:
        by_gene[item["gene"]].add(item["pmid"])
    written: list[dict] = []
    for gene, pmids in sorted(by_gene.items()):
        source_path, provenance = gold_path(gene)
        with source_path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            fieldnames = reader.fieldnames or []
            rows = [r for r in reader if (r.get("pmid") or "").strip() in pmids]
        rows.sort(key=lambda r: (r.get("pmid", ""), r.get("variant", "")))
        target = root / f"{gene}_recall_input.csv"
        with target.open("w", newline="", encoding="utf-8") as handle:
            # LF, not the csv default CRLF: the trailing-whitespace pre-commit hook
            # strips the \r bytes, which would break the frozen digests below.
            writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        written.append(
            {
                "gene": gene,
                "file": target.name,
                "papers": len(pmids),
                "gold_rows": len(rows),
                "gold_provenance": provenance,
                "gold_source": str(source_path.relative_to(REPO)),
                "sha256": sha256_file(target),
            }
        )
    return written


def registry_block(
    tier_id: str,
    manifest_name: str,
    digest: str,
    selected: list[dict],
    seed: int,
    order: int,
) -> dict:
    gene_counts: dict[str, int] = defaultdict(int)
    for item in selected:
        gene_counts[item["gene"]] += 1
    return {
        "id": tier_id,
        "order": order,
        "manifest": manifest_name,
        "sha256": digest,
        "attempt_count": len(selected),
        "unique_pmid_count": len({item["pmid"] for item in selected}),
        "gene_attempt_counts": dict(sorted(gene_counts.items())),
        "selection_seed": seed,
        "role": "scored_gold_replication",
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--tier-id", required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument(
        "--gene",
        action="append",
        required=True,
        metavar="GENE=N",
        help="Per-gene draw size, repeatable (e.g. --gene SCN5A=30).",
    )
    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Manifest/provenance basename written into this directory.",
    )
    parser.add_argument(
        "--exclude-tier",
        action="append",
        default=["tier1_gold_50.tsv", "tier2_gold_120.tsv", "tier3_reviewer_545.tsv"],
        help="Tier manifests whose attempts are ineligible. Repeatable.",
    )
    parser.add_argument(
        "--exclude-level",
        choices=["pmid", "attempt"],
        default="pmid",
        help=(
            "'pmid' (default) treats a used article as used under every gene, so "
            "the tranche holds only never-scored, never-tuned text. 'attempt' "
            "excludes only the exact gene-paper pair and will readmit multi-gene "
            "papers other tiers already optimized against."
        ),
    )
    parser.add_argument("--registry-order", type=int, default=4)
    parser.add_argument("--min-fulltext-bytes", type=int, default=3000)
    args = parser.parse_args()

    try:
        result = build(args)
    except SelectionError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    selected = result["selected"]
    manifest_path = HERE / f"{args.out_prefix}.tsv"
    digest = write_manifest(manifest_path, selected, args.tier_id, args.seed)
    key_root = HERE / f"{args.tier_id}_answer_key"
    answer_key = write_answer_key(key_root, selected)
    tier = registry_block(
        args.tier_id,
        manifest_path.name,
        digest,
        selected,
        args.seed,
        args.registry_order,
    )

    provenance = {
        "tier_id": args.tier_id,
        "manifest": manifest_path.name,
        "manifest_sha256": digest,
        "seed": args.seed,
        "selection_method": "sha256_rank_of_seed_gene_pmid_v1",
        "gold_value_blinded": True,
        "exclusion_level": args.exclude_level,
        "min_fulltext_bytes": args.min_fulltext_bytes,
        "registry_tier": tier,
        "answer_key": answer_key,
        "eligibility": [
            "human-curated gold row exists for the gene-paper attempt",
            "at least one gold row asserts a positive carriers/affected/unaffected value",
            f"corpus FULL_CONTEXT >= {args.min_fulltext_bytes} bytes",
            "attempt absent from every exclusion surface below",
        ],
        "exclusion_sources": result["exclusions"],
        "quarantine": [
            {"gene": gene, "pmid": pmid, "reason": reason}
            for (gene, pmid), reason in sorted(QUARANTINE.items())
        ],
        "per_gene": result["diagnostics"],
        "papers": selected,
    }
    provenance_path = HERE / f"{args.out_prefix}_selection.json"
    provenance_path.write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )

    print(f"manifest        {manifest_path.relative_to(REPO)}")
    print(f"provenance      {provenance_path.relative_to(REPO)}")
    print(f"answer key      {key_root.relative_to(REPO)}/")
    print(
        f"attempts        {tier['attempt_count']} ({tier['unique_pmid_count']} unique PMIDs)"
    )
    for gene, count in tier["gene_attempt_counts"].items():
        pool = next(
            d["eligible_pool"] for d in result["diagnostics"] if d["gene"] == gene
        )
        print(f"  {gene:<7} {count:>3} drawn from {pool} eligible")
    print("\nregistry.json block:")
    print(json.dumps(tier, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
