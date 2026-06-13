#!/usr/bin/env python3
"""Build a turnkey, blinded manual-curation packet for a cold-start gold standard.

A cold-start run (e.g. BRCA2) produces variants but no recall number — there is
no ground truth to score against. This script assembles a SELF-CONTAINED packet
a research assistant can curate on their own computer (no repo, no Python, no
API keys, no institutional access needed): it samples N full-text papers, copies
the already-fetched full text into the packet, and emits a curation spreadsheet
+ a protocol. The assistant fills the spreadsheet; ``score_curation_packet.py``
then scores the pipeline against their answers (recall / precision / F2).

It is BLINDED by construction: the packet contains the papers and a blank
template, never the pipeline's extracted variants, so the curator can't anchor
on them.

Usage:
    python scripts/build_curation_packet.py \
        --run-dir results/BRCA2/<ts>/BRCA2/<ts> --gene BRCA2 \
        --n 50 --seed 42 --out curation_packets/BRCA2_gold_50

Produces ``<out>/`` (and ``<out>.zip``) with:
    PROTOCOL.md            - the manual-extraction instructions
    manifest.csv           - the N papers (pmid, title, pubmed_url, file)
    curation_template.csv  - blank, one prefilled row per paper, gold schema
    papers/<pmid>.md       - the fetched full text to curate from
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

# Curation spreadsheet schema. The first five columns map onto the normalized
# gold recall input (variant, pmid, carriers, affected, unaffected); the extras
# capture the germline/somatic call (BRCA2 mixes both) and the evidence trail.
TEMPLATE_COLUMNS = [
    "pmid",
    "variant",
    "germline_or_somatic",
    "carriers",
    "affected",
    "unaffected",
    "evidence_note",
]


def _full_text_papers(run_dir: Path, min_bytes: int) -> list[tuple[str, Path, int]]:
    """Full-text papers = a FULL_CONTEXT.md with substantial body (not a stub)."""
    out: list[tuple[str, Path, int]] = []
    ft_dir = run_dir / "pmc_fulltext"
    for path in sorted(ft_dir.glob("*_FULL_CONTEXT.md")):
        pmid = path.name.split("_")[0]
        size = path.stat().st_size
        if pmid.isdigit() and size >= min_bytes:
            out.append((pmid, path, size))
    return out


def _title(run_dir: Path, pmid: str) -> str:
    aj = run_dir / "abstract_json" / f"{pmid}.json"
    if aj.exists():
        try:
            meta = json.loads(aj.read_text(encoding="utf-8")).get("metadata") or {}
            if isinstance(meta, dict):
                return str(meta.get("title") or "").strip()
        except (OSError, json.JSONDecodeError):
            pass
    return ""


def _protocol_md(gene: str, n: int) -> str:
    return f"""# Manual curation protocol — {gene} gold standard ({n} papers)

## What this is & why
You are building the **ground-truth answer key** for {n} papers so we can measure
how accurately an automated pipeline extracted {gene} variants from them. For each
paper you record every {gene} variant and its patient counts. We then compare the
pipeline against your answers. **Work independently** — you are NOT given the
pipeline's answers, and you should not try to find them. Honest, careful reading
is the whole point.

## What you need
Nothing but this folder and a spreadsheet program (Excel or Google Sheets). The
full text of every paper is included in `papers/<pmid>.md`, so you do **not** need
journal access or the internet. (A PubMed link is in `manifest.csv` if you want to
see figures/tables in the original — optional.)

## The task (per paper)
1. Open `papers/<pmid>.md` and read it (skim intro/methods, read results + tables).
2. Find every **{gene} variant** reported in patients/families in that paper.
3. For each distinct variant, add **one row** to `curation_template.csv`:

| column | what to put |
|---|---|
| `pmid` | the paper's PubMed ID (already filled in for you) |
| `variant` | the variant **as written** in the paper. Prefer HGVS: cDNA `c.5946delT` or protein `p.Ser1982fs`. Copy exactly; don't normalize. |
| `germline_or_somatic` | `germline` if inherited/constitutional (blood, saliva, family, carrier, proband, germline testing); `somatic` if tumor-only (tumor tissue, ctDNA, cell line, LOH); `unknown` if unclear. |
| `carriers` | number of **individuals** reported carrying this variant in this paper. If only families are given, put the family count and say so in the note. Blank if not reported. |
| `affected` | of those carriers, how many **had cancer** (breast/ovarian/the studied phenotype). Blank if not reported. |
| `unaffected` | of those carriers, how many were **cancer-free**. Blank if not reported. |
| `evidence_note` | where you found it (e.g. "Table 2") + a short quote. |

## Rules that matter
- **One row per distinct variant.** If a variant appears in two cohorts in the
  same paper, sum the carriers into one row (note it).
- **Germline is the target.** Still record `somatic` variants (we use them to
  calibrate), but mark them clearly — they are not heritable carriers.
- **Never guess a number.** If a count isn't stated, leave it blank. A blank
  means "not reported", which is different from 0.
- **No {gene} variant in the paper?** Keep the paper's row and put `NONE` in
  `variant`. (This matters — it tells us the pipeline shouldn't have found any.)
- **Don't normalize notation.** Copy the variant string verbatim; we handle
  matching. If the paper gives both cDNA and protein, put both (one in `variant`,
  the other in the note).

## Worked examples
| pmid | variant | germline_or_somatic | carriers | affected | unaffected | evidence_note |
|---|---|---|---|---|---|---|
| 12345678 | c.5946delT | germline | 14 | 9 | 5 | Table 2; "14 carriers, 9 with breast cancer" |
| 12345678 | p.Cys1365Tyr | germline | 1 | 1 |  | Case report, Results para 3 |
| 99999999 | NONE |  |  |  |  | Review article, no patient variants |

## When you're done
Save the spreadsheet as `curation_template_FILLED.csv` (keep CSV format) and send
it back. Budget ~3–4 days for {n} papers. Thank you!
"""


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--run-dir",
        required=True,
        type=Path,
        help="GVF run dir with pmc_fulltext/ + abstract_json/",
    )
    ap.add_argument("--gene", required=True)
    ap.add_argument("--n", type=int, default=50)
    ap.add_argument(
        "--seed", type=int, default=42, help="Deterministic sample seed (record it)."
    )
    ap.add_argument(
        "--min-fulltext-bytes", type=int, default=3000, help="Full-text proxy floor."
    )
    ap.add_argument(
        "--out", required=True, type=Path, help="Packet output dir (also zipped)."
    )
    ap.add_argument("--no-zip", action="store_true")
    args = ap.parse_args()

    run_dir = args.run_dir.expanduser().resolve()
    gene = args.gene.upper()
    pool = _full_text_papers(run_dir, args.min_fulltext_bytes)
    if len(pool) < args.n:
        print(
            f"warning: only {len(pool)} full-text papers >= {args.min_fulltext_bytes}B; "
            f"sampling all of them (< requested {args.n}).",
            file=sys.stderr,
        )
    sample = random.Random(args.seed).sample(pool, min(args.n, len(pool)))
    sample.sort(key=lambda r: r[0])

    out = args.out.expanduser().resolve()
    papers = out / "papers"
    if out.exists():
        shutil.rmtree(out)
    papers.mkdir(parents=True)

    (out / "PROTOCOL.md").write_text(_protocol_md(gene, len(sample)), encoding="utf-8")

    with (out / "manifest.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pmid", "title", "pubmed_url", "fulltext_file", "n_bytes"])
        for pmid, path, size in sample:
            shutil.copy2(path, papers / f"{pmid}.md")
            w.writerow(
                [
                    pmid,
                    _title(run_dir, pmid),
                    f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                    f"papers/{pmid}.md",
                    size,
                ]
            )

    with (out / "curation_template.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(TEMPLATE_COLUMNS)
        for pmid, _, _ in sample:  # one prefilled starter row per paper
            w.writerow([pmid] + [""] * (len(TEMPLATE_COLUMNS) - 1))

    meta = {
        "gene": gene,
        "run_dir": str(run_dir),
        "n_requested": args.n,
        "n_sampled": len(sample),
        "seed": args.seed,
        "min_fulltext_bytes": args.min_fulltext_bytes,
        "full_text_pool_size": len(pool),
        "blinded": True,
        "pmids": [pmid for pmid, _, _ in sample],
    }
    (out / "packet_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    archive = ""
    if not args.no_zip:
        archive = shutil.make_archive(
            str(out), "zip", root_dir=out.parent, base_dir=out.name
        )

    print(f"Packet: {out}  ({len(sample)} papers from a pool of {len(pool)} full-text)")
    print(f"  PROTOCOL.md, manifest.csv, curation_template.csv, papers/<pmid>.md")
    if archive:
        print(f"  zip -> {archive}")
    print("Hand off the zip; have the assistant return curation_template_FILLED.csv.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
