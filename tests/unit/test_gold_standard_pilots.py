"""Offline checks for gold-standard pilot generation."""

from __future__ import annotations

import csv
import subprocess
import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "build_gold_standard_pilots.py"
GOLD_DIR = REPO_ROOT / "gene_variant_fetcher_gold_standard"


class GoldStandardPilotTests(unittest.TestCase):
    def test_builds_recall_shaped_pilot_package(self) -> None:
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp) / "pilots"
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--gold-dir",
                    str(GOLD_DIR),
                    "--out-dir",
                    str(out_dir),
                    "--genes",
                    "KCNH2,KCNQ1,SCN5A",
                    "--pmids-per-gene",
                    "4",
                ],
                check=True,
                cwd=REPO_ROOT,
            )

            with (out_dir / "pilot_pmids.csv").open(newline="") as handle:
                pilot_pmids = list(csv.DictReader(handle))
            self.assertEqual(
                {row["gene"] for row in pilot_pmids}, {"KCNH2", "KCNQ1", "SCN5A"}
            )
            self.assertLessEqual(max(count_gene(pilot_pmids).values()), 4)

            for gene in ("KCNH2", "KCNQ1", "SCN5A"):
                recall_path = out_dir / "normalized" / f"{gene}_recall_input.csv"
                self.assertTrue(recall_path.exists())
                with recall_path.open(newline="") as handle:
                    reader = csv.DictReader(handle)
                    self.assertEqual(
                        reader.fieldnames,
                        ["variant", "pmid", "carriers", "affected", "unaffected"],
                    )
                    self.assertGreater(len(list(reader)), 0)


def count_gene(rows: list[dict[str, str]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["gene"]] = counts.get(row["gene"], 0) + 1
    return counts


if __name__ == "__main__":
    unittest.main()
