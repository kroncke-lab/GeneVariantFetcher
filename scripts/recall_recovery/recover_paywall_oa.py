"""Iterate recovery PMIDs through OARecoveryClient and write FULL_CONTEXT.md."""

import os
import sys
import json
import csv
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, "/Users/kronckbm/GitRepos/GeneVariantFetcher")

from dotenv import load_dotenv

load_dotenv()

from harvesting.oa_recovery import OARecoveryClient

PMC = Path("results/KCNH2/20260517_074737/pmc_fulltext")
RECOVERY = Path("/tmp/recovery_pmids.csv")
EMAIL = os.environ.get("NCBI_EMAIL") or "brett.kroncke@vanderbilt.edu"
NCBI_KEY = os.environ.get("NCBI_API_KEY")

client = OARecoveryClient(email=EMAIL, ncbi_api_key=NCBI_KEY)

import logging

logging.basicConfig(
    level=logging.WARNING, format="%(asctime)s %(name)s %(levelname)s %(message)s"
)
log = logging.getLogger("recover")
log.setLevel(logging.INFO)

successes = []
failures = []
with RECOVERY.open() as f:
    rdr = csv.DictReader(f)
    pmids = list(rdr)

print(f"Recovering {len(pmids)} PMIDs via OA paths...")
for i, row in enumerate(pmids):
    pmid = row["pmid"]
    print(
        f"[{i+1}/{len(pmids)}] PMID {pmid} (missing {row['missing_rows']} rows)... ",
        end="",
        flush=True,
    )
    try:
        result = client.recover(pmid)
    except Exception as e:
        print(f"ERROR {e}")
        failures.append((pmid, "exception", str(e)))
        continue
    if result.success:
        ctx_path = PMC / f"{pmid}_FULL_CONTEXT.md"
        backup = PMC / f"{pmid}_ABSTRACT_BACKUP.md"
        if ctx_path.exists() and not backup.exists():
            ctx_path.rename(backup)
        out = (
            f"# RECOVERED VIA {result.source}\n\n"
            f"PMID: {pmid}  Source: {result.source}  Chars: {result.n_chars}\n\n"
            + result.markdown
        )
        ctx_path.write_text(out, encoding="utf-8")
        cleaned_path = PMC / f"{pmid}_CLEANED.md"
        cleaned_path.write_text(result.markdown, encoding="utf-8")
        print(f"OK source={result.source} chars={result.n_chars}")
        successes.append((pmid, result.source, result.n_chars))
    else:
        attempts_brief = "; ".join(f"{a[0]}/{a[1]}" for a in result.attempts[:4])
        print(f"FAIL {result.error or '?'} attempts=[{attempts_brief}]")
        failures.append((pmid, result.error or "?", attempts_brief))

print()
print(f"=== Recovered {len(successes)} / {len(pmids)} ===")
with open("/tmp/recovery_results.csv", "w") as f:
    w = csv.writer(f)
    w.writerow(["pmid", "source", "chars"])
    for s in successes:
        w.writerow(s)
print("Top successes:")
for s in successes[:20]:
    print(f"  PMID {s[0]} via {s[1]}: {s[2]} chars")
print("\nTop failures:")
for fl in failures[:20]:
    print(f"  PMID {fl[0]}: {fl[1]} [{fl[2]}]")
