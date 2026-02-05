#!/usr/bin/env python3
"""Analyze what's blocking the failed downloads."""

import json
from pathlib import Path

results_file = (
    Path(__file__).parent.parent / "tests/fixtures/pmids/full_test_results.json"
)

with open(results_file) as f:
    data = json.load(f)

failed = [r for r in data if not r["success"]]

print(f"Failed papers: {len(failed)}")
print()

# Categorize failures
categories = {"timeout": [], "captcha": [], "no_content": [], "other": []}

for r in failed:
    error = r["error"] or "Unknown"

    if "Timeout" in error:
        categories["timeout"].append(r)
    elif (
        "captcha" in error.lower()
        or "robot" in error.lower()
        or "cloudflare" in error.lower()
    ):
        categories["captcha"].append(r)
    elif "No usable content" in error:
        categories["no_content"].append(r)
    else:
        categories["other"].append(r)

print("=" * 60)
print("FAILURE BREAKDOWN")
print("=" * 60)

print(f"\n🕐 TIMEOUT (page took >30s): {len(categories['timeout'])}")
for r in categories["timeout"]:
    print(f"   {r['pmid']} ({r['publisher']}) - {r['doi']}")

print(f"\n🤖 CAPTCHA/BOT DETECTION: {len(categories['captcha'])}")
for r in categories["captcha"]:
    print(f"   {r['pmid']} ({r['publisher']}) - {r['doi']}")

print(f"\n📄 NO USABLE CONTENT: {len(categories['no_content'])}")
for r in categories["no_content"]:
    print(f"   {r['pmid']} ({r['publisher']}) - {r['doi']}")

print(f"\n❓ OTHER: {len(categories['other'])}")
for r in categories["other"]:
    print(f"   {r['pmid']} ({r['publisher']})")
    print(f"      Error: {r['error'][:80] if r['error'] else 'Unknown'}...")

print("\n" + "=" * 60)
print("ANALYSIS")
print("=" * 60)

if categories["timeout"]:
    print("""
TIMEOUT issues (most common):
- The 30s timeout isn't enough for Cloudflare challenges
- FIX: Increase timeout to 60s and add retry logic
- These are likely recoverable with longer timeouts
""")

if categories["captcha"]:
    print("""
CAPTCHA issues:
- Explicit bot detection requiring human verification
- FIX: Manual download via fetch_manager.py
- Or use abstract-only extraction
""")

if categories["no_content"]:
    print("""
NO CONTENT issues:
- Page loaded but article text couldn't be extracted
- Might be paywall, might be unusual page structure
- FIX: Check page structure, may need custom scraper
""")
