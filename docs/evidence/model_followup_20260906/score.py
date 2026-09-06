"""Grade locked outputs against source records, never use integer matches as evidence."""

from collections import Counter, defaultdict
import json
import re
import time
import statistics
from client import HERE, sha, write

COLS = [
    "sources",
    "variants",
    "family",
    "person",
    "genotype",
    "phenotype",
    "index_case",
    "endpoint_status",
]
FIELDS = COLS[1:]


def norm(v):
    if isinstance(v, list):
        return sorted(norm(x) for x in v)
    if isinstance(v, str):
        return (
            re.sub(r"\s+", "", v)
            .replace("−", "-")
            .replace("–", "-")
            .replace(":", "-")
            .lower()
        )
    return v


def key(r):
    return (norm(r.get("family")) or "", norm(str(r.get("person", ""))))


def aggregate(rows):
    groups = defaultdict(dict)
    for r in rows:
        if r.get("genotype") != "Y":
            continue
        for v in r["variants"]:
            ident = key(r)
            group = groups[norm(v)]
            if ident in group and group[ident] != r["endpoint_status"]:
                group[ident] = "unknown"
            else:
                group[ident] = r["endpoint_status"]
    return {
        v: {
            "carriers": len(people),
            "positive": sum(e == "positive" for e in people.values()),
            "negative": sum(e == "negative" for e in people.values()),
            "unknown": sum(e == "unknown" for e in people.values()),
        }
        for v, people in groups.items()
    }


def grade(packet, expected, record):
    out = dict(
        packet=packet["name"],
        reference_people=len(expected),
        reference_counts=aggregate(expected),
        status="not_dispatched",
        exact_people=0,
        missing_people=len(expected),
        extra_people=0,
        field_differences=[],
        source_binding_errors=[],
        endpoint_errors=[],
        duplicate_people=[],
        counts={},
    )
    out["count_comparison"] = [
        dict(
            variant=v,
            field=f,
            actual=None,
            reference=n,
            exact=False,
            reference_positive=n > 0,
        )
        for v, cc in out["reference_counts"].items()
        for f, n in cc.items()
        if f != "unknown"
    ]
    if record.get("status") == "not_dispatched":
        return out
    path = HERE / record["path"]
    assert sha(path) == record["sha256"]
    d = json.loads(path.read_text())
    out.update(
        status=d["status"],
        seconds=d["seconds"],
        api_proxy_usd=d.get("api_proxy_usd"),
        accounted_usd=d["accounted_usd"],
    )
    if d["status"] != "returned":
        return out
    choice = d.get("response", {}).get("choices", [{}])[0]
    out["finish_reason"] = choice.get("finish_reason")
    out["reasoning_tokens"] = (
        d.get("usage", {}).get("completion_tokens_details") or {}
    ).get("reasoning_tokens", 0)
    text = choice.get("message", {}).get("content") or ""
    if choice.get("finish_reason") != "stop":
        out["status"] = (
            "length" if choice.get("finish_reason") == "length" else "incomplete"
        )
        return out
    try:
        parsed = json.loads(text)
        assert parsed["columns"] == COLS
        assert isinstance(parsed["rows"], list) and isinstance(
            parsed["limitations"], list
        )
        rows = []
        for a in parsed["rows"]:
            assert isinstance(a, list) and len(a) == len(COLS)
            r = dict(zip(COLS, a))
            assert (
                isinstance(r["sources"], list)
                and r["sources"]
                and all(isinstance(x, str) for x in r["sources"])
            )
            assert (
                isinstance(r["variants"], list)
                and r["variants"]
                and all(isinstance(x, str) for x in r["variants"])
            )
            assert isinstance(r["person"], str) and r["person"]
            assert r["family"] is None or isinstance(r["family"], str)
            assert r["genotype"] in ["Y", "N", "?"] and r["endpoint_status"] in [
                "positive",
                "negative",
                "unknown",
            ]
            assert r["index_case"] is None or type(r["index_case"]) is bool
            assert isinstance(r["phenotype"], str)
            rows.append(r)
    except (KeyError, TypeError, AssertionError, ValueError) as e:
        out.update(status="schema_or_json_failure", error=str(e))
        return out
    out.update(
        status="scored",
        returned_people=len(rows),
        limitations=parsed["limitations"],
        counts=aggregate(rows),
    )
    ek = {key(r): r for r in expected}
    seen = Counter(key(r) for r in rows)
    out["duplicate_people"] = [list(k) for k, n in seen.items() if n > 1]
    out["missing_people"] = len(set(ek) - set(seen))
    out["extra_people"] = sum(k not in ek for k in seen)
    valid = []
    for r in rows:
        k = key(r)
        e = ek.get(k)
        if e is None:
            continue
        differences = [f for f in FIELDS if norm(r[f]) != norm(e[f])]
        missing_sources = set(e["sources"]) - set(r["sources"])
        unknown_sources = set(r["sources"]) - set(packet["source_units"])
        if differences:
            out["field_differences"].append(
                dict(person=list(k), fields=differences, expected=e, observed=r)
            )
        if "endpoint_status" in differences:
            out["endpoint_errors"].append(list(k))
        if missing_sources or unknown_sources:
            out["source_binding_errors"].append(
                dict(
                    person=list(k),
                    missing=sorted(missing_sources),
                    unknown=sorted(unknown_sources),
                )
            )
        if (
            not differences
            and not missing_sources
            and not unknown_sources
            and seen[k] == 1
        ):
            valid.append(r)
    out["exact_people"] = len(valid)
    out["roster_values_exact"] = not (
        out["missing_people"]
        or out["extra_people"]
        or out["field_differences"]
        or out["duplicate_people"]
    )
    out["complete_source_bound_roster"] = (
        out["roster_values_exact"] and not out["source_binding_errors"]
    )
    out["correct_empty_roster"] = not rows if not expected else None
    # A source-adjudication audit, not a general-purpose production acceptance gate.
    # Emit a complete variant packet only if every reference person for that variant
    # is exact and source-bound, and there are no extra/incorrect carriers assigned.
    accepted = {}
    for v, counts in out["reference_counts"].items():
        expected_members = {
            key(r)
            for r in expected
            if r["genotype"] == "Y" and v in norm(r["variants"])
        }
        valid_members = {
            key(r) for r in valid if r["genotype"] == "Y" and v in norm(r["variants"])
        }
        raw_members = {
            key(r) for r in rows if r["genotype"] == "Y" and v in norm(r["variants"])
        }
        if expected_members == valid_members == raw_members:
            accepted[v] = counts
    out["source_adjudicated_complete_variant_counts"] = accepted
    fields = ["carriers", "positive", "negative"]
    pairs = []
    for v in sorted(set(out["reference_counts"]) | set(out["counts"])):
        for f in fields:
            actual = out["counts"].get(v, {}).get(f)
            reference = out["reference_counts"].get(v, {}).get(f)
            pairs.append(
                dict(
                    variant=v,
                    field=f,
                    actual=actual,
                    reference=reference,
                    exact=actual == reference,
                    reference_positive=reference is not None and reference > 0,
                )
            )
    out["count_comparison"] = pairs
    return out


def main():
    locks = [HERE / "outputs_locked.json", HERE / "grok_outputs_locked.json"]
    assert all(p.exists() for p in locks)
    prepared = json.loads((HERE / "prepared.json").read_text())
    assert sha(HERE / "source_reference.json") == prepared["reference_sha256"]
    refs = json.loads((HERE / "source_reference.json").read_text())["rows"]
    results = []
    for lock in locks:
        for c in json.loads(lock.read_text())["cells"]:
            prefix = (
                "grok_low_"
                if c["name"].startswith("grok_")
                else (
                    "astra_low_"
                    if c["name"].startswith("astra_low_")
                    else "astra_medium_"
                )
            )
            name = c["name"][len(prefix) :]
            p = HERE / "packets" / (name + ".json")
            assert sha(p) == prepared["packets"][name]
            result = grade(json.loads(p.read_text()), refs[name], c)
            result["arm"] = prefix[:-1]
            results.append(result)
    summary = {}
    for arm in ["astra_low", "astra_medium", "grok_low"]:
        rr = [r for r in results if r["arm"] == arm]
        scored = [r for r in rr if r["status"] == "scored"]
        pos = [r for r in rr if r["reference_people"] > 0]
        comparisons = [c for r in rr for c in r.get("count_comparison", [])]
        summary[arm] = dict(
            packets=len(rr),
            status=dict(Counter(r["status"] for r in rr)),
            exact_people=sum(r["exact_people"] for r in rr),
            reference_people=sum(r["reference_people"] for r in rr),
            positive_packets=len(pos),
            exact_roster_values=sum(r.get("roster_values_exact", False) for r in pos),
            complete_source_bound_rosters=sum(
                r.get("complete_source_bound_roster", False) for r in pos
            ),
            correct_empty_rosters=sum(
                r.get("correct_empty_roster") is True for r in rr
            ),
            source_adjudicated_complete_variants=sum(
                len(r.get("source_adjudicated_complete_variant_counts", {})) for r in rr
            ),
            exact_raw_count_fields=sum(c["exact"] for c in comparisons),
            compared_raw_count_fields=len(comparisons),
            positive_reference_exact=sum(
                c["exact"] and c["reference_positive"] for c in comparisons
            ),
            positive_reference_fields=sum(c["reference_positive"] for c in comparisons),
            known_api_proxy_usd=sum(r.get("api_proxy_usd") or 0 for r in rr),
            accounted_usd=sum(r.get("accounted_usd", 0) for r in rr),
            median_seconds=statistics.median(
                r["seconds"] for r in rr if "seconds" in r
            ),
        )
    write(
        HERE / "source_scores.json",
        dict(
            graded_unix=time.time(),
            lock_hashes={p.name: sha(p) for p in locks},
            scorer_sha256=sha(__import__("pathlib").Path(__file__)),
            classification="Source-adjudication diagnostic; counts use packet-specific endpoints/cohorts, no gold reference read",
            summary=summary,
            results=results,
        ),
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
