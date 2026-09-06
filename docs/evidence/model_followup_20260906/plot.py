"""Primary opened-source diagnostic only; the repeated DOC pair is excluded."""

import csv
import json
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from client import HERE, write, sha

s = json.loads((HERE / "source_scores.json").read_text())
data = []
for arm in ["astra_low", "astra_medium", "grok_low"]:
    rr = [
        r
        for r in s["results"]
        if r["arm"] == arm and r["packet"] != "myb204_structured"
    ]
    cs = [c for r in rr for c in r["count_comparison"]]
    data.append(
        dict(
            arm=arm,
            people=sum(r["reference_people"] for r in rr),
            source_bound=sum(r["exact_people"] for r in rr),
            count_fields=len(cs),
            count_exact=sum(c["exact"] for c in cs),
            positive_reference_fields=sum(c["reference_positive"] for c in cs),
            positive_reference_exact=sum(
                c["reference_positive"] and c["exact"] for c in cs
            ),
            api_proxy=sum(r["api_proxy_usd"] for r in rr),
        )
    )
labels = ["Astra low", "Astra medium", "Grok 4.6 low"]
colors = ["#3978ad", "#6a519f", "#de883f"]
fig, axes = plt.subplots(1, 3, figsize=(10.5, 4.2))
fig.subplots_adjust(left=0.07, right=0.98, bottom=0.23, top=0.75, wspace=0.4)
for ax, key, title, total in zip(
    axes,
    ["count_exact", "source_bound", "api_proxy"],
    [
        "Derived packet counts\nexact vs source",
        "People passing all\nstrict provenance checks",
        "API usage proxy\nfor seven packets ($)",
    ],
    [168, 93, None],
):
    vals = [d[key] for d in data]
    ax.bar(range(3), vals, color=colors, width=0.65)
    ax.set_xticks(range(3), labels, rotation=25, ha="right", fontsize=9)
    ax.set_title(title, fontsize=11, pad=12)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.18)
    ax.set_axisbelow(True)
    ax.set_ylim(0, (total if total else max(vals)) * 1.24)
    for i, v in enumerate(vals):
        ax.text(
            i,
            v + (total if total else max(vals)) * 0.035,
            f"{v}/{total}" if total else f"${v:.3f}",
            ha="center",
            fontsize=10,
        )
fig.suptitle("Bounded reading of opened source packets", fontsize=16, y=0.97)
fig.text(
    0.07,
    0.02,
    "Four patient-bearing papers (93 people), three roster-negative packets. Packet endpoints and cohorts differ from whole-paper gold.\nSource-literal totals are derived in code; integer agreement does not certify provenance. The repeated DOC representation is excluded.",
    fontsize=8,
    color="#444444",
)
for suffix in ["png", "svg", "pdf"]:
    fig.savefig(HERE / ("bounded_roster_comparison." + suffix), dpi=180)
with (HERE / "bounded_roster_comparison.csv").open("w") as f:
    w = csv.DictWriter(f, fieldnames=list(data[0]))
    w.writeheader()
    w.writerows(data)
write(
    HERE / "bounded_roster_comparison.json",
    dict(
        source_scores_sha256=sha(HERE / "source_scores.json"),
        data=data,
        classification="Primary source diagnostic, not population recall or benchmark A/U accuracy",
        duplicate_representation_excluded="myb204_structured",
    ),
)
