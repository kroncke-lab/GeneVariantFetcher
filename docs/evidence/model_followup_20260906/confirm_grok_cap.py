"""Counterbalanced parameter-only diagnosis; legacy cap uncertainty reserved first."""

from client import call, update, HERE, write
import time

write(
    HERE / "grok_confirmation_design.json",
    dict(
        prepared_unix=time.time(),
        order=["modern", "legacy", "modern"],
        prompt="Return exactly OK.",
        effort="low",
        cap=256,
        timeout=55,
        only_payload_change="max_completion_tokens versus max_tokens",
        legacy_unknown_output_reserve=128000,
    ),
)
body = dict(
    model="grok-4.6",
    messages=[dict(role="user", content="Return exactly OK.")],
    reasoning_effort="low",
)
call("grok_cap_aba_modern_before", {**body, "max_completion_tokens": 256})
extra = (128000 - 256) * 6 / 1e6
hold = "grok_legacy_confirmation_contingency"


def reserve(d):
    assert not any(c["name"] == hold for c in d["calls"])
    assert sum(c["accounted_usd"] for c in d["calls"]) + extra + 0.01 <= d["limit_usd"]
    d["calls"].append(
        dict(
            name=hold,
            kind="contingency_not_api_call",
            status="reserved",
            reserved_usd=extra,
            accounted_usd=extra,
            started_unix=time.time(),
        )
    )


update(reserve)
r = call("grok_cap_aba_legacy", {**body, "max_tokens": 256})


def finish(d):
    h = next(c for c in d["calls"] if c["name"] == hold)
    c = next(c for c in d["calls"] if c["name"] == "grok_cap_aba_legacy")
    h.update(
        status="transferred_or_released", accounted_usd=0, transferred_to=c["name"]
    )
    if not c.get("usage_known"):
        c["accounted_usd"] += extra
        c["legacy_cap_uncertainty"] = True


update(finish)
call("grok_cap_aba_modern_after", {**body, "max_completion_tokens": 256})
