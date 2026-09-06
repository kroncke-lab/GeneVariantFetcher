"""Bounded follow-up transport. No SDK retries; original campaign is immutable."""

import fcntl
import hashlib
import json
import subprocess
import time
from pathlib import Path

import httpx

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / "model_followup_20260906/budget.json"
ORIGINAL = HERE.parent / "model_routing_20260906/budget.json"
BASE = "https://magen-api-2-resource.services.ai.azure.com/openai/v1"
LIMIT = 2.17
RATES = {"gpt-6-astra": (10, 50), "grok-4.6": (2, 6)}


def write(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def update(fn):
    with (HERE / "budget.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        path = HERE / "budget.json"
        if path.exists():
            data = json.loads(path.read_text())
        else:
            old = json.loads(OLD.read_text())
            original = json.loads(ORIGINAL.read_text())
            assert old["prior_budget_sha256"] == sha(ORIGINAL)
            accounted = old["prior_accounted_usd"] + sum(
                c["accounted_usd"] for c in old["calls"]
            )
            assert accounted + LIMIT <= original["api_ceiling_usd"]
            data = dict(
                limit_usd=LIMIT,
                prior_accounted_usd=accounted,
                prior_budget_sha256=sha(OLD),
                original_budget_sha256=sha(ORIGINAL),
                calls=[],
            )
        assert data["prior_budget_sha256"] == sha(OLD)
        assert data["original_budget_sha256"] == sha(ORIGINAL)
        result = fn(data)
        write(path, data)
        return result


def key():
    return subprocess.check_output(
        [
            "az",
            "cognitiveservices",
            "account",
            "keys",
            "list",
            "-g",
            "MAGen",
            "-n",
            "magen-api-2-resource",
            "--query",
            "key1",
            "-o",
            "tsv",
        ],
        text=True,
    ).strip()


def call(name, body, path="/chat/completions", timeout=55):
    """Reserve before dispatch; keep full reservation if usage is unavailable."""
    target = HERE / "responses" / (name + ".json")
    target.parent.mkdir(exist_ok=True)
    assert not target.exists(), name
    model = body["model"]
    inp, out = RATES[model]
    cap = next(
        body[k]
        for k in ("max_completion_tokens", "max_output_tokens", "max_tokens")
        if k in body
    )
    bound = len(json.dumps(body).encode()) + 2048
    assert bound < 272000
    # Text-only byte bound including overhead; Astra cache-write premium covered.
    reserve = (bound * inp * 1.25 + cap * out) / 1e6
    secret = key()

    def begin(data):
        assert not any(c["name"] == name for c in data["calls"])
        used = sum(c["accounted_usd"] for c in data["calls"])
        if used + reserve > data["limit_usd"]:
            raise RuntimeError("Follow-up API budget exhausted before dispatch")
        data["calls"].append(
            dict(
                name=name,
                model=model,
                status="reserved",
                reserved_usd=reserve,
                accounted_usd=reserve,
                started_unix=time.time(),
                cap=cap,
            )
        )

    update(begin)
    started = time.monotonic()
    result = dict(
        name=name,
        request=body,
        endpoint=BASE + path,
        retries=0,
        timeout_seconds=timeout,
        reserved_usd=reserve,
    )
    try:
        with httpx.Client(
            timeout=timeout, transport=httpx.HTTPTransport(retries=0)
        ) as client:
            response = client.post(BASE + path, headers={"api-key": secret}, json=body)
        result["http_status"] = response.status_code
        result["headers"] = {
            k: v
            for k, v in response.headers.items()
            if any(s in k for s in ("request-id", "ratelimit", "region", "retry-after"))
        }
        result["response"] = response.json()
        result["status"] = "returned" if response.is_success else "http_error"
    except Exception as exc:
        result.update(
            status="failed", error_type=type(exc).__name__, error=str(exc)[:300]
        )
    result["seconds"] = time.monotonic() - started
    usage = result.get("response", {}).get("usage")
    if usage:
        ni = usage.get("prompt_tokens", usage.get("input_tokens", 0))
        no = max(
            usage.get("completion_tokens", usage.get("output_tokens", 0)),
            usage.get("total_tokens", 0) - ni,
        )
        result["usage"] = usage
        result["api_proxy_usd"] = (ni * inp + no * out) / 1e6
        result["accounted_usd"] = (ni * inp * 1.25 + no * out) / 1e6
    else:
        result["accounted_usd"] = reserve
    write(target, result)

    def finish(data):
        row = next(c for c in data["calls"] if c["name"] == name)
        row.update({k: result[k] for k in ("status", "seconds", "accounted_usd")})
        row["response_sha256"] = sha(target)
        row["usage_known"] = bool(usage)
        if usage:
            row["api_proxy_usd"] = result["api_proxy_usd"]

    update(finish)
    print(
        json.dumps(
            {k: v for k, v in result.items() if k not in ("request", "response")}
        ),
        flush=True,
    )
    return result
