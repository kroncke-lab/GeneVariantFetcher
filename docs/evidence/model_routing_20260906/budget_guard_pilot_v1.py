"""Campaign-only call accounting installed before importing the production CLI.

This is an experiment hook, not a default production policy. No prompts or keys
are persisted here; the ordinary tracer retains sanitized scientific evidence.
"""

from __future__ import annotations

import fcntl
import json
import os
import sys
import time
import uuid
from pathlib import Path

RATES = {
    "gpt-6-astra": (10.0, 50.0),
    "grok-4.6": (2.0, 6.0),
    "grok-4.3": (1.25, 2.5),
    "gpt-5.6-sol": (5.0, 30.0),
    "gpt-5.6-luna": (0.20, 1.20),
    "kimi-k2.6": (0.95, 4.0),
}


def update(fn):
    path = Path(os.environ["GVF_EXPERIMENT_BUDGET"])
    with path.with_suffix(".lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        data = json.loads(path.read_text())
        result = fn(data)
        temp = path.with_suffix(".tmp")
        temp.write_text(json.dumps(data, indent=2) + "\n")
        temp.replace(path)
        return result


def install():
    import litellm
    import utils.llm_trace as trace

    original = trace.capture_llm_call

    def guarded(**kwargs):
        model = kwargs["requested_model"]
        matches = [k for k in RATES if k in model.lower()]
        if len(matches) != 1 or "anthropic" in model.lower():
            raise RuntimeError(f"Campaign refuses unpriced model: {model}")
        rate_in, rate_out = RATES[matches[0]]
        request = kwargs["request"]
        # A byte bound deliberately overestimates normal text tokenization.
        # Payload bytes also overestimate base64 image tokens. Add overhead.
        bound_in = len(json.dumps(request, default=str).encode()) + 2048
        body = request.get("body", request)
        output_cap = next(
            (
                body[k]
                for k in ("max_output_tokens", "max_completion_tokens", "max_tokens")
                if body.get(k)
            ),
            None,
        )
        if not isinstance(output_cap, int) or output_cap <= 0:
            raise RuntimeError("Campaign call lacks an explicit output token cap")
        # Reserve long-context premium and room for two hidden transport retries.
        # The normal LiteLLM retry count is set to zero below; explicit caller
        # retries return through this wrapper and reserve independently.
        reserve = 3 * (2 * bound_in * rate_in + 2 * output_cap * rate_out) / 1e6
        identity = uuid.uuid4().hex
        context = trace.current_trace_context()

        def begin(data):
            committed = sum(
                x.get("accounted_usd", x["reserved_usd"]) for x in data["calls"]
            )
            if (
                reserve > data["per_call_reserve_limit_usd"]
                or committed + reserve + data["smoke_uncertainty_reserve_usd"]
                > data["api_ceiling_usd"]
            ):
                raise RuntimeError(
                    f"Campaign API budget prevents call to {model}; reservation ${reserve:.3f}"
                )
            row = {
                "id": identity,
                "model": model,
                "arm": os.environ.get("GVF_EXPERIMENT_ARM"),
                "gene": context.get("gene"),
                "pmid": context.get("pmid"),
                "stage": context.get("stage"),
                "status": "reserved",
                "reserved_usd": reserve,
                "output_cap": output_cap,
                "started_at_unix": time.time(),
            }
            data["calls"].append(row)

        update(begin)
        litellm.num_retries = 0
        try:
            response, ref = original(**kwargs)
        except BaseException:

            def failed(data):
                row = next(x for x in data["calls"] if x["id"] == identity)
                row["status"] = "failed_usage_unknown_reservation_retained"

            update(failed)
            raise
        payload = trace.json_safe(response)
        if not isinstance(payload, dict):
            payload = {}
        usage = payload.get("usage") or {}
        input_tokens = usage.get("input_tokens", usage.get("prompt_tokens"))
        output_tokens = usage.get("output_tokens", usage.get("completion_tokens"))

        def finish(data):
            row = next(x for x in data["calls"] if x["id"] == identity)
            row["usage"] = usage
            row["trace_id"] = ref.get("trace_id") if ref else None
            row["response_model"] = payload.get("model")
            if not isinstance(input_tokens, int) or not isinstance(output_tokens, int):
                row["status"] = "returned_usage_unknown_reservation_retained"
                return
            long = input_tokens > (272000 if "astra" in model else 200000)
            cost = (
                (2 if long else 1) * input_tokens * rate_in
                + (2 if long else 1) * output_tokens * rate_out
            ) / 1e6
            row.update(status="returned", api_proxy_usd=cost, accounted_usd=cost)
            if cost > reserve:
                raise RuntimeError(
                    "Returned usage exceeded campaign reservation; stop for accounting review"
                )

        update(finish)
        return response, ref

    # Importing utils initializes llm_utils before this hook. Replace those
    # already-bound aliases as well as future imports; otherwise primary calls
    # would bypass the ledger while vision calls were accounted.
    for module in tuple(sys.modules.values()):
        if module is not None and getattr(module, "capture_llm_call", None) is original:
            module.capture_llm_call = guarded
    litellm.num_retries = 0
