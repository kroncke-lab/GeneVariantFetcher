#!/usr/bin/env python3
"""Smoke-test configured Azure LLM deployments with LiteLLM.

This is intentionally live-network and opt-in. It makes tiny JSON calls through
the same LiteLLM `azure_ai/<deployment>` route used by GVF, then exits non-zero
if any configured deployment is missing, empty, or emits non-JSON.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from typing import Iterable

import litellm

from config.settings import get_settings
from utils.llm_utils import litellm_completion


REASONING_HINTS = ("kimi", "grok-4", "gpt-5", "gpt5", "deepseek", "o1", "o3")


def _unique(items: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for item in items:
        value = (item or "").strip()
        if not value or value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


def _configured_models(include_final: bool) -> list[str]:
    settings = get_settings()
    models: list[str] = [
        settings.get_tier2_model(),
        settings.get_table_router_model(),
        *settings.get_tier3_models(),
        *settings.get_tier3_adjudicator_models(),
        *settings.get_early_debate_models(),
    ]
    if include_final:
        models.extend(settings.get_final_adjudicator_models())
        models.append(settings.get_final_arbiter_model())
        models.append(settings.get_paper_final_check_model())
    return _unique(model for model in models if model.startswith("azure_ai/"))


def _max_tokens_for(model: str, requested: int) -> int:
    if any(hint in model.lower() for hint in REASONING_HINTS):
        return max(requested, 8192)
    return requested


def _smoke_one(
    model: str, timeout: int, max_tokens: int, attempts: int
) -> tuple[bool, str]:
    last_detail = ""
    for attempt in range(1, attempts + 1):
        start = time.time()
        try:
            # temperature=0 is dropped for gpt-5.6-* (only default 1 allowed).
            response = litellm_completion(
                model=model,
                messages=[
                    {"role": "system", "content": "Return strict JSON only."},
                    {
                        "role": "user",
                        "content": 'Return exactly {"ok": true, "answer": 7}.',
                    },
                ],
                temperature=0,
                max_tokens=_max_tokens_for(model, max_tokens),
                response_format={"type": "json_object"},
                timeout=timeout,
            )
        except Exception as exc:  # noqa: BLE001
            last_detail = f"{type(exc).__name__}: {str(exc)[:500]}"
            continue

        content = response.choices[0].message.content or ""
        finish_reason = getattr(response.choices[0], "finish_reason", None)
        elapsed = time.time() - start
        if not content.strip():
            last_detail = (
                f"empty content finish_reason={finish_reason} "
                f"elapsed={elapsed:.2f}s attempt={attempt}/{attempts}"
            )
            continue
        try:
            parsed = json.loads(content)
        except json.JSONDecodeError as exc:
            last_detail = f"non-JSON response: {exc}; content={content[:200]!r}"
            continue
        if parsed.get("ok") is not True:
            last_detail = f"JSON missing ok=true: {parsed!r}"
            continue
        suffix = f" after {attempt} attempts" if attempt > 1 else ""
        return True, f"finish_reason={finish_reason} elapsed={elapsed:.2f}s{suffix}"
    return False, last_detail


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "models",
        nargs="*",
        help=(
            "Optional explicit azure_ai/<deployment> model strings. Defaults to "
            "configured Azure routine/debate models."
        ),
    )
    parser.add_argument("--timeout", type=int, default=90)
    parser.add_argument("--max-tokens", type=int, default=1024)
    parser.add_argument("--attempts", type=int, default=3)
    parser.add_argument(
        "--include-final",
        action="store_true",
        help=(
            "Also include Azure-configured final-review models, including the "
            "independent per-paper final check."
        ),
    )
    args = parser.parse_args(argv)

    litellm.drop_params = True
    models = _unique(args.models or _configured_models(args.include_final))
    if not models:
        print("No Azure models configured to smoke-test.", file=sys.stderr)
        return 2

    failures = 0
    for model in models:
        ok, detail = _smoke_one(
            model,
            timeout=args.timeout,
            max_tokens=args.max_tokens,
            attempts=args.attempts,
        )
        status = "OK" if ok else "FAIL"
        print(f"{status} {model}: {detail}", flush=True)
        failures += 0 if ok else 1
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
