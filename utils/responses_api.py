"""Small, provider-aware client for bounded Responses API tool loops.

This module intentionally does not depend on an OpenAI SDK.  The extraction
experiment needs the wire-level Responses API contract (including
``previous_response_id`` continuations) while remaining easy to exercise with
an in-memory transport in unit tests.

The client does not estimate money.  It reports provider token usage verbatim
so callers can apply a versioned price table outside this module.
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict, dataclass, field
from typing import Any, Callable, Literal, Mapping, Sequence

import requests


REASONING_EFFORTS = ("none", "low", "medium", "high", "xhigh", "max")

Provider = Literal["azure", "openai"]
ToolHandler = Callable[[dict[str, Any]], Any]
Transport = Callable[..., Any]

_RETRYABLE_STATUS_CODES = frozenset({408, 409, 429})


@dataclass(frozen=True)
class ResponsesLoopLimits:
    """Hard limits applied to one tool loop."""

    max_rounds: int = 8
    max_tool_calls: int = 24
    max_tool_output_chars: int = 30_000
    max_total_tool_output_chars: int = 120_000

    def __post_init__(self) -> None:
        if self.max_rounds < 1:
            raise ValueError("max_rounds must be at least 1")
        if self.max_tool_calls < 0:
            raise ValueError("max_tool_calls must be non-negative")
        if self.max_tool_output_chars < 0:
            raise ValueError("max_tool_output_chars must be non-negative")
        if self.max_total_tool_output_chars < 0:
            raise ValueError("max_total_tool_output_chars must be non-negative")

    def to_dict(self) -> dict[str, int]:
        return asdict(self)


@dataclass
class ResponsesUsage:
    """Token usage accumulated across every response in a loop."""

    input_tokens: int = 0
    cached_input_tokens: int = 0
    cache_write_input_tokens: int = 0
    output_tokens: int = 0
    reasoning_tokens: int = 0
    total_tokens: int = 0

    def to_dict(self) -> dict[str, int]:
        return asdict(self)

    def add_response_usage(self, usage: Mapping[str, Any] | None) -> None:
        if not isinstance(usage, Mapping):
            return

        input_tokens = _nonnegative_int(
            usage.get("input_tokens", usage.get("prompt_tokens", 0))
        )
        output_tokens = _nonnegative_int(
            usage.get("output_tokens", usage.get("completion_tokens", 0))
        )
        input_details = usage.get("input_tokens_details")
        if not isinstance(input_details, Mapping):
            input_details = usage.get("prompt_tokens_details")
        if not isinstance(input_details, Mapping):
            input_details = {}
        output_details = usage.get("output_tokens_details")
        if not isinstance(output_details, Mapping):
            output_details = usage.get("completion_tokens_details")
        if not isinstance(output_details, Mapping):
            output_details = {}

        cached = _first_nonnegative_int(
            input_details,
            ("cached_tokens", "cached_input_tokens", "cache_read_input_tokens"),
        )
        if cached == 0:
            cached = _first_nonnegative_int(
                usage,
                ("cached_input_tokens", "cache_read_input_tokens"),
            )
        cache_write = _first_nonnegative_int(
            input_details,
            (
                "cache_write_tokens",
                "cache_write_input_tokens",
                "cache_creation_input_tokens",
            ),
        )
        if cache_write == 0:
            cache_write = _first_nonnegative_int(
                usage,
                ("cache_write_input_tokens", "cache_creation_input_tokens"),
            )
        reasoning = _first_nonnegative_int(
            output_details,
            ("reasoning_tokens",),
        )
        if reasoning == 0:
            reasoning = _nonnegative_int(usage.get("reasoning_tokens", 0))

        total = _nonnegative_int(usage.get("total_tokens", 0))
        if total == 0:
            total = input_tokens + output_tokens

        self.input_tokens += input_tokens
        self.cached_input_tokens += cached
        self.cache_write_input_tokens += cache_write
        self.output_tokens += output_tokens
        self.reasoning_tokens += reasoning
        self.total_tokens += total


@dataclass
class ResponsesTelemetry:
    """Operational telemetry accumulated across a tool loop."""

    # api_calls counts HTTP attempts; rounds counts logical Responses requests.
    api_calls: int = 0
    rounds: int = 0
    retries: int = 0
    # HTTP attempt latency; retry backoff is intentionally excluded.
    latency_seconds: float = 0.0
    tool_latency_seconds: float = 0.0
    tool_calls_requested: int = 0
    tool_calls_executed: int = 0
    tool_errors: int = 0
    tool_output_chars: int = 0
    truncated_tool_outputs: int = 0

    def to_dict(self) -> dict[str, int | float]:
        return asdict(self)


@dataclass
class ResponsesResult:
    """Final structured result or a bounded, inspectable failure."""

    ok: bool
    status: str
    value: Any | None
    output_text: str | None
    error: str | None
    response_id: str | None
    usage: ResponsesUsage
    telemetry: ResponsesTelemetry
    raw_responses: list[dict[str, Any]] = field(default_factory=list)

    def to_dict(self, *, include_raw_responses: bool = False) -> dict[str, Any]:
        result = {
            "ok": self.ok,
            "status": self.status,
            "value": self.value,
            "output_text": self.output_text,
            "error": self.error,
            "response_id": self.response_id,
            "usage": self.usage.to_dict(),
            "telemetry": self.telemetry.to_dict(),
        }
        if include_raw_responses:
            result["raw_responses"] = self.raw_responses
        return result


class _RequestFailure(Exception):
    def __init__(self, kind: str, message: str) -> None:
        super().__init__(message)
        self.kind = kind


class ResponsesAPIClient:
    """Execute Responses API requests and bounded local function tools.

    ``transport`` has the same call shape as ``requests.post``.  Supplying a
    fake transport makes the complete request/continuation loop testable
    without network access.
    """

    def __init__(
        self,
        provider: Provider,
        base_url: str,
        api_key: str,
        *,
        transport: Transport | None = None,
        timeout_s: float = 1_200,
        max_retries: int = 2,
        backoff_s: float = 1.0,
        sleep_fn: Callable[[float], None] = time.sleep,
        clock: Callable[[], float] = time.monotonic,
    ) -> None:
        if provider not in ("azure", "openai"):
            raise ValueError("provider must be 'azure' or 'openai'")
        if not base_url or not base_url.strip():
            raise ValueError("base_url is required")
        if not api_key or not api_key.strip():
            raise ValueError("api_key is required")
        if timeout_s <= 0:
            raise ValueError("timeout_s must be positive")
        if max_retries < 0:
            raise ValueError("max_retries must be non-negative")
        if backoff_s < 0:
            raise ValueError("backoff_s must be non-negative")

        self.provider = provider
        self.url = _responses_url(provider, base_url)
        self.api_key = api_key.strip()
        self.transport = transport or requests.post
        self.timeout_s = timeout_s
        self.max_retries = max_retries
        self.backoff_s = backoff_s
        self.sleep_fn = sleep_fn
        self.clock = clock

    @classmethod
    def azure(
        cls,
        base_url: str,
        api_key: str,
        **kwargs: Any,
    ) -> "ResponsesAPIClient":
        """Build a client for an Azure Foundry OpenAI v1 endpoint."""

        return cls("azure", base_url, api_key, **kwargs)

    @classmethod
    def openai(
        cls,
        api_key: str,
        base_url: str = "https://api.openai.com/v1",
        **kwargs: Any,
    ) -> "ResponsesAPIClient":
        """Build a client for the public OpenAI v1 endpoint."""

        return cls("openai", base_url, api_key, **kwargs)

    def run_tool_loop(
        self,
        *,
        model: str,
        reasoning_effort: str,
        instructions: str,
        initial_input: str | Sequence[Mapping[str, Any]],
        tools: Sequence[Mapping[str, Any]],
        tool_handlers: Mapping[str, ToolHandler],
        text: Mapping[str, Any] | None = None,
        max_output_tokens: int | None = None,
        limits: ResponsesLoopLimits | None = None,
        extra_body: Mapping[str, Any] | None = None,
    ) -> ResponsesResult:
        """Run a Responses function-tool loop and parse final ``output_text``.

        Every continuation repeats ``instructions``, ``tools``, and ``text``
        because Responses does not carry instructions forward automatically
        when ``previous_response_id`` is used.  Tool handler failures become
        bounded, matching ``function_call_output`` items instead of escaping
        the loop or invoking an unregistered function.
        """

        if reasoning_effort not in REASONING_EFFORTS:
            allowed = ", ".join(REASONING_EFFORTS)
            raise ValueError(f"reasoning_effort must be one of: {allowed}")
        if not model or not model.strip():
            raise ValueError("model is required")
        if max_output_tokens is not None and max_output_tokens < 1:
            raise ValueError("max_output_tokens must be positive when provided")

        active_limits = limits or ResponsesLoopLimits()
        usage = ResponsesUsage()
        telemetry = ResponsesTelemetry()
        raw_responses: list[dict[str, Any]] = []
        response_id: str | None = None
        seen_call_ids: set[str] = set()
        next_input: str | Sequence[Mapping[str, Any]] = initial_input

        common_body: dict[str, Any] = dict(extra_body or {})
        common_body.update(
            {
                "model": _bare_model_name(model),
                "instructions": instructions,
                "reasoning": {"effort": reasoning_effort},
                "tools": [dict(tool) for tool in tools],
                # Continuation by response id requires server-side response state.
                "store": True,
            }
        )
        if text is not None:
            common_body["text"] = dict(text)
        else:
            common_body.pop("text", None)
        if max_output_tokens is not None:
            common_body["max_output_tokens"] = max_output_tokens
        else:
            common_body.pop("max_output_tokens", None)

        for round_index in range(active_limits.max_rounds):
            body = dict(common_body)
            body["input"] = next_input
            if response_id is not None:
                body["previous_response_id"] = response_id
            else:
                body.pop("previous_response_id", None)

            telemetry.rounds += 1
            try:
                response = self._post_json(body, telemetry)
            except _RequestFailure as exc:
                return _failure(
                    status=exc.kind,
                    error=str(exc),
                    response_id=response_id,
                    usage=usage,
                    telemetry=telemetry,
                    raw_responses=raw_responses,
                )

            raw_responses.append(response)
            usage.add_response_usage(response.get("usage"))
            raw_id = response.get("id")
            if isinstance(raw_id, str) and raw_id:
                response_id = raw_id

            response_status = response.get("status")
            if response_status == "incomplete" or response.get("incomplete_details"):
                details = response.get("incomplete_details")
                return _failure(
                    status="incomplete",
                    error=f"Responses API output was incomplete: {details!r}",
                    response_id=response_id,
                    usage=usage,
                    telemetry=telemetry,
                    raw_responses=raw_responses,
                )
            if response_status in {"failed", "cancelled"}:
                return _failure(
                    status=str(response_status),
                    error=f"Responses API returned status {response_status!r}",
                    response_id=response_id,
                    usage=usage,
                    telemetry=telemetry,
                    raw_responses=raw_responses,
                )

            refusal = _response_refusal(response)
            if refusal is not None:
                return _failure(
                    status="refusal",
                    error=refusal,
                    response_id=response_id,
                    usage=usage,
                    telemetry=telemetry,
                    raw_responses=raw_responses,
                )

            function_calls = _response_function_calls(response)
            if not function_calls:
                output_text = _response_output_text(response)
                if output_text is None or not output_text.strip():
                    return _failure(
                        status="no_final_output",
                        error="Responses API returned no final output_text",
                        response_id=response_id,
                        usage=usage,
                        telemetry=telemetry,
                        raw_responses=raw_responses,
                    )
                try:
                    value = json.loads(output_text)
                except (TypeError, json.JSONDecodeError) as exc:
                    return ResponsesResult(
                        ok=False,
                        status="invalid_json",
                        value=None,
                        output_text=output_text,
                        error=f"Final output_text was not valid JSON: {exc}",
                        response_id=response_id,
                        usage=usage,
                        telemetry=telemetry,
                        raw_responses=raw_responses,
                    )
                return ResponsesResult(
                    ok=True,
                    status="completed",
                    value=value,
                    output_text=output_text,
                    error=None,
                    response_id=response_id,
                    usage=usage,
                    telemetry=telemetry,
                    raw_responses=raw_responses,
                )

            telemetry.tool_calls_requested += len(function_calls)
            batch_call_ids: set[str] = set()
            for function_call in function_calls:
                call_id = function_call.get("call_id")
                if not isinstance(call_id, str) or not call_id:
                    return _failure(
                        status="protocol_error",
                        error="Function call did not include a non-empty call_id",
                        response_id=response_id,
                        usage=usage,
                        telemetry=telemetry,
                        raw_responses=raw_responses,
                    )
                if call_id in seen_call_ids or call_id in batch_call_ids:
                    return _failure(
                        status="protocol_error",
                        error=f"Responses API repeated function call_id {call_id!r}",
                        response_id=response_id,
                        usage=usage,
                        telemetry=telemetry,
                        raw_responses=raw_responses,
                    )
                batch_call_ids.add(call_id)
            if telemetry.tool_calls_requested > active_limits.max_tool_calls:
                return _failure(
                    status="tool_call_limit_exceeded",
                    error=(
                        "Model requested more than "
                        f"{active_limits.max_tool_calls} function calls"
                    ),
                    response_id=response_id,
                    usage=usage,
                    telemetry=telemetry,
                    raw_responses=raw_responses,
                )
            if round_index + 1 >= active_limits.max_rounds:
                return _failure(
                    status="round_limit_exceeded",
                    error=(
                        "Model requested tools on the final permitted round "
                        f"({active_limits.max_rounds})"
                    ),
                    response_id=response_id,
                    usage=usage,
                    telemetry=telemetry,
                    raw_responses=raw_responses,
                )
            if response_id is None:
                return _failure(
                    status="protocol_error",
                    error="Tool-call response did not include an id for continuation",
                    response_id=response_id,
                    usage=usage,
                    telemetry=telemetry,
                    raw_responses=raw_responses,
                )

            next_items: list[dict[str, str]] = []
            for function_call in function_calls:
                call_id = function_call.get("call_id")
                # Validated for the complete batch before any handler executes.
                assert isinstance(call_id, str)
                raw_tool_output = self._execute_tool(
                    function_call,
                    tool_handlers,
                    telemetry,
                )
                remaining = max(
                    0,
                    active_limits.max_total_tool_output_chars
                    - telemetry.tool_output_chars,
                )
                bounded_output, truncated = _bound_text(
                    raw_tool_output,
                    min(active_limits.max_tool_output_chars, remaining),
                )
                if truncated:
                    telemetry.truncated_tool_outputs += 1
                telemetry.tool_output_chars += len(bounded_output)
                next_items.append(
                    {
                        "type": "function_call_output",
                        "call_id": call_id,
                        "output": bounded_output,
                    }
                )
            next_input = next_items
            seen_call_ids.update(batch_call_ids)

        # The loop always returns from inside, but keep a defensive terminal.
        return _failure(
            status="round_limit_exceeded",
            error=f"Tool loop exceeded {active_limits.max_rounds} rounds",
            response_id=response_id,
            usage=usage,
            telemetry=telemetry,
            raw_responses=raw_responses,
        )

    def _execute_tool(
        self,
        function_call: Mapping[str, Any],
        handlers: Mapping[str, ToolHandler],
        telemetry: ResponsesTelemetry,
    ) -> str:
        name = function_call.get("name")
        if not isinstance(name, str) or not name:
            telemetry.tool_errors += 1
            return _tool_error("malformed_tool_call", "Function name is missing")

        handler = handlers.get(name)
        if handler is None:
            telemetry.tool_errors += 1
            return _tool_error("unknown_tool", f"No handler registered for {name!r}")

        raw_arguments = function_call.get("arguments")
        if not isinstance(raw_arguments, str):
            telemetry.tool_errors += 1
            return _tool_error(
                "malformed_arguments",
                "Function arguments must be a JSON string",
            )
        try:
            arguments = json.loads(raw_arguments)
        except json.JSONDecodeError as exc:
            telemetry.tool_errors += 1
            return _tool_error(
                "malformed_arguments",
                f"Function arguments were not valid JSON: {exc.msg}",
            )
        if not isinstance(arguments, dict):
            telemetry.tool_errors += 1
            return _tool_error(
                "malformed_arguments",
                "Function arguments must decode to a JSON object",
            )

        tool_started = self.clock()
        try:
            result = handler(arguments)
        except Exception as exc:  # Tool failures are data, not loop failures.
            telemetry.tool_calls_executed += 1
            telemetry.tool_errors += 1
            return _tool_error(
                "tool_execution_error",
                f"{type(exc).__name__}: {exc}",
            )
        finally:
            telemetry.tool_latency_seconds += max(0.0, self.clock() - tool_started)
        telemetry.tool_calls_executed += 1
        try:
            return _serialize_tool_output(result)
        except (TypeError, ValueError) as exc:
            telemetry.tool_errors += 1
            return _tool_error(
                "unserializable_tool_output",
                f"Tool output was not JSON serializable: {type(exc).__name__}",
            )

    def _post_json(
        self,
        body: Mapping[str, Any],
        telemetry: ResponsesTelemetry,
    ) -> dict[str, Any]:
        headers = {"Content-Type": "application/json"}
        if self.provider == "azure":
            headers["api-key"] = self.api_key
        else:
            headers["Authorization"] = f"Bearer {self.api_key}"

        for attempt in range(self.max_retries + 1):
            started = self.clock()
            telemetry.api_calls += 1
            try:
                response = self.transport(
                    self.url,
                    headers=headers,
                    json=dict(body),
                    timeout=self.timeout_s,
                )
            except (
                requests.RequestException,
                TimeoutError,
                ConnectionError,
                OSError,
            ) as exc:
                telemetry.latency_seconds += max(0.0, self.clock() - started)
                if attempt >= self.max_retries:
                    raise _RequestFailure(
                        "transport_error",
                        f"Responses API transport failed: {type(exc).__name__}: {exc}",
                    ) from exc
                telemetry.retries += 1
                self._sleep_before_retry(attempt)
                continue

            telemetry.latency_seconds += max(0.0, self.clock() - started)
            status_code = getattr(response, "status_code", None)
            if not isinstance(status_code, int):
                raise _RequestFailure(
                    "protocol_error",
                    "Responses transport returned no integer status_code",
                )
            if status_code < 200 or status_code >= 300:
                retryable = status_code in _RETRYABLE_STATUS_CODES or (
                    500 <= status_code < 600
                )
                if retryable and attempt < self.max_retries:
                    telemetry.retries += 1
                    self._sleep_before_retry(attempt)
                    continue
                response_text = str(getattr(response, "text", ""))[:500]
                raise _RequestFailure(
                    "http_error",
                    f"Responses API returned HTTP {status_code}: {response_text}",
                )

            try:
                payload = response.json()
            except Exception as exc:
                raise _RequestFailure(
                    "protocol_error",
                    "Responses API returned a non-JSON success response",
                ) from exc
            if not isinstance(payload, dict):
                raise _RequestFailure(
                    "protocol_error",
                    "Responses API success payload was not a JSON object",
                )
            return payload

        raise AssertionError("retry loop exhausted without returning")

    def _sleep_before_retry(self, attempt: int) -> None:
        delay = self.backoff_s * (2**attempt)
        if delay > 0:
            self.sleep_fn(delay)


def _responses_url(provider: Provider, base_url: str) -> str:
    normalized = base_url.strip().rstrip("/")
    if normalized.endswith("/responses") or "/responses?" in normalized:
        return normalized
    if provider == "openai":
        if normalized.endswith("/v1"):
            return f"{normalized}/responses"
        return f"{normalized}/v1/responses"
    if normalized.endswith("/openai/v1"):
        return f"{normalized}/responses?api-version=v1"
    return f"{normalized}/openai/v1/responses?api-version=v1"


def _bare_model_name(model: str) -> str:
    normalized = model.strip()
    for prefix in ("azure_ai/", "openai/"):
        if normalized.startswith(prefix):
            return normalized[len(prefix) :]
    return normalized


def _nonnegative_int(value: Any) -> int:
    if isinstance(value, bool):
        return 0
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return 0
    return max(0, parsed)


def _first_nonnegative_int(values: Mapping[str, Any], keys: Sequence[str]) -> int:
    for key in keys:
        if key in values:
            return _nonnegative_int(values.get(key))
    return 0


def _response_function_calls(response: Mapping[str, Any]) -> list[dict[str, Any]]:
    output = response.get("output")
    if not isinstance(output, list):
        return []
    return [
        dict(item)
        for item in output
        if isinstance(item, Mapping) and item.get("type") == "function_call"
    ]


def _response_output_text(response: Mapping[str, Any]) -> str | None:
    direct = response.get("output_text")
    if isinstance(direct, str):
        return direct

    chunks: list[str] = []
    output = response.get("output")
    if not isinstance(output, list):
        return None
    for item in output:
        if not isinstance(item, Mapping) or item.get("type") != "message":
            continue
        content = item.get("content")
        if not isinstance(content, list):
            continue
        for block in content:
            if not isinstance(block, Mapping) or block.get("type") != "output_text":
                continue
            text = block.get("text")
            if isinstance(text, str):
                chunks.append(text)
    return "".join(chunks) if chunks else None


def _response_refusal(response: Mapping[str, Any]) -> str | None:
    output = response.get("output")
    if not isinstance(output, list):
        return None
    for item in output:
        if not isinstance(item, Mapping) or item.get("type") != "message":
            continue
        content = item.get("content")
        if not isinstance(content, list):
            continue
        for block in content:
            if not isinstance(block, Mapping) or block.get("type") != "refusal":
                continue
            refusal = block.get("refusal", block.get("text", "Model refused request"))
            return str(refusal)
    return None


def _serialize_tool_output(value: Any) -> str:
    if isinstance(value, str):
        return value
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def _tool_error(error_type: str, message: str) -> str:
    return json.dumps(
        {"ok": False, "error": {"type": error_type, "message": message}},
        ensure_ascii=False,
        separators=(",", ":"),
    )


def _bound_text(value: str, limit: int) -> tuple[str, bool]:
    if len(value) <= limit:
        return value, False
    if limit <= 0:
        return "", True
    suffix = "...[truncated]"
    if limit <= len(suffix):
        return suffix[:limit], True
    return f"{value[: limit - len(suffix)]}{suffix}", True


def _failure(
    *,
    status: str,
    error: str,
    response_id: str | None,
    usage: ResponsesUsage,
    telemetry: ResponsesTelemetry,
    raw_responses: list[dict[str, Any]],
) -> ResponsesResult:
    return ResponsesResult(
        ok=False,
        status=status,
        value=None,
        output_text=None,
        error=error,
        response_id=response_id,
        usage=usage,
        telemetry=telemetry,
        raw_responses=raw_responses,
    )


__all__ = [
    "REASONING_EFFORTS",
    "ResponsesAPIClient",
    "ResponsesLoopLimits",
    "ResponsesResult",
    "ResponsesTelemetry",
    "ResponsesUsage",
]
