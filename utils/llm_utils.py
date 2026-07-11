"""
LLM utilities for making structured calls to language models.

This module provides a base class and utilities for making LLM calls with
consistent error handling, retry logic, and JSON response parsing.
"""

import json
import logging
import os
import threading
import time
from typing import Any, Dict, List, Optional

import litellm
from litellm import completion

litellm.drop_params = True

from .retry_utils import llm_retry

# Configure LiteLLM retry behavior to reduce excessive retries
litellm.num_retries = (
    2  # Reduce from default 6 to 2 (our @llm_retry handles additional retries)
)


def _resolve_litellm_request_timeout() -> int:
    raw = os.getenv("LLM_REQUEST_TIMEOUT_S", "").strip()
    if raw:
        try:
            timeout_s = int(raw)
            if timeout_s > 0:
                return timeout_s
        except ValueError:
            pass
        logger = logging.getLogger(__name__)
        logger.warning(
            "Ignoring invalid LLM_REQUEST_TIMEOUT_S=%r; expected positive integer",
            raw,
        )
    return 1200  # 20 minute default for large extractions.


litellm.request_timeout = _resolve_litellm_request_timeout()

logger = logging.getLogger(__name__)


# =============================================================================
# MODEL TOKEN LIMITS REGISTRY
# =============================================================================
# Completion (output) token limits with safety margins for various models.
# Format: model_pattern -> (actual_limit, safe_limit)
# The safe_limit accounts for overhead and prevents edge-case failures.

MODEL_TOKEN_LIMITS = {
    # Azure AI Foundry deployments — match by deployment name suffix so they
    # work regardless of the "azure_ai/" prefix LiteLLM expects.
    "deepseek-v4-pro": (16384, 15000),  # DeepSeek V4 Pro via Azure
    "gpt-5.3-codex": (16384, 15000),  # OpenAI GPT-5.3 Codex via Azure
    "gpt-5.5": (16384, 15000),  # OpenAI GPT-5.5 via Azure
    "gpt-5.4": (16384, 15000),  # OpenAI GPT-5.4 via Azure
    # GPT-5.6 reasoning models spend hidden reasoning tokens against the output
    # budget before emitting JSON, so at xhigh a 15000 cap truncates to empty on
    # long inputs. Give ample headroom (verified on the source-grounded summary).
    "gpt-5.6-sol": (128000, 64000),  # OpenAI GPT-5.6 Sol via Azure Foundry
    "gpt-5.6": (128000, 64000),  # OpenAI GPT-5.6 family via Azure
    # Generic GPT-5 family fallback. Lookup is longest-pattern-first, so the
    # specific gpt-5.x entries above still win; this only catches gpt-5.x ids we
    # haven't enumerated yet. Without it an unrecognized gpt-5.x falls to
    # DEFAULT_TOKEN_LIMIT (4000) and gets truncated mid hidden-reasoning.
    "gpt-5": (16384, 15000),  # OpenAI GPT-5 family fallback
    "kimi-k2": (16384, 15000),  # Moonshot Kimi K2 via Azure
    "grok-4.3": (16384, 15000),  # xAI Grok-4.3 reasoning via Azure
    "grok-4": (16384, 15000),  # xAI Grok-4.x reasoning via Azure (also covers 4.3)
    # OpenAI models - gpt-4o has 16384 completion limit
    "gpt-4o-mini": (16384, 15000),
    "gpt-4o": (16384, 15000),
    "gpt-4-turbo": (4096, 4000),
    "gpt-4-": (8192, 8000),  # Standard GPT-4 variants
    "gpt-3.5-turbo": (4096, 4000),
    "gpt-3.5": (4096, 4000),
    # Anthropic Claude models. Patterns are matched as substrings (case-insensitive)
    # against the LiteLLM model string, so "anthropic/claude-sonnet-4-6" matches
    # "claude-sonnet-4" and "claude-4-".
    # Claude 4.x family — 64k output tokens for Sonnet/Opus, 32k for Haiku per
    # Anthropic's published limits. Safe limits leave headroom for system overhead.
    "claude-opus-4-7": (64000, 60000),
    "claude-opus-4-6": (64000, 60000),
    "claude-opus-4": (64000, 60000),
    "claude-sonnet-4-6": (64000, 60000),
    "claude-sonnet-4": (64000, 60000),
    "claude-haiku-4-5": (32000, 30000),
    "claude-haiku-4": (32000, 30000),
    "claude-4": (64000, 60000),
    # Claude 3.x family
    "claude-3-opus": (4096, 4000),
    "claude-3-sonnet": (4096, 4000),
    "claude-3-haiku": (4096, 4000),
    "claude-3.5-sonnet": (8192, 8000),
    "claude-3-5-sonnet": (8192, 8000),
    "claude-3.5": (8192, 8000),
    "claude-2": (4096, 4000),
    # Google models
    "gemini-1.5-pro": (8192, 8000),
    "gemini-1.5-flash": (8192, 8000),
    "gemini-2.0-flash": (8192, 8000),
    "gemini-2.0": (8192, 8000),
    "gemini-pro": (8192, 8000),
    # Mistral models
    "mistral-large": (8192, 8000),
    "mistral-medium": (8192, 8000),
    "mistral-small": (8192, 8000),
}

# Default for unknown models (conservative)
DEFAULT_TOKEN_LIMIT = (4096, 4000)


# =============================================================================
# REASONING EFFORT
# =============================================================================
# OpenAI-style reasoning models accept a `reasoning_effort` knob.
# GPT-5 family historically: minimal|low|medium|high.
# GPT-5.6 on Azure/OpenAI also accepts none|xhigh (and OpenAI docs mention max
# as the deepest single-agent setting; Azure currently rejects "max" and takes
# "xhigh" as the deepest available). We treat "max" as an alias of "xhigh".
# Anthropic exposes the same capability through extended `thinking`, which
# additionally requires temperature=1 — so we do NOT route Claude through this
# helper yet. Grok-4-class models reason by default and ignore an effort param.
# litellm.drop_params=True would silently drop an unsupported value, so this
# allow-list exists to make the no-op explicit (and logged) rather than leaving
# callers believing effort was applied.
REASONING_EFFORT_LEVELS = ("none", "minimal", "low", "medium", "high", "xhigh")
REASONING_EFFORT_ALIASES = {"max": "xhigh"}
_REASONING_EFFORT_MODEL_HINTS = ("gpt-5", "gpt5", "o1", "o3", "o4-mini")


# =============================================================================
# AZURE FOUNDRY OPENAI v1 ENDPOINT
# =============================================================================
# Foundry portal "OpenAI-compatible" endpoints look like:
#   https://<resource>.services.ai.azure.com/openai/v1
# with deployment names passed as the model id (e.g. gpt-5.6-sol).
#
# LiteLLM's azure_ai/* provider rewrites services.ai.azure.com hosts to
# /models/chat/completions, which is the wrong path for this endpoint. When
# AZURE_AI_API_BASE ends with /openai/v1 we rewrite azure_ai/<deployment> to
# openai/<deployment> and pass api_base/api_key explicitly.


def normalize_azure_ai_api_base(base: Optional[str] = None) -> str:
    """Return AZURE_AI_API_BASE with trailing slash stripped."""
    raw = base if base is not None else os.environ.get("AZURE_AI_API_BASE", "")
    return (raw or "").strip().rstrip("/")


def is_azure_openai_v1_base(base: Optional[str] = None) -> bool:
    """True when the configured Azure base is the Foundry OpenAI v1 endpoint."""
    normalized = normalize_azure_ai_api_base(base)
    return normalized.endswith("/openai/v1")


def azure_responses_api_url(base: Optional[str] = None) -> str:
    """Build the Azure Foundry Responses API URL for the configured base.

    Accepts either the resource root or a base that already ends in
    ``/openai/v1`` so callers never double the path segment.
    """
    normalized = normalize_azure_ai_api_base(base)
    if not normalized:
        return ""
    if is_azure_openai_v1_base(normalized):
        return f"{normalized}/responses?api-version=v1"
    return f"{normalized}/openai/v1/responses?api-version=v1"


def _model_rejects_nondefault_temperature(model: str) -> bool:
    """GPT-5.6 (and some GPT-5) endpoints only accept the default temperature=1."""
    m = (model or "").lower()
    # Strip provider prefixes for matching.
    if "/" in m:
        m = m.rsplit("/", 1)[-1]
    return m.startswith("gpt-5.6") or m.startswith("gpt-5.6-")


def resolve_litellm_model_and_kwargs(
    model: str, **kwargs: Any
) -> tuple[str, Dict[str, Any]]:
    """Rewrite azure_ai/* for Foundry OpenAI v1 bases; otherwise pass through."""
    out = dict(kwargs)
    base = normalize_azure_ai_api_base()
    key = (os.environ.get("AZURE_AI_API_KEY") or "").strip()
    resolved = model
    if model.startswith("azure_ai/") and is_azure_openai_v1_base(base):
        deployment = model[len("azure_ai/") :]
        out.setdefault("api_base", base)
        if key:
            out.setdefault("api_key", key)
        resolved = f"openai/{deployment}"
    # gpt-5.6-sol rejects temperature != 1 (and temperature=0 is common in GVF).
    # Omit the param so the provider default applies.
    if _model_rejects_nondefault_temperature(resolved) and "temperature" in out:
        temp = out.get("temperature")
        if temp is not None and float(temp) != 1.0:
            logger.debug(
                "Dropping temperature=%r for model %r (only default 1 supported)",
                temp,
                resolved,
            )
            out.pop("temperature", None)
    return resolved, out


def litellm_completion(*, model: str, **kwargs: Any) -> Any:
    """LiteLLM completion with Azure Foundry OpenAI v1 routing when configured."""
    resolved_model, resolved_kwargs = resolve_litellm_model_and_kwargs(model, **kwargs)
    return completion(model=resolved_model, **resolved_kwargs)


def normalize_reasoning_effort(effort: Optional[str]) -> Optional[str]:
    """Normalize effort aliases (e.g. max→xhigh); return None when unset."""
    if effort is None:
        return None
    v = str(effort).strip().lower()
    if not v:
        return None
    return REASONING_EFFORT_ALIASES.get(v, v)


def build_reasoning_effort_kwargs(
    model: Optional[str], effort: Optional[str]
) -> Dict[str, Any]:
    """Return completion() kwargs requesting a reasoning effort, or {} for default.

    Safe no-op when ``effort`` is falsy or the model's provider does not accept an
    OpenAI-style ``reasoning_effort`` parameter. Keeping the allow-list explicit
    (instead of relying on litellm.drop_params) means an effort set on an
    unsupported model is logged as ignored rather than silently dropped.
    """
    effort = normalize_reasoning_effort(effort)
    if not effort:
        return {}
    m = (model or "").lower()
    if any(hint in m for hint in _REASONING_EFFORT_MODEL_HINTS):
        return {"reasoning_effort": effort}
    logger.debug(
        "reasoning_effort=%r ignored for model %r (no OpenAI-style effort knob)",
        effort,
        model,
    )
    return {}


def build_responses_reasoning_param(
    model: Optional[str], effort: Optional[str]
) -> Dict[str, Any]:
    """Return a Responses-API ``reasoning`` body fragment, or {} for default.

    The OpenAI / Azure AI Foundry Responses API expresses reasoning effort as a
    nested ``{"reasoning": {"effort": ...}}`` request-body field rather than the
    chat-completions ``reasoning_effort`` kwarg. Same allow-list and no-op
    semantics as :func:`build_reasoning_effort_kwargs`, so vision/figure code
    that POSTs to ``/responses`` can spread it straight into the JSON body.
    """
    effort = normalize_reasoning_effort(effort)
    if not effort:
        return {}
    m = (model or "").lower()
    if any(hint in m for hint in _REASONING_EFFORT_MODEL_HINTS):
        return {"reasoning": {"effort": effort}}
    logger.debug(
        "reasoning effort=%r ignored for Responses model %r (no effort support)",
        effort,
        model,
    )
    return {}


def get_model_token_limit(model: str, use_safe_limit: bool = True) -> int:
    """
    Get the token limit for a model.

    Args:
        model: Model name/identifier
        use_safe_limit: If True, return the safe limit with margin;
                       if False, return the actual API limit

    Returns:
        Token limit for the model
    """
    m = model.lower() if model else ""

    # Check patterns from most specific to least specific
    for pattern in sorted(MODEL_TOKEN_LIMITS.keys(), key=len, reverse=True):
        if pattern in m:
            actual, safe = MODEL_TOKEN_LIMITS[pattern]
            return safe if use_safe_limit else actual

    # Return default
    actual, safe = DEFAULT_TOKEN_LIMIT
    return safe if use_safe_limit else actual


def clamp_max_tokens(model: str, requested: int, warn: bool = True) -> int:
    """
    Clamp max_tokens to the safe limit for a model.

    Args:
        model: Model name/identifier
        requested: Requested max_tokens value
        warn: If True, log when clamping occurs

    Returns:
        Clamped max_tokens value
    """
    limit = get_model_token_limit(model, use_safe_limit=True)

    if requested > limit:
        if warn:
            logger.info(f"Clamping max_tokens for '{model}': {requested} -> {limit}")
        return limit

    return requested


# Rate limiter: Simple token bucket to avoid hitting OpenAI rate limits
class RateLimiter:
    """Thread-safe rate limiter using token bucket algorithm."""

    def __init__(self, requests_per_minute: int = 50):
        self.requests_per_minute = requests_per_minute
        self.min_interval = 60.0 / requests_per_minute  # seconds between requests
        self.last_request_time = 0.0
        self._lock = threading.Lock()

    def wait_if_needed(self):
        """Sleep if necessary to respect rate limit.

        Claim the next slot under the lock (so concurrent workers each get a
        distinct slot), then release the lock before sleeping. The prior
        implementation slept *inside* the lock, which serialized N workers
        through ~min_interval each and collapsed effective parallelism — a
        10-worker filter run with RPM=200 ran at ~21 RPM instead of ~200.
        """
        with self._lock:
            current_time = time.time()
            next_slot = max(current_time, self.last_request_time + self.min_interval)
            sleep_time = next_slot - current_time
            self.last_request_time = next_slot
        if sleep_time > 0:
            logger.debug(f"Rate limiting: sleeping {sleep_time:.2f}s")
            time.sleep(sleep_time)


# Rate limiter resolution order:
#   1. LLM_REQUESTS_PER_MINUTE env var (explicit override)
#   2. Model-provider default from Settings (200 for Anthropic, 50 for Azure)
#   3. MODEL_PROVIDER-aware default from Settings when no model is available
#   4. Hard fallback of 50 if Settings can't be loaded yet (e.g. import-time
#      circular issues during early bootstrap or tests that stub the env).
# The previous flat 50 RPM was the default for OpenAI tiers and bottlenecked
# the pipeline once Anthropic became the default provider — its tier-4 RPM
# budget is several hundred RPM for Sonnet/Haiku.
def _model_provider_hint(model: Optional[str]) -> Optional[str]:
    m = (model or "").strip().lower()
    if m.startswith("anthropic/") or "claude" in m:
        return "anthropic"
    if "deepseek" in m:
        return "deepseek"
    if m.startswith("azure_ai/"):
        return "azure"
    return None


def _resolve_rpm_limit(model: Optional[str] = None) -> int:
    raw = os.getenv("LLM_REQUESTS_PER_MINUTE", "").strip()
    if raw:
        try:
            value = int(raw)
            if value > 0:
                return value
        except ValueError:
            logger.warning(
                "LLM_REQUESTS_PER_MINUTE=%r is not a positive integer; falling back to provider default",
                raw,
            )

    # Defer import to avoid circular dependency at module load (utils.llm_utils
    # is imported widely; config.settings imports nothing from utils).
    try:
        from config.settings import get_settings

        settings = get_settings()
        provider = _model_provider_hint(model)
        if provider == "anthropic":
            return settings.anthropic_rpm
        if provider == "deepseek":
            return settings.deepseek_rpm
        if provider == "azure":
            return settings.azure_rpm
        return settings.get_requests_per_minute()
    except Exception as exc:  # pragma: no cover — defensive fallback
        logger.debug("Could not load Settings for RPM (%s); using 50", exc)
        return 50


_rate_limiters: Dict[int, RateLimiter] = {}
_rate_limiters_lock = threading.Lock()


def _get_rate_limiter(model: Optional[str] = None) -> RateLimiter:
    rpm = _resolve_rpm_limit(model)
    with _rate_limiters_lock:
        limiter = _rate_limiters.get(rpm)
        if limiter is None:
            limiter = RateLimiter(requests_per_minute=rpm)
            _rate_limiters[rpm] = limiter
        return limiter


def wait_for_llm_rate_limit(model: Optional[str] = None) -> None:
    """Wait according to the throttle appropriate for the specific model."""
    _get_rate_limiter(model).wait_if_needed()


class BaseLLMCaller:
    """
    Base class for components that make LLM API calls.

    This class encapsulates common LLM calling patterns including:
    - Consistent parameter configuration (model, temperature, max_tokens)
    - Retry logic for transient failures
    - JSON response parsing and validation
    - Error handling and logging

    Attributes:
        model: The LLM model identifier (e.g., "gpt-4o", "claude-3-opus")
        temperature: Sampling temperature for response generation (0.0-1.0)
        max_tokens: Maximum number of tokens in the response

    Example:
        >>> class MyFilter(BaseLLMCaller):
        >>>     def __init__(self):
        >>>         super().__init__(model="gpt-4o", temperature=0.3)
        >>>
        >>>     def analyze(self, text: str) -> dict:
        >>>         prompt = f"Analyze this text: {text}"
        >>>         return self.call_llm_json(prompt)
    """

    def __init__(
        self,
        model: str = "gpt-4o",
        temperature: float = 0.3,
        max_tokens: float = 4000,
        reasoning_effort: Optional[str] = None,
    ):
        """
        Initialize the LLM caller with model parameters.

        Args:
            model: The LLM model to use (default: "gpt-4o")
            temperature: Sampling temperature (default: 0.3)
            max_tokens: Maximum response length (default: 4000)
            reasoning_effort: Optional OpenAI-style reasoning effort
                ("minimal"|"low"|"medium"|"high"). None leaves the provider
                default. Silently ignored for models without an effort knob
                (see build_reasoning_effort_kwargs).
        """
        self.model = model
        self.temperature = temperature
        self.max_tokens = max_tokens
        self.reasoning_effort = reasoning_effort
        logger.info(
            f"Initialized {self.__class__.__name__} with model={model}, "
            f"temperature={temperature}, max_tokens={max_tokens}, "
            f"reasoning_effort={reasoning_effort}"
        )

    def _attempt_json_repair(self, raw_text: str) -> Optional[Dict[str, Any]]:
        """
        Best-effort repair for malformed/truncated JSON emitted by the LLM.

        Uses the same model with a constrained prompt to fix the structure while
        trimming incomplete trailing items. Returns None on failure.

        For large variant extractions, uses higher token limits to preserve data.
        """
        # Use a larger char limit for repair - variant tables can be very large
        # Each variant is ~300-500 chars of JSON, so 150 variants needs ~75K chars
        max_repair_chars = 80000
        repair_chunk = raw_text[:max_repair_chars]

        repair_prompt = (
            "The following text is intended to be a JSON object but is malformed or truncated. "
            "Return a valid JSON object that keeps ALL intact content and discards only incomplete tail fragments. "
            "CRITICAL: Preserve all special characters in protein notation, especially asterisks (*) for stop codons "
            "(e.g., 'p.Arg412*' or 'p.Gly24fs*58'). "
            "Do not add new information. Respond with JSON only.\n\n"
            f"{repair_chunk}"
        )

        # Use higher token limit for repair to handle large variant tables
        # Default repair needs more tokens for 100+ variant extractions
        repair_max_tokens = min(self.max_tokens, 16000)

        try:
            wait_for_llm_rate_limit(self.model)
            response = litellm_completion(
                model=self.model,
                messages=[
                    {
                        "role": "system",
                        "content": "You fix malformed JSON. Preserve all data including special characters like asterisks (*) in protein notation. Respond with JSON only.",
                    },
                    {"role": "user", "content": repair_prompt},
                ],
                temperature=0,
                max_tokens=repair_max_tokens,
                response_format={"type": "json_object"},
            )
            repaired_text = response.choices[0].message.content
            repaired = parse_llm_json_response(repaired_text)
            # Reasoning models (e.g. Kimi-K2.6) sometimes truncate the original
            # output mid-prefix; repair then "succeeds" by emitting an empty
            # `{}`. That looks like a successful parse to callers, but it's
            # missing the entire schema. Treat empty objects (or empty
            # lists/strings) as a failed repair so callers fall through to
            # their JSONDecodeError handler / fail-open path instead of
            # silently using default values for every key.
            if isinstance(repaired, dict) and not repaired:
                logger.warning(
                    "JSON repair returned empty object — treating as failure."
                )
                return None
            if isinstance(repaired, list) and not repaired:
                logger.warning("JSON repair returned empty list — treating as failure.")
                return None
            return repaired
        except Exception as repair_exc:
            logger.error(f"JSON repair attempt failed: {repair_exc}")
            return None

    @llm_retry
    def call_llm_json_with_status(
        self,
        prompt: str,
        system_message: Optional[str] = None,
        response_format: Optional[Dict[str, Any]] = None,
    ) -> tuple:
        """
        Call the LLM and parse JSON, returning (data, was_truncated, raw_text).

        This variant returns additional metadata needed for continuation extraction.
        """
        messages = []
        if system_message:
            messages.append({"role": "system", "content": system_message})
        messages.append({"role": "user", "content": prompt})

        if response_format is None:
            response_format = {"type": "json_object"}

        def _make_call():
            wait_for_llm_rate_limit(self.model)
            return litellm_completion(
                model=self.model,
                messages=messages,
                temperature=self.temperature,
                max_tokens=self.max_tokens,
                response_format=response_format,
                **build_reasoning_effort_kwargs(self.model, self.reasoning_effort),
            )

        response = _make_call()

        result_text = response.choices[0].message.content
        finish_reason = getattr(response.choices[0], "finish_reason", None)
        if not (result_text or "").strip():
            logger.warning(
                "LLM returned empty JSON content for model=%s finish_reason=%s; retrying once",
                self.model,
                finish_reason,
            )
            response = _make_call()
            result_text = response.choices[0].message.content
            finish_reason = getattr(response.choices[0], "finish_reason", None)
        was_truncated = finish_reason == "length"

        try:
            result_data = parse_llm_json_response(result_text)
            return result_data, was_truncated, result_text
        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse LLM JSON response: {e}")
            if was_truncated:
                logger.warning(
                    "LLM response was cut off due to max_tokens; attempting repair."
                )
            repaired = self._attempt_json_repair(result_text)
            if repaired is not None:
                logger.info("JSON repair succeeded after initial parse failure.")
                return repaired, was_truncated, result_text
            raise

    @llm_retry
    def call_llm_json(
        self,
        prompt: str,
        system_message: Optional[str] = None,
        response_format: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """
        Call the LLM with a prompt and parse the JSON response.

        This method handles the full lifecycle of an LLM call:
        1. Constructs the message array with optional system message
        2. Calls the LLM API with retry logic
        3. Extracts and parses the JSON response
        4. Validates the response and handles errors

        Args:
            prompt: The user prompt to send to the LLM
            system_message: Optional system message to set context
            response_format: Optional response format specification
                           (default: {"type": "json_object"})

        Returns:
            Parsed JSON response as a dictionary

        Raises:
            json.JSONDecodeError: If the response is not valid JSON
            Exception: For other LLM API errors

        Example:
            >>> caller = BaseLLMCaller()
            >>> result = caller.call_llm_json("What is 2+2?")
            >>> print(result)
        """
        # Build messages array
        messages = []
        if system_message:
            messages.append({"role": "system", "content": system_message})
        messages.append({"role": "user", "content": prompt})

        # Use default JSON response format if not specified
        if response_format is None:
            response_format = {"type": "json_object"}

        logger.debug(
            f"Calling LLM with model={self.model}, prompt length={len(prompt)}"
        )

        try:

            def _make_call():
                wait_for_llm_rate_limit(self.model)
                return litellm_completion(
                    model=self.model,
                    messages=messages,
                    temperature=self.temperature,
                    max_tokens=self.max_tokens,
                    response_format=response_format,
                    **build_reasoning_effort_kwargs(self.model, self.reasoning_effort),
                )

            # Make the LLM API call
            response = _make_call()

            # Extract response text
            result_text = response.choices[0].message.content
            finish_reason = getattr(response.choices[0], "finish_reason", None)
            if not (result_text or "").strip():
                logger.warning(
                    "LLM returned empty JSON content for model=%s finish_reason=%s; retrying once",
                    self.model,
                    finish_reason,
                )
                response = _make_call()
                result_text = response.choices[0].message.content
                finish_reason = getattr(response.choices[0], "finish_reason", None)

            # Parse JSON response
            result_data = parse_llm_json_response(result_text)

            logger.debug(
                f"LLM call successful, response keys: {list(result_data.keys())}"
            )
            return result_data

        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse LLM JSON response: {e}")
            logger.error(f"Response text: {result_text[:500]}")
            if finish_reason == "length":
                logger.warning(
                    "LLM response was cut off due to max_tokens; attempting repair."
                )
            repaired = self._attempt_json_repair(result_text)
            if repaired is not None:
                logger.info("JSON repair succeeded after initial parse failure.")
                return repaired
            raise
        except Exception as e:
            logger.error(f"LLM API call failed: {e}")
            raise

    def call_llm_text(self, prompt: str, system_message: Optional[str] = None) -> str:
        """
        Call the LLM and return the raw text response.

        Use this method when you need free-form text rather than structured JSON.

        Args:
            prompt: The user prompt to send to the LLM
            system_message: Optional system message to set context

        Returns:
            Raw text response from the LLM

        Example:
            >>> caller = BaseLLMCaller()
            >>> text = caller.call_llm_text("Write a haiku about code")
            >>> print(text)
        """
        # Build messages array
        messages = []
        if system_message:
            messages.append({"role": "system", "content": system_message})
        messages.append({"role": "user", "content": prompt})

        logger.debug(f"Calling LLM for text response, prompt length={len(prompt)}")

        try:
            wait_for_llm_rate_limit(self.model)
            response = litellm_completion(
                model=self.model,
                messages=messages,
                temperature=self.temperature,
                max_tokens=self.max_tokens,
                **build_reasoning_effort_kwargs(self.model, self.reasoning_effort),
            )

            result_text = response.choices[0].message.content
            logger.debug(f"LLM text response received, length={len(result_text)}")
            return result_text

        except Exception as e:
            logger.error(f"LLM API call failed: {e}")
            raise


def parse_llm_json_response(response_text: str) -> Dict[str, Any]:
    """
    Parse a JSON response from an LLM.

    This function handles common edge cases in LLM JSON responses:
    - Strips markdown code blocks (```json ... ```)
    - Handles extra whitespace
    - Provides clear error messages on parse failures

    Args:
        response_text: Raw text response from the LLM

    Returns:
        Parsed JSON as a dictionary

    Raises:
        json.JSONDecodeError: If the response cannot be parsed as JSON

    Example:
        >>> response = '```json\\n{"result": "success"}\\n```'
        >>> data = parse_llm_json_response(response)
        >>> print(data)
        {'result': 'success'}
    """
    # Strip markdown code blocks if present
    cleaned_text = response_text.strip()
    if cleaned_text.startswith("```json"):
        cleaned_text = cleaned_text[7:]  # Remove ```json
    if cleaned_text.startswith("```"):
        cleaned_text = cleaned_text[3:]  # Remove ```
    if cleaned_text.endswith("```"):
        cleaned_text = cleaned_text[:-3]  # Remove trailing ```

    cleaned_text = cleaned_text.strip()

    # Try strict parse first — fast path when the model returns clean JSON.
    try:
        return json.loads(cleaned_text)
    except json.JSONDecodeError as strict_err:
        # Common failure mode: model emits valid JSON followed by reasoning
        # prose ("Extra data: line 6 column 1"). Claude Haiku 4.5 does this
        # frequently. raw_decode parses the first valid JSON value and
        # returns the index where it stopped, so we can ignore trailing prose.
        decoder = json.JSONDecoder()
        candidate_starts = [
            idx for idx, char in enumerate(cleaned_text) if char in "{["
        ]
        for start in candidate_starts:
            try:
                obj, _end = decoder.raw_decode(cleaned_text[start:])
            except json.JSONDecodeError:
                continue
            if isinstance(obj, dict):
                return obj
            # Tolerate top-level lists by wrapping; callers expect dicts.
            if isinstance(obj, list):
                return {"items": obj}

        logger.error(f"Failed to parse JSON: {strict_err}")
        logger.error(f"Response text (first 500 chars): {cleaned_text[:500]}")
        raise strict_err
