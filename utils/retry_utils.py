"""
Retry utilities for handling transient failures in network operations.

This module provides standardized retry configurations used across the project
for API calls, LLM requests, and web scraping operations.
"""

import os
from typing import Tuple, Type

import requests
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
    wait_exponential_jitter,
)

# Pull LiteLLM's transient error types so llm_retry can include them. These
# are the standard exceptions LiteLLM raises across providers; we want to
# retry rate-limit, timeout, and service-unavailable errors but NOT auth or
# bad-request errors (those are programmer-side bugs that won't fix on retry).
try:
    from litellm import RateLimitError as _LiteLLMRateLimitError
    from litellm import Timeout as _LiteLLMTimeout
    from litellm import APIConnectionError as _LiteLLMAPIConnectionError
    from litellm import ServiceUnavailableError as _LiteLLMServiceUnavailable
    from litellm import InternalServerError as _LiteLLMInternalServerError

    LITELLM_TRANSIENT_ERRORS: Tuple[Type[Exception], ...] = (
        _LiteLLMRateLimitError,
        _LiteLLMTimeout,
        _LiteLLMAPIConnectionError,
        _LiteLLMServiceUnavailable,
        _LiteLLMInternalServerError,
    )
except ImportError:
    LITELLM_TRANSIENT_ERRORS = ()

# Standard retry configuration used throughout the project
# Retries up to 3 times with exponential backoff (1s, 2s, 4s)
standard_retry = retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=1, max=10),
    reraise=True,
)


def get_standard_retry_decorator(
    max_attempts: int = 3,
    multiplier: float = 1.0,
    min_wait: float = 1.0,
    max_wait: float = 10.0,
    retry_exceptions: Tuple[Type[Exception], ...] = (
        requests.exceptions.RequestException,
        ConnectionError,
        TimeoutError,
    ),
):
    """
    Create a retry decorator with custom configuration.

    This factory function allows creating retry decorators with custom parameters
    while maintaining consistency across the codebase.

    Args:
        max_attempts: Maximum number of retry attempts (default: 3)
        multiplier: Exponential backoff multiplier (default: 1.0)
        min_wait: Minimum wait time in seconds (default: 1.0)
        max_wait: Maximum wait time in seconds (default: 10.0)
        retry_exceptions: Tuple of exception types to retry on

    Returns:
        A configured retry decorator

    Example:
        >>> @get_standard_retry_decorator(max_attempts=5)
        >>> def fetch_data():
        >>>     return requests.get("https://api.example.com")
    """
    return retry(
        stop=stop_after_attempt(max_attempts),
        wait=wait_exponential(multiplier=multiplier, min=min_wait, max=max_wait),
        retry=retry_if_exception_type(retry_exceptions),
        reraise=True,
    )


# Pre-configured retry decorators for common use cases

# For API calls that may have transient network issues
api_retry = get_standard_retry_decorator(
    max_attempts=3,
    retry_exceptions=(
        requests.exceptions.RequestException,
        ConnectionError,
        TimeoutError,
    ),
)

# For LLM calls that may hit rate limits or transient errors. Includes
# LiteLLM's rate-limit/timeout/connection/server-error types so a Tier 2/3
# call doesn't immediately fail when an Azure deployment quota briefly
# saturates from parallel workers.
#
# Bumped to 8 attempts with jittered exponential backoff because:
#   - Anthropic returns HTTP 429 (rate limit) and 529 (overloaded) under
#     bursty load; both are LiteLLM RateLimitError / InternalServerError.
#   - 529 (overloaded) means *all* Anthropic users are saturated, so it can
#     persist for tens of seconds — needs longer total budget than 4-60s.
#   - With 10 parallel workers, all retries fire on the same wall-clock
#     boundary without jitter, producing thundering-herd retries that hit
#     the rate limit again. wait_exponential_jitter spreads them out.
_LLM_RETRY_EXCEPTIONS: Tuple[Type[Exception], ...] = (
    requests.exceptions.RequestException,
    ConnectionError,
    TimeoutError,
    *LITELLM_TRANSIENT_ERRORS,
)


def _positive_int_env(name: str, default: int) -> int:
    raw = os.getenv(name, "").strip()
    if raw:
        try:
            value = int(raw)
            if value > 0:
                return value
        except ValueError:
            pass
    return default


llm_retry = retry(
    stop=stop_after_attempt(_positive_int_env("LLM_RETRY_ATTEMPTS", 8)),
    # Jittered exponential: nominal 4s, 8s, 16s, 32s, 60s, 60s, 60s with up to
    # 10s additive jitter per attempt. Total budget ~5 minutes, enough to
    # outlast a typical Anthropic 529 spike without burning all attempts in
    # the first 30 seconds.
    wait=wait_exponential_jitter(initial=4.0, max=60.0, exp_base=2.0, jitter=10.0),
    retry=retry_if_exception_type(_LLM_RETRY_EXCEPTIONS),
    reraise=True,
)

# For web scraping operations
scraping_retry = get_standard_retry_decorator(
    max_attempts=2,  # Fewer retries for scraping
    multiplier=1.0,
    min_wait=1.0,
    max_wait=5.0,
    retry_exceptions=(
        requests.exceptions.RequestException,
        ConnectionError,
        TimeoutError,
    ),
)
