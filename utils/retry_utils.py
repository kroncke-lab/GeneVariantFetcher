"""
Retry utilities for handling transient failures in network operations.

This module provides standardized retry configurations used across the project
for API calls, LLM requests, and web scraping operations.
"""

from typing import Tuple, Type

import requests
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

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

# For LLM calls that may hit rate limits or transient errors
llm_retry = get_standard_retry_decorator(
    max_attempts=3,
    multiplier=2.0,  # Longer backoff for rate limits
    min_wait=2.0,
    max_wait=30.0,
    retry_exceptions=(
        requests.exceptions.RequestException,
        ConnectionError,
        TimeoutError,
    ),
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
