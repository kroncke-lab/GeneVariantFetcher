"""
Retry Manager - Exponential backoff retry logic for GVF downloads.

Handles transient failures with intelligent retry strategies.
"""

import logging
import time
import random
from typing import Callable, TypeVar, Optional, Any
from functools import wraps
from dataclasses import dataclass, field
from datetime import datetime
import json
from pathlib import Path

logger = logging.getLogger(__name__)

T = TypeVar("T")


@dataclass
class RetryConfig:
    """Configuration for retry behavior."""

    max_retries: int = 3
    base_delay: float = 1.0  # seconds
    max_delay: float = 60.0  # seconds
    exponential_base: float = 2.0
    jitter: bool = True  # Add randomness to prevent thundering herd
    retry_on_exceptions: tuple = (Exception,)


@dataclass
class RetryAttempt:
    """Record of a single retry attempt."""

    attempt_number: int
    timestamp: str
    success: bool
    error: Optional[str] = None
    delay_used: float = 0.0


@dataclass
class RetryResult:
    """Result of a retry operation."""

    success: bool
    result: Any = None
    total_attempts: int = 0
    total_time: float = 0.0
    attempts: list = field(default_factory=list)
    final_error: Optional[str] = None


class RetryManager:
    """Manages retry logic with exponential backoff."""

    def __init__(
        self, config: Optional[RetryConfig] = None, log_file: Optional[Path] = None
    ):
        """
        Initialize retry manager.

        Args:
            config: Retry configuration
            log_file: Optional path to log retry attempts
        """
        self.config = config or RetryConfig()
        self.log_file = log_file
        self._retry_history = []

    def calculate_delay(self, attempt: int) -> float:
        """Calculate delay for given attempt number using exponential backoff."""
        delay = self.config.base_delay * (self.config.exponential_base**attempt)
        delay = min(delay, self.config.max_delay)

        if self.config.jitter:
            # Add up to 25% jitter
            jitter = delay * random.uniform(0, 0.25)
            delay += jitter

        return delay

    def retry(self, func: Callable[[], T], context: str = "") -> RetryResult:
        """
        Execute function with retry logic.

        Args:
            func: Function to execute (should take no arguments)
            context: Description for logging

        Returns:
            RetryResult with success status and result/error
        """
        attempts = []
        start_time = time.time()
        last_error = None

        for attempt in range(self.config.max_retries + 1):
            attempt_start = datetime.now().isoformat()
            delay_used = 0.0

            if attempt > 0:
                delay_used = self.calculate_delay(attempt - 1)
                logger.info(
                    f"Retry {attempt}/{self.config.max_retries} for {context}, waiting {delay_used:.1f}s"
                )
                time.sleep(delay_used)

            try:
                result = func()
                attempts.append(
                    RetryAttempt(
                        attempt_number=attempt + 1,
                        timestamp=attempt_start,
                        success=True,
                        delay_used=delay_used,
                    )
                )

                total_time = time.time() - start_time
                retry_result = RetryResult(
                    success=True,
                    result=result,
                    total_attempts=attempt + 1,
                    total_time=total_time,
                    attempts=attempts,
                )

                self._log_result(context, retry_result)
                return retry_result

            except self.config.retry_on_exceptions as e:
                last_error = str(e)
                logger.warning(
                    f"Attempt {attempt + 1} failed for {context}: {last_error}"
                )

                attempts.append(
                    RetryAttempt(
                        attempt_number=attempt + 1,
                        timestamp=attempt_start,
                        success=False,
                        error=last_error,
                        delay_used=delay_used,
                    )
                )

        # All retries exhausted
        total_time = time.time() - start_time
        retry_result = RetryResult(
            success=False,
            total_attempts=self.config.max_retries + 1,
            total_time=total_time,
            attempts=attempts,
            final_error=last_error,
        )

        self._log_result(context, retry_result)
        return retry_result

    def _log_result(self, context: str, result: RetryResult):
        """Log retry result to file if configured."""
        if not self.log_file:
            return

        log_entry = {
            "context": context,
            "timestamp": datetime.now().isoformat(),
            "success": result.success,
            "total_attempts": result.total_attempts,
            "total_time": result.total_time,
            "final_error": result.final_error,
        }

        self._retry_history.append(log_entry)

        try:
            # Append to log file
            with open(self.log_file, "a") as f:
                f.write(json.dumps(log_entry) + "\n")
        except Exception as e:
            logger.warning(f"Failed to write retry log: {e}")

    def get_stats(self) -> dict:
        """Get retry statistics."""
        if not self._retry_history:
            return {"total": 0, "success_rate": 0.0}

        total = len(self._retry_history)
        successes = sum(1 for r in self._retry_history if r["success"])

        return {
            "total": total,
            "successes": successes,
            "failures": total - successes,
            "success_rate": successes / total if total > 0 else 0.0,
            "avg_attempts": sum(r["total_attempts"] for r in self._retry_history)
            / total,
            "avg_time": sum(r["total_time"] for r in self._retry_history) / total,
        }


def with_retry(config: Optional[RetryConfig] = None, context: str = ""):
    """
    Decorator to add retry logic to a function.

    Usage:
        @with_retry(RetryConfig(max_retries=3), context="fetch paper")
        def fetch_paper(doi):
            ...
    """

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            manager = RetryManager(config)
            result = manager.retry(
                lambda: func(*args, **kwargs), context or func.__name__
            )
            if result.success:
                return result.result
            raise Exception(f"All retries failed: {result.final_error}")

        return wrapper

    return decorator


# Convenience function for quick retries
def retry_call(
    func: Callable[[], T], max_retries: int = 3, context: str = ""
) -> RetryResult:
    """Quick retry helper with default settings."""
    config = RetryConfig(max_retries=max_retries)
    manager = RetryManager(config)
    return manager.retry(func, context)
