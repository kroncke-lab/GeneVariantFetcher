"""
Utilities for retry logic with exponential backoff.
"""

import time
import logging
from functools import wraps
from typing import Callable, Any

logger = logging.getLogger(__name__)

def llm_retry(func: Callable) -> Callable:
    """
    A decorator to retry a function call with exponential backoff.

    This decorator will retry a function up to 5 times if it raises an
    exception. The delay between retries increases exponentially.

    Args:
        func: The function to decorate.

    Returns:
        The wrapped function with retry logic.
    """
    @wraps(func)
    def wrapper(*args, **kwargs) -> Any:
        max_retries = 5
        base_delay = 1  # in seconds

        for attempt in range(max_retries):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                if attempt == max_retries - 1:
                    logger.error(
                        f"Function {func.__name__} failed after {max_retries} attempts."
                    )
                    raise

                delay = base_delay * (2 ** attempt)
                logger.warning(
                    f"Attempt {attempt + 1}/{max_retries} failed with error: {e}. "
                    f"Retrying in {delay} seconds..."
                )
                time.sleep(delay)
    return wrapper