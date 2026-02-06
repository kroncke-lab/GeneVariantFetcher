"""
Resilience utilities for external API calls with circuit breaker pattern.
"""

import time
import logging
from typing import Any, Callable, Optional
from datetime import datetime, timedelta

logger = logging.getLogger(__name__)


class CircuitBreaker:
    """Circuit breaker for external API calls to prevent cascading failures."""
    
    def __init__(self, name: str, max_failures: int = 5, reset_timeout: int = 60):
        """
        Initialize circuit breaker.
        
        Args:
            name: Identifier for the circuit breaker (used in logging)
            max_failures: Number of failures before circuit opens
            reset_timeout: Seconds to wait before attempting half-open state
        """
        self.name = name
        self.failure_count = 0
        self.max_failures = max_failures
        self.reset_timeout = reset_timeout
        self.state = "closed"  # closed, open, half-open
        self.last_failure_time = None
        self.consecutive_successes = 0
        
    def is_open(self) -> bool:
        """Check if circuit is currently open."""
        if self.state == "closed":
            return False
        elif self.state == "open":
            if self.last_failure_time and \
               datetime.now() - self.last_failure_time >= timedelta(seconds=self.reset_timeout):
                # Time to try half-open
                self.state = "half-open"
                logger.info(f"Circuit breaker '{self.name}' transitioning to HALF-OPEN")
                return False
            return True
        else:  # half-open
            return False
    
    def call(self, func: Callable, *args, **kwargs) -> Any:
        """
        Execute a function with circuit breaker protection.
        
        Args:
            func: Function to execute
            *args: Positional arguments for the function
            **kwargs: Keyword arguments for the function
            
        Returns:
            Function return value
            
        Raises:
            Exception: Any exception from called function or circuit breaker
        """
        if self.is_open():
            raise CircuitBreakerOpenError(
                f"Circuit breaker '{self.name}' is OPEN ({self.state})"
            )
        
        try:
            result = func(*args, **kwargs)
            self.record_success()
            return result
        except Exception as e:
            self.record_failure()
            logger.warning(f"Circuit breaker '{self.name}' recorded failure: {e}")
            raise
    
    def record_success(self) -> None:
        """Record a successful call."""
        if self.state == "half-open":
            self.consecutive_successes += 1
            if self.consecutive_successes >= 3:  # Need 3 successes in half-open to close
                self.state = "closed"
                self.failure_count = 0
                self.consecutive_successes = 0
                logger.info(f"Circuit breaker '{self.name}' transitioning to CLOSED")
        elif self.state == "closed":
            # Already closed, just reset failure count
            self.failure_count = 0
            self.consecutive_successes = 0
    
    def record_failure(self) -> None:
        """Record a failed call."""
        self.failure_count += 1
        self.last_failure_time = datetime.now()
        self.consecutive_successes = 0
        
        if self.failure_count >= self.max_failures and self.state == "closed":
            self.state = "open"
            logger.error(
                f"Circuit breaker '{self.name}' OPENED after {self.failure_count} failures"
            )
    
    def get_state(self) -> dict:
        """Get current state for monitoring."""
        return {
            "name": self.name,
            "state": self.state,
            "failure_count": self.failure_count,
            "max_failures": self.max_failures,
            "last_failure_time": self.last_failure_time.isoformat() if self.last_failure_time else None,
            "reset_timeout": self.reset_timeout,
            "is_open": self.is_open()
        }


class CircuitBreakerOpenError(Exception):
    """Raised when attempting to call a function with an open circuit breaker."""
    pass


class ResilientAPIClient:
    """Wrapper class that adds circuit breaker protection to API clients."""
    
    def __init__(self, base_client: Any, circuit_breaker: CircuitBreaker):
        """
        Initialize resilient API client.
        
        Args:
            base_client: The underlying API client
            circuit_breaker: Circuit breaker for protection
        """
        self.base_client = base_client
        self.circuit_breaker = circuit_breaker
    
    def __getattr__(self, name: str) -> Any:
        """Proxy all method calls through the circuit breaker."""
        original_method = getattr(self.base_client, name)
        
        if callable(original_method):
            def wrapped_method(*args, **kwargs):
                return self.circuit_breaker.call(original_method, *args, **kwargs)
            return wrapped_method
        return original_method