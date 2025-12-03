"""
Centralized logging configuration for the Gene Variant Fetcher pipeline.

This module provides standardized logging setup to ensure consistent
log formatting and behavior across all pipeline components.
"""

import logging
import sys
from pathlib import Path
from typing import Optional


# Standard log format used across the project
DEFAULT_LOG_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

# Compact format for CLI tools
CLI_LOG_FORMAT = '%(levelname)s: %(message)s'


def setup_logging(
    level: str | int = logging.INFO,
    log_file: Optional[str | Path] = None,
    format_string: Optional[str] = None,
    use_cli_format: bool = False
) -> None:
    """
    Configure logging for the application.

    This function should be called once at application startup to configure
    the root logger. All modules using `logging.getLogger(__name__)` will
    inherit this configuration.

    Args:
        level: Logging level (e.g., 'INFO', 'DEBUG', or logging.INFO)
        log_file: Optional path to a log file. If provided, logs will be
                  written to both console and file.
        format_string: Custom format string. If None, uses DEFAULT_LOG_FORMAT.
        use_cli_format: If True, uses compact CLI_LOG_FORMAT instead of default.

    Example:
        >>> from utils.logging_utils import setup_logging
        >>> setup_logging(level='DEBUG', log_file='pipeline.log')
    """
    # Determine format
    if format_string:
        log_format = format_string
    elif use_cli_format:
        log_format = CLI_LOG_FORMAT
    else:
        log_format = DEFAULT_LOG_FORMAT

    # Convert string level to int if needed
    if isinstance(level, str):
        level = getattr(logging, level.upper(), logging.INFO)

    # Configure handlers
    handlers = []

    # Console handler (always present)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(logging.Formatter(log_format))
    handlers.append(console_handler)

    # File handler (if log_file specified)
    if log_file:
        file_path = Path(log_file)
        file_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(file_path, encoding='utf-8')
        file_handler.setFormatter(logging.Formatter(DEFAULT_LOG_FORMAT))
        handlers.append(file_handler)

    # Configure root logger
    logging.basicConfig(
        level=level,
        format=log_format,
        handlers=handlers,
        force=True  # Override any existing configuration
    )


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger with the specified name.

    This is a convenience wrapper around logging.getLogger() that ensures
    the logger follows project conventions.

    Args:
        name: Logger name, typically __name__ from the calling module.

    Returns:
        Configured Logger instance.
    """
    return logging.getLogger(name)
