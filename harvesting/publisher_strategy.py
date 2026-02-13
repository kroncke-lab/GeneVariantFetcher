"""Publisher fallback strategy helpers for DOI-based full-text retrieval."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List


@dataclass(frozen=True)
class PublisherAttempt:
    """One provider attempt in fallback order."""

    name: str
    should_try: bool
    try_fetch: Callable[[str, str], tuple]


def build_publisher_attempt_plan(
    doi: str,
    attempts: List[PublisherAttempt],
) -> List[PublisherAttempt]:
    """Return ordered provider attempts: matched DOI providers first, then others."""
    matched = [attempt for attempt in attempts if attempt.should_try]
    unmatched = [attempt for attempt in attempts if not attempt.should_try]
    return matched + unmatched
