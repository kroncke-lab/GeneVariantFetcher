"""Offline integration checks across resilience, manifests, and harvesting."""

import json


def test_imports():
    """Test that all new modules can be imported."""
    from utils.resilience import (
        CircuitBreaker,
        CircuitBreakerOpenError,
        ResilientAPIClient,
    )
    from utils.run_manifest import RunManifest

    assert all(
        symbol is not None
        for symbol in (
            CircuitBreaker,
            CircuitBreakerOpenError,
            ResilientAPIClient,
            RunManifest,
        )
    )


def test_circuit_breaker():
    """Test circuit breaker functionality."""
    from utils.resilience import CircuitBreaker

    cb = CircuitBreaker("test_api", max_failures=3, reset_timeout=1)

    # Should start closed
    assert cb.state == "closed", f"Expected closed, got {cb.state}"
    # Record failures to trip the breaker
    for _ in range(3):
        cb.record_failure()

    assert cb.state == "open", f"Expected open after 3 failures, got {cb.state}"

    # Should block calls when open
    assert cb.is_open() is True


def test_run_manifest(tmp_path):
    """Test run manifest generation."""
    from utils.run_manifest import RunManifest

    manifest = RunManifest(output_dir=tmp_path, gene_symbol="KCNH2")
    manifest.set_config(test=True)
    assert manifest.run_id is not None

    manifest.save()
    manifest_path = tmp_path / "run_manifest.json"
    assert manifest_path.exists()

    manifest.finalize()
    data = json.loads(manifest_path.read_text())
    assert data["status"] == "completed"


def test_orchestrator_imports():
    """Test that orchestrator can import with new circuit breakers."""
    from harvesting.orchestrator import PMCHarvester

    assert PMCHarvester is not None
