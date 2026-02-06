#!/usr/bin/env python3
"""Quick integration test for new GVF modules."""

import sys
import tempfile
import json
from pathlib import Path

def test_imports():
    """Test that all new modules can be imported."""
    print("Testing imports...")
    
    try:
        from utils.resilience import CircuitBreaker, ResilientAPIClient, CircuitBreakerOpenError
        print("  ✓ CircuitBreaker, ResilientAPIClient")
    except ImportError as e:
        print(f"  ✗ resilience: {e}")
        return False
    
    try:
        from utils.pmid_status import get_pmid_status, get_failed_pmids, get_stats_summary
        print("  ✓ pmid_status utilities")
    except ImportError as e:
        print(f"  ✗ pmid_status: {e}")
        return False
    
    try:
        from utils.run_manifest import RunManifest
        print("  ✓ RunManifest")
    except ImportError as e:
        print(f"  ✗ run_manifest: {e}")
        return False
    
    return True

def test_circuit_breaker():
    """Test circuit breaker functionality."""
    print("\nTesting CircuitBreaker...")
    from utils.resilience import CircuitBreaker, CircuitBreakerOpenError
    
    cb = CircuitBreaker("test_api", max_failures=3, reset_timeout=1)
    
    # Should start closed
    assert cb.state == "closed", f"Expected closed, got {cb.state}"
    print("  ✓ Starts in closed state")
    
    # Record failures to trip the breaker
    for i in range(3):
        cb.record_failure()
    
    assert cb.state == "open", f"Expected open after 3 failures, got {cb.state}"
    print("  ✓ Opens after max_failures")
    
    # Should block calls when open
    assert cb.is_open() == True
    print("  ✓ is_open() returns True when open")
    
    return True

def test_pmid_status():
    """Test PMID status tracking."""
    print("\nTesting PMID status...")
    from utils.pmid_status import get_pmid_status, get_failed_pmids, get_stats_summary
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir)
        status_dir = output_dir / "pmid_status"
        status_dir.mkdir()
        
        # Write a test status file
        test_status = {
            "pmid": "12345678",
            "status": "extracted",
            "download_timestamp": "2026-02-06T14:00:00",
            "variant_count": 5,
            "failure_reason": None
        }
        with open(status_dir / "12345678.json", "w") as f:
            json.dump(test_status, f)
        
        # Write a failed status
        failed_status = {
            "pmid": "87654321",
            "status": "failed",
            "failure_reason": "Paywall"
        }
        with open(status_dir / "87654321.json", "w") as f:
            json.dump(failed_status, f)
        
        # Test get_pmid_status
        status = get_pmid_status(output_dir, "12345678")
        assert status is not None, "Should find status for 12345678"
        assert status["status"] == "extracted"
        print("  ✓ get_pmid_status works")
        
        # Test get_failed_pmids
        failed = get_failed_pmids(output_dir)
        assert "87654321" in failed, "Should find failed PMID"
        print("  ✓ get_failed_pmids works")
        
        # Test get_stats_summary
        stats = get_stats_summary(output_dir)
        assert stats["extracted"] == 1
        assert stats["failed"] == 1
        print("  ✓ get_stats_summary works")
    
    return True

def test_run_manifest():
    """Test run manifest generation."""
    print("\nTesting RunManifest...")
    from utils.run_manifest import RunManifest
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir)
        
        # Create manifest
        manifest = RunManifest(
            output_dir=output_dir,
            gene_symbol="KCNH2"
        )
        manifest.set_config(test=True)
        
        assert manifest.run_id is not None
        print(f"  ✓ Created manifest with run_id: {manifest.run_id[:8]}...")
        
        # Save it
        manifest.save()
        manifest_path = output_dir / "run_manifest.json"
        assert manifest_path.exists(), "Manifest file should exist"
        print("  ✓ Manifest saved to file")
        
        # Finalize it
        manifest.finalize()
        
        # Reload and check
        with open(manifest_path) as f:
            data = json.load(f)
        assert data["status"] == "completed"
        print("  ✓ Manifest finalized")
    
    return True

def test_orchestrator_imports():
    """Test that orchestrator can import with new circuit breakers."""
    print("\nTesting orchestrator integration...")
    try:
        from harvesting.orchestrator import PMCHarvester
        print("  ✓ PMCHarvester imports successfully")
        return True
    except ImportError as e:
        print(f"  ✗ Failed to import PMCHarvester: {e}")
        return False

def main():
    print("=" * 50)
    print("GVF Integration Tests")
    print("=" * 50)
    
    all_passed = True
    
    if not test_imports():
        all_passed = False
    
    if not test_circuit_breaker():
        all_passed = False
    
    if not test_pmid_status():
        all_passed = False
    
    if not test_run_manifest():
        all_passed = False
    
    if not test_orchestrator_imports():
        all_passed = False
    
    print("\n" + "=" * 50)
    if all_passed:
        print("✅ ALL TESTS PASSED")
        return 0
    else:
        print("❌ SOME TESTS FAILED")
        return 1

if __name__ == "__main__":
    sys.exit(main())
