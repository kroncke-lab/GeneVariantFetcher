#!/usr/bin/env python3
"""
Vanderbilt Institutional Access Validator

Validates all institutional API connections and generates access report.
"""

import os
import sys
import json
from pathlib import Path
from datetime import datetime

def add_project_root():
    """Add project root to Python path."""
    project_root = Path(__file__).parent.parent
    sys.path.append(str(project_root))
    return project_root

project_root = add_project_root()

from harvesting.vandy_access import InstitutionalAccessManager

def main():
    """Main validation script."""
    
    print("=" * 60)
    print("VANDERBILT INSTITUTIONAL ACCESS VALIDATION")
    print("=" * 60)
    print("Testing at: %s" % datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    print("Project root: %s" % str(project_root))
    print()
    
    # Initialize manager
    manager = InstitutionalAccessManager()
    
    # Display current configuration
    print("🔍 CURRENT CONFIGURATION")
    print("-" * 30)
    config = manager.get_access_summary()
    for key, value in config.items():
        status = "✅" if "VALID" in value or "ENABLED" in value else "❌"
        print(f"{status} {key.upper()}: {value}")
    print()
    
    # Validate API access
    print("🧪 VALIDATING ACCESS")
    print("-" * 30)
    access_results = manager.validate_access()
    
    # Create validation report
    validation_path = manager.output_dir / "validation_report.json"
    with open(validation_path, 'w') as f:
        json.dump(access_results, f, indent=2, sort_keys=True)
    
    # Display results
    all_valid = True
    for provider, results in access_results.items():
        status = results.get('status', 'unknown')
        key_valid = results.get('key_valid', False)
        
        if status == 'VALID' and key_valid:
            symbol = "✅"
        elif status in ['NO_KEY', 'NEED_KEY']:
            symbol = "❌"
            all_valid = False
        elif status in ['ERROR']:
            symbol = "⚠️"
            all_valid = False
        else:
            symbol = "🤔"
            all_valid = False
            
        print(f"{symbol} {provider.title()}:")
        print(f"   Status: {status}")
        print(f"   Key: {'Valid' if key_valid else 'Invalid/Missing'}")
        if results.get('error'):
            print(f"   Error: {results['error']}")
    
    print()
    
    # Overall assessment
    print("📊 OVERALL ASSESSMENT")
    print("-" * 30)
    
    required_providers = ['elsevier', 'wiley', 'springer']
    available = [p for p in required_providers 
                if p in access_results 
                and access_results[p].get('status') == 'VALID' 
                and access_results[p].get('key_valid')]
    
    print("Available providers: %d/%d" % (len(available), len(required_providers)))
    print("Achievable coverage: %d/%d%%" % (len(available), len(required_providers) * 85))
    
    if len(available) >= 2:
        print("✅ **READY FOR 224-PAPER DOWNLOAD**")
    elif len(available) == 1:
        print("⚠️ **PROCEED WITH PARTIAL COVERAGE**")
        missing = [p for p in required_providers if p not in available]
        print("  Need: %s" % str(missing))
    else:
        print("❌ **BLOCKED - NONE AVAILABLE**")
    
    print()
    
    # Next steps
    print("🎯 NEXT STEPS")
    print("-" * 20)
    
    if 'springer' in [p for p in required_providers if p not in available]:
        print("1. Register for Springer API key:")
        print("   - Go to: https://dev.springernature.com/")
        print("   - Use institutional email: brett.kroncke@vanderbilt.edu")
        print("   - Wait for email confirmation (check spam)")
    
    for provider in ['elsevier', 'wiley']:
        if provider in [p for p in required_providers if p not in available]:
            print("1. Update %s_API_KEY in .env file" % provider.upper())
    
    print("2. Run test validation:")
    print("   python scripts/validate_vandy_access.py")
    print("3. Once 2+ APIs are working, test with 10 papers")
    print("4. Scale to full 224-paper queue")
    
    print()
    print("Report saved: %s" % str(validation_path))

if __name__ == "__main__":
    main()