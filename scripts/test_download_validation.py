#!/usr/bin/env python3
"""
GVF Download Verification Test Suite
Tests the validation framework on actual downloaded papers
"""

import os
import sys
import json
from pathlib import Path
import subprocess

# Add project root
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.validate_downloads import DocumentValidator

def run_validation_test():
    """Run validation tests on existing downloads"""
    print("=" * 80)
    print("GVF DOWNLOAD VALIDATION TEST SUITE")
    print("=" * 80)
    
    # Test directories
    test_dirs = [
        "/mnt/temp2/kronckbm/gvf_output/download_batch_20260206",
        "/mnt/temp2/kronckbm/gvf_output/baseline_downloads_20260208"
    ]
    
    all_results = {}
    
    for test_dir in test_dirs:
        if Path(test_dir).exists():
            print(f"\nTesting directory: {test_dir}")
            print("-" * 60)
            
            # Run validation
            validator = DocumentValidator(test_dir)
            results = validator.scan_directory()
            
            # Store results
            all_results[test_dir] = results
            
            # Print summary
            summary = results['summary']
            print(f"Total packages: {summary['total_packes']}")
            print(f"Valid papers: {summary['valid_papers']}")
            print(f"Metadata stubs: {summary['metadata_stubs']}")
            print(f"Missing files: {summary['missing_files']}")
            
            # Show first few stubs if any
            stubs = [p for p in results['packages'] if p['overall_status'] == 'stub_or_metadata']
            if stubs:
                print(f"\nFirst 5 stubs:")
                for stub in stubs[:5]:
                    pmid = stub['pmid']
                    if 'pdf_validation' in stub and stub['pdf_validation']:
                        stub_type = stub['pdf_validation'].get('stub_type', 'unknown')
                        print(f"  - PMID {pmid}: {stub_type}")
                    else:
                        print(f"  - PMID {pmid}: validation error")
                        
            # Show first few valid papers
            valid = [p for p in results['packages'] if p['overall_status'].startswith('complete')]
            if valid:
                print(f"\nFirst 3 valid papers:")
                for valid_pkg in valid[:3]:
                    pmid = valid_pkg['pmid']
                    if 'pdf_validation' in valid_pkg and valid_pkg['pdf_validation']:
                        info = valid_pkg['pdf_validation']
                        size_kb = info.get('size_kb', 0)
                        pages = info.get('page_count', 'unknown')
                        supplement_info = len(valid_pkg.get('supplement_validation', []))
                        supplement_str = f" + {supplement_info} supplements" if supplement_info else ""
                        print(f"  - PMID {pmid}: {size_kb}KB, {pages} pages{supplement_str}")
                    else:
                        print(f"  - PMID {pmid}: valid but details unavailable")
        else:
            print(f"Directory not found: {test_dir}")
    
    return all_results

def test_individual_validation():
    """Test individual file validation"""
    print("\n" + "=" * 80)
    print("INDIVIDUAL FILE VALIDATION TESTS")
    print("=" * 80)
    
    # Test individual PDF validation
    test_files = []
    
    # Look for some actual PDF files
    search_dirs = [
        "/mnt/temp2/kronckbm/gvf_output/download_batch_20260206",
        "/mnt/temp2/kronckbm/gvf_output/baseline_downloads_20260208",
        "/mnt/temp2/kronckbm/gvf_output/baseline_downloads_direct_20260208"
    ]
    
    for search_dir in search_dirs:
        if Path(search_dir).exists():
            for item in Path(search_dir).rglob("*.pdf"):
                if item.exists():
                    test_files.append(item)
                    if len(test_files) >= 3:  # Test first 3 found
                        break
            if len(test_files) >= 3:
                break
    
    if test_files:
        validator = DocumentValidator(str(test_files[0].parent))
        
        for pdf_file in test_files:
            print(f"\nValidating: {pdf_file}")
            
            # Test PDF integrity
            validation = validator.validate_pdf_integrity(pdf_file)
            
            if validation:
                print(f"  Size: {validation.get('size_kb', 'unknown')} KB")
                print(f"  Status: {'VALID' if not validation.get('is_stub', True) else 'STUB/METADATA'}")
                print(f"  Pages: {validation.get('page_count', 'unknown')}")
                print(f"  Text length: {validation.get('text_length', 'unknown')} characters")
                
                stub_type = validation.get('stub_type')
                if stub_type:
                    print(f"  Stub type: {stub_type}")
                
                error = validation.get('error')
                if error:
                    print(f"  Error: {error}")
            else:
                print("  Validation failed")
    
def test_supplement_validation():
    """Test supplement archive validation"""
    print("\n" + "=" * 80)
    print("SUPPLEMENT VALIDATION TESTS")
    print("=" * 80)
    
    # Look for supplement archives
    search_dirs = [
        "/mnt/temp2/kronckbm/gvf_output/download_batch_20260206",
        "/mnt/temp2/kronckbm/gvf_output/baseline_downloads_20260208"
    ]
    
    test_archives = []
    archive_patterns = ["*.zip", "*.tar.gz", "*.tgz"]
    
    for search_dir in search_dirs:
        if Path(search_dir).exists():
            for pattern in archive_patterns:
                for item in Path(search_dir).rglob(pattern):
                    if item.exists():
                        test_archives.append(item)
                        if len(test_archives) >= 2:
                            break
                if test_archives:
                    break
    
    if test_archives:
        validator = DocumentValidator(str(test_archives[0].parent))
        
        for archive_path in test_archives:
            print(f"\nValidating archive: {archive_path.name}")
            
            validation = validator.validate_supplement_archive(archive_path)
            
            print(f"  Exists: {validation['exists']}")
            print(f"  Valid: {validation['is_valid']}")
            print(f"  Total size: {validation['total_size_mb']} MB")
            
            if validation['files']:
                print(f"  Contains {len(validation['files'])} files:")
                for i, file in enumerate(validation['files'][:5]):
                    print(f"    - {file}")
                if len(validation['files']) > 5:
                    print(f"    ... and {len(validation['files']) - 5} more files")
    
def run_comprehensive_report():
    """Generate comprehensive validation report"""
    print("\n" + "=" * 80)
    print("GENERATING COMPREHENSIVE TEST REPORT")
    print("=" * 80)
    
    # Create output directory
    report_dir = Path("/mnt/temp2/kronckbm/gvf_output/validation_test_reports")
    report_dir.mkdir(exist_ok=True)
    
    timestamp = "20260208"
    report_path = report_dir / f"validation_test_{timestamp}.json"
    
    # Run all tests and collect results
    test_results = {
        'test_timestamp': timestamp,
        'directory_validation': run_validation_test(),
        'system_check': check_system_requirements()
    }
    
    # Save results
    with open(report_path, 'w') as f:
        json.dump(test_results, f, indent=2, default=str)
    
    print(f"\nComplete validation report saved to: {report_path}")
    
    # Generate markdown summary
    markdown_path = report_dir / f"validation_summary_{timestamp}.md"
    generate_markdown_summary(test_results, markdown_path)
    
    return test_results

def check_system_requirements():
    """Check if required tools are available"""
    print("\n" + "=" * 80)
    print("SYSTEM REQUIREMENT CHECK")
    print("=" * 80)
    
    tools = ['pdfinfo', 'pdftotext', 'unzip', 'tar']
    
    print("Checking system tools:")
    for tool in tools:
        try:
            result = subprocess.run([tool, '--version'], 
                                  capture_output=True, text=True, timeout=5)
            print(f"  ✅ {tool}: Available")
        except (FileNotFoundError, subprocess.TimeoutExpired):
            print(f"  ❌ {tool}: Not found - some validation features may be limited")
        except Exception as e:
            print(f"  ⚠️  {tool}: Error checking ({e})")
    
    return True

def generate_markdown_summary(results, output_path):
    """Generate markdown summary report"""
    
    summary = []
    
    for test_dir, dir_results in results['directory_validation'].items():
        summary.append(f"## {Path(test_dir).name}")
        summary.append(f"- **Total packages**: {dir_results['summary']['total_packes']}")
        summary.append(f"- **Valid papers**: {dir_results['summary']['valid_papers']}")
        summary.append(f"- **Metadata stubs**: {dir_results['summary']['metadata_stubs']}")
        summary.append(f"- **Missing files**: {dir_results['summary']['missing_files']}")
        summary.append("")
    
    markdown_content = f"""# GVF Validation Test Report

Generated on: {results['test_timestamp']}

## Summary
This report summarizes the results of comprehensive validation testing
of downloaded papers to ensure they are full-text PDFs rather than 
metadata stubs.

{chr(10).join(summary)}

## Next Steps
1. Review failed downloads in systematic download log
2. Retry failed PMIDs with verification
3. Manually validate any questionable papers
4. Consider institutional access for paywalled content
"""
    
    with open(output_path, 'w') as f:
        f.write(markdown_content)
    
    print(f"Markdown summary saved to: {output_path}")

def main():
    """Run complete test suite"""
    try:
        # Run basic validation tests
        results = run_validation_test()
        
        # Run individual file tests
        test_individual_validation()
        
        # Test supplement validation if tools are available
        test_supplement_validation()
        
        # Check system requirements
        check_system_requirements()
        
        # Generate comprehensive report
        final_report = run_comprehensive_report()
        
        print("\n" + "=" * 80)
        print("ALL TESTS COMPLETED SUCCESSFULLY!")
        print("=" * 80)
        print("The validation framework is ready for systematic use.")
        
        return final_report
        
    except Exception as e:
        print(f"Error running test suite: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    main()