#!/usr/bin/env python3
"""
GVF Verification Framework - Full-text PDF Validation
Validates that downloads contain actual papers vs metadata stubs
"""

import os
import re
import hashlib
import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import subprocess

class DocumentValidator:
    """Validates downloaded documents are real papers, not metadata stubs"""
    
    def __init__(self, base_path: str):
        self.base_path = Path(base_path)
        self.validation_cache = {}
        
    def get_file_statistics(self, file_path: Path) -> Dict:
        """Get basic file statistics"""
        try:
            stat = file_path.stat()
            return {
                'size_bytes': stat.st_size,
                'size_kb': round(stat.st_size / 1024, 2),
                'exists': file_path.exists(),
                'extension': file_path.suffix.lower()
            }
        except Exception as e:
            return {'error': str(e), 'exists': False}
    
    def validate_pdf_integrity(self, pdf_path: Path) -> Dict:
        """Check if PDF is a real document vs stub"""
        is_valid = False
        content_check = {
            'is_stub': True,
            'stub_type': None,
            'page_count': 0,
            'text_length': 0,
            'has_metadata': False,
            'validation_hash': None
        }
        
        if not pdf_path.exists():
            return content_check
            
        try:
            # Check file size - stubs are usually <10KB
            file_size = pdf_path.stat().st_size
            content_check['size_bytes'] = file_size
            content_check['size_kb'] = round(file_size / 1024, 2)
            
            if file_size < 5000:  # Less than 5KB
                content_check['stub_type'] = 'tiny_file'
                return content_check
                
            # Try to get page count using pdfinfo
            try:
                result = subprocess.run(['pdfinfo', str(pdf_path)], 
                                      capture_output=True, text=True, timeout=10)
                if result.returncode == 0:
                    pages_match = re.search(r'Pages:\s*(\d+)', result.stdout)
                    if pages_match:
                        pages = int(pages_match.group(1))
                        content_check['page_count'] = pages
                        
                    # Check for actual content
                    if 'Title:' in result.stdout and 'Author:' in result.stdout:
                        content_check['has_metadata'] = True
                    
            except (subprocess.TimeoutExpired, FileNotFoundError):
                # Fallback to pdftotext
                pass
            
            # Extract text to check actual content
            try:
                text_result = subprocess.run(['pdftotext', str(pdf_path), '-'], 
                                           capture_output=True, text=True, timeout=15)
                if text_result.returncode == 0:
                    text = text_result.stdout
                    content_check['text_length'] = len(text)
                    
                    # Check for stub indicators
                    stub_indicators = [
                        'metadata not available',
                        'this article is not available',
                        'paywall',
                        'subscription required',
                        'login required',
                        'access denied',
                        'server not found',
                        'error 404',
                        'error 403',
                        'doi not found'
                    ]
                    
                    text_lower = text.lower()
                    for stub_marker in stub_indicators:
                        if stub_marker in text_lower:
                            content_check['stub_type'] = stub_marker
                            break
                    else:
                        # Not a stub if we have substantial text
                        if len(text) > 1000:  # At least 1000 characters
                            content_check['is_stub'] = False
                            
            except (subprocess.TimeoutExpired, FileNotFoundError):
                content_check['stub_type'] = 'extraction_error'
                
        except Exception as e:
            content_check['error'] = str(e)
            
        return content_check
    
    def validate_supplement_archive(self, archive_path: Path) -> Dict:
        """Validate supplement archive files"""
        result = {
            'exists': archive_path.exists(),
            'is_valid': False,
            'files': [],
            'total_size_mb': 0
        }
        
        if not archive_path.exists():
            return result
            
        try:
            # Check if it's a valid archive
            if archive_path.suffix.lower() in ['.zip', '.gz']:
                # Unzip to check contents
                import zipfile
                import tarfile
                
                if archive_path.suffix == '.zip':
                    with zipfile.ZipFile(archive_path, 'r') as zip_file:
                        files = zip_file.namelist()
                        total_size = sum(item.file_size for item in zip_file.filelist)
                        result.update({
                            'files': files,
                            'total_size_mb': round(total_size / (1024*1024), 2),
                            'is_valid': len(files) > 0
                        })
                        
                elif archive_path.suffix in ['.tar.gz', '.tgz']:
                    with tarfile.open(archive_path, 'r:gz') as tar_file:
                        files = tar_file.getnames()
                        total_size = sum(item.size for item in tar_file.getmembers())
                        result.update({
                            'files': files,
                            'total_size_mb': round(total_size / (1024*1024), 2),
                            'is_valid': len(files) > 0
                        })
        except Exception as e:
            result['error'] = str(e)
            result['is_valid'] = False
            
        return result
    
    def calculate_hash(self, file_path: Path, algorithm: str = 'md5') -> str:
        """Calculate file hash for duplicate detection"""
        try:
            hash_func = hashlib.md5() if algorithm == 'md5' else hashlib.sha256()
            with open(file_path, 'rb') as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    hash_func.update(chunk)
            return hash_func.hexdigest()
        except Exception:
            return None
    
    def validate_pmid_package(self, pmid: str, download_dir: Path) -> Dict:
        """Validate a complete PMID package"""
        package_result = {
            'pmid': pmid,
            'directory': str(download_dir),
            'pdf_validation': None,
            'supplement_validation': [],
            'overall_status': 'missing'
        }
        
        # Look for PDF
        pdf_candidates = list(download_dir.glob(f"{pmid}*.pdf")) + list(download_dir.glob("*.pdf"))
        
        if pdf_candidates:
            main_pdf = pdf_candidates[0]  # Take first PDF found
            package_result['pdf_validation'] = self.validate_pdf_integrity(main_pdf)
            
        # Look for supplements
        supplement_dirs = [d for d in download_dir.iterdir() if 'supplement' in str(d).lower()]
        supplement_archives = (list(download_dir.glob(f"{pmid}*.tar.gz")) + 
                              list(download_dir.glob(f"{pmid}*.zip")) + 
                              list(download_dir.glob("*supplement*.tar.gz")) + 
                              list(download_dir.glob("*supplement*.zip")))
        
        supplement_validations = []
        
        # Check supplement directories
        for supp_dir in supplement_dirs:
            if supp_dir.is_dir():
                archive_files = list(supp_dir.glob("*"))
                for archive_path in archive_files:
                    if archive_path.is_file():
                        validation = self.validate_supplement_archive(archive_path)
                        validation['archive_path'] = str(archive_path)
                        supplement_validations.append(validation)
        
        # Check supplement archives directly
        for archive_path in supplement_archives:
            if archive_path.exists():
                validation = self.validate_supplement_archive(archive_path)
                validation['archive_path'] = str(archive_path)
                supplement_validations.append(validation)
        
        package_result['supplement_validation'] = supplement_validations
        
        # Determine overall status
        pdf_validation = package_result.get('pdf_validation', {})
        if pdf_validation:
            pdf_valid = not pdf_validation.get('is_stub', True)
            if pdf_valid:
                has_supplements = len(supplement_validations) > 0
                if has_supplements:
                    package_result['overall_status'] = 'complete_plus_supplements'
                else:
                    package_result['overall_status'] = 'complete'
            elif pdf_validation.get('exists', False):
                package_result['overall_status'] = 'stub_or_metadata'
            else:
                package_result['overall_status'] = 'missing'
        else:
            # Check if there are any valid supplement archives as fallback
            valid_supplements = [v for v in supplement_validations if v.get('is_valid', False)]
            if valid_supplements:
                package_result['overall_status'] = 'supplements_only'
            else:
                package_result['overall_status'] = 'missing'
            
        return package_result
    
    def scan_directory(self, base_download_path: str = None) -> Dict:
        """Scan entire download directory and validate all packages"""
        if base_download_path is None:
            base_download_path = self.base_path
            
        results = {
            'packages': [],
            'summary': {
                'total_packes': 0,
                'valid_papers': 0,
                'metadata_stubs': 0,
                'missing_files': 0,
                'with_supplements': 0
            }
        }
        
        # Find all download directories
        download_dirs = [d for d in Path(base_download_path).glob("*/") if d.is_dir()]
        
        for download_dir in download_dirs:
            # Look for PMID identifiers
            pmid_match = re.search(r'(\d+)', str(download_dir.name))
            if pmid_match:
                pmid = pmid_match.group(1)
                validation_result = self.validate_pmid_package(pmid, download_dir)
                results['packages'].append(validation_result)
                
                # Update summary
                results['summary']['total_packes'] += 1
                
                status = validation_result['overall_status']
                if status == 'complete_plus_supplements':
                    results['summary']['valid_papers'] += 1
                    results['summary']['with_supplements'] += 1
                elif status == 'complete':
                    results['summary']['valid_papers'] += 1
                elif status == 'stub_or_metadata':
                    results['summary']['metadata_stubs'] += 1
                elif status == 'missing':
                    results['summary']['missing_files'] += 1
        
        return results
    
    def generate_report(self, validation_results: Dict, output_path: str) -> None:
        """Generate a detailed PDF validation report"""
        report_lines = [
            "# GVF Full-Text Paper Validation Report",
            "",
            "## Summary Statistics",
            f"- Total packages: {validation_results['summary']['total_packes']}",
            f"- Valid papers: {validation_results['summary']['valid_papers']}",
            f"- Metadata stubs: {validation_results['summary']['metadata_stubs']}",
            f"- Missing files: {validation_results['summary']['missing_files']}",
            f"- With supplements: {validation_results['summary']['with_supplements']}",
            "",
            "## Detail Analysis",
            ""
        ]
        
        stubs = [p for p in validation_results['packages'] if p['overall_status'] == 'stub_or_metadata']
        missing = [p for p in validation_results['packages'] if p['overall_status'] == 'missing']
        valid = [p for p in validation_results['packages'] if p['overall_status'].startswith('complete')]
        
        if stubs:
            report_lines.append("### Metadata Stubs (Need Re-download)")
            for stub in stubs:
                report_lines.append(f"- PMID {stub['pmid']}: {stub['pdf_validation'].get('stub_type', 'unknown stub')}")
            report_lines.append("")
        
        if missing:
            report_lines.append("### Missing Files")
            for missing_pkg in missing:
                report_lines.append(f"- PMID {missing_pkg['pmid']}: No download found")
            report_lines.append("")
        
        if valid:
            report_lines.append("### Valid Papers")
            for valid_pkg in valid[:20]:  # Show first 20
                report_lines.append(f"- PMID {valid_pkg['pmid']}: Verified full-text")  
            if len(valid) > 20:
                report_lines.append(f"... and {len(valid) - 20} more")
        
        # Save report
        with open(output_path, 'w') as f:
            f.write('\n'.join(report_lines))
        print(f"Validation report saved to: {output_path}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        download_path = sys.argv[1]
    else:
        download_path = "/mnt/temp2/kronckbm/gvf_output/"
        
    validator = DocumentValidator(download_path)
    results = validator.scan_directory()
    
    # Generate report
    report_path = str(Path(download_path) / "validation_report.md")
    validator.generate_report(results, report_path)
    
    # Save detailed results
    json_path = str(Path(download_path) / "validation_results.json")
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"Validation complete. Found:")
    print(f"  Valid papers: {results['summary']['valid_papers']}")
    print(f"  Metadata stubs: {results['summary']['metadata_stubs']}")
    print(f"  Missing files: {results['summary']['missing_files']}")
    print(f"  With supplements: {results['summary']['with_supplements']}")