#!/usr/bin/env python3
"""
Rebuild FULL_CONTEXT.md files using intelligent source selection.

This script fixes the #1 recall problem identified in the architecture reassessment:
- FULL_CONTEXT.md files were being generated from PDFs that produced garbage
- Better sources (fulltext.txt, BioC XML, publisher XML) were ignored

Source Priority (highest to lowest quality):
1. Publisher XML/HTML (Elsevier, Wiley, Springer) - structured, complete, has tables
2. BioC JSON/XML from NCBI - well-structured text with annotations
3. fulltext.txt from NCBI - plain text, usually clean
4. PMC XML - structured but may lack tables
5. PDF-extracted text - often garbled or truncated

Quality Checks:
- Reject sources with avg line length < 10 chars (garbled)
- Reject sources with total length < 500 chars (truncated)
- Prefer sources with more variant-like patterns

Created: 2026-02-10 for GVF pipeline improvement
"""

import sys
import os
import re
import json
import logging
import shutil
from pathlib import Path
from typing import Optional, Dict, List, Tuple
from dataclasses import dataclass
from datetime import datetime
import xml.etree.ElementTree as ET

# Add GVF to path
gvf_root = Path("/mnt/temp2/kronckbm/gitrepos/GeneVariantFetcher")
sys.path.insert(0, str(gvf_root))

from dotenv import load_dotenv
load_dotenv(gvf_root / ".env")

from harvesting.format_converters import FormatConverter

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


@dataclass
class SourceScore:
    """Score for a potential source file."""
    path: Path
    source_type: str  # elsevier_xml, bioc_xml, fulltext_txt, pmc_xml, pdf, etc.
    priority: int  # Lower is better (1-6)
    length: int
    avg_line_length: float
    variant_mentions: int  # Count of variant-like patterns
    quality_score: float  # Combined score
    is_valid: bool
    rejection_reason: Optional[str] = None
    content: Optional[str] = None


class SourceSelector:
    """Intelligently select the best source for FULL_CONTEXT generation."""
    
    # Patterns that indicate variant data
    VARIANT_PATTERNS = [
        r'\b[A-Z]\d{1,4}[A-Z]\b',  # e.g., W1001X, R850Q
        r'\bp\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}\b',  # p.Trp1001Ter
        r'\bc\.\d+[ACGT]>[ACGT]\b',  # c.123A>G
        r'\brs\d{6,}\b',  # rsIDs
        r'\bmutation\b',
        r'\bvariant\b',
        r'\bpathogenic\b',
        r'\bbenign\b',
    ]
    
    def __init__(self, converter: Optional[FormatConverter] = None):
        self.converter = converter or FormatConverter()
        
    def find_all_sources(self, pmid: str, search_dirs: List[Path]) -> List[Path]:
        """Find all potential source files for a PMID across directories."""
        sources = []
        patterns = [
            f"{pmid}_*.txt",
            f"{pmid}_*.xml",
            f"{pmid}_*.json",
            f"{pmid}.pdf",
            f"*_{pmid}_*.xml",
            f"*_{pmid}.xml",
            f"pmc_{pmid}.xml",
            f"PMID_{pmid}_*.*",
            f"PMID{pmid}_*.*",
            f"*_PMID_{pmid}*.*",
        ]
        
        for search_dir in search_dirs:
            if not search_dir.exists():
                continue
            
            # Check direct patterns
            for pattern in patterns:
                sources.extend(search_dir.glob(pattern))
            
            # Check in PMID subdirectories
            pmid_dir = search_dir / f"PMID_{pmid}"
            if pmid_dir.exists():
                for pattern in patterns:
                    sources.extend(pmid_dir.glob(pattern))
            
            pmid_dir2 = search_dir / pmid
            if pmid_dir2.exists():
                for pattern in patterns:
                    sources.extend(pmid_dir2.glob(pattern))
        
        # Remove duplicates and exclude non-content files
        seen = set()
        unique = []
        exclude_patterns = [
            "FULL_CONTEXT",
            "DATA_ZONES",
            "_metadata.json",  # Metadata files, not content
            "_download_result.json",
            "pmid_status",
            ".json",  # Generally exclude JSON unless it's BioC
        ]
        
        for p in sources:
            name = p.name.lower()
            
            # Skip if already seen
            if p in seen:
                continue
                
            # Skip excluded patterns (but allow bioc.json)
            if "bioc" not in name:
                if any(excl in name for excl in exclude_patterns):
                    continue
            
            seen.add(p)
            unique.append(p)
        
        return unique
    
    def classify_source(self, path: Path) -> Tuple[str, int]:
        """Classify a source file and assign priority (1=best, 6=worst)."""
        name = path.name.lower()
        
        # Skip annotation/metadata files entirely
        if "pubtator" in name:
            return "pubtator", 99  # Will be filtered as too low priority
        if "metadata" in name:
            return "metadata", 99
        
        # Publisher XML - highest priority
        if "elsevier" in name and name.endswith(".xml"):
            return "elsevier_xml", 1
        if "wiley" in name and name.endswith(".xml"):
            return "wiley_xml", 1
        if "springer" in name and name.endswith(".xml"):
            return "springer_xml", 1
        
        # BioC - second priority (well-structured text from NCBI)
        if "_bioc." in name or name.endswith("_bioc.xml") or name.endswith("_bioc.json"):
            return "bioc", 2
        
        # Fulltext.txt from NCBI - usually clean plain text
        if name.endswith("_fulltext.txt") or "_fulltext." in name:
            return "fulltext_txt", 3
        
        # PMC XML - structured but may be incomplete
        if "pmc" in name and name.endswith(".xml"):
            return "pmc_xml", 4
        
        # Other text files
        if name.endswith(".txt"):
            return "text", 4
        
        # PDF - often garbled, use as last resort
        if name.endswith(".pdf"):
            return "pdf", 6
        
        # Other XML files
        if name.endswith(".xml"):
            return "other_xml", 5
        
        # JSON files (rare, usually metadata)
        if name.endswith(".json"):
            return "json", 99  # Skip most JSON
        
        return "unknown", 99
    
    def count_variants(self, text: str) -> int:
        """Count variant-like patterns in text."""
        count = 0
        for pattern in self.VARIANT_PATTERNS:
            count += len(re.findall(pattern, text, re.IGNORECASE))
        return count
    
    def read_source_content(self, path: Path, source_type: str) -> Optional[str]:
        """Read and convert source to markdown/text."""
        try:
            if source_type == "elsevier_xml":
                return self._convert_elsevier_xml(path)
            elif source_type == "wiley_xml":
                return self._convert_wiley_xml(path)
            elif source_type == "springer_xml":
                return self._convert_springer_xml(path)
            elif source_type == "bioc":
                return self._convert_bioc(path)
            elif source_type in ("fulltext_txt", "text"):
                return path.read_text(encoding="utf-8", errors="ignore")
            elif source_type == "pmc_xml":
                xml_content = path.read_text(encoding="utf-8", errors="ignore")
                return self.converter.xml_to_markdown(xml_content)
            elif source_type == "pdf":
                return self.converter.pdf_to_markdown(path)
            elif source_type == "other_xml":
                xml_content = path.read_text(encoding="utf-8", errors="ignore")
                return self.converter.xml_to_markdown(xml_content)
            elif source_type == "json":
                data = json.loads(path.read_text(encoding="utf-8", errors="ignore"))
                return self._json_to_text(data)
            else:
                return path.read_text(encoding="utf-8", errors="ignore")
        except Exception as e:
            logger.warning(f"Error reading {path}: {e}")
            return None
    
    def _convert_elsevier_xml(self, path: Path) -> str:
        """Convert Elsevier XML to markdown, preserving tables."""
        try:
            content = path.read_text(encoding="utf-8", errors="ignore")
            
            # Parse with namespace handling
            root = ET.fromstring(content)
            
            # Find all namespaces
            nsmap = {}
            for elem in root.iter():
                for ns in re.findall(r'xmlns(?::\w+)?="([^"]+)"', ET.tostring(elem, encoding='unicode')):
                    pass
            
            markdown_parts = ["# MAIN TEXT\n\n"]
            
            # Extract title
            title_elem = root.find('.//{http://purl.org/dc/elements/1.1/}title')
            if title_elem is not None and title_elem.text:
                markdown_parts.append(f"## {title_elem.text.strip()}\n\n")
            
            # Extract all text content
            def extract_text(elem):
                """Recursively extract text from element."""
                text_parts = []
                if elem.text:
                    text_parts.append(elem.text.strip())
                for child in elem:
                    child_text = extract_text(child)
                    if child_text:
                        text_parts.append(child_text)
                    if child.tail:
                        text_parts.append(child.tail.strip())
                return " ".join(text_parts)
            
            # Find all sections and paragraphs
            for elem in root.iter():
                tag = elem.tag.split('}')[-1] if '}' in elem.tag else elem.tag
                
                if tag == 'section':
                    # Section title
                    title = elem.find('.//{http://www.elsevier.com/xml/common/dtd}section-title')
                    if title is not None:
                        title_text = extract_text(title)
                        if title_text:
                            markdown_parts.append(f"\n### {title_text}\n\n")
                
                elif tag == 'para' or tag == 'simple-para':
                    para_text = extract_text(elem)
                    if para_text and len(para_text) > 20:
                        markdown_parts.append(f"{para_text}\n\n")
                
                elif tag == 'table':
                    # Extract table content
                    table_md = self._extract_elsevier_table(elem)
                    if table_md:
                        markdown_parts.append(f"\n{table_md}\n\n")
            
            result = "".join(markdown_parts)
            
            # If structured extraction failed, fall back to text extraction
            if len(result.strip()) < 1000:
                from bs4 import BeautifulSoup
                soup = BeautifulSoup(content, 'xml')
                text = soup.get_text(separator='\n', strip=True)
                if len(text) > len(result):
                    result = "# MAIN TEXT\n\n" + text
            
            return result
            
        except Exception as e:
            logger.warning(f"Error parsing Elsevier XML {path}: {e}")
            # Fallback: extract raw text
            try:
                from bs4 import BeautifulSoup
                content = path.read_text(encoding="utf-8", errors="ignore")
                soup = BeautifulSoup(content, 'xml')
                return "# MAIN TEXT\n\n" + soup.get_text(separator='\n', strip=True)
            except:
                return None
    
    def _extract_elsevier_table(self, table_elem) -> str:
        """Extract table from Elsevier XML as markdown."""
        try:
            rows = []
            
            # Find all rows
            for row in table_elem.iter():
                tag = row.tag.split('}')[-1] if '}' in row.tag else row.tag
                if tag == 'row':
                    cells = []
                    for entry in row.iter():
                        entry_tag = entry.tag.split('}')[-1] if '}' in entry.tag else entry.tag
                        if entry_tag == 'entry':
                            text = "".join(entry.itertext()).strip()
                            cells.append(text.replace('|', '/'))
                    if cells:
                        rows.append(cells)
            
            if not rows:
                return ""
            
            # Build markdown table
            max_cols = max(len(r) for r in rows)
            md_lines = []
            for i, row in enumerate(rows):
                padded = row + [''] * (max_cols - len(row))
                md_lines.append('| ' + ' | '.join(padded) + ' |')
                if i == 0:
                    md_lines.append('|' + '---|' * max_cols)
            
            return '\n'.join(md_lines)
        except:
            return ""
    
    def _convert_wiley_xml(self, path: Path) -> str:
        """Convert Wiley XML to markdown."""
        try:
            from bs4 import BeautifulSoup
            content = path.read_text(encoding="utf-8", errors="ignore")
            soup = BeautifulSoup(content, 'xml')
            return "# MAIN TEXT\n\n" + soup.get_text(separator='\n', strip=True)
        except Exception as e:
            logger.warning(f"Error parsing Wiley XML {path}: {e}")
            return None
    
    def _convert_springer_xml(self, path: Path) -> str:
        """Convert Springer XML to markdown."""
        try:
            from bs4 import BeautifulSoup
            content = path.read_text(encoding="utf-8", errors="ignore")
            soup = BeautifulSoup(content, 'xml')
            return "# MAIN TEXT\n\n" + soup.get_text(separator='\n', strip=True)
        except Exception as e:
            logger.warning(f"Error parsing Springer XML {path}: {e}")
            return None
    
    def _convert_bioc(self, path: Path) -> str:
        """Convert BioC JSON/XML to text."""
        try:
            content = path.read_text(encoding="utf-8", errors="ignore")
            
            if path.suffix == ".json":
                data = json.loads(content)
                return self._bioc_json_to_text(data)
            else:
                # BioC XML
                from bs4 import BeautifulSoup
                soup = BeautifulSoup(content, 'xml')
                
                text_parts = []
                for passage in soup.find_all('passage'):
                    text_elem = passage.find('text')
                    if text_elem and text_elem.string:
                        text_parts.append(text_elem.string.strip())
                
                if text_parts:
                    return "# MAIN TEXT\n\n" + "\n\n".join(text_parts)
                
                # Fallback
                return "# MAIN TEXT\n\n" + soup.get_text(separator='\n', strip=True)
                
        except Exception as e:
            logger.warning(f"Error parsing BioC {path}: {e}")
            return None
    
    def _bioc_json_to_text(self, data: dict) -> str:
        """Convert BioC JSON to text."""
        text_parts = ["# MAIN TEXT\n\n"]
        
        if isinstance(data, dict):
            if 'documents' in data:
                for doc in data['documents']:
                    if 'passages' in doc:
                        for passage in doc['passages']:
                            if 'text' in passage:
                                text_parts.append(passage['text'])
            elif 'passages' in data:
                for passage in data['passages']:
                    if 'text' in passage:
                        text_parts.append(passage['text'])
        
        return "\n\n".join(text_parts)
    
    def _json_to_text(self, data) -> str:
        """Convert JSON data to readable text."""
        if isinstance(data, str):
            return data
        elif isinstance(data, dict):
            if 'text' in data:
                return data['text']
            elif 'content' in data:
                return data['content']
            elif 'fulltext' in data:
                return data['fulltext']
            else:
                return json.dumps(data, indent=2)
        else:
            return str(data)
    
    def score_source(self, path: Path) -> SourceScore:
        """Score a source file for quality."""
        source_type, priority = self.classify_source(path)
        
        # Skip metadata/annotation files entirely
        if priority >= 99:
            return SourceScore(
                path=path,
                source_type=source_type,
                priority=priority,
                length=0,
                avg_line_length=0,
                variant_mentions=0,
                quality_score=0,
                is_valid=False,
                rejection_reason=f"Skipped non-content file type: {source_type}"
            )
        
        # Read content
        content = self.read_source_content(path, source_type)
        
        if content is None:
            return SourceScore(
                path=path,
                source_type=source_type,
                priority=priority,
                length=0,
                avg_line_length=0,
                variant_mentions=0,
                quality_score=0,
                is_valid=False,
                rejection_reason="Failed to read content"
            )
        
        # Calculate metrics
        length = len(content)
        lines = [line for line in content.split('\n') if line.strip()]
        avg_line_length = sum(len(line) for line in lines) / max(len(lines), 1)
        variant_mentions = self.count_variants(content)
        
        # Quality checks
        is_valid = True
        rejection_reason = None
        
        if length < 500:
            is_valid = False
            rejection_reason = f"Too short ({length} chars < 500 min)"
        elif avg_line_length < 10:
            is_valid = False
            rejection_reason = f"Garbled text (avg line {avg_line_length:.1f} chars < 10)"
        elif source_type == "pdf" and avg_line_length < 20:
            # PDFs are more prone to garbled text
            is_valid = False
            rejection_reason = f"Likely garbled PDF (avg line {avg_line_length:.1f} chars < 20)"
        
        # Calculate composite quality score
        # Higher is better
        quality_score = 0
        if is_valid:
            # Base score from priority (inverse - lower priority number = higher score)
            quality_score = (7 - priority) * 1000
            # Add length bonus (log scale to not over-weight length)
            quality_score += min(length / 100, 500)
            # Add variant mentions bonus
            quality_score += variant_mentions * 10
            # Add avg line length bonus (indicates readable text)
            quality_score += min(avg_line_length, 100)
        
        return SourceScore(
            path=path,
            source_type=source_type,
            priority=priority,
            length=length,
            avg_line_length=avg_line_length,
            variant_mentions=variant_mentions,
            quality_score=quality_score,
            is_valid=is_valid,
            rejection_reason=rejection_reason,
            content=content if is_valid else None
        )
    
    def select_best_source(self, pmid: str, search_dirs: List[Path]) -> Optional[SourceScore]:
        """Find and select the best source for a PMID."""
        sources = self.find_all_sources(pmid, search_dirs)
        
        if not sources:
            logger.warning(f"PMID {pmid}: No sources found")
            return None
        
        # Score all sources
        scores = []
        for path in sources:
            score = self.score_source(path)
            scores.append(score)
            if score.is_valid:
                logger.debug(f"  {path.name}: {score.source_type}, score={score.quality_score:.0f}, {score.length} chars, {score.variant_mentions} variants")
            else:
                logger.debug(f"  {path.name}: REJECTED - {score.rejection_reason}")
        
        # Filter to valid sources
        valid_scores = [s for s in scores if s.is_valid]
        
        if not valid_scores:
            logger.warning(f"PMID {pmid}: All {len(sources)} sources failed quality checks")
            for s in scores:
                logger.warning(f"  - {s.path.name}: {s.rejection_reason}")
            return None
        
        # Sort by quality score (highest first)
        valid_scores.sort(key=lambda s: s.quality_score, reverse=True)
        
        best = valid_scores[0]
        logger.info(f"PMID {pmid}: Selected {best.path.name} ({best.source_type}, {best.length} chars, {best.variant_mentions} variants)")
        
        return best


def rebuild_full_context(
    pmid: str,
    output_dir: Path,
    search_dirs: List[Path],
    selector: SourceSelector,
    backup: bool = True
) -> Tuple[bool, str]:
    """
    Rebuild FULL_CONTEXT.md for a single PMID.
    
    Returns:
        Tuple of (success, message)
    """
    output_file = output_dir / f"{pmid}_FULL_CONTEXT.md"
    
    # Backup existing file
    if backup and output_file.exists():
        backup_file = output_dir / f"{pmid}_FULL_CONTEXT_old.md"
        if not backup_file.exists():  # Don't overwrite existing backup
            shutil.copy2(output_file, backup_file)
    
    # Find best source
    best_source = selector.select_best_source(pmid, search_dirs)
    
    if best_source is None:
        return False, "No valid sources found"
    
    # Write new FULL_CONTEXT.md
    content = best_source.content
    
    # Add header with source info
    header = f"""<!-- 
FULL_CONTEXT.md rebuilt by rebuild_full_context.py
Date: {datetime.now().isoformat()}
Source: {best_source.path}
Source type: {best_source.source_type}
Quality score: {best_source.quality_score:.0f}
Length: {best_source.length} chars
Variant mentions: {best_source.variant_mentions}
-->

"""
    
    output_file.write_text(header + content, encoding="utf-8")
    
    return True, f"Rebuilt from {best_source.source_type} ({best_source.length} chars)"


def main():
    """Main function to rebuild all FULL_CONTEXT.md files."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Rebuild FULL_CONTEXT.md files with intelligent source selection")
    parser.add_argument("--pmid", help="Process a single PMID")
    parser.add_argument("--output-dir", type=Path, required=True, help="Output directory for FULL_CONTEXT.md files")
    parser.add_argument("--search-dirs", type=Path, nargs="+", help="Directories to search for sources")
    parser.add_argument("--no-backup", action="store_true", help="Don't backup existing files")
    parser.add_argument("--test-only", action="store_true", help="Only test known-bad PMIDs")
    args = parser.parse_args()
    
    # Default search directories
    if args.search_dirs:
        search_dirs = args.search_dirs
    else:
        base = Path("/mnt/temp2/kronckbm/gvf_output")
        search_dirs = [
            base / "KCNH2" / "fulltext_downloads",
            base / "KCNH2" / "europepmc_downloads",
            base / "KCNH2" / "tier2_downloads",
            base / "KCNH2" / "tier3_downloads" / "direct_api_20260207_053427",
            base / "KCNH2" / "browser_downloads",
            base / "verified_downloads_20260208",
            base / "papers",
        ]
    
    # Ensure output dir exists
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    selector = SourceSelector()
    
    # Test PMIDs first
    test_pmids = ["24218437", "15840476", "19038855"]
    
    if args.pmid:
        pmids_to_process = [args.pmid]
    elif args.test_only:
        pmids_to_process = test_pmids
    else:
        # Find all PMIDs with sources
        all_pmids = set()
        for search_dir in search_dirs:
            if not search_dir.exists():
                continue
            for f in search_dir.rglob("*"):
                if f.is_file():
                    # Extract PMID from filename
                    match = re.search(r'(\d{7,8})', f.name)
                    if match:
                        all_pmids.add(match.group(1))
        pmids_to_process = sorted(all_pmids)
    
    logger.info(f"Processing {len(pmids_to_process)} PMIDs")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Search directories: {[str(d) for d in search_dirs]}")
    
    results = {"success": 0, "failed": 0, "details": []}
    
    for i, pmid in enumerate(pmids_to_process, 1):
        logger.info(f"\n[{i}/{len(pmids_to_process)}] Processing PMID {pmid}...")
        
        success, message = rebuild_full_context(
            pmid=pmid,
            output_dir=args.output_dir,
            search_dirs=search_dirs,
            selector=selector,
            backup=not args.no_backup
        )
        
        if success:
            results["success"] += 1
            logger.info(f"  ✓ {message}")
        else:
            results["failed"] += 1
            logger.warning(f"  ✗ {message}")
        
        results["details"].append({
            "pmid": pmid,
            "success": success,
            "message": message
        })
    
    # Summary
    logger.info("\n" + "=" * 60)
    logger.info("SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Total PMIDs: {len(pmids_to_process)}")
    logger.info(f"Successful: {results['success']}")
    logger.info(f"Failed: {results['failed']}")
    
    # Save results
    results_file = args.output_dir / "rebuild_results.json"
    results_file.write_text(json.dumps(results, indent=2), encoding="utf-8")
    logger.info(f"\nResults saved to: {results_file}")


if __name__ == "__main__":
    main()
