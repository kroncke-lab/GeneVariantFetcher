#!/usr/bin/env python3
"""
Deterministic pre-processor for FULL_CONTEXT.md files.
Runs BEFORE any LLM extraction. Cheap, fast, no API calls.

Purpose:
1. Strip excessive XML/HTML tags
2. Remove reference lists (bibliographies)
3. Remove figure/image markup that adds no variant data
4. Remove boilerplate (copyright, author bios, funding, COI)
5. Inject PubMed abstract as guaranteed baseline
6. Flag file quality (empty, abstract-only, full-text, garbled)

Output: CLEANED_CONTEXT.md — ready for LLM extraction
"""

import os
import re
import sys
import json
from pathlib import Path


class PaperPreprocessor:
    """Clean and categorize downloaded paper content."""
    
    # Sections to strip entirely
    STRIP_SECTIONS = [
        r'(?i)##?\s*(?:references|bibliography|works cited|literature cited)\s*\n.*',
        r'(?i)##?\s*(?:acknowledgments?|acknowledgements?)\s*\n.*?(?=\n##?\s|\Z)',
        r'(?i)##?\s*(?:conflict[s]? of interest|competing interests?|declaration[s]?)\s*\n.*?(?=\n##?\s|\Z)',
        r'(?i)##?\s*(?:author contributions?|author information|affiliations?)\s*\n.*?(?=\n##?\s|\Z)',
        r'(?i)##?\s*(?:funding|financial support|grant support)\s*\n.*?(?=\n##?\s|\Z)',
        r'(?i)##?\s*(?:supplementary|supporting) (?:information|material|data)\s*\n.*?(?=\n##?\s|\Z)',
    ]
    
    # XML/HTML patterns to clean
    XML_PATTERNS = [
        (r'<ce:para[^>]*>', ''),
        (r'</ce:para>', '\n'),
        (r'<ce:section[^>]*>', ''),
        (r'</ce:section>', ''),
        (r'<ce:section-title[^>]*>(.*?)</ce:section-title>', r'## \1'),
        (r'<ce:italic>(.*?)</ce:italic>', r'\1'),
        (r'<ce:bold>(.*?)</ce:bold>', r'**\1**'),
        (r'<ce:sup>(.*?)</ce:sup>', r'\1'),
        (r'<ce:sub>(.*?)</ce:sub>', r'\1'),
        (r'<ce:cross-ref[^>]*>(.*?)</ce:cross-ref>', r'\1'),
        (r'<ce:inter-ref[^>]*>(.*?)</ce:inter-ref>', r'\1'),
        (r'<xocs:[^>]*/?>', ''),
        (r'</xocs:[^>]*>', ''),
        (r'<ja:[^>]*/?>', ''),
        (r'</ja:[^>]*>', ''),
        (r'<sb:[^>]*>(.*?)</sb:[^>]*>', r'\1'),
        (r'<mml:[^>]*>(.*?)</mml:[^>]*>', r'\1'),
        # Generic cleanup
        (r'<[^>]{1,200}>', ''),  # Remove remaining tags (cap at 200 chars to avoid catastrophic backtracking)
        (r'&lt;', '<'),
        (r'&gt;', '>'),
        (r'&amp;', '&'),
        (r'&nbsp;', ' '),
        (r'&#\d+;', ''),
    ]
    
    # Lines that are just noise
    NOISE_PATTERNS = [
        r'^\s*https?://api\.(elsevier|springer|wiley)\.com/.*$',
        r'^\s*doi:.*$',
        r'^\s*\d+-s\d+\.\d+-\S+$',  # Elsevier PII
        r'^\s*text/plain\s*$',
        r'^\s*Copyright ©.*$',
        r'^\s*All rights reserved\.?\s*$',
        r'^\s*Published by .*$',
        r'^\s*\[Wiley PDF content for.*\]$',
        r'^\s*Received \d+.*Accepted \d+.*$',
        r'^\s*Available online.*$',
    ]
    
    # Figure/table captions - keep the text but mark it
    FIGURE_PATTERNS = [
        r'(?i)(?:figure|fig\.?)\s*\d+\.\s*(?:\.eps|\.tif|\.jpg|\.gif|\.png|image)',
    ]

    def __init__(self, abstracts_dir=None):
        self.abstracts_dir = abstracts_dir
    
    def classify(self, text):
        """Classify file content quality."""
        size = len(text)
        if size < 200:
            return "empty"
        if size < 3000:
            return "abstract_only"
        
        # Check for garbled/XML content
        tag_count = len(re.findall(r'<[^>]+>', text))
        text_chars = len(re.sub(r'<[^>]+>', '', text).strip())
        if tag_count > 0 and tag_count / max(1, text_chars) > 0.1:
            return "garbled_xml"
        
        # Check for section headers
        sections = re.findall(r'(?i)(introduction|methods|results|discussion|conclusion|background|patients)', text)
        if len(sections) >= 3:
            return "full_text"
        
        if size > 10000:
            return "content_no_sections"  # Probably a real paper, just no standard headers
        
        return "partial"
    
    def clean(self, text):
        """Apply deterministic cleaning to paper text."""
        # Step 1: Strip XML/HTML tags
        for pattern, replacement in self.XML_PATTERNS:
            text = re.sub(pattern, replacement, text, flags=re.DOTALL)
        
        # Step 2: Remove noise lines
        lines = text.split('\n')
        cleaned_lines = []
        for line in lines:
            skip = False
            for pattern in self.NOISE_PATTERNS:
                if re.match(pattern, line):
                    skip = True
                    break
            if not skip:
                cleaned_lines.append(line)
        text = '\n'.join(cleaned_lines)
        
        # Step 3: Strip reference sections (everything after "References" header)
        # Be careful - only strip if it looks like a real reference section
        ref_match = re.search(r'\n##?\s*(?:References|Bibliography|REFERENCES)\s*\n', text)
        if ref_match:
            # Check if there's substantial content after (numbered refs)
            after = text[ref_match.start():]
            ref_count = len(re.findall(r'^\s*\d+[\.\)]\s+', after, re.MULTILINE))
            if ref_count > 3:
                text = text[:ref_match.start()] + '\n\n[References section removed]\n'
        
        # Step 4: Collapse excessive whitespace
        text = re.sub(r'\n{4,}', '\n\n\n', text)
        text = re.sub(r' {3,}', ' ', text)
        
        # Step 5: Remove rebuild_full_context.py metadata comment (keep it but trim)
        # Keep the PMID/source info, strip the rest
        
        return text.strip()
    
    def inject_abstract(self, text, pmid):
        """Inject PubMed abstract at the top if available and not already present.
        
        Deduplication: Check if the abstract text is already in the document
        (full-text papers often include their own abstract). Only inject if
        the abstract content is NOT found in the existing text.
        """
        if not self.abstracts_dir:
            return text
        
        abs_file = os.path.join(self.abstracts_dir, f"{pmid}.txt")
        if not os.path.exists(abs_file):
            return text
        
        with open(abs_file, 'r') as f:
            abstract = f.read().strip()
        
        if not abstract:
            return text
        
        # Deduplication check: if the first 200 chars of abstract are already
        # present in the text (fuzzy match), skip injection to avoid doubling
        abstract_sample = abstract[:200].lower().replace('\n', ' ').replace('  ', ' ')
        text_lower = text[:10000].lower().replace('\n', ' ').replace('  ', ' ')  # Check first 10K chars
        
        if abstract_sample in text_lower:
            # Abstract already present in document
            return text
        
        header = "## PUBMED ABSTRACT (Authoritative Source)\n\n"
        return header + abstract + "\n\n---\n\n## FULL TEXT CONTENT\n\n" + text
    
    def process_file(self, filepath, pmid=None, output_dir=None):
        """Process a single FULL_CONTEXT.md file."""
        with open(filepath, 'r', errors='replace') as f:
            text = f.read()
        
        # Auto-detect PMID from filename if not provided
        if not pmid:
            basename = os.path.basename(filepath)
            pmid = basename.replace('_FULL_CONTEXT.md', '').replace('KCNH2_PMID_', '')
        
        # Classify
        classification = self.classify(text)
        
        # Clean
        cleaned = self.clean(text)
        
        # Inject abstract
        cleaned = self.inject_abstract(cleaned, pmid)
        
        # Stats
        stats = {
            "pmid": pmid,
            "original_size": len(text),
            "cleaned_size": len(cleaned),
            "classification": classification,
            "reduction_pct": round((1 - len(cleaned) / max(1, len(text))) * 100, 1),
        }
        
        # Write output
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            out_path = os.path.join(output_dir, f"{pmid}_CLEANED.md")
            with open(out_path, 'w') as f:
                f.write(cleaned)
        
        return cleaned, stats
    
    def process_directory(self, input_dir, output_dir):
        """Process all FULL_CONTEXT.md files in a directory."""
        os.makedirs(output_dir, exist_ok=True)
        
        files = [f for f in os.listdir(input_dir) if f.endswith('_FULL_CONTEXT.md')]
        print(f"Processing {len(files)} files...")
        
        all_stats = []
        classifications = {}
        
        for i, filename in enumerate(sorted(files)):
            filepath = os.path.join(input_dir, filename)
            pmid = filename.replace('_FULL_CONTEXT.md', '').replace('KCNH2_PMID_', '')
            
            _, stats = self.process_file(filepath, pmid, output_dir)
            all_stats.append(stats)
            
            cls = stats['classification']
            classifications[cls] = classifications.get(cls, 0) + 1
            
            if (i + 1) % 100 == 0:
                print(f"  {i+1}/{len(files)} processed...")
        
        # Summary
        print(f"\n=== PREPROCESSING SUMMARY ===")
        print(f"Total files: {len(files)}")
        for cls, count in sorted(classifications.items(), key=lambda x: -x[1]):
            print(f"  {cls}: {count} ({count*100//len(files)}%)")
        
        avg_reduction = sum(s['reduction_pct'] for s in all_stats) / max(1, len(all_stats))
        print(f"Average size reduction: {avg_reduction:.1f}%")
        
        # Save stats
        with open(os.path.join(output_dir, '_preprocessing_stats.json'), 'w') as f:
            json.dump({"summary": classifications, "files": all_stats}, f, indent=2)
        
        return all_stats


if __name__ == "__main__":
    input_dir = sys.argv[1] if len(sys.argv) > 1 else "/mnt/temp2/kronckbm/gvf_output/KCNH2/bulk_download"
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "/mnt/temp2/kronckbm/gvf_output/KCNH2/cleaned_context"
    abstracts_dir = sys.argv[3] if len(sys.argv) > 3 else "/mnt/temp2/kronckbm/gvf_output/KCNH2/pubmed_abstracts"
    
    preprocessor = PaperPreprocessor(abstracts_dir=abstracts_dir)
    preprocessor.process_directory(input_dir, output_dir)
