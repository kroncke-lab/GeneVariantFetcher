#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
KCNH2 Tier 3 Variant Extraction Script
Extracts KCNH2/HERG variants from PMC XML files using pattern matching
"""

import re
import json
import os
import glob
from pathlib import Path
from xml.dom import minidom

class KCNH2VariantExtractor:
    def __init__(self, gold_standard_path):
        """Initialize with gold standard data"""
        self.gold_standard_path = gold_standard_path
        self.gold_standard_variants = self.load_gold_standard()
        
        # Patterns for variant detection
        self.protein_pattern = re.compile(r'[A-Z][0-9]+[A-Z]')
        self.frameshift_pattern = re.compile(r'[A-Z][0-9]+fs')
        
        # KCNH2 specific terms
        self.kcnh2_terms = {'KCNH2', 'HERG', 'hERG', 'Kv11.1', 'erg'}
        
    def load_gold_standard(self):
        """Load variants from gold standard JSON"""
        try:
            with open(self.gold_standard_path, 'r') as f:
                data = json.load(f)
                return set(str(v) for v in data.get('matched_variants', []))
        except Exception as e:
            print("Warning: Could not load gold standard: %s" % e)
            return set()
    
    def extract_text_from_xml(self, xml_path):
        """Extract text content from PMC XML file"""
        try:
            # Read XML and extract text content
            with open(xml_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            # Simple text extraction by removing XML tags
            text = re.sub(r'<[^>]+>', ' ', content)
            text = re.sub(r'\s+', ' ', text)
            
            return text.strip()
            
        except Exception as e:
            print("Error parsing %s: %s" % (xml_path, e))
            return ""
    
    def find_variants(self, text, paper_id):
        """Find variants in text using patterns"""
        variants = []
        
        if not text:
            return variants
            
        text = re.sub(r'\s+', ' ', text)
        text_lower = text.lower()
        
        # Check if paper mentions KCNH2 or related terms (case-insensitive)
        kcnh2_found = any(term.lower() in text_lower for term in self.kcnh2_terms)
        
        if not kcnh2_found:
            return variants
        
        # Find protein variants using regex
        # Pattern: One letter AA, number(s), one letter AA
        for match in re.finditer(r'(?<!\w)([A-Z])([0-9]{1,4})([A-Z])(?!\w)', text):
            ref_aa, pos, var_aa = match.groups()
            variant = ref_aa + pos + var_aa
            
            # Validate position makes sense (reasonable protein length)
            if int(pos) > 2000:
                continue
                
            # Get context
            context_start = max(0, match.start() - 15)
            context_end = min(len(text), match.end() + 15)
            context = text[context_start:context_end]
            
            variants.append({
                'variant': variant,
                'type': 'protein_change',
                'position': int(pos),
                'context': context.strip(),
                'paper_id': paper_id
            })
        
        # Find frameshift variants
        for match in re.finditer(r'(?<!\w)([A-Z])([0-9]{1,4})fs(?!\w)', text):
            ref_aa, pos = match.groups()
            variant = ref_aa + pos + 'fs'
            
            context_start = max(0, match.start() - 15)
            context_end = min(len(text), match.end() + 15)
            context = text[context_start:context_end]
            
            variants.append({
                'variant': variant,
                'type': 'frameshift',
                'position': int(pos),
                'context': context.strip(),
                'paper_id': paper_id
            })
        
        return variants
    
    def process_papers(self, xml_dir):
        """Process all XML files in directory"""
        results = {
            'total_papers': 0,
            'papers_with_variants': 0,
            'variants_found': [],
            'unique_variants': set(),
            'new_matches': set(),
            'paper_summary': {}
        }
        
        # Find all XML files
        xml_files = []
        for subdir in ['direct_api_20260207_053427', 'remaining']:
            search_path = os.path.join(xml_dir, subdir, "*.xml")
            xml_files.extend(glob.glob(search_path))
        
        print("Found %d XML files in subdirectories" % len(xml_files))
        
        # Process first 15 papers for Tier 3 analysis
        processed = 0
        for xml_path in xml_files[:15]:
            paper_id = Path(xml_path).stem
            print("Processing %s..." % paper_id)
            
            results['total_papers'] += 1
            processed += 1
            
            # Extract text
            text = self.extract_text_from_xml(xml_path)
            if not text:
                print("  No text extracted")
                continue
            
            # Find variants
            variants = self.find_variants(text, paper_id)
            
            if variants:
                results['papers_with_variants'] += 1
                
                # Process each variant
                seen_in_paper = set()
                for var_data in variants:
                    variant = var_data['variant']
                    
                    # Avoid duplicates in same paper
                    if variant in seen_in_paper:
                        continue
                    seen_in_paper.add(variant)
                    
                    # Check against gold standard
                    normalized = self.normalize_variant(variant)
                    var_data['normalized'] = normalized
                    
                    if normalized in self.gold_standard_variants:
                        results['new_matches'].add(normalized)
                    
                    results['variants_found'].append(var_data)
                    results['unique_variants'].add(normalized)
                
                results['paper_summary'][paper_id] = {
                    'variant_count': len(seen_in_paper),
                    'variants': list(seen_in_paper)
                }
                
                print("  Found %d variants: %s" % (len(seen_in_paper), ', '.join(seen_in_paper)))
            else:
                print("  No variants found")
            
            if processed >= 15:
                break
        
        return results
    
    def normalize_variant(self, variant):
        """Normalize variant format"""
        return str(variant).strip().upper()

def main():
    """Main execution"""
    # Paths
    gold_standard_path = "/mnt/temp2/kronckbm/gvf_output/final_recall_20260206_details.json"
    xml_dir = "/mnt/temp2/kronckbm/gvf_output/KCNH2/tier3_downloads"
    
    # Initialize extractor
    extractor = KCNH2VariantExtractor(gold_standard_path)
    
    print("Starting KCNH2 Tier 3 variant extraction...")
    print("Gold standard loaded: %d variants" % len(extractor.gold_standard_variants))
    
    # Process papers
    results = extractor.process_papers(xml_dir)
    
    # Calculate final statistics
    total_variants = len(results['variants_found'])
    unique_variants = len(results['unique_variants'])
    new_matches = len(results['new_matches'])
    baseline_recall = 14.6
    
    # Estimated recall improvement
    gold_standard_missed = 743  # From final_recall_20260206_details.json
    recall_improvement = (new_matches / gold_standard_missed) * 100 if gold_standard_missed > 0 else 0
    
    # Generate report
    report = "# KCNH2 Tier 3 Variant Extraction Results\n\n"
    report += "## Executive Summary\n\n"
    report += "- **Papers processed**: %d\n" % results['total_papers']
    report += "- **Papers with identifiable variants**: %d\n" % results['papers_with_variants']
    report += "- **Total variants found**: %d\n" % total_variants
    report += "- **Unique variants**: %d\n" % unique_variants
    report += "- **New matches against gold standard**: %d\n" % new_matches
    report += "- **Baseline recall (Tier 2)**: %.1f%%\n" % baseline_recall
    report += "- **Estimated recall improvement**: %.2f%%\n" % recall_improvement
    
    report += "\n## Detailed Analysis\n\n"
    report += "### Variant Types\n"
    
    # Break down by type
    protein_changes = [v for v in results['variants_found'] if v['type'] == 'protein_change']
    frameshifts = [v for v in results['variants_found'] if v['type'] == 'frameshift']
    
    report += "- **Protein changes**: %d variants\n" % len(protein_changes)
    report += "- **Frameshift mutations**: %d variants\n" % len(frameshifts)
    
    report += "\n### Unique New Variants Discovered\n\n"
    unique_list = sorted(results['unique_variants'])
    if unique_list:
        for variant in unique_list:
            count = sum(1 for v in results['variants_found'] if v['variant'] == variant)
            report += "- **%s**: Found in %d papers\n" % (variant, count)
    else:
        report += "- No new variants meeting quality criteria\n"
    
    report += "\n### Paper-by-Paper Results\n\n"
    
    for paper_id in sorted(results['paper_summary'].keys()):
        summary = results['paper_summary'][paper_id]
        report += "#### **%s**\n\n" % paper_id
        report += "- **Variants**: %s\n" % ", ".join(summary['variants'])
        report += "- **Count**: %d\n\n" % summary['variant_count']
    
    if not results['paper_summary']:
        report += "#### No successful variant extraction from processed papers\n\n"
    
    report += "### Technical Implementation Notes\n\n"
    report += "**Pattern Matching Used:**\n"
    report += "- Protein changes: `[A-Z][0-9]{1,4}[A-Z]`\n"
    report += "- Frameshift: `[A-Z][0-9]{1,4}fs`\n"
    report += "- Position validation: max 2000 residues\n"
    report += "- Context extraction: ±15 characters\n\n"
    
    report += "**Quality Filters:**\n"
    report += "- Paper must mention KCNH2/HERG/hERG\n"
    report += "- Duplicate removal within papers\n"
    report += "- Validate reasonable amino acid positions\n"
    report += "- Cross-reference with 861 gold standard variants\n"
    
    # Save markdown report
    output_path = "/mnt/temp2/kronckbm/gvf_output/tier3_extraction_results.md"
    with open(output_path, 'w') as f:
        f.write(report)
    
    # Create JSON for programmatic use
    json_output = {
        'summary': {
            'papers_processed': results['total_papers'],
            'papers_with_variants': results['papers_with_variants'],
            'total_variants': total_variants,
            'unique_variants': unique_variants,
            'new_matches': new_matches,
            'recall_improvement_pct': recall_improvement,
            'protein_changes': len(protein_changes),
            'frameshifts': len(frameshifts)
        },
        'all_variants': results['variants_found'],
        'unique_variants': list(results['unique_variants']),
        'new_matches': list(results['new_matches']),
        'paper_summary': results['paper_summary']
    }
    
    json_path = "/mnt/temp2/kronckbm/gvf_output/tier3_extraction_details.json"
    with open(json_path, 'w') as f:
        json.dump(json_output, f, indent=2, sort_keys=True)
    
    print("=" * 60)
    print("KCNH2 TIER 3 VARIANT EXTRACTION COMPLETE")
    print("=" * 60)
    print("Results saved to:")
    print("  - %s (human-readable)" % output_path)
    print("  - %s (machine-readable)" % json_path)
    print()
    print("Key findings:")
    print("  Papers processed: %d" % results['total_papers'])
    print("  Unique variants: %d" % unique_variants)
    print("  New gold standard matches: %d" % new_matches)
    print("  Estimated recall improvement: %.2f%%" % recall_improvement)

if __name__ == "__main__":
    main()