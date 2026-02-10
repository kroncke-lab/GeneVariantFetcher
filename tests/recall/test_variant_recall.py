"""Variant extraction recall tests for GVF."""
import pytest
import re


class TestVariantRecall:
    """Variant extraction recall tests."""
    
    def test_overall_variant_recall(self, gold_standard, extracted_variants):
        """Test overall variant extraction recall."""
        gold = gold_standard['variants']
        
        matches = gold & extracted_variants
        recall = len(matches) / len(gold) if gold else 0
        
        print(f"\nVariant Recall: {len(matches)}/{len(gold)} = {recall:.1%}")
        print(f"Extracted: {len(extracted_variants)}")
        print(f"Novel (not in gold): {len(extracted_variants - gold)}")
        
        # Floor test - should not regress below 15%
        assert recall >= 0.15, f"Variant recall {recall:.1%} regressed below 15% floor"
    
    def test_variant_type_breakdown(self, gold_standard, extracted_variants):
        """Break down recall by variant type."""
        gold = gold_standard['variants']
        
        categories = {
            'missense': [v for v in gold if re.match(r'^[A-Z]\d+[A-Z]$', v) and v[-1] != 'X'],
            'nonsense': [v for v in gold if re.match(r'^[A-Z]\d+X$', v)],
            'frameshift': [v for v in gold if 'FS' in v],
        }
        
        print("\nRecall by variant type:")
        for cat, variants in categories.items():
            if variants:
                matched = sum(1 for v in variants if v in extracted_variants)
                recall = matched / len(variants)
                print(f"  {cat}: {matched}/{len(variants)} = {recall:.1%}")
    
    def test_high_frequency_variants(self, gold_standard, extracted_variants):
        """Check for most commonly reported variants."""
        df = gold_standard['df']
        
        # Count variant frequency
        var_counts = df['Variant_norm'].value_counts()
        top_10 = var_counts.head(10)
        
        print("\nTop 10 gold standard variants (by paper count):")
        found = 0
        for var, count in top_10.items():
            if pd.notna(var):
                status = "✓" if var in extracted_variants else "✗"
                print(f"  {status} {var}: {count} papers")
                if var in extracted_variants:
                    found += 1
        
        # At least 50% of top-10 should be found
        assert found >= 5, f"Only {found}/10 high-frequency variants found"
    
    def test_position_coverage(self, gold_standard, extracted_variants):
        """Check coverage across protein positions."""
        def get_position(v):
            if pd.isna(v):
                return None
            m = re.search(r'(\d+)', str(v))
            return int(m.group(1)) if m else None
        
        gold = gold_standard['variants']
        gold_positions = {get_position(v) for v in gold if get_position(v)}
        extracted_positions = {get_position(v) for v in extracted_variants if get_position(v)}
        
        coverage = len(gold_positions & extracted_positions) / len(gold_positions) if gold_positions else 0
        
        print(f"\nPosition coverage: {len(gold_positions & extracted_positions)}/{len(gold_positions)} = {coverage:.1%}")
        
        # At least 30% of positions should be covered
        assert coverage >= 0.30, f"Position coverage {coverage:.1%} below 30% floor"


class TestSampleMissedVariants:
    """Diagnostic tests for missed variants."""
    
    def test_show_missed_variants(self, gold_standard, extracted_variants):
        """Display sample of missed variants for debugging."""
        gold = gold_standard['variants']
        missed = gold - extracted_variants
        
        print(f"\nSample missed variants ({len(missed)} total):")
        for v in sorted(missed)[:20]:
            print(f"  {v}")


# Import pandas for value_counts
import pandas as pd
