"""Carrier count recall tests for GVF."""
import pytest


class TestCarrierRecall:
    """Carrier count recall tests."""
    
    def test_carrier_recall(self, gold_standard, downloaded_pmids):
        """Test carrier count recall based on paper coverage."""
        df = gold_standard['df']
        
        # Count carriers in papers we have downloaded
        carriers_captured = len(df[df['PMID_int'].isin(downloaded_pmids)])
        total = gold_standard['total_carriers']
        
        recall = carriers_captured / total if total > 0 else 0
        
        # Calculate theoretical max (carriers with PMIDs)
        theoretical_max = gold_standard['carriers_with_pmid'] / total
        
        print(f"\nCarrier Recall: {carriers_captured}/{total} = {recall:.1%}")
        print(f"Theoretical max (carriers with PMIDs): {theoretical_max:.1%}")
        print(f"Carriers with no PMID: {total - gold_standard['carriers_with_pmid']}")
        
        # Floor test
        assert recall >= 0.28, f"Carrier recall {recall:.1%} regressed below 28% floor"
    
    def test_carrier_distribution(self, gold_standard):
        """Analyze carrier distribution across papers."""
        df = gold_standard['df']
        
        # Carriers per PMID
        carriers_per_pmid = df.groupby('PMID_int').size().sort_values(ascending=False)
        
        print("\nTop 10 papers by carrier count:")
        for pmid, count in carriers_per_pmid.head(10).items():
            if pmid:
                print(f"  PMID {int(pmid)}: {count} carriers")
    
    def test_carriers_without_pmid(self, gold_standard):
        """Check carriers that have no PMID (unpublished cohorts)."""
        df = gold_standard['df']
        
        no_pmid = df[df['PMID_int'].isna()]
        print(f"\nCarriers with no PMID: {len(no_pmid)} ({len(no_pmid)/len(df):.1%})")
        
        # These represent unpublished cohorts - informational only
        if len(no_pmid) > 0:
            variants_no_pmid = no_pmid['Variant_norm'].dropna().unique()
            print(f"Unique variants in unpublished cohorts: {len(variants_no_pmid)}")
