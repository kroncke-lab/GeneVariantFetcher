"""LQT-affected individual recall tests for GVF."""
import pytest


class TestLQTRecall:
    """LQT-affected individual recall tests."""
    
    def test_lqt_affected_recall(self, gold_standard, downloaded_pmids):
        """Test LQT-affected individual recall."""
        df = gold_standard['df']
        
        # LQT-affected in papers we have
        captured = len(df[
            (df['PMID_int'].isin(downloaded_pmids)) & 
            (df['is_LQT_affected'])
        ])
        
        total = gold_standard['total_lqt_affected']
        recall = captured / total if total > 0 else 0
        
        print(f"\nLQT-Affected Recall: {captured}/{total} = {recall:.1%}")
        
        # Floor test
        assert recall >= 0.27, f"LQT recall {recall:.1%} regressed below 27% floor"
    
    def test_lqt_phenotype_distribution(self, gold_standard):
        """Analyze LQT phenotype distribution."""
        df = gold_standard['df']
        
        if 'LQT' not in df.columns:
            pytest.skip("LQT column not found")
        
        lqt_counts = df['LQT'].value_counts(dropna=False)
        
        print("\nLQT phenotype distribution:")
        for val, count in lqt_counts.head(10).items():
            print(f"  LQT={val}: {count} carriers")
    
    def test_lqt_by_paper_source(self, gold_standard, downloaded_pmids):
        """Break down LQT capture by paper source."""
        df = gold_standard['df']
        
        # In downloaded papers
        in_downloaded = df[df['PMID_int'].isin(downloaded_pmids)]['is_LQT_affected'].sum()
        
        # In missing papers  
        missing_pmids = gold_standard['pmids'] - downloaded_pmids
        in_missing = df[df['PMID_int'].isin(missing_pmids)]['is_LQT_affected'].sum()
        
        # No PMID
        no_pmid = df[df['PMID_int'].isna()]['is_LQT_affected'].sum()
        
        print(f"\nLQT-affected breakdown:")
        print(f"  In downloaded papers: {int(in_downloaded)}")
        print(f"  In missing papers: {int(in_missing)}")
        print(f"  No PMID (unpublished): {int(no_pmid)}")
