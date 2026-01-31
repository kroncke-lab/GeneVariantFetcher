#!/usr/bin/env python3
"""Test circuit breaker for low-quality extraction input."""

import pytest
from pipeline.extraction import ExpertExtractor
from utils.models import Paper


class TestCircuitBreaker:
    """Test the extraction circuit breaker functionality."""

    @pytest.fixture
    def extractor(self):
        """Create an ExpertExtractor instance for testing."""
        return ExpertExtractor(models=["gpt-4"])

    def test_short_text_triggers_circuit_breaker(self, extractor):
        """Text under 500 chars should be skipped."""
        paper = Paper(
            pmid="12345678",
            gene_symbol="BRCA1",
            full_text="Short text",  # ~10 chars
        )
        result = extractor.extract(paper)
        assert not result.success
        assert "SKIPPED" in result.error
        assert "too short" in result.error.lower()

    def test_empty_text_triggers_circuit_breaker(self, extractor):
        """Empty text should be skipped."""
        paper = Paper(
            pmid="12345678",
            gene_symbol="BRCA1",
            full_text="",
        )
        result = extractor.extract(paper)
        assert not result.success
        assert "SKIPPED" in result.error

    def test_failed_extraction_markers_trigger_circuit_breaker(self, extractor):
        """Text with multiple failed extraction markers should be skipped."""
        garbage_text = """
        [PDF file available at: http://example.com/paper.pdf]
        [Error converting PDF to text]
        [Legacy .doc file available at: http://example.com/paper.doc]
        Some minimal text content here that won't be enough.
        """ + "x" * 500  # Pad to meet length threshold
        
        paper = Paper(
            pmid="12345678",
            gene_symbol="BRCA1",
            full_text=garbage_text,
        )
        result = extractor.extract(paper)
        assert not result.success
        assert "SKIPPED" in result.error
        assert "failed extraction" in result.error.lower() or "marker" in result.error.lower()

    def test_html_garbage_triggers_circuit_breaker(self, extractor):
        """Text with HTML markup garbage should be skipped."""
        html_garbage = """
        <div class="article-content">
        <span style="font-size: 12px">
        <script type="text/javascript">
        var x = {"key": "value"};
        </script>
        </span>
        </div>
        """ + "x" * 500
        
        paper = Paper(
            pmid="12345678",
            gene_symbol="BRCA1",
            full_text=html_garbage,
        )
        result = extractor.extract(paper)
        assert not result.success
        assert "SKIPPED" in result.error

    def test_low_alphanumeric_ratio_triggers_circuit_breaker(self, extractor):
        """Text with low alphanumeric ratio should be skipped."""
        punctuation_garbage = "!@#$%^&*()_+{}|:<>?[];',./" * 50 + "\n" * 100
        
        paper = Paper(
            pmid="12345678",
            gene_symbol="BRCA1",
            full_text=punctuation_garbage,
        )
        result = extractor.extract(paper)
        assert not result.success
        assert "SKIPPED" in result.error

    def test_no_variant_content_triggers_circuit_breaker(self, extractor):
        """Text without variant content or gene mentions should be skipped."""
        irrelevant_text = """
        This is a long article about cooking recipes. It discusses
        various ingredients and methods for preparing delicious meals.
        The author shares their experience with different cuisines from
        around the world. Topics include Italian pasta, Japanese sushi,
        and Mexican tacos. No scientific or medical content here at all.
        """ * 10  # Make it long enough
        
        paper = Paper(
            pmid="12345678",
            gene_symbol="BRCA1",
            full_text=irrelevant_text,
        )
        result = extractor.extract(paper)
        assert not result.success
        assert "SKIPPED" in result.error
        assert "variant" in result.error.lower() or "gene" in result.error.lower()

    def test_valid_text_passes_circuit_breaker(self, extractor):
        """Valid scientific text with variant content should pass."""
        valid_text = """
        We identified a novel BRCA1 mutation c.5266dupC (p.Gln1756fs)
        in a patient with early-onset breast cancer. The variant was
        confirmed by Sanger sequencing. Additional mutations included
        p.R1699Q and c.181T>G. A total of 15 carriers were identified
        in the study cohort, with penetrance estimated at 85%.
        """ * 10  # Make it long enough
        
        paper = Paper(
            pmid="12345678",
            gene_symbol="BRCA1",
            full_text=valid_text,
        )
        # Note: This will still likely fail because we're not mocking the LLM
        # But it should NOT fail with "SKIPPED" - it should try to extract
        is_usable, reason = extractor._assess_input_quality(valid_text, "BRCA1")
        assert is_usable, f"Valid text should pass quality check: {reason}"

    def test_gene_mention_alone_passes_circuit_breaker(self, extractor):
        """Text mentioning the gene but without variant notation should pass (with warning)."""
        gene_only_text = """
        This study examines the role of BRCA1 in DNA repair pathways.
        We performed functional assays to characterize BRCA1 activity
        in cell lines. The BRCA1 protein was found to interact with
        multiple partners in the homologous recombination pathway.
        """ * 10
        
        paper = Paper(
            pmid="12345678",
            gene_symbol="BRCA1",
            full_text=gene_only_text,
        )
        is_usable, reason = extractor._assess_input_quality(gene_only_text, "BRCA1")
        assert is_usable, f"Gene-mentioning text should pass: {reason}"
        assert "no variant patterns" in reason.lower()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
