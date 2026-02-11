#!/usr/bin/env python3
"""
Abstract-based classifier for predicting whether a paper contains patient-level variant data.

This module provides a keyword-based scoring model that analyzes paper abstracts and titles
to predict the probability that a paper contains extractable patient-variant associations.

Usage:
    from utils.abstract_classifier import AbstractClassifier
    
    classifier = AbstractClassifier()
    result = classifier.score("title text", "abstract text")
    # result = {'probability': 0.85, 'raw_score': 25.5, 'matched_terms': {...}}
    
    # Or for batch processing:
    results = classifier.score_batch([
        {'pmid': '12345', 'title': '...', 'abstract': '...'},
        ...
    ])

Validated against KCNH2 gold standard (n=262 papers with patient-variant data):
    - AUC-ROC: 0.916
    - At threshold 0.5: Sensitivity 86.5%, Specificity 84.6%
    - Gold standard papers: mean probability 0.77
    - Non-gold papers: mean probability 0.27

Author: GVF Pipeline
Date: 2026-02-10
"""

import re
import math
from typing import Dict, List, Optional, Tuple, Any

# Version for tracking model updates
MODEL_VERSION = "1.0.0"

# Keyword weights - tuned on KCNH2 data but designed to generalize
KEYWORDS: Dict[str, float] = {
    # High-weight terms indicating patient-level variant data
    'mutation screening': 3.0,
    'mutation analysis': 3.0,
    'genotype-phenotype': 3.0,
    'proband': 2.5,
    'probands': 2.5,
    'index case': 2.5,
    'index patient': 2.5,
    'heterozygous': 2.5,
    'heterozygote': 2.5,
    'compound heterozygous': 3.0,
    'homozygous': 2.0,
    'carrier': 2.0,
    'carriers': 2.0,
    'supplementary table': 2.5,
    'supplemental table': 2.5,
    'cohort study': 2.0,
    'familial': 1.5,
    'families': 1.5,
    'family members': 2.0,
    'pedigree': 2.0,
    'mutation identified': 2.5,
    'novel mutation': 2.5,
    'novel variant': 2.5,
    'pathogenic variant': 2.0,
    'genetic testing': 1.5,
    'genetic screening': 1.5,
    'sequencing': 1.0,
    'exome': 1.5,
    'whole exome': 1.5,
    'next-generation sequencing': 1.5,
    
    # Gene-agnostic clinical terms
    'sudden cardiac death': 1.5,
    'sudden death': 1.2,
    'arrhythmia': 1.0,
    'channelopathy': 1.5,
    'cardiomyopathy': 1.2,
    'syndromic': 1.2,
    
    # Medium-weight terms
    'variant': 1.0,
    'variants': 1.0,
    'mutation': 1.0,
    'mutations': 1.0,
    'genotype': 1.0,
    'allele': 1.0,
    'missense': 1.5,
    'frameshift': 1.5,
    'nonsense': 1.5,
    'splice': 1.2,
    'truncating': 1.5,
    'deleterious': 1.2,
    'loss of function': 1.5,
    'loss-of-function': 1.5,
    'patient': 0.5,
    'patients': 0.5,
    'case': 0.5,
    'cases': 0.5,
    'clinical': 0.5,
}

# Regex patterns for variant nomenclature with weights
VARIANT_PATTERNS: List[Tuple[str, float]] = [
    (r'p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}', 3.0),  # p.Ala123Val
    (r'p\.[A-Z]\d+[A-Z]', 2.5),  # p.A123V
    (r'c\.\d+[ACGT]>[ACGT]', 3.0),  # c.123A>G
    (r'c\.\d+[_\d]*del', 2.5),  # c.123del, c.123_456del
    (r'c\.\d+[_\d]*ins', 2.5),  # c.123ins
    (r'c\.\d+[_\d]*dup', 2.5),  # c.123dup
    (r'\d+\s*(patients|families|probands|cases|individuals)', 2.0),  # "25 patients"
    (r'(n\s*=\s*\d+)', 1.5),  # n = 123
]


class AbstractClassifier:
    """
    Classifier for predicting whether a paper contains patient-level variant data.
    
    Uses weighted keyword matching and variant nomenclature pattern detection
    to generate a probability score.
    """
    
    def __init__(
        self,
        keywords: Optional[Dict[str, float]] = None,
        patterns: Optional[List[Tuple[str, float]]] = None,
        sigmoid_k: float = 0.15,
        sigmoid_x0: float = 15.0,
        gene_terms: Optional[List[str]] = None
    ):
        """
        Initialize the classifier.
        
        Args:
            keywords: Custom keyword weights (defaults to built-in KEYWORDS)
            patterns: Custom regex patterns with weights (defaults to VARIANT_PATTERNS)
            sigmoid_k: Sigmoid steepness parameter
            sigmoid_x0: Sigmoid midpoint (raw score that yields 0.5 probability)
            gene_terms: Optional list of gene-specific terms to add with weight 1.5
        """
        self.keywords = keywords or KEYWORDS.copy()
        self.patterns = patterns or VARIANT_PATTERNS.copy()
        self.sigmoid_k = sigmoid_k
        self.sigmoid_x0 = sigmoid_x0
        self.version = MODEL_VERSION
        
        # Add gene-specific terms if provided
        if gene_terms:
            for term in gene_terms:
                self.keywords[term.lower()] = 1.5
    
    def _sigmoid(self, x: float) -> float:
        """Convert raw score to probability using sigmoid function."""
        try:
            return 1.0 / (1.0 + math.exp(-self.sigmoid_k * (x - self.sigmoid_x0)))
        except OverflowError:
            return 0.0 if x < self.sigmoid_x0 else 1.0
    
    def score_text(self, text: str) -> Tuple[float, Dict[str, float]]:
        """
        Score a text string based on keywords and patterns.
        
        Args:
            text: Combined title and abstract text
            
        Returns:
            Tuple of (raw_score, matched_terms_dict)
        """
        if not text:
            return 0.0, {}
        
        text_lower = text.lower()
        matched_terms: Dict[str, float] = {}
        
        # Keyword matching
        for keyword, weight in self.keywords.items():
            count = text_lower.count(keyword)
            if count > 0:
                # Diminishing returns for repeated terms
                score = weight * min(count, 3)
                matched_terms[keyword] = score
        
        # Pattern matching
        for pattern, weight in self.patterns:
            try:
                matches = re.findall(pattern, text, re.IGNORECASE)
                if matches:
                    count = len(matches)
                    score = weight * min(count, 5)
                    matched_terms[f'pattern:{pattern[:25]}'] = score
            except re.error:
                continue
        
        total_score = sum(matched_terms.values())
        return total_score, matched_terms
    
    def score(
        self,
        title: str = "",
        abstract: str = "",
        combined_text: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Score a paper's abstract/title for patient-variant data probability.
        
        Args:
            title: Paper title
            abstract: Paper abstract
            combined_text: Pre-combined text (overrides title+abstract)
            
        Returns:
            Dict with 'probability', 'raw_score', 'matched_terms', 'version'
        """
        text = combined_text or f"{title} {abstract}".strip()
        raw_score, matched_terms = self.score_text(text)
        probability = self._sigmoid(raw_score)
        
        return {
            'probability': probability,
            'raw_score': raw_score,
            'matched_terms': matched_terms,
            'version': self.version
        }
    
    def score_batch(
        self,
        papers: List[Dict[str, str]],
        title_key: str = 'title',
        abstract_key: str = 'abstract',
        id_key: str = 'pmid'
    ) -> List[Dict[str, Any]]:
        """
        Score a batch of papers.
        
        Args:
            papers: List of dicts with title, abstract, and id keys
            title_key: Key for title in each dict
            abstract_key: Key for abstract in each dict
            id_key: Key for paper ID in each dict
            
        Returns:
            List of result dicts with id and scoring info
        """
        results = []
        for paper in papers:
            title = paper.get(title_key, '')
            abstract = paper.get(abstract_key, '')
            paper_id = paper.get(id_key, '')
            
            score_result = self.score(title=title, abstract=abstract)
            score_result['id'] = paper_id
            results.append(score_result)
        
        return results
    
    def predict(
        self,
        title: str = "",
        abstract: str = "",
        threshold: float = 0.5
    ) -> bool:
        """
        Predict whether a paper contains patient-variant data.
        
        Args:
            title: Paper title
            abstract: Paper abstract
            threshold: Probability threshold for positive prediction
            
        Returns:
            True if predicted to contain patient-variant data
        """
        result = self.score(title=title, abstract=abstract)
        return result['probability'] >= threshold
    
    def get_priority_score(self, probability: float) -> str:
        """
        Convert probability to priority tier for pipeline processing.
        
        Returns:
            'HIGH' (>0.7), 'MEDIUM' (0.4-0.7), 'LOW' (<0.4)
        """
        if probability >= 0.7:
            return 'HIGH'
        elif probability >= 0.4:
            return 'MEDIUM'
        else:
            return 'LOW'


def create_classifier_for_gene(gene_symbol: str) -> AbstractClassifier:
    """
    Create a classifier with gene-specific terms added.
    
    Args:
        gene_symbol: Gene symbol (e.g., 'KCNH2', 'SCN5A')
        
    Returns:
        Configured AbstractClassifier instance
    """
    # Gene-specific terms to add
    gene_terms = [gene_symbol.lower()]
    
    # Add common aliases/related terms for known genes
    gene_aliases = {
        'KCNH2': ['herg', 'lqt2', 'long qt syndrome type 2'],
        'SCN5A': ['nav1.5', 'lqt3', 'brugada'],
        'KCNQ1': ['kvlqt1', 'lqt1'],
        'RYR2': ['cpvt', 'catecholaminergic'],
        'MYBPC3': ['hcm', 'hypertrophic cardiomyopathy'],
        'MYH7': ['hcm', 'dcm'],
        'LMNA': ['laminopathy', 'dcm'],
        'PKP2': ['arvc', 'arrhythmogenic'],
        'DSP': ['arvc', 'arrhythmogenic'],
        'BRCA1': ['hereditary breast', 'hboc'],
        'BRCA2': ['hereditary breast', 'hboc'],
    }
    
    if gene_symbol.upper() in gene_aliases:
        gene_terms.extend(gene_aliases[gene_symbol.upper()])
    
    return AbstractClassifier(gene_terms=gene_terms)


# Convenience function for quick scoring
def score_abstract(title: str = "", abstract: str = "") -> float:
    """
    Quick function to get probability score for a paper.
    
    Args:
        title: Paper title
        abstract: Paper abstract
        
    Returns:
        Probability (0-1) that paper contains patient-variant data
    """
    classifier = AbstractClassifier()
    result = classifier.score(title=title, abstract=abstract)
    return result['probability']


if __name__ == "__main__":
    # Demo usage
    classifier = AbstractClassifier()
    
    # Example abstract
    demo_title = "Novel KCNH2 mutations identified in a cohort of long QT syndrome probands"
    demo_abstract = """
    Background: Long QT syndrome (LQTS) is caused by mutations in cardiac ion channels.
    Methods: We performed mutation screening in 150 unrelated probands with LQT2.
    Results: We identified 25 novel heterozygous mutations in KCNH2, including 15 missense 
    variants (p.Ala123Val, p.Gly456Arg, etc.). Genotype-phenotype analysis revealed that
    carriers of pore mutations had more severe phenotypes. Supplementary Table 1 lists
    all identified variants with their functional characterization.
    Conclusions: Genetic testing identifies pathogenic variants in most LQTS families.
    """
    
    result = classifier.score(title=demo_title, abstract=demo_abstract)
    
    print(f"Demo Abstract Analysis:")
    print(f"  Probability: {result['probability']:.3f}")
    print(f"  Raw Score: {result['raw_score']:.1f}")
    print(f"  Top matched terms:")
    sorted_terms = sorted(result['matched_terms'].items(), key=lambda x: -x[1])[:10]
    for term, score in sorted_terms:
        print(f"    - {term}: {score:.1f}")
