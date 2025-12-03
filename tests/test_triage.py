#!/usr/bin/env python3
"""
Test cases for the Clinical Data Triage Filter.

These test cases demonstrate various scenarios to ensure the triage filter
correctly identifies papers with original clinical data vs. reviews, animal studies, etc.
"""

import json
from pipeline.filters import ClinicalDataTriageFilter


def test_case_report():
    """Test Case 1: Should KEEP - Case report with patient data."""
    print("\n" + "="*80)
    print("TEST 1: Case Report with Patient Data")
    print("="*80)

    triage = ClinicalDataTriageFilter()

    title = "A novel SCN5A mutation causing Brugada syndrome in a young patient"
    abstract = """
    We report a 28-year-old male patient presenting with syncope and a family history
    of sudden cardiac death. Electrocardiogram showed type 1 Brugada pattern.
    Genetic testing revealed a novel heterozygous mutation in SCN5A (c.2345G>A, p.Arg782His).
    The patient was treated with an implantable cardioverter-defibrillator (ICD).
    Functional studies confirmed the mutation reduced sodium current by 60%.
    This case expands the spectrum of SCN5A mutations associated with Brugada syndrome.
    """

    result = triage.triage(
        title=title,
        abstract=abstract,
        gene="SCN5A",
        pmid="12345678"
    )

    print(json.dumps(result, indent=2))
    print(f"\nExpected: KEEP (case report with patient phenotype)")
    print(f"Actual: {result['decision']}")
    assert result['decision'] == 'KEEP', "Failed: Should KEEP case reports"
    print("✓ PASSED")


def test_review_article():
    """Test Case 2: Should DROP - Review article without new data."""
    print("\n" + "="*80)
    print("TEST 2: Review Article")
    print("="*80)

    triage = ClinicalDataTriageFilter()

    title = "SCN5A variants in cardiac arrhythmias: A comprehensive review"
    abstract = """
    This review summarizes the current knowledge on SCN5A mutations in cardiac
    arrhythmias including Long QT syndrome and Brugada syndrome. We discuss the
    molecular mechanisms, genotype-phenotype correlations, and clinical management
    strategies. Over 500 mutations have been reported in the literature. We analyze
    patterns of inheritance, penetrance, and therapeutic approaches based on
    published studies. This comprehensive review provides insights for clinicians
    and researchers working with SCN5A-related disorders.
    """

    result = triage.triage(
        title=title,
        abstract=abstract,
        gene="SCN5A",
        pmid="12345679"
    )

    print(json.dumps(result, indent=2))
    print(f"\nExpected: DROP (review article)")
    print(f"Actual: {result['decision']}")
    assert result['decision'] == 'DROP', "Failed: Should DROP review articles"
    print("✓ PASSED")


def test_animal_study():
    """Test Case 3: Should DROP - Animal study without patient data."""
    print("\n" + "="*80)
    print("TEST 3: Animal Study (Mouse Model)")
    print("="*80)

    triage = ClinicalDataTriageFilter()

    title = "SCN5A knockout mice exhibit arrhythmia phenotypes"
    abstract = """
    We generated Scn5a knockout mice to study the role of sodium channels in
    cardiac function. Homozygous knockout mice died in utero, while heterozygous
    mice exhibited prolonged QRS duration and increased susceptibility to
    ventricular arrhythmias. Electrophysiological studies showed reduced sodium
    current in isolated cardiomyocytes. These findings provide insights into the
    mechanisms of SCN5A-related arrhythmias in mouse models.
    """

    result = triage.triage(
        title=title,
        abstract=abstract,
        gene="SCN5A",
        pmid="12345680"
    )

    print(json.dumps(result, indent=2))
    print(f"\nExpected: DROP (animal study)")
    print(f"Actual: {result['decision']}")
    assert result['decision'] == 'DROP', "Failed: Should DROP animal studies"
    print("✓ PASSED")


def test_cell_study():
    """Test Case 4: Should DROP - Cell study only."""
    print("\n" + "="*80)
    print("TEST 4: Cell Study Only (No Patients)")
    print("="*80)

    triage = ClinicalDataTriageFilter()

    title = "Functional characterization of SCN5A variants in HEK293 cells"
    abstract = """
    We expressed wild-type and mutant SCN5A channels in HEK293 cells and performed
    patch-clamp recordings. Five variants showed reduced sodium current, three showed
    enhanced persistent current, and two had altered voltage-dependence of activation.
    These in vitro findings provide mechanistic insights into channel dysfunction
    but require validation in clinical contexts.
    """

    result = triage.triage(
        title=title,
        abstract=abstract,
        gene="SCN5A",
        pmid="12345681"
    )

    print(json.dumps(result, indent=2))
    print(f"\nExpected: DROP (cell study only)")
    print(f"Actual: {result['decision']}")
    assert result['decision'] == 'DROP', "Failed: Should DROP cell-only studies"
    print("✓ PASSED")


def test_cohort_study():
    """Test Case 5: Should KEEP - Clinical cohort study."""
    print("\n" + "="*80)
    print("TEST 5: Clinical Cohort Study")
    print("="*80)

    triage = ClinicalDataTriageFilter()

    title = "Clinical outcomes in 150 patients with SCN5A mutations"
    abstract = """
    We enrolled 150 patients carrying pathogenic SCN5A mutations from 8 medical centers.
    Mean age was 42 ± 15 years, 60% were male. Phenotypes included Brugada syndrome (n=85),
    Long QT syndrome (n=42), and mixed phenotypes (n=23). During a median follow-up of
    7.2 years, 18 patients experienced cardiac events. Multivariate analysis identified
    mutation location and baseline QTc as independent predictors of outcomes. This large
    cohort provides prognostic insights for SCN5A mutation carriers.
    """

    result = triage.triage(
        title=title,
        abstract=abstract,
        gene="SCN5A",
        pmid="12345682"
    )

    print(json.dumps(result, indent=2))
    print(f"\nExpected: KEEP (clinical cohort)")
    print(f"Actual: {result['decision']}")
    assert result['decision'] == 'KEEP', "Failed: Should KEEP cohort studies"
    print("✓ PASSED")


def test_functional_with_phenotype():
    """Test Case 6: Should KEEP - Functional study that also describes phenotype."""
    print("\n" + "="*80)
    print("TEST 6: Functional Study + Patient Phenotype")
    print("="*80)

    triage = ClinicalDataTriageFilter()

    title = "Novel SCN5A mutation causes gain-of-function in patient with Long QT syndrome"
    abstract = """
    A 35-year-old female presented with recurrent syncope and QTc of 520 ms. Genetic
    testing identified a novel SCN5A mutation (c.1234G>T, p.Gly412Val). The patient's
    daughter also carried the mutation and had borderline QT prolongation. We expressed
    the mutant channel in CHO cells and performed patch-clamp analysis. The mutation
    caused a significant increase in persistent sodium current (3.2% vs 0.5% for WT),
    explaining the Long QT phenotype. Both patients were treated with beta-blockers.
    This study links genotype to phenotype through functional characterization.
    """

    result = triage.triage(
        title=title,
        abstract=abstract,
        gene="SCN5A",
        pmid="12345683"
    )

    print(json.dumps(result, indent=2))
    print(f"\nExpected: KEEP (functional + phenotype)")
    print(f"Actual: {result['decision']}")
    assert result['decision'] == 'KEEP', "Failed: Should KEEP functional studies with phenotype"
    print("✓ PASSED")


def test_meta_analysis():
    """Test Case 7: Should DROP - Meta-analysis without individual data."""
    print("\n" + "="*80)
    print("TEST 7: Meta-Analysis")
    print("="*80)

    triage = ClinicalDataTriageFilter()

    title = "Meta-analysis of SCN5A mutations in Brugada syndrome: pooled analysis"
    abstract = """
    We performed a systematic review and meta-analysis of studies reporting SCN5A
    mutations in Brugada syndrome. Searching PubMed and Embase yielded 87 eligible
    studies including 3,245 patients. The pooled prevalence of SCN5A mutations was
    23% (95% CI: 19-27%). Mutation-positive patients had higher rates of cardiac
    events (OR=2.1, p<0.001). Subgroup analyses revealed geographic variation in
    mutation frequency. This meta-analysis summarizes the literature on SCN5A mutations.
    """

    result = triage.triage(
        title=title,
        abstract=abstract,
        gene="SCN5A",
        pmid="12345684"
    )

    print(json.dumps(result, indent=2))
    print(f"\nExpected: DROP (meta-analysis)")
    print(f"Actual: {result['decision']}")
    assert result['decision'] == 'DROP', "Failed: Should DROP meta-analyses"
    print("✓ PASSED")


def test_guidelines():
    """Test Case 8: Should DROP - Clinical guidelines without new cases."""
    print("\n" + "="*80)
    print("TEST 8: Clinical Guidelines")
    print("="*80)

    triage = ClinicalDataTriageFilter()

    title = "ACMG guidelines for interpretation of SCN5A sequence variants"
    abstract = """
    The American College of Medical Genetics and Genomics (ACMG) provides recommendations
    for classifying SCN5A sequence variants. We propose a framework for variant
    interpretation incorporating functional evidence, segregation data, population
    frequency, and computational predictions. Specific criteria are provided for
    pathogenic, likely pathogenic, uncertain significance, likely benign, and benign
    classifications. These guidelines aim to standardize variant interpretation in
    clinical laboratories.
    """

    result = triage.triage(
        title=title,
        abstract=abstract,
        gene="SCN5A",
        pmid="12345685"
    )

    print(json.dumps(result, indent=2))
    print(f"\nExpected: DROP (guidelines)")
    print(f"Actual: {result['decision']}")
    assert result['decision'] == 'DROP', "Failed: Should DROP guidelines"
    print("✓ PASSED")


def test_case_series():
    """Test Case 9: Should KEEP - Case series with multiple patients."""
    print("\n" + "="*80)
    print("TEST 9: Case Series")
    print("="*80)

    triage = ClinicalDataTriageFilter()

    title = "Five patients with de novo SCN5A mutations and severe phenotypes"
    abstract = """
    We describe five unrelated patients with de novo SCN5A mutations. Patient 1
    (c.1234G>A) presented with Brugada syndrome at age 19. Patient 2 (c.2345C>T)
    had Long QT syndrome with cardiac arrest at age 22. Patient 3 (c.3456del)
    showed conduction disease and atrial fibrillation. Patient 4 (c.4567T>G) had
    mixed Brugada-Long QT phenotype. Patient 5 (c.5678A>C) presented with dilated
    cardiomyopathy. All mutations were absent in parents and arose de novo.
    These cases highlight the phenotypic diversity of SCN5A mutations.
    """

    result = triage.triage(
        title=title,
        abstract=abstract,
        gene="SCN5A",
        pmid="12345686"
    )

    print(json.dumps(result, indent=2))
    print(f"\nExpected: KEEP (case series)")
    print(f"Actual: {result['decision']}")
    assert result['decision'] == 'KEEP', "Failed: Should KEEP case series"
    print("✓ PASSED")


def run_all_tests():
    """Run all test cases."""
    print("\n" + "#"*80)
    print("# CLINICAL DATA TRIAGE FILTER - TEST SUITE")
    print("#"*80)

    tests = [
        test_case_report,
        test_review_article,
        test_animal_study,
        test_cell_study,
        test_cohort_study,
        test_functional_with_phenotype,
        test_meta_analysis,
        test_guidelines,
        test_case_series,
    ]

    passed = 0
    failed = 0

    for test in tests:
        try:
            test()
            passed += 1
        except AssertionError as e:
            failed += 1
            print(f"✗ FAILED: {e}")
        except Exception as e:
            failed += 1
            print(f"✗ ERROR: {e}")

    print("\n" + "="*80)
    print(f"TEST SUMMARY: {passed} passed, {failed} failed out of {len(tests)} tests")
    print("="*80 + "\n")

    return failed == 0


if __name__ == '__main__':
    import sys

    # Note: These tests require API access and will make real LLM calls
    print("\n⚠️  WARNING: These tests will make real API calls to the LLM service.")
    print("    This will consume API credits (approximately $0.01-0.02 total).")
    response = input("\nProceed with tests? (y/n): ")

    if response.lower() != 'y':
        print("Tests cancelled.")
        sys.exit(0)

    success = run_all_tests()
    sys.exit(0 if success else 1)
