#!/usr/bin/env python3
"""
Simple example demonstrating the Clinical Data Triage Filter.

This script shows basic usage of the triage functionality to classify
papers as containing original clinical data or not.
"""

from pipeline.filters import ClinicalDataTriageFilter
import json
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()


def main():
    """Run simple triage examples."""

    print("\n" + "="*80)
    print("Clinical Data Triage Filter - Example Usage")
    print("="*80 + "\n")

    # Initialize the triage filter
    triage = ClinicalDataTriageFilter()

    # Example 1: Case report (should KEEP)
    print("Example 1: Case Report")
    print("-" * 80)

    result1 = triage.triage(
        title="Novel SCN5A Mutation Causing Brugada Syndrome in a Young Patient",
        abstract="""
        We report a 28-year-old male patient who presented with syncope and a family
        history of sudden cardiac death. ECG showed type 1 Brugada pattern. Genetic
        testing revealed a novel heterozygous SCN5A mutation (c.2345G>A, p.Arg782His).
        The patient was successfully treated with an ICD implant.
        """,
        gene="SCN5A",
        pmid="EX001"
    )

    print(json.dumps(result1, indent=2))
    print()

    # Example 2: Review article (should DROP)
    print("Example 2: Review Article")
    print("-" * 80)

    result2 = triage.triage(
        title="SCN5A Variants in Cardiac Arrhythmias: A Comprehensive Review",
        abstract="""
        This review summarizes current knowledge on SCN5A mutations in cardiac
        arrhythmias including Long QT syndrome and Brugada syndrome. We discuss
        molecular mechanisms, genotype-phenotype correlations, and clinical management
        based on analysis of published literature. Over 500 mutations have been
        reported to date.
        """,
        gene="SCN5A",
        pmid="EX002"
    )

    print(json.dumps(result2, indent=2))
    print()

    # Example 3: Animal study (should DROP)
    print("Example 3: Animal Study")
    print("-" * 80)

    result3 = triage.triage(
        title="Scn5a Knockout Mice Exhibit Arrhythmia Phenotypes",
        abstract="""
        We generated Scn5a knockout mice to study sodium channel function. Homozygous
        knockout mice died in utero, while heterozygous mice showed prolonged QRS
        duration and increased ventricular arrhythmia susceptibility. These findings
        provide insights into SCN5A function in mouse cardiac tissue.
        """,
        gene="SCN5A",
        pmid="EX003"
    )

    print(json.dumps(result3, indent=2))
    print()

    # Example 4: Cohort study (should KEEP)
    print("Example 4: Clinical Cohort Study")
    print("-" * 80)

    result4 = triage.triage(
        title="Clinical Outcomes in 150 Patients with SCN5A Mutations",
        abstract="""
        We enrolled 150 patients carrying pathogenic SCN5A mutations from 8 medical
        centers. Mean age was 42±15 years, 60% male. Phenotypes included Brugada
        syndrome (n=85), Long QT syndrome (n=42), and mixed phenotypes (n=23).
        During median follow-up of 7.2 years, 18 patients experienced cardiac events.
        Mutation location predicted outcomes.
        """,
        gene="SCN5A",
        pmid="EX004"
    )

    print(json.dumps(result4, indent=2))
    print()

    # Summary
    print("="*80)
    print("Summary")
    print("="*80)

    results = [result1, result2, result3, result4]
    keep_count = sum(1 for r in results if r['decision'] == 'KEEP')
    drop_count = sum(1 for r in results if r['decision'] == 'DROP')

    print(f"Total papers triaged: {len(results)}")
    print(f"KEEP decisions: {keep_count}")
    print(f"DROP decisions: {drop_count}")
    print()

    # Display which papers were kept
    print("Papers to KEEP (original clinical data):")
    for i, result in enumerate(results, 1):
        if result['decision'] == 'KEEP':
            print(f"  - Example {i} (PMID: {result['pmid']}): {result['reason']}")

    print()
    print("Papers to DROP (no original clinical data):")
    for i, result in enumerate(results, 1):
        if result['decision'] == 'DROP':
            print(f"  - Example {i} (PMID: {result['pmid']}): {result['reason']}")

    print("\n" + "="*80 + "\n")


if __name__ == '__main__':
    # Check for API key
    import os
    if not os.getenv('AI_INTEGRATIONS_OPENAI_API_KEY'):
        print("\n⚠️  Warning: AI_INTEGRATIONS_OPENAI_API_KEY not set!")
        print("   Set it with: export AI_INTEGRATIONS_OPENAI_API_KEY='your-key'\n")
        print("   This example requires API access to run.\n")
        exit(1)

    main()
