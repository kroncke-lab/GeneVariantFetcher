"""
Gene Configuration for Variant Validation

Contains gene-specific metadata for validating extracted variants:
- Protein length (amino acids) for p. notation validation
- Transcript length (nucleotides) for c. notation validation
- Common aliases for gene symbol matching

Data sources:
- Protein lengths: UniProt canonical sequences
- Transcript lengths: RefSeq canonical transcripts (coding sequence only)

Usage:
    from config.gene_config import GENE_CONFIG, validate_variant_position
    
    if validate_variant_position("p.G628S", "KCNH2"):
        # Valid position
    
    if validate_variant_position("c.1882G>A", "KCNH2"):
        # Valid cDNA position
"""

import re
import logging
from typing import Optional, Dict, Any

logger = logging.getLogger(__name__)

# Gene configuration: protein_length (aa), transcript_length (nt coding), aliases
GENE_CONFIG: Dict[str, Dict[str, Any]] = {
    # Long QT / Cardiac channelopathies
    "KCNH2": {
        "protein_length": 1159,      # UniProt Q12809
        "transcript_length": 3480,   # NM_000238.4 (1159 * 3 = 3477, +3 for stop)
        "aliases": ["HERG", "LQT2", "ERG", "ERG1", "HERG1", "Kv11.1", "SQT1"],
        "chromosome": "7",
        "description": "Potassium voltage-gated channel subfamily H member 2",
    },
    "KCNQ1": {
        "protein_length": 676,       # UniProt P51787
        "transcript_length": 2031,   # NM_000218.3
        "aliases": ["LQT1", "KVLQT1", "Kv7.1", "KCNA9"],
        "chromosome": "11",
        "description": "Potassium voltage-gated channel subfamily Q member 1",
    },
    "SCN5A": {
        "protein_length": 2016,      # UniProt Q14524
        "transcript_length": 6051,   # NM_198056.3
        "aliases": ["LQT3", "Nav1.5", "HB1", "HB2", "HBBD", "HH1", "ICCD", "IVF", "SSS1", "VF1"],
        "chromosome": "3",
        "description": "Sodium voltage-gated channel alpha subunit 5",
    },
    "KCNE1": {
        "protein_length": 129,       # UniProt P15382
        "transcript_length": 390,    # NM_000219.6
        "aliases": ["LQT5", "JLNS2", "MinK", "ISK"],
        "chromosome": "21",
        "description": "Potassium voltage-gated channel subfamily E member 1",
    },
    "KCNE2": {
        "protein_length": 123,       # UniProt Q9Y6J6
        "transcript_length": 372,    # NM_172201.2
        "aliases": ["LQT6", "MIRP1", "MiRP1"],
        "chromosome": "21",
        "description": "Potassium voltage-gated channel subfamily E member 2",
    },
    "KCNJ2": {
        "protein_length": 427,       # UniProt P63252
        "transcript_length": 1284,   # NM_000891.3
        "aliases": ["LQT7", "Kir2.1", "HHBIRK1", "IRK1"],
        "chromosome": "17",
        "description": "Potassium inwardly rectifying channel subfamily J member 2",
    },
    "CACNA1C": {
        "protein_length": 2221,      # UniProt Q13936
        "transcript_length": 6666,   # NM_000719.7
        "aliases": ["LQT8", "Cav1.2", "TS", "CACN2", "CCHL1A1"],
        "chromosome": "12",
        "description": "Calcium voltage-gated channel subunit alpha1 C",
    },
    # Brugada / other arrhythmias
    "SCN1B": {
        "protein_length": 218,       # UniProt Q07699
        "transcript_length": 657,    # NM_001037.5
        "aliases": ["BRGDA5"],
        "chromosome": "19",
        "description": "Sodium voltage-gated channel beta subunit 1",
    },
    "SCN2B": {
        "protein_length": 215,       # UniProt O60939
        "transcript_length": 648,    # NM_004588.5
        "aliases": [],
        "chromosome": "11",
        "description": "Sodium voltage-gated channel beta subunit 2",
    },
    "SCN3B": {
        "protein_length": 215,       # UniProt Q9NY72
        "transcript_length": 648,    # NM_018400.4
        "aliases": ["BRGDA7"],
        "chromosome": "11",
        "description": "Sodium voltage-gated channel beta subunit 3",
    },
    # RyR2 - very large protein
    "RYR2": {
        "protein_length": 4967,      # UniProt Q92736
        "transcript_length": 14904,  # NM_001035.3
        "aliases": ["CPVT1", "ARVC2", "ARVD2"],
        "chromosome": "1",
        "description": "Ryanodine receptor 2",
    },
    # CALM genes
    "CALM1": {
        "protein_length": 149,       # UniProt P0DP23
        "transcript_length": 450,    # NM_006888.6
        "aliases": ["CPVT4", "LQT14"],
        "chromosome": "14",
        "description": "Calmodulin 1",
    },
    "CALM2": {
        "protein_length": 149,       # UniProt P0DP24
        "transcript_length": 450,    # NM_001743.6
        "aliases": ["LQT15"],
        "chromosome": "2",
        "description": "Calmodulin 2",
    },
    "CALM3": {
        "protein_length": 149,       # UniProt P0DP25
        "transcript_length": 450,    # NM_005184.4
        "aliases": ["LQT16"],
        "chromosome": "19",
        "description": "Calmodulin 3",
    },
}


def get_gene_config(gene_symbol: str) -> Optional[Dict[str, Any]]:
    """
    Get configuration for a gene, checking aliases.
    
    Args:
        gene_symbol: Gene symbol or alias (case-insensitive)
        
    Returns:
        Gene config dict or None if not found
    """
    gene_upper = gene_symbol.upper()
    
    # Direct match
    if gene_upper in GENE_CONFIG:
        return GENE_CONFIG[gene_upper]
    
    # Check aliases
    for gene, config in GENE_CONFIG.items():
        if gene_upper in [a.upper() for a in config.get("aliases", [])]:
            return config
    
    return None


def validate_variant_position(variant_notation: str, gene_symbol: str) -> bool:
    """
    Validate that a variant position is within the gene's coding region.
    
    Handles both protein (p.) and cDNA (c.) notations.
    
    Args:
        variant_notation: Variant string (e.g., "p.G628S", "c.1882G>A")
        gene_symbol: Gene symbol to validate against
        
    Returns:
        True if position is valid, False if definitely invalid,
        True if gene not in config (permissive fallback)
    """
    config = get_gene_config(gene_symbol)
    if not config:
        # Gene not in our config - be permissive
        logger.debug(f"Gene {gene_symbol} not in config, skipping position validation")
        return True
    
    # Extract position from protein notation (p.G628S, p.Gly628Ser, p.V2016M, etc.)
    # Handles both 1-letter (p.G628S) and 3-letter (p.Gly628Ser) amino acid codes
    protein_match = re.search(r"p\.(?:[A-Z][a-z]{2}|[A-Z])(\d+)", variant_notation, re.IGNORECASE)
    if protein_match:
        position = int(protein_match.group(1))
        max_pos = config["protein_length"]
        if position > max_pos:
            logger.warning(
                f"Invalid protein position: {variant_notation} "
                f"(pos {position} > {max_pos} for {gene_symbol})"
            )
            return False
        return True
    
    # Extract position from cDNA notation (c.1882G>A, c.1234+5G>A, etc.)
    # Handle simple positions and intronic positions
    cdna_match = re.search(r"c\.(\d+)", variant_notation, re.IGNORECASE)
    if cdna_match:
        position = int(cdna_match.group(1))
        max_pos = config["transcript_length"]
        if position > max_pos:
            logger.warning(
                f"Invalid cDNA position: {variant_notation} "
                f"(pos {position} > {max_pos} for {gene_symbol})"
            )
            return False
        return True
    
    # Could not parse position - be permissive
    return True


def gene_symbol_in_text(gene_symbol: str, text: str) -> bool:
    """
    Check if gene symbol or any of its aliases appear in text.
    
    Args:
        gene_symbol: Primary gene symbol
        text: Text to search (case-insensitive)
        
    Returns:
        True if gene or alias found in text
    """
    text_upper = text.upper()
    
    # Check primary symbol
    if gene_symbol.upper() in text_upper:
        return True
    
    # Check aliases
    config = get_gene_config(gene_symbol)
    if config:
        for alias in config.get("aliases", []):
            if alias.upper() in text_upper:
                return True
    
    return False


# For backwards compatibility / quick access
def get_protein_length(gene_symbol: str) -> Optional[int]:
    """Get protein length for a gene."""
    config = get_gene_config(gene_symbol)
    return config["protein_length"] if config else None


def get_transcript_length(gene_symbol: str) -> Optional[int]:
    """Get transcript (CDS) length for a gene."""
    config = get_gene_config(gene_symbol)
    return config["transcript_length"] if config else None
