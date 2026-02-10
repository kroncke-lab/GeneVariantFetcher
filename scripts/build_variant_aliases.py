#!/usr/bin/env python3
"""
Build comprehensive KCNH2 variant alias dictionary.

Generates all possible naming forms for each gold standard variant:
- Missense: R176W, p.R176W, p.Arg176Trp, Arg176Trp
- Frameshift: L987fsX, L987fs, p.Leu987fs, p.Leu987fsTer, p.Leu987Profs*, L987fsX10
- Nonsense: W1001X, W1001*, p.Trp1001Ter, p.Trp1001*, Trp1001X
- Splice: IVS9+1G>A, c.2398+1G>A
- Deletions: L552del, p.Leu552del, del552
"""

import pandas as pd
import json
import re
from pathlib import Path

# Amino acid mappings
AA_3_TO_1 = {
    'Ala': 'A', 'Cys': 'C', 'Asp': 'D', 'Glu': 'E', 'Phe': 'F',
    'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Lys': 'K', 'Leu': 'L',
    'Met': 'M', 'Asn': 'N', 'Pro': 'P', 'Gln': 'Q', 'Arg': 'R',
    'Ser': 'S', 'Thr': 'T', 'Val': 'V', 'Trp': 'W', 'Tyr': 'Y',
    'Ter': '*', 'Stop': '*', 'Xaa': 'X'
}

AA_1_TO_3 = {v: k for k, v in AA_3_TO_1.items()}
AA_1_TO_3['X'] = 'Ter'  # X represents stop codon
AA_1_TO_3['*'] = 'Ter'

# IVS to cDNA position mapping for KCNH2 (intron end positions)
# Based on KCNH2 transcript NM_000238.4
KCNH2_IVS_MAP = {
    'IVS1': 'c.175',   # Intron 1 starts after position 175
    'IVS2': 'c.453',   # Intron 2
    'IVS3': 'c.575',   # Intron 3
    'IVS4': 'c.787',   # Intron 4
    'IVS5': 'c.1129',  # Intron 5
    'IVS6': 'c.1351',  # Intron 6
    'IVS7': 'c.1592',  # Intron 7
    'IVS8': 'c.1759',  # Intron 8
    'IVS9': 'c.2006',  # Intron 9
    'IVS10': 'c.2398', # Intron 10
    'IVS11': 'c.2599', # Intron 11
    'IVS12': 'c.2769', # Intron 12
    'IVS13': 'c.3040', # Intron 13
    'IVS14': 'c.3152', # Intron 14
}


def parse_variant(variant_str):
    """Parse a variant string and return its components."""
    if not variant_str or not isinstance(variant_str, str):
        return None
    
    v = variant_str.strip()
    
    # Skip chromosomal deletions, complex rearrangements
    if any(x in v.lower() for x in ['del(7)', 'q34', 'q36', 'exon', 'ex']):
        return {'type': 'chromosomal', 'original': v}
    
    # cDNA notation: c.1234A>G, c.1234+1G>A
    if v.startswith('c.'):
        return {'type': 'cdna', 'original': v, 'cdna': v}
    
    # IVS notation: IVS9+1G>A, IVS9-28A/G
    ivs_match = re.match(r'^IVS(\d+)([\+\-]\d+)([ACGT])[>/]([ACGT])$', v, re.IGNORECASE)
    if ivs_match:
        ivs_num, offset, ref, alt = ivs_match.groups()
        return {
            'type': 'splice',
            'original': v,
            'ivs_num': ivs_num,
            'offset': offset,
            'ref': ref.upper(),
            'alt': alt.upper()
        }
    
    # Missense: A561V, p.Ala561Val, Ala561Val
    # Single-letter with p. prefix
    m = re.match(r'^(?:p\.)?([A-Z])(\d+)([A-Z])$', v, re.IGNORECASE)
    if m:
        return {
            'type': 'missense',
            'original': v,
            'ref': m.group(1).upper(),
            'pos': int(m.group(2)),
            'alt': m.group(3).upper()
        }
    
    # Three-letter: Ala561Val, p.Ala561Val
    m = re.match(r'^(?:p\.)?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$', v, re.IGNORECASE)
    if m:
        ref_3 = m.group(1).capitalize()
        alt_3 = m.group(3).capitalize()
        return {
            'type': 'missense',
            'original': v,
            'ref': AA_3_TO_1.get(ref_3, ref_3[0].upper()),
            'pos': int(m.group(2)),
            'alt': AA_3_TO_1.get(alt_3, alt_3[0].upper()),
            'ref_3': ref_3,
            'alt_3': alt_3
        }
    
    # Nonsense/Stop: R864X, R864*, p.Arg864Ter, Arg864*
    m = re.match(r'^(?:p\.)?([A-Z])(\d+)(X|\*|Ter|stop)$', v, re.IGNORECASE)
    if m:
        return {
            'type': 'nonsense',
            'original': v,
            'ref': m.group(1).upper(),
            'pos': int(m.group(2))
        }
    
    # Three-letter nonsense: Arg864Ter, p.Arg864*
    m = re.match(r'^(?:p\.)?([A-Z][a-z]{2})(\d+)(Ter|\*|stop|X)$', v, re.IGNORECASE)
    if m:
        ref_3 = m.group(1).capitalize()
        return {
            'type': 'nonsense',
            'original': v,
            'ref': AA_3_TO_1.get(ref_3, ref_3[0].upper()),
            'pos': int(m.group(2)),
            'ref_3': ref_3
        }
    
    # Frameshift: A193fsX, A193fs, A193fsX10, A193fs*10
    m = re.match(r'^(?:p\.)?([A-Z])(\d+)(?:[A-Za-z]*)?(fs)(?:X|\*|Ter)?(\d*)$', v, re.IGNORECASE)
    if m:
        return {
            'type': 'frameshift',
            'original': v,
            'ref': m.group(1).upper(),
            'pos': int(m.group(2)),
            'ter_pos': m.group(4) if m.group(4) else None
        }
    
    # Three-letter frameshift: Ala193fs, Ala193Profs*10
    m = re.match(r'^(?:p\.)?([A-Z][a-z]{2})(\d+)(?:[A-Z][a-z]{2})?(fs)(?:X|\*|Ter)?(\d*)$', v, re.IGNORECASE)
    if m:
        ref_3 = m.group(1).capitalize()
        return {
            'type': 'frameshift',
            'original': v,
            'ref': AA_3_TO_1.get(ref_3, ref_3[0].upper()),
            'pos': int(m.group(2)),
            'ref_3': ref_3,
            'ter_pos': m.group(4) if m.group(4) else None
        }
    
    # Deletion: I30del, I30Del, p.Ile30del
    m = re.match(r'^(?:p\.)?([A-Z])(\d+)(del|dup|ins)$', v, re.IGNORECASE)
    if m:
        return {
            'type': 'indel',
            'subtype': m.group(3).lower(),
            'original': v,
            'ref': m.group(1).upper(),
            'pos': int(m.group(2))
        }
    
    # Three-letter indel: Ile30del
    m = re.match(r'^(?:p\.)?([A-Z][a-z]{2})(\d+)(del|dup|ins)$', v, re.IGNORECASE)
    if m:
        ref_3 = m.group(1).capitalize()
        return {
            'type': 'indel',
            'subtype': m.group(3).lower(),
            'original': v,
            'ref': AA_3_TO_1.get(ref_3, ref_3[0].upper()),
            'pos': int(m.group(2)),
            'ref_3': ref_3
        }
    
    return {'type': 'unknown', 'original': v}


def generate_all_forms(parsed):
    """Generate all possible naming forms for a parsed variant."""
    forms = set()
    
    if not parsed or parsed.get('type') in ['unknown', 'chromosomal']:
        if parsed:
            forms.add(parsed.get('original', '').upper())
        return list(forms)
    
    if parsed['type'] == 'missense':
        ref = parsed['ref']
        pos = parsed['pos']
        alt = parsed['alt']
        ref_3 = AA_1_TO_3.get(ref, 'Xaa')
        alt_3 = AA_1_TO_3.get(alt, 'Xaa')
        
        # All possible forms
        forms.update([
            f"{ref}{pos}{alt}",                      # A561V
            f"p.{ref}{pos}{alt}",                    # p.A561V
            f"{ref_3}{pos}{alt_3}",                  # Ala561Val
            f"p.{ref_3}{pos}{alt_3}",                # p.Ala561Val
            f"{ref.lower()}{pos}{alt.lower()}",     # a561v (lowercase)
        ])
    
    elif parsed['type'] == 'nonsense':
        ref = parsed['ref']
        pos = parsed['pos']
        ref_3 = AA_1_TO_3.get(ref, 'Xaa')
        
        forms.update([
            f"{ref}{pos}X",                          # R864X
            f"{ref}{pos}*",                          # R864*
            f"p.{ref}{pos}X",                        # p.R864X
            f"p.{ref}{pos}*",                        # p.R864*
            f"p.{ref}{pos}Ter",                      # p.R864Ter
            f"{ref_3}{pos}Ter",                      # Arg864Ter
            f"{ref_3}{pos}X",                        # Arg864X
            f"{ref_3}{pos}*",                        # Arg864*
            f"p.{ref_3}{pos}Ter",                    # p.Arg864Ter
            f"p.{ref_3}{pos}*",                      # p.Arg864*
            f"{ref}{pos}STOP",                       # R864STOP
            f"{ref.lower()}{pos}x",                  # r864x
        ])
    
    elif parsed['type'] == 'frameshift':
        ref = parsed['ref']
        pos = parsed['pos']
        ref_3 = AA_1_TO_3.get(ref, 'Xaa')
        ter = parsed.get('ter_pos')
        
        # Base forms
        base_forms = [
            f"{ref}{pos}fsX",                        # A193fsX (canonical)
            f"{ref}{pos}fs",                         # A193fs
            f"{ref}{pos}fs*",                        # A193fs*
            f"p.{ref}{pos}fs",                       # p.A193fs
            f"p.{ref}{pos}fsX",                      # p.A193fsX
            f"p.{ref}{pos}fs*",                      # p.A193fs*
            f"{ref_3}{pos}fs",                       # Ala193fs
            f"{ref_3}{pos}fsX",                      # Ala193fsX
            f"{ref_3}{pos}fs*",                      # Ala193fs*
            f"p.{ref_3}{pos}fs",                     # p.Ala193fs
            f"p.{ref_3}{pos}fs*",                    # p.Ala193fs*
            f"p.{ref_3}{pos}fsTer",                  # p.Ala193fsTer
            f"{ref}{pos}FS",                         # A193FS (uppercase)
            f"{ref.lower()}{pos}fsx",                # a193fsx
        ]
        forms.update(base_forms)
        
        # With termination position
        if ter:
            forms.update([
                f"{ref}{pos}fsX{ter}",               # A193fsX10
                f"{ref}{pos}fs*{ter}",               # A193fs*10
                f"p.{ref}{pos}fsX{ter}",             # p.A193fsX10
                f"p.{ref}{pos}fs*{ter}",             # p.A193fs*10
                f"{ref_3}{pos}fs*{ter}",             # Ala193fs*10
                f"p.{ref_3}{pos}fsTer{ter}",         # p.Ala193fsTer10
            ])
    
    elif parsed['type'] == 'indel':
        ref = parsed['ref']
        pos = parsed['pos']
        subtype = parsed['subtype']
        ref_3 = AA_1_TO_3.get(ref, 'Xaa')
        
        forms.update([
            f"{ref}{pos}{subtype}",                  # I30del
            f"p.{ref}{pos}{subtype}",                # p.I30del
            f"{ref_3}{pos}{subtype}",                # Ile30del
            f"p.{ref_3}{pos}{subtype}",              # p.Ile30del
            f"{ref}{pos}{subtype.upper()}",          # I30DEL
            f"{subtype}{pos}",                       # del30 (less common)
        ])
    
    elif parsed['type'] == 'splice':
        ivs_num = parsed['ivs_num']
        offset = parsed['offset']
        ref = parsed['ref']
        alt = parsed['alt']
        
        # IVS forms
        forms.update([
            f"IVS{ivs_num}{offset}{ref}>{alt}",      # IVS9+1G>A
            f"IVS{ivs_num}{offset}{ref}/{alt}",      # IVS9+1G/A
            f"IVS{ivs_num}{offset}{ref}->{alt}",     # IVS9+1G->A
        ])
        
        # cDNA forms if we have mapping
        ivs_key = f"IVS{ivs_num}"
        if ivs_key in KCNH2_IVS_MAP:
            base = KCNH2_IVS_MAP[ivs_key]
            forms.update([
                f"{base}{offset}{ref}>{alt}",        # c.2006+1G>A
                f"{base}{offset}{ref}>{alt}".lower(),# c.2006+1g>a
            ])
    
    elif parsed['type'] == 'cdna':
        cdna = parsed.get('cdna', parsed['original'])
        forms.add(cdna)
        forms.add(cdna.upper())
        forms.add(cdna.lower())
        # Add without c. prefix
        if cdna.startswith('c.'):
            forms.add(cdna[2:])
    
    return list(forms)


def main():
    # Read the gold standard Excel
    xls_path = "comparison_results/KCNH2 HeterozygoteDatabase-4-clinical_w_paris_japan_mayo_and_italycohort.xls"
    df = pd.read_excel(xls_path)
    
    print(f"=== BUILDING KCNH2 VARIANT ALIAS DICTIONARY ===")
    print(f"Loaded {len(df)} rows from gold standard")
    
    # Get EXCLUDE column (filter to included only)
    if 'EXCLUDE?' in df.columns:
        include_mask = df['EXCLUDE?'] == 0
        df_include = df[include_mask]
        print(f"Included variants (EXCLUDE?=0): {len(df_include)}")
    else:
        df_include = df
    
    # Extract unique variants from Variant column
    variants = df['Variant'].dropna().unique()
    print(f"Unique variants in gold standard: {len(variants)}")
    
    # Build alias dictionary: alias -> canonical form
    # Canonical = single-letter, no p. prefix (e.g., A561V, A193fsX, R864X)
    alias_dict = {}
    variant_info = {}
    parse_failures = []
    
    for v in variants:
        if not isinstance(v, str):
            continue
        v = v.strip()
        if not v:
            continue
        
        parsed = parse_variant(v)
        
        if not parsed or parsed.get('type') == 'unknown':
            parse_failures.append(v)
            # Still add the original
            alias_dict[v.upper()] = v.upper()
            continue
        
        # Determine canonical form (single-letter, no prefix)
        if parsed['type'] == 'missense':
            canonical = f"{parsed['ref']}{parsed['pos']}{parsed['alt']}"
        elif parsed['type'] == 'nonsense':
            canonical = f"{parsed['ref']}{parsed['pos']}X"
        elif parsed['type'] == 'frameshift':
            canonical = f"{parsed['ref']}{parsed['pos']}fsX"
        elif parsed['type'] == 'indel':
            canonical = f"{parsed['ref']}{parsed['pos']}{parsed['subtype']}"
        elif parsed['type'] == 'splice':
            # Use IVS notation as canonical for splice
            canonical = f"IVS{parsed['ivs_num']}{parsed['offset']}{parsed['ref']}>{parsed['alt']}"
        elif parsed['type'] == 'cdna':
            canonical = parsed.get('cdna', v)
        elif parsed['type'] == 'chromosomal':
            canonical = v.upper()
        else:
            canonical = v.upper()
        
        # Generate all forms and map to canonical
        all_forms = generate_all_forms(parsed)
        for form in all_forms:
            alias_dict[form.upper()] = canonical
        
        # Also add original form
        alias_dict[v.upper()] = canonical
        
        # Store info
        variant_info[canonical] = {
            'original': v,
            'parsed': parsed,
            'all_forms': all_forms
        }
    
    print(f"\nParsed {len(variant_info)} unique canonical variants")
    print(f"Parse failures: {len(parse_failures)}")
    if parse_failures:
        print("  Examples:", parse_failures[:10])
    
    print(f"Total aliases: {len(alias_dict)}")
    
    # Save alias dictionary
    output_path = Path("utils/kcnh2_variant_aliases.json")
    with open(output_path, "w") as f:
        json.dump({
            "metadata": {
                "description": "KCNH2 variant alias dictionary: maps any form to canonical single-letter format",
                "canonical_format": "Single-letter, no prefix (A561V, R864X, A193fsX)",
                "total_variants": len(variant_info),
                "total_aliases": len(alias_dict),
                "source": "Gold standard Excel + generated forms"
            },
            "aliases": alias_dict,
            "variant_info": {k: {"original": v["original"], "all_forms": v["all_forms"]} 
                           for k, v in variant_info.items()}
        }, f, indent=2)
    
    print(f"\n=== SAVED TO {output_path} ===")
    
    # Print some examples
    print("\n=== SAMPLE ALIASES ===")
    sample_canonicals = list(variant_info.keys())[:10]
    for canonical in sample_canonicals:
        info = variant_info[canonical]
        forms = info['all_forms'][:5]
        print(f"  {canonical}: {forms}...")
    
    return alias_dict, variant_info


if __name__ == '__main__':
    main()
