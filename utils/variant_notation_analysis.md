Analysis of Variant Notation Issues - KCNH2 Phase 3

KEY FINDINGS FROM DISCREPANCIES ANALYSIS:

1. MOST COMMON NOTATION FORMATS:
   - Single-letter amino acid codes: "A561V", "R134W", etc.
   - Frame-shift notation: "A1017fsX" → normalized to "A1017fs*"
   - Three-letter codes in Excel: "p.Ala561Val" (when normalized)

2. CRITICAL MISMATCH PATTERNS:
   
   A. SINGLE vs THREE-LETTER CODES:
   - Excel: "A561T" but SQLite: "p.Ala561Thr" (no match)
   - Excel: "R752W" but SQLite: "p.Arg752Trp" (no match)
   
   B. MISSING 'p.' PREFIX:
   - Excel variants lack 'p.' prefix but SQLite expects it
   - Classic case: "A561V" vs "p.Ala561Val"
   
   C. FRAMESHIFT NOTATION:
   - Excel format: "A193fsX" 
   - Expected format: "p.Ala193fs*" (A193 → Ala193, fsX → fs*)
   
   D. NON-STANDARD SUFFIXES:
   - "sp" endings (e.g., "A715sp", "D864sp")
   - "Del" vs "del"

3. INSERTION/DELETION NOTATION:
   - "G189Ins" vs "p.Gly189ins"
   - "A671Del" vs "p.Ala671del"

4. POSITION CONSISTENCY:
   - All positions seem to be 1-based amino acid positions
   - Gene length (KCNH2) = 1159 amino acids

5. WORKING EXAMPLES FROM DATA:
   
   FROM DISCREPANCIES:
   - ''31020140,A121fsX,A121fs*,p.Gly121fs*,p.Gly121fs*,fuzzy,0.8888888888888888''
     → Shows A121fsX → p.Gly121fs* conversion needed
   
   - ''29752375,A190T,A190T,p.Ala190Thr,p.Ala190Thr,exact,1.0''
     → Shows correct A190T → p.Ala190Thr conversion

COMPARISON DATA:
- Total variants in missing_in_sqlite: 920 variants
- All show "none" match_type, indicating zero normalization/recall
- Most frequent patterns: single-letter codes that need three-letter + 'p.' prefix

PRIORITY FIXES NEEDED:
1. Single-letter → three-letter conversion (A561V → p.Ala561Val)
2. Handle fsX → fs* frameshift notation
3. Normalize case: 'Ins' → 'ins', 'Del' → 'del'
4. Add 'p.' prefix where missing
5. Handle 'sp' suffix endings
6. Support insertion/deletion notation normalization