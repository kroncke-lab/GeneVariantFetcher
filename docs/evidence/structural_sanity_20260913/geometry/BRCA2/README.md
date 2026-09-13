# BRCA2 ordered geometry and full-sequence polymer correction

The previous BRCA2 run analyzed only 72 resolved canonical COM positions in two
experimental assemblies. Its builder explicitly supplied **no polymer outside
those constructs**. The AlphaFold API returned no full BRCA2 model, but twelve
AlphaFold v6 fragments were already cached in PPA. That cache was overlooked.
Consequently, the previous 127 supported variants described small bound peptide
regions; the result did not evaluate BRCA2's extensive candidate disorder or its
larger predicted ordered regions. Earlier evidence is preserved unchanged.

This correction maps all **3,418 canonical positions** and supplies 3,248 positions
with usable primary COM/polymer geometry. It is not a complete native biological
unit or a global BRCA2 structure. Experimental assemblies, predicted fragment
frames, and the canonical polymer frame retain separate identities.

| Geometry policy | Experimental + AF structured positions | Polymer positions | Usable COM positions | Unavailable COM positions |
| --- | ---: | ---: | ---: | ---: |
| Primary: experimental, local AF, conservative disorder | 984 | 2,264 | 3,248 | 170 |
| Experimental + conservative disorder comparator | 72 | 2,264 | 2,336 | 1,082 |
| Primary structure + all unresolved assumed polymer | 984 | 2,433 | 3,417 | 1 |

These are residue geometry counts, not supported variant counts. A polymer target
still needs at least one other eligible variant in its contiguous segment.

## Exclusive source routing

Routing applies equally to targets and donors. A canonical position never mixes
experimental, AF and polymer contexts in one density:

1. **Experimental precedence:** retain the earlier 7LDG and 8PBC biological
   assembly coordinates, chain IDs and local frames exactly. Their 73 observed
   canonical Cα positions are excluded from every AF and polymer donor route.
   Native methionine COM remains withheld at chemically modified MSE2322;
   Cα sensitivity remains available there. Experimental COM covers 72 positions.
2. **Predicted ordered positions:** among the other positions, require maximum
   pLDDT across covering fragments ≥70. Each fragment contributes coordinates
   only when its own residue pLDDT is ≥70. This adds 912 unique ordered positions.
   No coordinates are merged across fragments and no distances cross frames.
3. **Conservative candidate disorder:** outside experimental positions, require
   *every* covering fragment's pLDDT to be below 50. These 2,264 positions use one
   canonical polymer frame and one chain. Contiguous segments are defined over
   the full canonical sequence, independently of experimental construct bounds.
4. **Ambiguous/unresolved:** 169 remaining positions stay unavailable in primary
   geometry. Together with MSE2322's native COM, that makes 170 unavailable COM
   positions. The broader sensitivity explicitly assumes those 169 positions
   are polymer, retaining all primary ordered geometry. Adjacent assumed and
   consensus polymer residues form continuous segments in that sensitivity.

The polymer distance is `3.8 * sqrt(abs(i-j))` Å, within the same segment and
chain only. The runner retains the requested normalized sigmoid with half
weight at 3 Å and positive tails beyond 20 Å. Missing or differently assigned
geometry does not create a cross-domain distance. Target-only variant exclusion
still permits other substitutions at the same residue to donate.

## Fragment and disorder evidence

All twelve cached PDB fragments independently have a **unique, exact sequence
match** to canonical human BRCA2, P51587. Each uses local numbering beginning at
1. Verified starts are 1, 201, …, 2201; fragments 1–11 have 1,400 residues and
fragment 12 ends at canonical 3418. Every observed WT residue matches, every
canonical position is covered, and overlapping fragments agree on identity.
There are 16,618 fragment-residue records and 1–7 covering fragments per position.

Across the full canonical sequence, maximum fragment pLDDT is below 50 at 2,315
positions, 50–70 at 191 positions, and at least 70 at 912 positions. There are
137 positions where covering fragments disagree across the 50 threshold. Taking
the maximum for the disorder gate requires unanimity that a residue is below 50,
so a single poor fragment cannot label it disordered. The 51 low-confidence
positions that are experimentally ordered retain their experimental geometry.

Low pLDDT is a candidate-disorder indicator, not proof of intrinsic disorder.
The frozen UniProt/MobiDB-lite annotations are retained in the confidence ledger
as independent evidence. Seven annotated positions, 37–43, disagree with the
max-pLDDT gate: positions 37–40 have confident predicted geometry and 41–43 remain
ambiguous. This disagreement is explicit; the annotation is not silently used
to overwrite the declared confidence policy.

The coordinate metric is the frozen complete side-chain heavy-atom
mass-weighted COM, with Cα for glycine. Every AF sidechain passes completeness
checks. Cα columns are retained for sensitivity. Source fragment models and
pLDDTs remain inspectable in the cache and compressed confidence ledger.

## Interpretation limits

The 12 predicted frames overlap; they are alternative local model contexts,
not independent experiments or biological copies. The existing PPA calculation
normalizes donors within each available context, then averages those contexts.
It does not multiply observations by the number of fragments. Nevertheless,
fragment choice and overlap can affect neighborhood estimates. No per-fragment
PAE files are available in this cache. Even with confident local residues,
relative domain placement within an AF fragment can be uncertain. The
experimental-plus-IDR comparator isolates the effect of adding AF geometry.

The experimental 7LDG unit contains small BRCA2 fragments with MEILB2 and concerns
a meiosis-specific interaction; 8PBC contains TR2 peptides in a finite RAD51/DNA
filament. Those contexts have the same biological limitations as before. Neither
provides a full cancer-relevant BRCA2 unit. Partners remain in frozen raw
assemblies; variant donors are BRCA2 positions only.

Sources: [AlphaFold proteome downloads](https://alphafold.ebi.ac.uk/download),
[versioned human archive](https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v6.tar),
[UniProt P51587](https://rest.uniprot.org/uniprotkb/P51587.json),
[7LDG](https://www.rcsb.org/structure/7LDG),
[8PBC](https://www.rcsb.org/structure/8PBC).
The original assembly mapping and chemical-modification provenance are preserved
in the earlier `missense_structural_extension_20260912/geometry/BRCA2/` report;
the new manifest records its hash along with the original geometry hash.

## Files and reproduction

- `primary_geometry.csv.gz`: primary source-exclusive mapping.
- `experimental_idr_geometry.csv.gz`: experimental and conservative IDR comparator.
- `all_unresolved_polymer_geometry.csv.gz`: primary structure with explicit broader polymer assumption.
- `canonical_confidence_map.csv.gz`: one canonical row per residue, fragment
  confidence range and overlaps, UniProt disorder flags, final source assignment.
- `fragment_confidence.csv.gz`: every mapped per-fragment confidence value.
- `geometry_manifest.json`: source hashes, exact mapping, policies and file hashes.
- `geometry_audit.json`: independent parsed-data and selected PPA reference checks.

```sh
../ProteinProximityAnalysis/.venv/bin/python docs/evidence/structural_sanity_20260913/geometry/BRCA2/build_geometry.py
../ProteinProximityAnalysis/.venv/bin/python docs/evidence/structural_sanity_20260913/geometry/BRCA2/audit_geometry.py
```

No download, paid model call, count refit or scientific outcome fit occurs in
these builders/checks. Inputs are read from the pinned prior evidence and local
AF fragment cache. Every output CSV uses LF, including deterministic gzip
payloads, and remains below 1.15 MB.
