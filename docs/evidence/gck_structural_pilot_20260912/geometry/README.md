# GCK biological monomer geometry

The pilot uses separately mapped coordinate frames for the closed,
glucose/activator-bound **1V4S**, the glucose-free super-open **1V4T**, and the
full canonical **AF-P35557-F1 model v6**. Both experimental files are the actual
author-assigned biological assembly 1, with one protein chain; resolution is
2.3 and 3.4 Å, respectively. These are monomeric units with explicit unresolved
residues, rather than complete experimental coverage of every amino acid.
Sources: [RCSB 1V4S](https://www.rcsb.org/structure/1V4S),
[RCSB 1V4T](https://www.rcsb.org/structure/1V4T),
[UniProt GCK](https://rest.uniprot.org/uniprotkb/P35557.json).

| Frame | Canonical rows | Structured | Polymer eligible | Ambiguous | Missing |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1V4S | 465 | 446 | 1 | 3 | 15 |
| 1V4T | 465 | 423 | 24 | 3 | 15 |
| AlphaFold v6 | 465 | 455 | 1 | 9 | 0 |

“Polymer eligible” describes a residue's geometry class, not guaranteed donor
support. All three frames include the single low-confidence terminal Q465;
this alone supplies no neighbors at other residues. A target still needs
eligible other variants in the same segment.

## Numbering and coordinate checks

The displayed UniProt P35557-1 sequence has 465 amino acids. Isoform 2 changes
the first 15 canonical residues into a different 16-residue sequence, as
specified by UniProt feature VSP_002074. Consequently, matching author numbers
or isolated amino-acid letters at the N terminus is insufficient.

RCSB's SIFTS mapping identifies the shared tail: 1V4S entity positions 6–455
and 1V4T entity positions 2–451 both map to canonical positions 16–465. The
builder independently checks that this entire 450-residue sequence occurs
exactly once in the canonical sequence and validates every observed WT
residue. It excludes two observed isoform-specific residues in 1V4S and one
construct residue in 1V4T, retaining all canonical rows with missingness flags.
The AlphaFold sequence and all 465 coordinate residues match the full
canonical sequence exactly.

The primary metric is the mass-weighted center of complete side-chain heavy
atoms, with Cα for glycine. All observed sidechains in the downloaded models
are complete under this check; an incomplete sidechain makes the builder fail.
Separate Cα coordinates permit metric sensitivity. Experimental B-factors stay
in their own column and never determine disorder. The `plddt` column always
comes from the separately mapped AlphaFold model.

## State-specific disorder

Canonical E157–N179 is unresolved and experimentally disordered in the
super-open 1V4T structure. The source describes exactly these 23 residues and
their disorder; each position was checked against the SIFTS mapping, mmCIF
sequence scheme, and canonical WT identity. This segment therefore uses the
same-chain polymer rule in 1V4T even though the AlphaFold model has high
confidence there. The corresponding resolved 1V4S loop retains its measured
coordinates. See the primary study's results and modeling methods:
[Molnes et al., ATP binding and GCK conformational changes](https://pmc.ncbi.nlm.nih.gov/articles/PMC3531626/).

`1V4T_missing_loop_control.csv` preserves the alternative that leaves this loop
unavailable. The three main manifests use the biological state-specific rule.
No inferred or AlphaFold coordinates are spliced into either experimental
frame. Polymer distances apply only when both endpoints belong to the same
disordered segment; mixed structured/disordered pairs remain unavailable.

For AlphaFold geometry, pLDDT ≥70 permits coordinates; <50 marks candidate
disorder; 50–70 remains ambiguous. Low confidence is an operational proxy,
not independent experimental proof of intrinsic disorder.

## Reproduce

From the GVF root, using a Python environment with Biopython, NumPy and requests:

```sh
../ProteinProximityAnalysis/.venv/bin/python docs/evidence/gck_structural_pilot_20260912/geometry/build_geometry.py
```

Raw public files are cached in the ignored `results/gck_structural_pilot_20260912/raw/`
directory. `geometry_identity_report.json` records exact URLs, source hashes,
sequence checks, assembly metadata, and the disorder source. Compact canonical
coordinate manifests are the durable inputs to the density runner. The
downloaded PAE file is retained as provenance; no PAE threshold was used in this
pilot's geometry eligibility.
