# LDLR canonical geometry

The primary full-sequence reference is **AF-P01130-F1-model_v6**, consistent with
PPA's explicit LDLR monomer policy. Its sequence exactly matches the **860-residue
P01130-1** reference and **ENST00000558518.6** used for the frozen missense data.
This is a predicted monomer; it does not establish every relative domain position
or reproduce the receptor's ligand-bound, membrane, or endosomal states.

| Geometry file | Structured residues | Use |
| --- | --- | --- |
| `AF_P01130_canonical_geometry.csv` | 625 | Full AF monomer, COM or Cα |
| `AF_P01130_domain_separated_geometry.csv` | 625 | Conservative domain-placement sensitivity, COM or Cα |
| `1N7D_canonical_geometry.csv` | 580 | Partial experimental endosomal ectodomain, complete-sidechain COM or Cα |
| `1N7D_CA_canonical_geometry.csv` | 637 | Partial experimental ectodomain, **Cα metric only** |

Every manifest includes all 860 canonical positions. The primary coordinate is
the mass-weighted sidechain heavy-atom center of mass, with Cα for glycine. AF
Cartesian coordinates require pLDDT ≥70. The cleaved signal peptide, positions
1–21, remains unavailable for mature-receptor geometry. Low-confidence residues
in the annotated membrane helix 789–810 are not labeled IDR. Other unavailable
positions become candidate IDRs only with AF pLDDT <50 or UniProt's MobiDB-lite
predicted disorder annotation; this is not experimental proof of disorder.
Only contiguous final-state IDR segments on the same chain support the runner's
sequence-polymer distance. No experimental/AF coordinate hybrid is constructed.

**Local confidence does not validate domain packing.** Among structured AF Cα
pairs within 20 Å but separated by more than 20 residues, 1,975 of 16,113 have
maximum directional PAE >10 Å; 1,055 exceed 20 Å. These are geometry diagnostics,
not variant counts. The uncertain pairs are listed in
`AF_P01130_uncertain_close_pairs.csv`. The domain-separated sensitivity keeps the
same coordinates but places UniProt-defined domains in independent frames,
preventing comparisons across those frames. The six class-B repeats 397–658 are
retained as one beta propeller; high-confidence unassigned stretches form separate
local blocks. IDR segments remain polymer-only. This conservative comparison also
removes some real interdomain neighbors; it is not a claim that domains never touch.

**1N7D assembly 1** is a monomeric ectodomain structure at pH 5.3 and 3.7 Å
resolution. Its construct maps to canonical positions 22–720; 639 residues have
coordinates. Engineered substitutions **N515Q and N657Q** fail the canonical WT
check and are excluded from both experimental manifests. Another 57 WT residues
have incomplete sidechains, so the complete-COM manifest withholds their COM
geometry. The separate Cα-only manifest retains their observed Cα positions.
It must not be run with metric `com`. The two construct substitutions can affect
the surrounding structure, even though their own coordinates are excluded; this
frame is consequently a sensitivity comparison, not a pristine WT reference.
Glycans and other nonprotein components are retained in raw assembly metadata but
are not penetrance donors. SIFTS placement is independently checked against the
canonical sequence; author-number offsets are never applied blindly.

Recent ApoB-bound candidates 9BD8, 9BDE, and 9COO are inventoried in provenance.
Their deposited full target sequence is not evidence of full observed coverage.
Their target coordinates were not acquired in this bounded extension; no claim
is made that 1N7D is the latest available receptor complex.

Sources: [P01130](https://rest.uniprot.org/uniprotkb/P01130.json),
[AlphaFold metadata](https://alphafold.ebi.ac.uk/api/prediction/P01130),
[1N7D and primary citation](https://www.rcsb.org/structure/1N7D), and
[9BDE](https://www.rcsb.org/structure/9BDE).
`geometry_provenance.json` contains source URLs and hashes, assembly details,
mapping checks, candidates, exclusions, and the PAE audit.

Rebuild from the GVF root using the shared builder:

```sh
/Users/kronckbm/GitRepos/ProteinProximityAnalysis/.venv/bin/python docs/evidence/missense_structural_extension_20260912/geometry/HNF1A/build_geometry.py --gene LDLR
```

Add `--download` to restore absent raw files. Existing bytes are preserved and
validated against provenance hashes. Raw files and the full PAE matrix are under
ignored `results/missense_structural_extension_20260912/raw/LDLR/`.
