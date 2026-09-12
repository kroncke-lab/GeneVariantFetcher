# KCNQ1 cardiac biological-assembly geometry

Primary geometry is the actual **9U7F biological assembly 1**: four human
KCNQ1 chains (A/D/G/J), four KCNE1 chains, and four calmodulin chains. **9UC8
assembly 1** is a separate PIP2-added state for sensitivity. They describe
closed and open channel conformations, respectively. All KCNQ1 coordinate WT
residues and the deposited canonical reference sequence match P51787 exactly;
RCSB and PDBe SIFTS agree on the 1:1 canonical mapping. Sources:
[9U7F](https://www.rcsb.org/structure/9U7F),
[9UC8](https://www.rcsb.org/structure/9UC8), and
[the primary study](https://www.nature.com/articles/s41422-025-01182-9).

These are **partial resolved biological tetramers, not complete full-length
experimental structures**. The primary study's construct contains KCNQ1
residues 76–620, fused through a linker to KCNE1. PDB's 676-residue reference
SEQRES does not establish that the experiment contained the full-length
protein. The mapping retains all 676 canonical rows for each KCNQ1 copy and
explicitly records whether a position lies in that experimental construct.

| Frame | Per-copy canonical rows | Observed Cα | Complete COM | Candidate IDR | Ambiguous | Missing COM geometry |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 9U7F | 676 | 360 | 359 | 256 | 32 | 29 |
| 9UC8 | 676 | 339 | 339 | 252 | 38 | 47 |

The four copies have the same state counts. Main `geometry_state` follows
complete side-chain heavy-atom mass-weighted COM, with Cα for glycine. Lys569
in 9U7F lacks CG/CD/CE/NZ in every copy: its COM is withheld, while observed Cα remains
available with `ca_geometry_state=structured` for an explicitly selected
metric sensitivity. No incomplete side-chain centroid is substituted.

At absent experimental positions, independently mapped AlphaFold model v6
pLDDT below 50 identifies candidate IDR; 50–70 remains ambiguous, and higher
confidence remains unavailable. Low confidence is an operational proxy, not
experimental proof of disorder. This rule is independently supported by the
canonical AF sequence, not inferred from missing coordinates or experimental
B-factors. Candidate segments are split at every non-IDR position and qualified
by frame and chain. The same-chain, same-segment polymer rule applies; missing
tails receive no invented global coordinates. No AlphaFold coordinates are
spliced into either experimental frame.

The raw assembly retains KCNE1, CaM, ligands, and the complete experimental
coordinate context. `geometry_identity_report.json` records partner identities,
stoichiometry, observed counts, mapping, and raw source hashes. Geometry CSVs
contain only KCNQ1 chains so partner residues cannot become extra variant
donors. The downstream density runner must exclude the target variant from
all four copies and collapse donor copies according to the declared rule.

Reproduce from the GVF root:

```sh
../ProteinProximityAnalysis/.venv/bin/python docs/evidence/missense_structural_extension_20260912/geometry/KCNQ1/build_geometry.py
```

Public raw bytes are cached in the ignored
`results/missense_structural_extension_20260912/raw/KCNQ1/` directory. A rerun
checks source hashes against the frozen geometry report. This acquisition
does not run density or fit penetrance models.
