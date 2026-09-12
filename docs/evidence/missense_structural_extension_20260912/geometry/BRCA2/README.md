# BRCA2 domain-local experimental geometry

This is an explicitly exploratory analysis of two small bound BRCA2 peptides
in their actual downloaded biological assemblies. It is **not a full-gene
structure**, a merged fragment model, or a successful strict full-unit BRCA2
pipeline. The combined manifest preserves separate `frame_id` values, and
distances between these frames are unavailable.

| Frame | Complex | Construct interval | Resolved canonical positions | Primary COM positions |
| --- | --- | --- | --- | ---: |
| 7LDG assembly 1 | MEILB2–BRCA2 | 2271–2335 | 2276–2282 and 2286–2335 | 56 |
| 8PBC assembly 1 | RAD51/ssDNA–BRCA2 TR2 | 3260–3308 | 3289–3304 | 16 |

Combined primary COM coverage is **72 of 3,418 positions (2.11%)**. The
reference/construct intervals exceed the observed-coordinate intervals; they
must not be reported as resolved coverage. Every frame contains all 3,418
canonical rows for each target chain, with unobserved positions explicitly
missing. No polymer fallback, AlphaFold coordinates, or cross-domain distances
are supplied. The per-frame and combined CSVs are compressed deterministically.
Sources: [7LDG](https://www.rcsb.org/structure/7LDG),
[8PBC](https://www.rcsb.org/structure/8PBC),
[UniProt BRCA2](https://rest.uniprot.org/uniprotkb/P51587.json).

## Chain and residue identity

RCSB and PDBe SIFTS map 7LDG entity residues 2–66 to canonical 2271–2335,
and 8PBC residues 1–49 to 3260–3308. Each mapping also passes an independent
unique exact-sequence match and checks every observed WT residue. The 7LDG
extra construct residue outside the SIFTS mapping has no observed coordinate.

The downloaded **7LDG assembly has B/D/B-2/D-2 BRCA2 fragment chain IDs**.
B and B-2 each resolve 50 positions; D and D-2 each resolve seven. These
fragment IDs are preserved rather than inventing connections between them.
The study describes a 4:2 MEILB2:BRCA2 architecture, while the deposited
assembly counts eight protein fragment chains. Neither representation means
four complete BRCA2 proteins. This structure concerns a meiosis-specific
interaction, which further limits interpretation as a cancer-specific
functional unit. [Primary 7LDG study](https://doi.org/10.1038/s41594-021-00635-0).

**8PBC** contains ten BRCA2 peptide copies, L–U, against eleven RAD51 chains
and a synthetic ssDNA chain. It is a finite filament model; target positions
near the model ends need not have equivalent geometric contexts. Partner
coordinates stay in the raw assembly; only BRCA2 chains appear as candidate
variant donors. [Primary 8PBC study](https://doi.org/10.1038/s41467-023-42830-1).

## COM and modified residues

Primary coordinates require a complete side-chain heavy-atom mass-weighted
center, with Cα for glycine. The only chemical modification in the target
coordinates is selenomethionine at canonical Met2322 in 7LDG B/B-2. Its
canonical parent identity is verified against the official MSE chemical
component record, but its COM is withheld because selenium changes the
native methionine center. Observed Cα is retained with
`ca_geometry_state=structured`. Consequently Cα sensitivity covers 73 unique
canonical positions. Every other observed target sidechain is complete.

Individual `*_ca_canonical_geometry.csv.gz` and combined
`BRCA2_ca_canonical_geometry.csv.gz` files expose the Cα eligibility state
explicitly. They must be paired with the Cα metric; they are not replacement
COM manifests. `coordinate_exclusions.json` records the modified-residue
exclusions. Missingness does not imply intrinsic disorder.

## Reproduce

```sh
../ProteinProximityAnalysis/.venv/bin/python docs/evidence/missense_structural_extension_20260912/geometry/BRCA2/build_geometry.py
```

The builder reuses the sibling **new evidence folder's** KCNQ1 coordinate and
identity helpers, with its source hash recorded. Raw public bytes remain in
the ignored `results/missense_structural_extension_20260912/raw/BRCA2/` cache.
`geometry_identity_report.json` records URLs, hashes, mapping, assembly
metadata, copy counts, source studies, and acceptance checks. This acquisition
does not fit penetrance or run disease density.
