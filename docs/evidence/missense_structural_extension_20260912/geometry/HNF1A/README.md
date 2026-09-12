# HNF1A canonical geometry

The primary frame is **8PI8 biological assembly 1**, containing two wild-type
HNF1A DNA-binding domains bound to the wild-type P2-HNF4A promoter. The independent
**1IC8 assembly 1** comparison contains two HNF1A domains bound to another promoter.
These are **partial experimental biological units**, not full-length HNF1A dimers.
They do not establish distances from the DNA-binding domain to the omitted
N-terminal dimerization domain or C-terminal transactivation region.

| Geometry file | Construct mapped to P20823-1 | Observed structured residues, C / D | Rows |
| --- | --- | --- | --- |
| `8PI8_canonical_geometry.csv` | 83–279 | 174 / 164 | 1,262 |
| `1IC8_canonical_geometry.csv` | 85–278 | 170 / 173 | 1,262 |

Each file retains all **631 canonical positions for each target chain**. Chain
identifiers are mmCIF label chains **C/D**, corresponding to author chains **A/B**.
The two DNA strands are present in the downloaded assembly and partner provenance;
they are not variant-bearing donor chains. In the assembly metadata, the total
polymer count is four because it includes those DNA strands. The protein count is
two. The 8PI8 construct's extra N-terminal tag residue is excluded from canonical
mapping. Every structured coordinate matches its canonical wild-type residue.

Canonical sequence P20823-1 matches the frozen missense annotation reference and
UniProt's link to **ENST00000257555.11**. Mapping uses SIFTS label-sequence segments,
an independently unique sequence placement in the 631-residue reference, and
per-coordinate residue checks. Author numbering is recorded, never assumed to be
canonical. Both structures have complete observed sidechains. The primary point
is the mass-weighted sidechain heavy-atom center of mass, with Cα for glycine;
the same files also support the Cα distance sensitivity.

Missing coordinates alone do not indicate disorder. Unresolved positions use a
candidate IDR assignment from AlphaFold v6 pLDDT <50 or UniProt's **MobiDB-lite
predicted** disorder annotations. These are prediction-based assignments, not
experimental proof of disorder. Positions 1–32 stay unavailable because they form
a known dimerization domain omitted by these constructs. Observed experimental
coordinates take precedence over disorder predictions. Each final IDR segment is
contiguous within one chain; a resolved residue breaks a segment. The runner may
use its sequence-polymer model within such segments, with no cross-chain polymer
distances. No AlphaFold coordinates are inserted into an experimental frame.

Primary 8PI8 has 338 structured, 822 candidate-IDR, 26 ambiguous, and 76 missing
rows across two copies. Thus much of the full-sequence coverage is **polymer-only
prediction support**, not experimental Cartesian geometry. The full AlphaFold
model is retained in ignored raw storage for its confidence values; it is not
substituted for the DNA-bound dimer. The other UniProt-indexed dimerization candidate,
2GYP, contains a D-Ala20 analog and is excluded as a wild-type reference.

Sources: [8PI8](https://www.rcsb.org/structure/8PI8),
[Kind et al., 2024](https://insight.jci.org/articles/view/175278),
[1IC8](https://www.rcsb.org/structure/1IC8),
[P20823](https://rest.uniprot.org/uniprotkb/P20823.json), and
[AlphaFold model metadata](https://alphafold.ebi.ac.uk/api/prediction/P20823).
`geometry_provenance.json` records all URLs, raw-file SHA-256 hashes, assembly
metadata, transcript identity, disorder evidence, and candidate-selection reasons.

Rebuild both genes from frozen raw sources, from the GVF repository root:

```sh
/Users/kronckbm/GitRepos/ProteinProximityAnalysis/.venv/bin/python docs/evidence/missense_structural_extension_20260912/geometry/HNF1A/build_geometry.py
```

Add `--download` only to restore absent raw files; existing source bytes are
preserved and checked against the provenance hashes. Add `--gene HNF1A` to rebuild
only this gene. Raw files live under ignored
`results/missense_structural_extension_20260912/raw/HNF1A/`.
