# Historical kernel and biological-unit preflight — 2026-09-12

Final execution choices are in [PLAN.md](PLAN.md). In particular, the final
primary polymer path requires both endpoints in the same IDR on the same chain;
mixed ordered/IDR geometry stays unavailable. Earlier suggestions below to
route a pair when either endpoint is disordered are superseded by that rule.
The user subsequently selected a positive-tail sigmoid with half weight at 3 Å
and strong downweighting toward 20 Å. The compact sine and N<=2 statements
below describe the recovered historical method, not the current primary plan.

Read-only inspection of ProteinProximityAnalysis (PPA), its saved outputs, BayesianPenetranceEstimator (BPE), and the author's public historical code. No sibling code, data, or Git state changed. No structure prediction or full analysis run performed.

## Recovered original implementation

Author source: [Bayes_BrS1_Penetrance/func_dist_seq.R](https://github.com/kroncke-lab/Bayes_BrS1_Penetrance/blob/master/func_dist_seq.R), retrieved directly from raw GitHub. The function explicitly excludes `deval$var != var`: only the target variant is omitted, and other substitutions at the same residue remain donors. Each donor variant has one kernel weight; it is not multiplied by its carrier count.

The original sine kernel is:

```
K(d; a) = 1                         d <= a - 3.14
          0.5 - 0.5*sin((d-a)/2)     a - 3.14 < d <= a + 3.14
          0                         d > a + 3.14
```

Therefore `a=3 Å` is a literal historical-form implementation of the requested half weight at 3 Å. Its weight at distance zero is approximately 0.99875, at 3 Å exactly 0.5, and it reaches zero at approximately 6.14 Å. Use exact pi instead of the historical rounded 3.14 only if explicitly recorded as a numerical correction. A clean reparameterized version can set `K(0)=1`, `K(h)=0.5`, and `K(2h)=0` using `0.5*(1+cos(pi*d/(2*h)))` on `[0,2h]`; this is a rescaled sine taper, not byte-for-byte historical replication. Default midpoint h=3; sensitivity h=2, 4, 5 Å.

Historical sigmoid is `K(d)=1/(1+exp(d/2))`. It has K(0)=0.5; its **relative** half weight is at `2*ln(3)=2.197 Å`, so the argument `7` in a sigmoid funcdist call does not set a 7 Å bandwidth. A normalized tunable version is `2/(1+exp(log(3)*d/h))`, which has exactly unit weight at zero and half weight at h. Use h=3 as the sigmoid sensitivity rather than mislabeling the historical argument 7.

Historical polymer fallback is `d=3.8*sqrt(abs(i-j))`, same sequence, within ±30 residues when no structure-distance list is available for the target. With the requested h≈3 sine cutoff, only sequence separations N=0,1,2 survive: distances are 0, 3.8, 5.37, then 6.58 Å for N=3. This is appropriately very local but can leave zero donors after target exclusion. Return an explicit missing density plus zero support; do not manufacture confidence or broaden the kernel silently. The sigmoid remains nonzero farther away.

The [2020 SCN5A report](https://github.com/kroncke-lab/Bayes_BrS1_Penetrance/blob/master/SCN5A-BrS1-penetrance-report.Rmd), line 392, actually called sigmoid on empirical posterior means `p_mean_w`. The [repository README](https://github.com/kroncke-lab/Bayes_BrS1_Penetrance) describes `distance_file` as distances between centroid atoms. The exact centroid construction is not documented in the inspected function. Thus a side-chain center of mass is available and reasonable, but should not be claimed to exactly reproduce the historical centroid definition until the distance-file producer is recovered.

The later [KCNH2 analysis](https://github.com/kroncke-lab/LQTS-Penetrance-APC-MAVE-Events/blob/main/predict_penetrance_kcnh2-v7.Rmd) loads an external `func_dist_seq.R` from `Bayes_KCNH2_LQT2_Penetrance` and `5va1.dists.txt`. Its shown call passes the literal string `"var"` as the excluded variant while iterating residues (lines 419–421); that call is not an implementation of the user's requested variant-specific LOO. Do not reuse its precomputed density as the new result.

## Existing PPA versus required changes

Local source root: `/Users/kronckbm/GitRepos/ProteinProximityAnalysis`.

- `src/alphafold_rin/penetrance_density.py:163` and all density paths currently use `exp(-d/sigma)`, default sigma=7 Å. No sine/sigmoid selection exists.
- `prepare_variants` at line 66 retains resnum, chain, penetrance, carrier weight and mechanism annotations; it drops the canonical key. Preserve a stable variant ID and normalized allele identity through every observation expansion. Target-ID exclusion must remove **all target chain copies**, while retaining every distinct donor variant at its residue.
- `_expand_variants_to_residue_indices` at line 202 expands sequence-level observations to every matching target-protein chain. Keep variant identity attached and do not count homomer copies as independent clinical donors. Report unique donor count/effective support after grouping copies. Inter-chain geometry can affect the kernel, but copy multiplicity must be specified rather than accidentally multiplying statistical precision.
- Current helpers multiply kernel weight by `total_carriers`. To replicate the historical requested empirical-posterior feature, set donor response to `posterior_empirical_mean` and use equal per-variant donor weighting as primary. Counts already determine each donor's shrinkage. An explicitly labeled precision-weighted sensitivity can be evaluated later, not quietly retained as the default.
- `src/alphafold_rin/disorder.py:72` implements tunable `b*N^nu`, currently b=5.5 Å and nu=.55. Requested historical baseline is b=3.8 Å, nu=.5.
- Current disorder routing at density lines 471, 617, and 811 is **target-only**. A structured target still sees an unstructured donor's unreliable coordinates. New routing must inspect both endpoints: same-chain pair uses polymer if either endpoint is disordered; cross-chain IDR geometry is unavailable without a justified ensemble/tether model. Keep structured-to-structured geometry from the biological assembly.
- A same-chain endpoint rule is a conservative first implementation. Structured residues separated by an IDR need confidence in relative domain placement; PAE/assembly context should flag unreliable interdomain geometry rather than assuming high local pLDDT establishes the entire arrangement.
- `src/alphafold_rin/structure.py:212` already supports `ca`, `cb`, and `com`; `com` is the **mass-weighted side-chain heavy-atom center**, with glycine fallback to Cα. `distances.coordinate_array(...use_ca=False)` uses the chosen coordinate. However, the pipeline currently extracts mode `ca`, and coordinate/fragment density paths hardcode `ca_coord`. Therefore choosing `com` at the parser alone is insufficient: add an explicit density-distance metric propagated identically through dense, coordinate and fragment paths. Proposed primary is documented side-chain centroid/COM; Cα is a sensitivity. Exact legacy centroid reproduction remains unconfirmed.

## Biological-unit readiness and first pilot

| Gene | Readiness | Next structural step |
|---|---|---|
| GCK | Best initial target: true monomer, 465-aa P35557. No cached GCK PPA output found. | Use full-length correctly mapped monomer; compare an experimental biological monomer assembly. 1V4S assembly1 (448 modeled residues, no mutations, 2.3 Å) and 1V4T assembly1 (424, no mutations, 3.4 Å) are active/inactive-conformation candidates, not silently interchangeable. Both are labeled isoform2 by RCSB, so map the shared sequence explicitly to the chosen clinical transcript/P35557 isoform before calculation. |
| LDLR | PPA has an explicit full-length AF monomer policy for P01130, but no cached LDLR structure found. | Good second target after canonical numbering and domain/IDR geometry checks. The local override states there is no curated oligomeric unit for this workflow; do not overstate that as proof of all functional contexts. |
| HNF1A | UniProt P20823 says DNA-binding dimer, with PCBD1 interactions/dimer-of-dimers context. No local output found. | Select the intended DNA-bound/dimer biological unit and partner context; no silent monomer shortcut. |
| KCNQ1 | Saved 3BJ4 assembly1 has four chains but only 150 total residues, 5.5% canonical coverage. Also a full AF monomer and unrun A4 prediction inputs exist. | Reject saved 3BJ4 as a full-gene unit. Obtain an appropriate numbered channel assembly/model including the chosen partner context; local A4 inputs are useful but not a completed biological unit. |
| BRCA2 | Twelve local AF fragments; strict PPA policy calls for full-length prediction. | No cross-fragment distance invention. Defer full-gene assembly result; fragment-domain analyses must remain labeled and IDR polymer can provide sequence-local information. |

Local KCNQ1 evidence:
`output_biounit_20A/biounit_20A_summary.csv`,
`output_biounit_20A/KCNQ1/structure/KCNQ1_assembly.json`,
`output_biounit_20A/multimer_jobs/KCNQ1_P51787_homooligomer_A4_job.json`.

Local BRCA2 evidence: `output/BRCA2/structure/fragments/`; policies in `src/alphafold_rin/biounit.py:127` and `:150`.

Public primary preflight sources read 2026-09-12:

- [UniProt GCK P35557](https://rest.uniprot.org/uniprotkb/P35557.txt): SUBUNIT monomer; canonical sequence and alternative isoforms.
- [GCK 1V4S](https://www.rcsb.org/structure/1V4S) and [GCK 1V4T](https://www.rcsb.org/structure/1V4T): authors' biological assembly1 is monomer A1; modeled coverage and construct annotations.
- [UniProt HNF1A P20823](https://rest.uniprot.org/uniprotkb/P20823.txt): dimer/PCBD1 subunit context.

## Numbering acceptance requirements

Keep accession plus isoform/transcript, reference sequence hash, assembly ID/state, entity IDs, author chain/residue/insertion code, canonical residue, reference amino acid, and mapping status. Apply biological-assembly symmetry mates, map target and partner chains to accessions, retain chain-qualified IDs, and check expected target copy count. Primary variants anchor only to target-protein chains.

PPA `sifts.py:108` provides best-effort author-to-UniProt mapping and can silently retain author numbering when mapping is missing/low coverage. For this analysis, require a positive mapping/sequence-alignment record for every scored residue; no silent assumed offset. Missing experimental coordinates do not automatically mean disorder: use independent disorder/confidence evidence and distinguish unresolved ordered positions. Do not invent a global 3D frame by appending monomer coordinates to an experimental assembly.

## Focused acceptance tests for execution

1. Perturb target counts radically: its LOO feature cannot change through its own donor row; a different substitution at the same residue still changes it. Target aliases and every expanded chain copy remain excluded.
2. Equivalent relabeling/duplication of assembly chains cannot create extra independent clinical support; all real interface contributions remain represented under the declared geometry aggregation.
3. Structured/IDR pair reversal chooses the same distance mode; cross-chain IDR pairs never acquire a fabricated polymer distance from matching sequence numbers.
4. Sine/sigmoid boundaries and half-distance values are checked numerically. IDR N=0,1,2,3 tests establish expected very short support.
5. Dense/coordinate/fragment implementations agree when they share the same frame and metric; explicitly selected COM cannot silently revert to Cα.
6. Zero donors returns missing density/zero support plus declared empirical-baseline fallback at the subsequent prediction stage, never a manufactured structural feature.
