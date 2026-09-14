# Independent five-gene refresh validation

The refreshed inputs, type-specific priors and primary structural calculations
pass the numerical audit. This is validation of the declared descriptive
features and source-correction arithmetic, not certification of every remaining
patient observation or disease-risk calibration. Source-specific decisions and
their remaining endpoint/identity limits are documented in the parent report.

[audit_refresh.py](audit_refresh.py) independently checks each saved complete
union against the frozen eligible population inventory, reconstructs unit counts
from member alleles, and recomputes clinical-key changes from observation-level
decision files. Every eligible genomic allele appears exactly once. Existing
identity quarantines remain explicit in the clinical join ledgers; they are not
mistaken for lost source records or zero counts. The per-gene `*_audit.json`
receipts report each exclusion reason and all-consequence carrier totals.

| Gene | Missense units | Affected | Literature unaffected | gnomAD unaffected | Prior mean | A=0/U=1 own posterior | Supported targets | Supported residues |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| HNF1A | 781 | 400 | 42 | 1,721,270 | 10.2801% | 7.0925% | 722 | 426 |
| GCK | 614 | 603 | 32 | 4,369 | 34.4230% | 25.7804% | 594 | 339 |
| LDLR | 1,361 | 1,346 | 95 | 108,024 | 14.8131% | 10.5940% | 1,192 | 637 |
| BRCA2 | 4,444 | 1,858 | 1,482 | 1,516,774 | 3.8159% | 3.0941% | 4,235 | 2,481 |
| KCNQ1 | 1,061 | 5,845 | 1,594 | 22,640 | 29.4672% | 20.1604% | 975 | 508 |

All ten missense/nonsense fits were recomputed with the exact historical
weighted mean and MSE divisor equal to the number of variant units. Every saved
posterior has **α+A** and **β+literature U+gnomAD U**. Nondegenerate A/U reversal
swaps α and β; adding an affected observation raises the posterior mean and
adding an unaffected observation lowers it, with the prior held fixed. All
actual fits have finite positive α and β. The existing fit rejects boundary
means, zero variance and invalid Beta moments with `ValueError`; no clipping or
fallback prior was introduced. Exact values are in the five `*_prior_checks.csv`
files. The singleton values above are examples of shrinkage, not lower bounds
for variants with larger unaffected counts.

BRCA2's four-carrier change within nonsense is fully explained by three released
canonical frameshift alleles: N1023X (one carrier), D1547X (one) and S2201X (two).
Their old clinical X keys assigned them to nonsense; canonical DNA annotation
and non-triplet indel lengths place them in frameshift. All remain in the
complete union. [Exact identities and transcript evidence](BRCA2_released_frameshift_alleles.csv)
are preserved.

All **67 primary weight shards** were checked across **8,261 missense units**.
Supported rows sum to one, unsupported rows have zero weights and missing final
scores, and every target's own identity has exactly zero weight. Distinct units
at the same residue remain: 495 HNF1A, 433 GCK, 924 LDLR, 3,002 BRCA2 and 764
KCNQ1 targets have such donors. The saved posterior averages, raw-fraction
averages, prior retention/component, affected-weight share, Kish donor count and
conditional donor variance reproduce from the matrices. Maximum posterior-score
reconstruction error is 2.23×10⁻¹⁶.

Fresh reference enumeration checks **55 supported contexts at 34 selected
targets**, including separately selected unsupported targets. The five
`*_context_checks.csv` files preserve direct calculations. Canonical amino-acid
identity matches the geometry manifest for every donor. Selected pair distances
are independently recomputed from COM coordinates within one frame, or from
3.8√|Δresidue| only within the same chain and IDR segment. Mixed/unavailable
geometry does not become a distance. The h=3 sigmoid stays positive; repeated
copies contribute at most one donor vote per context, and supported contexts
are averaged equally. Raw-count diagnostics use the absolute, unnormalized K,
with one prior per context, before context averaging.

Changing a target's posterior α while holding the shared prior, geometry and
eligibility fixed leaves its own excluded score and weights unchanged. A
constant donor-label control remains constant. All five final runs use the frozen
h=3 geometry and the same numerical engine. All five clinical sources now have
bounded corrections; reproducing their old scores is therefore not an acceptance
criterion. The source deltas, rather than an estimator change, explain the
changed inputs.

Zero observed affected donor evidence and missing support remain distinct.
HNF1A has one supported target with no affected donor, LDLR twenty, BRCA2 forty
and KCNQ1 twenty-seven; GCK has none. Their observed-fraction features are
exactly zero, while none of the primary posterior features is zero. These
counts mean no affected donor anywhere in a geometry-permitted context, not
absence of cases within an arbitrary distance cutoff. The `*_zero_case_idr_segments.csv`
files separately list zero-case candidate-IDR contexts; repeated chain contexts
are identified and must not be summed as independent carrier evidence.

The [residue-table audit](plot_table_audit.json) checks all 6,050 canonical
positions, every numeric field, variant multiplicity, support and the GCK/KCNQ1
endpoint sensitivity comparisons. The [layout audit](plot_layout_audit.json)
checks actual plotted primary/own-posterior arrays, visible text bounds and that
all positive diagnostic values lie within the log axis. All six final PNGs were
visually inspected; the [visual receipt](visual_audit.json) pins their hashes.
The two detected rendering issues were repaired: clipped log-axis labels and
positive HNF1A diagnostic values below a fixed axis limit. Exact-zero markers
now sit below every positive value and their explanation sits in the footer.
The same Matplotlib figures were exported to PDF; PDFs were not separately
rasterized for this audit.

[Source-to-fit mappings](clinical_type_scope_changes.csv) and their
[summary](clinical_type_scope_summary.csv) explain why source class totals can
differ from the canonical union. In KCNQ1, incompatible shared-allele keys A344D,
E261L and R555Q remove 8 A/1 U, while canonical-cDNA resolution of A287T and
c.1795G>A (V599M) adds 5 A/0 U to missense: net −3 A/−1 U. The existing identity
rules explain the difference. The [scope audit](scope_mapping_audit.json) also
independently verifies the KCNQ1 cardiac-event sensitivity source arithmetic.

Reproduce from the repository root with the scientific Python environment:

```bash
OPENBLAS_NUM_THREADS=1 ../BayesianPenetranceEstimator/.venv/bin/python \
  docs/evidence/residue_density_refresh_20260914/validation/audit_refresh.py
OPENBLAS_NUM_THREADS=1 ../BayesianPenetranceEstimator/.venv/bin/python \
  docs/evidence/residue_density_refresh_20260914/validation/audit_plot_tables.py
OPENBLAS_NUM_THREADS=1 ../BayesianPenetranceEstimator/.venv/bin/python \
  docs/evidence/residue_density_refresh_20260914/validation/audit_plot_layout.py
OPENBLAS_NUM_THREADS=1 ../BayesianPenetranceEstimator/.venv/bin/python \
  docs/evidence/residue_density_refresh_20260914/validation/audit_scope_mappings.py
```

The audit checks each finalized run input manifest and saved weight-file hash.
[audit_source_hashes.json](audit_source_hashes.json) pins the audit, executed
refresh runner and both PPA numerical modules. No source DB or core module was
modified by this validation.
