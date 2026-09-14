**Local visual and scientific review — 14 September 2026**

**Ready within the reviewed scope, with the scientific caveats below.** All six
final PNGs were opened with the local image viewer: the five-gene overview and
the HNF1A, GCK, LDLR, BRCA2 and KCNQ1 individual figures. No material clipping,
title/legend collision, missing axis label, or misleading connection across an
unestimated primary-score residue was found. The final HNF1A/GCK footers were
reopened after the small-number formatting correction. This is a local review;
no figure, source record, count or internal fitted parameter was sent to an
external model. The separate CLI reviews are conceptual methods reviews.

The review read `plot_residues.py`, all five source READMEs, the final residue
tables, `score_summary.csv` and `prior_summary.csv`. A read-only Python check
confirmed complete canonical position sequences, missing primary values where
there are no supported targets, missing own means where there is no observed
variant, all plotted positive log values inside the final bounds, and agreement
between the first/highest supported residue means and their original variant
rows for all five genes. Those checks passed. The independent full numerical
audit is documented in [validation](../validation/README.md); this review does
not replace its weight-matrix or source-accounting checks.

| Gene | Supported / canonical residues | Geometry at supported residues | Fitted missense prior | Primary residue-mean range |
|---|---:|---|---:|---:|
| HNF1A | 426 / 631 | 115 structured, 309 polymer, 2 mixed | 10.2801% | 2.9990–51.5284% |
| GCK | 339 / 465 | 339 structured | 34.4230% | 15.4722–49.0540% |
| LDLR | 637 / 860 | 533 structured, 104 polymer | 14.8131% | 5.2810–34.9220% |
| BRCA2 | 2,481 / 3,418 | 749 structured, 1,732 polymer | 3.8159% | 0.1615–18.3315% |
| KCNQ1 | 508 / 676 | 292 structured, 216 polymer | 29.4672% | 9.3606–60.0522% |

The ranges above are **means across supported variants at each residue**.
They are not the variant-level ranges in `score_summary.csv`, clinical risk
intervals, or confidence bounds.

**What the figures communicate correctly**

- The overview uses the same 0–100% primary axis for each gene, clearly names
  each intended endpoint and shows its actual fitted prior and residue coverage.
  BRCA2 variation is compressed at this scale; its separate log panel makes the
  small scores inspectable without changing the comparison axis.
- Each individual figure separates the posterior neighborhood, the variant's
  own count posterior, and observed counts. The primary line averages only
  supported target-excluded variant scores; the own-posterior line averages all
  observed eligible variants at that residue, including variants without
  geometry support. That is an intentional difference in support. The pale band
  is explicitly a between-variant range, not uncertainty.
- Primary NaNs break the lines. For example, KCNQ1 residue 589 has three observed
  units and residue totals A=777/U=370, yet no supported neighborhood score. Its
  own evidence remains visible in the lower panels while the primary remains a
  gap. GCK residues 99 and 456 have no retained eligible observed units in this
  input and remain unestimated; removing an incompatible source did not create
  negative examples at those residues.
- The support strip is outside the score axis, with structured, polymer, mixed
  and unestimated states named. It shows geometry availability, not confidence
  in a physical contact or a measured global protein conformation. Same-IDR
  polymer distances and separate local frames require the accompanying methods
  and numerical audit; they cannot be established by looking at a curve.
- Log panels use points rather than connections across selected nonzero rows.
  Every positive value fits above the dynamic lower limit. Downward triangles
  identify exact zeros at a documented artificial plotting position, below
  every positive value. The 0.1% reference line is correctly expressed as a
  percentage; it is not a fitted population-risk baseline. Log count panels
  display positive observed totals only, not hypothetical unobserved alleles.
- Distinct same-residue alternatives are retained by the stated variant-only
  leave-one-out rule. These are descriptive features, as the footers state;
  none of the panels establishes disease-probability calibration.

**Selected biological and numerical interpretation checks**

HNF1A's central elevation is compatible with retained MODY evidence, but the
picture does not independently validate that evidence. The [source review](../source/HNF1A/README.md)
separates type-2 diabetes and liver-tumor/adenoma endpoints and corrects one
R272H family to four diabetes-affected and one unaffected at follow-up. The
residue-272 table pools all retained variants and papers at that site (A=29/U=1),
so its 51.5284% neighborhood mean must not be attributed to that single family.

GCK remains elevated after [source-specific endpoint repairs](../source/GCK/README.md).
That is not evidence that every residue causes MODY. At residue 455, the retained
units have A=0/U=75, mean own posterior 10.9791%, and neighborhood score 28.6218%.
Its direct prior component is 18.0556 percentage points. The neighborhood
therefore answers a different question from those variants' own counts. The
reviewed hypoglycemia observations and diagnosis-only proxy remain separately
handled; removing them was not justified by a desired lower curve.

LDLR's [reviewed myocardial-infarction counts](../source/LDLR/README.md) do not
establish FH status; the Greek relative counts instead required restoration to
the affected group. The final plot uses the corrected 14.8131% prior. Residue
775 illustrates the residual prior effect: A=0/U=51, own mean 1.4304%, primary
6.3430%, raw-fraction neighborhood 0.00530%, and direct prior component 6.3378
percentage points. These quantities are not interchangeable penetrance estimates.

BRCA2 is no longer the earlier approximately 10% profile. At residues 865 and
866, A=0 and U=12/8; primary means are 2.4562% and 2.5632%, respectively. Each
equals its direct prior component, while the observed-fraction feature is
exactly zero. These are clear counterexamples to explaining every nonzero
baseline by a distant affected donor. They also do not justify forcing a score
below 0.1% without changing its statistical meaning. The [source review](../source/BRCA2/README.md)
quarantines unsupported hereditary-endpoint assignments while preserving
unknown status and valid observations; its scope is bounded, not a full
adjudication of every remaining paper. Whole-gene own-posterior points are dense,
so exact allele comparisons require the tables rather than reading individual
points from the PNG.

KCNQ1's corrected curve looks very similar to the preceding curve. The final
table confirms the changed A/U input rather than a stale image. The empirical
prior moves only from 29.4760% to 29.4672%, and the largest reviewed count change
is at an unsupported neighborhood site. The [source review](../source/KCNQ1/README.md)
distinguishes genotype-positive carriers, cardiac events and clinical/ECG LQTS.
A broad distribution therefore does not itself demonstrate clinical calibration
or make KCNQ1 an independently validated control for the other genes.

No figure change is required by this review. Source heterogeneity, repeated
people, endpoint adjudication outside the reviewed rows, the shared empirical
prior and the difference between normalized posterior averages and pooled-count
diagnostics remain substantive interpretation limits. No bandwidth or prior was
retuned as part of this review.

**Final artifact receipt (SHA-256)**

The overview is 1,680×1,960 pixels; each individual PNG is 1,800×1,740 pixels.
The renderer is `14fcbf296d9ede45bea24c401b4a3e9f9093098ebab7b12d9531f84aeafcbcd5`;
`plot_input_hashes.json` is
`fb950141354b5135ec45d05926964c02786b9e18d43337d2e12da57515d520bc`.

| PNG under `plots/` | SHA-256 |
|---|---|
| ALL_GENES_RESIDUE_DENSITY.png | `189424cab38057b3d676aeccae2949e2877027a0ae32d3bd989a1b4152eda7d5` |
| HNF1A_RESIDUE_DENSITY.png | `d546a31bfaf7a6e82eaf955c02c51eb4159e94c8fc312f2dd20a23fca45f744f` |
| GCK_RESIDUE_DENSITY.png | `cc98ebecc54d4d816022ebc907d1131e25ce52acb4af3fe8cdd7d7dc1cf4e0ed` |
| LDLR_RESIDUE_DENSITY.png | `ceb484990e219418391ca07126b37238fe5bd25994abc593f05695d346018893` |
| BRCA2_RESIDUE_DENSITY.png | `e338feaa6551b265d29022f69c6520acbfe251deb4915990db99ace0ce62ca00` |
| KCNQ1_RESIDUE_DENSITY.png | `a48295d2217cd886ea598b60b5389fb70bdb5c6e2f67c550004a28cd0c19c008` |

| Residue table under `tables/` | SHA-256 |
|---|---|
| HNF1A_residue_density.csv | `7964f19bc9ab268f051275a86145e4be04e2a8a42bff8312c3970b880752376a` |
| GCK_residue_density.csv | `223ae699937f01e444a33f0e1e1696ad11abffee7235a0b2528ca5ad8fad323d` |
| LDLR_residue_density.csv | `c24b4efb1259551bdaf03484725020a52858793a9498312b5466f930d69cace5` |
| BRCA2_residue_density.csv | `c3bae778cb537eeb94937ed892497fca60340aabbdbf93aabd703ff7ae1fb577` |
| KCNQ1_residue_density.csv | `5ee75eba60249a56bcd40e27b20fddad680b8adbf7ca5a78b0ffe14e3821b563` |
