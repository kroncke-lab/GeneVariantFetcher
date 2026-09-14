**Methodological Review of Residue-Neighborhood Plotting Framework**

This independent review evaluates the proposed generic residue-neighborhood plotting framework. It analyzes the statistical methodology, identifies structural failure modes, and proposes falsifiable diagnostics to ensure the framework’s validity as an exploratory descriptive tool. It strictly avoids treating descriptive features as independently calibrated clinical probabilities.

### 1. Implementation Failure Modes and Falsifiable Checks

The framework relies on explicit mathematical conventions and data mapping steps that are susceptible to silent implementation errors. These must be tested against intentional model assumptions.

*   **Fitting Formula & Prior Orientation:** The historical weighted mean/MSE formula for the Beta prior parameters ($\alpha, \beta$) relies on $y=A/n$ and $w=1-1/(n+0.01)$. A catastrophic failure mode is swapped orientation ($\alpha$ tracking unaffected, $\beta$ tracking affected), which inverses the target signal.
    *   *Falsifiable Check:* Input a synthetic dataset where every variant has $A=100, U=0$. The resulting $\alpha$ must be strictly, massively greater than $\beta$. If $\alpha < \beta$, the implementation is inverted. Furthermore, verify that the variance calculation ($v$) strictly divides by $M$ (the number of eligible units), not the sum of weights, to ensure adherence to the stated historical convention rather than standard weighted variance.
*   **Target Exclusion & Same-Residue Alternatives:** The neighborhood feature $D_i = \sum W_{ij} p_j$ must exclude the *target variant unit* $i$, but retain distinct variants at the same residue.
    *   *Falsifiable Check:* Create a synthetic dataset containing exactly two identical variant units at the same residue (same $A, U$). The neighborhood score $D_i$ for one must exactly equal the own posterior $p_j$ of the other. If $D_i$ is 0 or missing, same-residue alternatives are incorrectly excluded. If $D_i > p_j$, the target was erroneously included in its own neighborhood.
*   **Coordinate Mapping:** Mapping 1D sequence to 3D geometry invites off-by-one errors or structural misalignment.
    *   *Falsifiable Check:* Introduce a synthetic point mutation with a massive $A$ count at a known catalytic residue. Its 3D neighbors (via Angstrom distance) should overwhelmingly weight other structurally adjacent active-site residues. If the highest weights map to a distant surface loop, the coordinate frames or sequence-to-structure indices are misaligned.

### 2. Normalization, Flat Surfaces, and Spatial Influence

The framework uses a smooth kernel $K(d)$ and normalized weights $W_{ij}$. This design can unintentionally mask local data sparsity and create artificially "flat" risk surfaces that appear confident but are merely empty.

Because weights are normalized ($\sum W_{ij} = 1$), even if the raw kernel mass $\sum K$ is infinitesimally small (e.g., all neighbors are far away), the distant neighbors' posteriors are still averaged at full weight. When combined with repeated sparse-variant shrinkage (many variants having $A=0, U=1$ reverting heavily to the Beta prior), distant, sparse neighborhoods will invariably average to a flat value near the prior mean. The steep tail decay of $K(d)$ does *not* prevent this; normalization specifically circumvents raw tail decay in empty regions.

*   **Diagnostics to Distinguish Flatness from Consensus:**
    *   *Raw Kernel Mass & Nearest Donor:* Track $\sum K$ and the distance to the nearest donor with $n_j > 0$. A flat surface with high raw mass indicates genuine regional consensus. A flat surface with near-zero mass indicates a void reverting to the prior.
    *   *Normalized-Weight Radii:* Calculate the physical distance enclosing 80% of the normalized weight. If this radius exceeds physically meaningful bounds (e.g., 20+ Angstroms), the local score is a statistical artifact of distant shrinkage, not a local neighborhood effect.
    *   *Diagnostic Sums:* Compare the primary feature $D_i = \sum (W_{ij} p_j)$ against the diagnostic one-prior kernel pool: $(\alpha + \sum K \cdot A) / (\alpha + \beta + \sum K \cdot n)$. The latter is sensitive to absolute $K$ scale. In an empty region, the primary feature force-averages distant shrunk posteriors, whereas the diagnostic pool elegantly reverts to exactly $\alpha/(\alpha+\beta)$ due to $\sum K \approx 0$.

### 3. Biological Plausibility and Confounding

The framework treats variants broadly, which risks obscuring fundamental biological realities.

*   **Variant Effects:** The model assumes a generic spatial neighborhood effect. However, not every substitution at a structurally sensitive residue causes disease. Tolerated variants exist even in critical domains. Furthermore, the framework does not currently distinguish between loss-of-function (often broadly destructive, matching neighborhood trends) and gain-of-function (highly specific, often defying local averages).
*   **Confounding Variables:** Mixing germline and somatic records, or variants associated with opposite phenotype directions, will scramble the spatial signal into a meaningless average.
*   **Count Interpretations:** Treating categorical catalogue classifications (e.g., a "Pathogenic" database flag) as observed patient counts ($A, n$) is statistically invalid. A catalogue entry is a curated endpoint judgment, not an empirical exposure. Assigning a flag as $A=1, U=0$ arbitrarily defines the denominator and creates fake statistical confidence.

### 4. Truthful Plot Semantics

Visualizations must strictly reflect what the model computes, avoiding implied certainty or false continuity.

*   **Plotting Rules & Pass/Fail Checks:**
    *   *Unsupported vs. Zero:* Regions lacking geometry or nearby variant data must be plotted as explicitly missing, not zero. *Pass/Fail Check:* A region with entirely missing data must not render as a continuous line at $y=0$ or $y=\text{prior\_mean}$.
    *   *Interpolation:* *Pass/Fail Check:* Long gaps in 1D sequence or breaks between 3D domains must not be bridged by line interpolation. Lines must break.
    *   *Multiplicity:* Multiple variant units map to single residues. A plotted residue summary must declare its aggregation method. *Pass/Fail Check:* A residue containing one pathogenic variant ($p=0.9$) and one benign variant ($p=0.1$) should not simply plot as $0.5$ without a visual indication of the between-variant spread.
    *   *No Invented Uncertainty:* The framework provides point estimates. Do not invent error bars, confidence bands, or credible intervals for the spatial feature, as the statistical model does not generate a posterior distribution for the neighborhood sum itself.

### 5. Source Validation and Independent Evaluation

Internal consistency—such as the spatial feature correlating well with the own-posterior labels used to construct it—is mathematically guaranteed by construction and does *not* constitute clinical calibration or valid risk prediction.

*   **Source Integrity:** The input matrix must be rigorously validated before tuning target thresholds. Are duplicate individuals in overlapping cohorts deduplicated? Genuine distinct relatives represent observations, but their statistical dependence must be modeled or acknowledged.
*   **Independent Evaluation:** True evaluation requires data external to the fitting process.
    *   *Functional Assays:* Deep mutational scanning (DMS) can independently validate if the spatial neighborhood correctly identifies biochemically sensitive regions, but DMS scores are not human disease probabilities.
    *   *Clinical Cohorts:* True clinical calibration requires a temporally or geographically distinct hold-out cohort. One must compare the model's $D_i$ scores against the actual, out-of-sample incidence of disease in newly sequenced carriers.
    *   *Status of Output:* Until validated against strictly independent clinical endpoints, all outputs—both the primary normalized feature and the raw-fraction diagnostics—remain descriptive, explicitly labeled exploratory features, not calibrated probabilities of human disease.
