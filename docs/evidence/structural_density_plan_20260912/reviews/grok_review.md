# Review of the requested structural penetrance-density plan

This is a methodology review only. Primary analysis is the four user-specified choices. Leave-one-residue-out, dropping other alleles at the target residue, demanding an external cohort before starting, and re-litigating gnomAD-as-unaffected are out of scope. AlphaMissense is a comparator, not a donor prior.

Throughout: **fact** = identity or sampling mathematics; **recommendation** = default that implements the request without silently changing it; **sensitivity** = labeled alternate that must not overwrite the primary run.

---

## 1. Target-variant leave-one-out (keep other substitutions at the same residue)

**What must be left out.** The statistical and structural unit to exclude is the *substitution*, not the residue. If the target is `GENE p.Arg123Cys` in this gene–disease set, every structural image of that same substitution is removed from the donor pool: every chain, every biological-assembly copy, every repeated domain instance that maps to the same canonical variant key. Other alleles at the same canonical residue (`p.Arg123His`, `p.Arg123Gly`, …) stay in as donors.

**Why “every structural copy.”** In a homomer or a duplicated domain, the same variant occupies several coordinates. Leaving one image in is direct leakage: the target’s own posterior would enter the kernel sum at a short distance, often at \(d=0\) after superposition. Unique key: `(gene, canonical_protein_change, disease)`. Spatial query may return several hits; all of them are ineligible when that key is the target.

**Direct vs indirect influence when hyperparameters are fit on the full set.**

- **Direct (forbidden):** target used as a 3D/polymer donor at any copy.
- **Indirect (allowed in the requested primary run):** target included in the gene–disease moment calculation that produces a single \((\alpha_{\mathrm{emp}},\beta_{\mathrm{emp}})\). That is a \(O(1/J)\) effect on the common prior, not a local spike at the target site.

**Recommendation.** Primary run = **fixed full-dataset** \((\alpha_{\mathrm{emp}},\beta_{\mathrm{emp}})\) + **hard variant-key exclusion** of the target from all donor lists. Do not re-estimate the Beta for every left-out variant in the primary run.

**Sensitivity (optional, not primary).** Refit \((\alpha,\beta)\) after dropping the target (and only the target). Report \(\Delta\alpha,\Delta\beta\) and \(\Delta\rho\) at the target site. This quantifies indirect leakage; it is not a reason to switch the primary definition.

**Implementation invariant.** For target residue \(r\) and target amino-acid change \(a\to a^*\): donors at \(r\) with \(a\to a'\) for \(a'\neq a^*\) have positive kernel eligibility; donors with key equal to the target have weight exactly 0, including assembly copies.

---

## 2. Common empirical Beta prior from full-dataset mean and MAE/MSE

Notation. Variant \(j=1,\ldots,J\) in the **eligible gene–disease set** (not a pool of unrelated diseases). Counts \(A_j\) affected, \(U_j\) unaffected (gnomAD counted as unaffected, as specified). \(n_j=A_j+U_j>0\), \(p_j=A_j/n_j\). Variants with \(n_j=0\) are ineligible, not \(p=0\).

Posterior donor mean after the common prior (this is the quantity that builds the density, not raw \(p_j\) and not AlphaMissense):

\[
\tilde p_j=\frac{\alpha_{\mathrm{emp}}+A_j}{\alpha_{\mathrm{emp}}+\beta_{\mathrm{emp}}+n_j}.
\]

### 2.1 Which mean and which second moment are coherent

A \(\mathrm{Beta}(\alpha,\beta)\) on variant penetrance has

\[
\mu=\frac{\alpha}{\alpha+\beta},\qquad
\sigma^2=\frac{\mu(1-\mu)}{s+1},\qquad s=\alpha+\beta.
\]

Inversion (moment matching) is

\[
s=\frac{\mu(1-\mu)}{\sigma^2}-1,\qquad
\alpha=\mu s,\qquad \beta=(1-\mu)s.
\]

**Fact — three different functionals.**

| Estimator | Definition | What it is the moment of |
|---|---|---|
| Unweighted variant mean | \(\hat\mu_{\mathrm{unw}}=J^{-1}\sum_j p_j\) | A randomly drawn **variant** |
| Carrier-weighted mean | \(\hat\mu_{\mathrm{w}}=\sum_j n_j p_j/\sum_j n_j=\sum A_j/\sum n_j\) | A randomly drawn **carrier** |
| Binomial sampling noise | \(\mathbb{E}[(p_j-\pi_j)^2\mid \pi_j]=\pi_j(1-\pi_j)/n_j\) | Measurement error around the latent penetrance \(\pi_j\) |

If \(\pi_j\sim\mathrm{Beta}(\alpha,\beta)\) i.i.d. and \(n_j\perp\pi_j\), both means are unbiased for \(\mu\), but the **MSEs are not interchangeable**.

Unweighted MSE of the **observed** frequencies:

\[
\widehat{\mathrm{MSE}}_{\mathrm{unw}}=J^{-1}\sum_j(p_j-\hat\mu_{\mathrm{unw}})^2
=\underbrace{\mathrm{Var}(\pi)}_{\text{Beta process}}+\underbrace{\text{average of }\pi(1-\pi)/n}_{\text{binomial noise}}+\text{finite-}J\text{ error}.
\]

Using \(\widehat{\mathrm{MSE}}_{\mathrm{unw}}\) as \(\sigma^2\) **overstates** process variance and **understates** \(s\) whenever some \(n_j\) are small. A first-order noise correction (still method of moments) is

\[
\hat\sigma^2_{\mathrm{proc}}=\widehat{\mathrm{MSE}}_{\mathrm{unw}}-J^{-1}\sum_j\frac{p_j(1-p_j)}{n_j}.
\]

The subtracted term is slightly biased because \(p(1-p)\) is biased for \(\pi(1-\pi)\); it is adequate as a documented MOM correction, not as MLE.

Carrier-weighted MSE \(\sum n_j(p_j-\hat\mu_{\mathrm{w}})^2/\sum n_j\) is **not** \(\mathrm{Var}(\pi)\). It is a size-biased second moment. It matches a model in which a typical **carrier** is drawn, not a typical **variant**. That would silently re-define a variant-level Beta prior.

**Recommendation (primary, coherent with variant-level donors).**

1. Restrict to the eligible **gene–disease** matrix.
2. \(\hat\mu=\hat\mu_{\mathrm{unw}}\).
3. Invert with **MSE**, using \(\hat\sigma^2_{\mathrm{proc}}\) if a non-trivial fraction of \(n_j\) is small; if all \(n_j\) are large, raw unweighted MSE is acceptable and should be labeled as “noise-uncorrected.”
4. Build every donor \(\tilde p_j\) from that single \((\alpha_{\mathrm{emp}},\beta_{\mathrm{emp}})\).

**Sensitivity (labeled).** Carrier-weighted \(\hat\mu_{\mathrm{w}}\) and/or carrier-weighted MSE; full beta-binomial MLE of \((\alpha,\beta)\). Do not replace the primary MOM run.

### 2.2 MAE is not a variance

**Fact.** MAE \(=J^{-1}\sum_j|p_j-\hat\mu|\) has units of probability. Variance has units of probability squared. For any centered random variable, \(\mathrm{MAE}\le\mathrm{RMSE}\), with equality only in degenerate two-point cases of a specific form. MAE does **not** determine \(\sigma^2\), hence does **not** determine \(s\).

What MAE **can** do:

- Describe typical absolute deviation.
- Bound RMSE from below.
- After a Beta is fit from **MSE**, compare empirical MAE to the Beta’s theoretical MAE \(\int_0^1|x-\mu|f_{\alpha,\beta}(x)\,dx\) (incomplete-beta split at \(\mu\)) as a **shape diagnostic**. Large disagreement means the Beta is a bad model (mass at 0/1, multimodality, extreme skew).

What MAE **cannot** do:

- Be plugged into \(s=\mu(1-\mu)/\sigma^2-1\) in place of \(\sigma^2\) or in place of \(\sigma\).
- Separate process variance from binomial noise.
- Uniquely identify \((\alpha,\beta)\) without an extra family assumption.

Gaussian bridge (optional sensitivity only): \(\mathrm{MAE}\approx\sigma\sqrt{2/\pi}\) implies \(\sigma\approx\mathrm{MAE}\sqrt{\pi/2}\approx 1.253\,\mathrm{MAE}\). Beta is bounded and often skewed; this is not the primary inversion.

**Recommendation.** Primary inversion uses **MSE** (noise-corrected as above). Report MAE as a companion statistic and as a posterior predictive check. If legacy notes say “mean and MAE/MSE,” treat MSE as the moment and MAE as diagnostics, and write that down in the decision log.

### 2.3 Valid domain and honest boundaries

**Fact.** For \(\mu\in(0,1)\), a Beta variance cannot exceed the Bernoulli bound \(\mu(1-\mu)\). Equivalently \(s>0\) iff \(0<\sigma^2<\mu(1-\mu)\). Also \(\mu\in\{0,1\}\) is not interior: it would require \(\alpha=0\) or \(\beta=0\), which is not a proper density on \((0,1)\).

**Honest handling (recommendation).** Do **not** silently clip \(\hat\sigma^2\) or \(s\).

- If \(\hat\mu\in\{0,1\}\): stop and report “empirical mean on the boundary; Beta prior not identified from moments.” Surface the variant mix (all \(p_j=0\), or gnomAD-dominated zeros).
- If \(\hat\sigma^2_{\mathrm{proc}}\le 0\): observed scatter is \(\le\) binomial noise. Report “no extra-binomial heterogeneity.” Document a **named** fallback, e.g. \(s\to\infty\) (common-\(p\) model, \(\tilde p_j\approx\hat\mu\)) or a weakly informative floor such as \(\mathrm{Beta}(1,1)\) or \(\mathrm{Beta}(1/2,1/2)\), with a flag on every downstream density.
- If \(\hat\sigma^2_{\mathrm{proc}}\ge\hat\mu(1-\hat\mu)\): implied \(s\le 0\). Report “dispersion at or beyond Bernoulli; MOM Beta impossible.” Named fallback: \(s\downarrow s_{\min}\) (e.g. \(s_{\min}=2\), Uniform) **with a warning**, or switch that gene–disease to MLE/sensitivity, but do not pretend MOM succeeded.

Write \((\hat\mu,\widehat{\mathrm{MSE}},\mathrm{MAE},\hat\sigma^2_{\mathrm{proc}},s,\alpha,\beta,\texttt{boundary_flag})\) to an artifact. That is part of the analysis, not a failure.

### 2.4 Full dataset and gnomAD

**Recommendation.** “Full dataset” = all eligible variants in **this gene–disease** analysis set, including the target, for the single shared prior. Do not pool unrelated diseases into one Beta.

gnomAD-as-unaffected is taken as specified: it increases \(U_j\), decreases \(p_j\) and \(\tilde p_j\), more so for variants with non-trivial population frequency. Consequence (not a veto): the empirical prior and all donor posteriors are shifted toward lower penetrance relative to a clinic-only \(U_j\). Record that \(U_j\) definition in the donor table.

Do not substitute AlphaMissense scores for \(\tilde p_j\), and do not train the density on raw \(p_j\) in the primary run.

---

## 3. Sine/sigmoid distance decay with half-weight near 3 Å

Need \(K(0)=1\), \(K(3\,\text{Å})=1/2\), rapid decay, optional compact support.

### 3.1 Explicit kernels

**Primary compact raised-cosine (sine/cosine family).**

\[
K_{\cos}(d)=\begin{cases}
\cos^2\!\left(\dfrac{\pi d}{12\,\text{Å}}\right)=\dfrac{1+\cos(\pi d/6\,\text{Å})}{2}, & 0\le d\le 6\,\text{Å},\\
0, & d>6\,\text{Å}.
\end{cases}
\]

Checks: \(K(0)=1\), \(K(3\,\text{Å})=1/2\), \(K(6\,\text{Å})=0\), \(K'(0)=0\) (smooth at the origin), compact first-shell support. This is the Hann / raised-cosine window.

**Alternate normalized logistic (non-compact).**

\[
K_{\mathrm{sig}}(d)=\frac{2}{1+\exp(d\ln 3/3\,\text{Å})}.
\]

Checks: \(K(0)=1\), \(K(3\,\text{Å})=1/2\), \(K(d)\to 0\) as \(d\to\infty\), never exactly 0, \(K'(0)\neq 0\). Faster or slower decay is a single rate parameter \(k=\ln 3/d_{\mathrm{half}}\).

**Optional even sigmoid with zero slope at 0:** \(K(d)=\mathrm{sech}^2(d/\lambda)\) with \(\lambda=3/\mathrm{acosh}(\sqrt{2})\approx 3.404\,\text{Å}\) also has \(K(0)=1\), \(K(3\,\text{Å})=1/2\). Sensitivity only.

Do not use \(\mathrm{sinc}\) (oscillatory negative lobes). Negative kernel weight would mix high- and low-penetrance donors with the wrong sign.

### 3.2 Cα vs centroid vs heavy-atom: why 3 Å is not portable

**Fact.** For a typical contacting pair,

\[
d_{\min\text{-heavy}} \;<\; d_{\text{side-chain centroid}} \;<\; d_{\mathrm{C}\alpha}.
\]

A 3 Å **minimum heavy-atom** distance is a real steric contact. A 3 Å **Cα–Cα** distance is inside backbone excluded volume (Cα–Cα of consecutive residues is already \(\sim 3.8\,\text{Å}\); contacting non-neighbors are usually 5–10 Å Cα–Cα). Therefore \(K_{\cos}\) with \(d_{\mathrm{half}}=3\,\text{Å}\) on Cα–Cα is nearly a self-weight: almost no legitimate neighbor receives \(K>1/2\).

**Recommendation.** State the metric in the primary run. For missense side-chain effects, **minimum heavy-atom distance** (fallback: side-chain centroid) is the metric for which “half-weight at 3 Å” is physically on-scale. Use Cα only when side-chain atoms are missing, and **recalibrate** \(d_{\mathrm{half}}\) (do not keep 3 Å blindly).

### 3.3 Prespecified sensitivity grid (small, not a search)

| Axis | Primary | Grid |
|---|---|---|
| Family | \(K_{\cos}\), compact at \(2 d_{\mathrm{half}}\) | \(K_{\cos}\) vs \(K_{\mathrm{sig}}\) |
| Metric | min heavy-atom | + centroid; + Cα |
| \(d_{\mathrm{half}}\) | 3 Å (heavy-atom/centroid) | \{2, 3, 4\} Å |
| Cα-only \(d_{\mathrm{half}}\) | 6 Å if Cα is forced | \{4, 6, 8\} Å |

Lock this grid before looking at gene-level summaries. Primary remains \(K_{\cos}\), 3 Å, heavy-atom.

---

## 4. Disordered regions, polymer \(\sqrt{N}\), biological unit, low-confidence coordinates

### 4.1 What \(N\) is

**Recommendation.** Polymer distances apply only to a pair of residues that:

- lie on the **same polypeptide chain**,
- both sit inside the **same contiguous disordered segment**,
- are numbered on the **canonical biological sequence** of that chain.

Then \(N=|i-j|\) in canonical residue indices (contour separation inside that coil).

**Not \(N\):** separation that jumps a folded domain; separation that jumps to another chain; author/PDB numbering with insertions; UniProt index on chain B minus index on chain A.

A structured spacer between two loops means those loops are **not** one ideal chain. Using \(\sqrt{|i-j|}\) across a domain invents a coil distance that does not exist. Using Euclidean distances from low-pLDDT loop coordinates invents 3D contacts that do not exist. Mixed ordered–disordered pairs: **no polymer shortcut and no low-confidence Euclidean edge**; contribution 0 unless both endpoints pass the coordinate-quality gate.

### 4.2 Scale of \(d=\lambda\sqrt{N}\)

User request is distance **proportional** to \(\sqrt{N}\) (ideal-chain / Gaussian-chain, Flory \(\nu=1/2\)). Units must match the Euclidean kernel.

Physically, RMS Cα end-to-end for an unfolded peptide is on the order of \(\lambda\sqrt{N}\) with \(\lambda\sim 5\)–\(6\,\text{Å}\) (virtual-bond \(3.8\,\text{Å}\) times \(\sqrt{C_\infty}\)). Under the primary \(K_{\cos}\) that dies at 6 Å, a physically scaled polymer channel is **very local** (only tiny \(N\) gets non-zero weight). That is a feature of combining a 3 Å contact kernel with coil statistics, not a bug.

**Recommendation.** Primary: \(d_{ij}=\lambda\sqrt{|i-j|}\) with a **stated** \(\lambda\) (e.g. \(5.5\,\text{Å}\)), same \(K\) as Euclidean. If that channel is effectively off for typical loop lengths, that goes in the support summary (`polymer_N_eff`), not an unannounced retune.

**Sensitivity.** \(\lambda\) such that \(K=1/2\) at \(N\in\{4,9,16\}\); optional \(\nu=0.588\) (good solvent) as a second-order check. Do not mix Euclidean \(d\) and dimensionless \(\sqrt{N}\) in one sum without converting units.

### 4.3 Do not invent 3D contacts from low-confidence coordinates

Gate **donor and target** atoms:

- Drop Euclidean edges if either residue is below the pLDDT (or equivalent) cutoff, or is in an unmodeled gap.
- If the **target** site fails the gate: no Euclidean density; polymer-within-segment only if the target is in a defined coil; otherwise \(\rho\) is undefined and reported as missing, not as 0.
- If a **donor** fails the gate: that donor is ineligible for Euclidean \(K\); it may still donate via polymer if both ends are in the same coil and numbering is canonical.
- Never impute a fake Cα for a missing loop and then apply \(K_{\cos}\).

Prespecify the cutoff (common choices: pLDDT 70 “confident,” or 50 “not very low”). Primary: 70 for Euclidean eligibility.

### 4.4 Biological functional unit and duplicate carriers

Use the **biological assembly** that matches function (true monomer \(\Rightarrow\) monomer pilot is correct; homodimer \(\Rightarrow\) assembly, not the asymmetric monomer). Numbering: map every chain to the canonical sequence; do not mix author residue numbers.

**Fact.** Two chains in a homodimer carrying the same heterozygous variant are **one statistical carrier**, two spatial images.

**Recommendation.**

- **Counts \(A_j,U_j\):** one row per variant key, never multiplied by chain count.
- **When the variant is the target:** exclude **all** images of that key (Section 1).
- **When it is a donor:** \(\tilde p_j\) enters **once**. If several images are geometrically close to the query site, collapse kernels with a stated rule, e.g. \(K_j=\max_c K(d_c)\) (primary, conservative against double-counting) or \(1-\prod_c(1-K(d_c))\) (union of “contact” events). Do **not** add \(\tilde p_j\) twice.

---

## Donor weighting, support, uncertainty

The spatial estimator at a query site \(x\) (target CA/centroid/heavy-atom representative) is a kernel-weighted mean of **posterior means**:

\[
\rho(x)=\frac{\sum_j w_j K(d_{xj})\,\tilde p_j}{\sum_j w_j K(d_{xj})}.
\]

**Primary \(w_j=1\)** (equal variants). That matches a variant-level Beta and the unweighted moment construction. It does **not** silently turn the density into a carrier-prevalence map.

**Sensitivities (label, do not swap into primary):**

- Carrier weights \(w_j=n_j\).
- Posterior precision \(w_j=\alpha_{\mathrm{emp}}+\beta_{\mathrm{emp}}+n_j\) (or \(1/\mathrm{Var}(\tilde p_j)\)). This is the statistically natural alternate if unequal \(n_j\) are severe.

If the denominator is 0, \(\rho\) is missing, not \(\hat\mu\). Optionally report a **separate** shrinkage display \(\rho_{\mathrm{shrunk}}=(s_0\hat\mu+\sum wK\tilde p)/(s_0+\sum wK)\) with prespecified \(s_0\), as a sensitivity.

**Support summaries (per target, required artifacts).**

- \(N_{\mathrm{eff}}=\sum K_j\) and Kish \(N_{\mathrm{eff}}^*=(\sum K_j)^2/\sum K_j^2\)
- \(C_{\mathrm{eff}}=\sum K_j n_j\)
- Distance to nearest **retained** donor; count of same-residue other alleles used
- Euclidean vs polymer split of \(\sum K\)
- Coordinate-quality flags (target and donors)
- Assembly multiplicity and the kernel-collapse rule used

**Uncertainty (conditional on fixed \((\alpha_{\mathrm{emp}},\beta_{\mathrm{emp}})\)).** Donor posteriors are independent given the prior, with

\[
\mathrm{Var}(\tilde p_j)=\frac{\alpha_j'\beta_j'}{(\alpha_j'+\beta_j')^2(\alpha_j'+\beta_j'+1)},\quad
\alpha_j'=\alpha_{\mathrm{emp}}+A_j,\;\beta_j'=\beta_{\mathrm{emp}}+U_j.
\]

Delta-method:

\[
\widehat{\mathrm{Var}}(\rho)=\frac{\sum_j (w_j K_j)^2\mathrm{Var}(\tilde p_j)}{\bigl(\sum_j w_j K_j\bigr)^2}.
\]

This ignores kernel-distance error and hyperparameter uncertainty (the latter is the \(O(1/J)\) indirect term). Report \(\rho\pm 1.96\sqrt{\widehat{\mathrm{Var}}}\) as **conditional** uncertainty, plus \(N_{\mathrm{eff}}\) so empty neighborhoods are not over-interpreted. Do not treat \(\rho\) as a known penetrance.

---

## Fair AlphaMissense / density / combined comparison

**Fact.** Variant-LOO density is a **descriptive internal** smoother of the same gene–disease labels. It is not independent external validation. Nearby variants share ascertainment, phenotype definitions, and the common empirical prior. AlphaMissense does not use this cohort’s \(A_j,U_j\), but is also not a clean “external truth,” and must not replace the empirical Beta donors.

**Fair protocol (recommendation).**

1. Scores compared on the same target set: AM score, LOO density \(\rho_{-j}\), and a **combined** score.
2. Labels: observed \(p_j\) or the binomial likelihood \((A_j,n_j)\), **not** \(\tilde p_j\) as both feature and label.
3. Combined model: prespecified (e.g. logistic or beta-binomial regression of \((A,n)\) on \(\rho_{-j}\) and AM) with coefficients fit **without** the target (simple outer variant-LOO, or a locked 50/50 weight as a non-fitted baseline).
4. Metrics: (i) association of scores with \(p_j\) (Pearson/Spearman, and \(n_j\)-weighted); (ii) beta-binomial or binomial calibration (reliability slope, intercept); (iii) **stratify by \(N_{\mathrm{eff}}\)** — density is only expected to help when local support exists; (iv) same-residue-other-allele present vs absent.
5. Language in outputs: “internal LOO descriptive association,” never “independent validation” or “held-out cohort performance.”

AM remains a baseline head-to-head, not a prior on \(\pi_j\).

---

## Concrete tests and artifacts

**Tests (fail the build if violated).**

1. Target key absent from donor list; at least one other substitution at the same canonical residue, if it exists in the set, remains.
2. Homomer fixture: two chains, one variant key \(\Rightarrow n_j\) not doubled; LOO of that variant zeros **both** images.
3. Synthetic \(\pi_j\sim\mathrm{Beta}(2,8)\), large equal \(n_j\): recovered \((\hat\alpha,\hat\beta)\) within MOM tolerance.
4. Same synthetic, small heterogeneous \(n_j\): noise-uncorrected \(s\) smaller than noise-corrected \(s\); both written out.
5. Boundary fixture: \(\hat\mu=0.3\), \(\mathrm{MSE}=0.25>\mu(1-\mu)\Rightarrow\) `boundary_flag`, no silent clip.
6. Kernel unit test: \(K(0)=1\), \(K(3\,\text{Å})=0.5\), \(K(6\,\text{Å})=0\) for \(K_{\cos}\).
7. Metric fixture: contacting pair has \(d_{\mathrm{heavy}}<d_{\mathrm{centroid}}<d_{\mathrm{C}\alpha}\).
8. Coil fixture: same disordered segment \(\Rightarrow\) polymer edge; pair spanning a helix \(\Rightarrow\) no polymer edge; different chains \(\Rightarrow\) no polymer edge.
9. Low pLDDT fixture: no Euclidean edge; target-low \(\Rightarrow\) \(\rho\) missing unless polymer path exists.
10. Self-leak fixture: including the target as donor moves \(\rho\) at the target site; a distant variant does not, within tolerance.
11. Numbering fixture: canonical \(|i-j|=1\) maps to a peptide neighbor, not an author-number gap.

**Artifacts.**

- Hyperparameter row: gene–disease, \(J\), \(\hat\mu\), MSE, MAE, \(\hat\sigma^2_{\mathrm{proc}}\), \(s\), \(\alpha_{\mathrm{emp}}\), \(\beta_{\mathrm{emp}}\), `boundary_flag`, \(U_j\) rule.
- Donor table: variant key, \(A,U,n,p,\tilde p\), chain(s), canonical index, pLDDT, ordered/disordered segment id, assembly copy count, LOO-excluded flag.
- Kernel curves for the locked grid.
- Per-target output: \(\rho_{-j}\), \(\widehat{\mathrm{se}}\), \(N_{\mathrm{eff}}\), \(N_{\mathrm{eff}}^*\), \(C_{\mathrm{eff}}\), nearest donor, Euclidean/polymer split, AM, combined, \(N_{\mathrm{eff}}\) stratum.
- Decision log: primary vs each labeled sensitivity.

---

## Ordered implementation plan (primary path)

1. **Define the eligible gene–disease table.** Canonical variant keys, \(A_j\), \(U_j\) with gnomAD as unaffected, \(n_j>0\). Map each key to residue(s) on the **biological unit** (monomer if the functional unit is a monomer). Record chain copies without multiplying counts.

2. **Fit one empirical Beta on that full table.** Unweighted \(\hat\mu\), MSE with documented noise correction, invert to \((\alpha_{\mathrm{emp}},\beta_{\mathrm{emp}})\). Apply the boundary protocol; write the hyperparameter artifact. Do not refit per target in the primary run.

3. **Compute donor posteriors** \(\tilde p_j=(\alpha_{\mathrm{emp}}+A_j)/(\alpha_{\mathrm{emp}}+\beta_{\mathrm{emp}}+n_j)\) for every eligible variant. These, not AM and not raw \(p_j\), are density donors.

4. **Coordinate module.** Biological-unit structure; canonical numbering; pLDDT/quality mask; DSSP or pLDDT run-length disordered segments. Distances: primary min heavy-atom (fallback centroid); polymer \(d=\lambda\sqrt{|i-j|}\) only inside a common disordered segment on the same chain.

5. **Kernel module.** Primary \(K_{\cos}\) with \(d_{\mathrm{half}}=3\,\text{Å}\), compact at 6 Å. Unit tests in §3.1 and the test list.

6. **LOO density.** For each target key: drop **all spatial copies** of that key; keep other alleles at the same residue; collapse remaining multi-chain images with \(\max_c K\); equal-variant Nadaraya–Watson mean of \(\tilde p_j\); write support and conditional SE. Missing denominator \(\Rightarrow\) missing \(\rho\).

7. **Locked sensitivity grid only after the primary is frozen:** \(d_{\mathrm{half}\) and metric; \(K_{\mathrm{sig}}\); \(\lambda\) / \(N_{\mathrm{half}}\); precision and carrier \(w_j\); optional hyperparameter refit without the target.

8. **Comparators.** AM vs \(\rho_{-j}\) vs combined under the fair protocol; stratify by \(N_{\mathrm{eff}}\); label as internal LOO association.

9. **Ship artifacts and the decision log.** Primary numbers come only from steps 1–6.

---

## Decisions still needing a named choice

These do not block starting the primary path; they need an explicit line in the decision log rather than an implicit default.

1. **Noise correction on MSE:** on (recommended when \(n_j\) vary) vs raw unweighted MSE (acceptable if all \(n_j\) large).
2. **Boundary fallback** when \(s\le 0\) or \(\hat\sigma^2_{\mathrm{proc}}\le 0\): common-\(p\) (\(s\to\infty\)) vs named weak Beta vs \(s_{\min}=2\).
3. **Distance metric** if heavy atoms are incomplete: centroid vs recalibrated Cα (not 3 Å Cα).
4. **pLDDT cutoff** for Euclidean eligibility (70 vs 50) and disordered-segment definition (pLDDT run vs DSSP vs external IDP call).
5. **Polymer prefactor \(\lambda\)** (physical \(\sim 5.5\,\text{Å}\) vs \(N_{\mathrm{half}}\) calibration).
6. **Multi-image kernel collapse:** \(\max_c K\) (recommended) vs noisy-OR vs mean.
7. **Precision-weighted \(w_j\)** as the first labeled sensitivity, with equal-variant remaining primary.
8. **Combined AM rule:** locked 50/50 vs beta-binomial regression under variant-LOO of the combiner only.

Mathematical facts in this review are the Beta/Bernoulli domain, the split between process variance and binomial noise, MAE vs MSE, and the metric-dependence of a 3 Å half-distance. Everything else is either the primary implementation of the four requested choices or a labeled sensitivity.
