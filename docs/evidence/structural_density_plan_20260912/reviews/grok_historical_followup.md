**Methodological review (recovered public historical pipeline).** Scope: empirical Beta prior, variant-level LOO, compact sine/sigmoid spatial kernel (midpoint \(a\approx 3\)Å), and \(3.8\sqrt{N}\)Å polymer in disorder. gnomAD is outside these assumptions. Historical conventions are the primary baseline; this note flags mathematics without substituting them.

### 1. Empirical prior (weighted intercept-only MOM)
Let \(y_i=A_i/n_i\), \(w_i=1-1/(n_i+0.01)\), \(\mu=\sum w_i y_i/\sum w_i\). Primary \(y\) keeps exact \(0/1\); historical \(.0005/.9995\) is a sensitivity only. The intended second moment is
\[
\mathrm{MSE}=\frac{1}{M}\sum_i w_i(y_i-\mu)^2,
\]
\(M=\#\) count-bearing variants, **not** \(\sum w_i\). Then \(\kappa=\mu(1-\mu)/\mathrm{MSE}-1\), \(\alpha=\mu\kappa\), \(\beta=(1-\mu)\kappa\), posterior \(\mathrm{Beta}(\alpha+A,\beta+U)\). Fit once on all count-bearing classes for the gene–disease, **before** donor filtering.

**Matching vs estimation.** \(\mathrm{MSE}\) is a **variance-matching convention** for the Beta identity \(\mathrm{Var}(\theta)=\mu(1-\mu)/(\kappa+1)\). It is **not** a latent-variance estimator of \(\mathrm{Var}(\theta)\): \(y\) still carries binomial noise, so \(\mathbb{E}[(y-\mu)^2]=\mathrm{Var}(\theta)+\mathbb{E}[\theta(1-\theta)/n]\). Plugging raw \(\mathrm{MSE}\) into \(\kappa\) therefore attributes sampling noise to the prior (typically **underconcentrated** \(\kappa\)). It is also not \(\widehat{\mathrm{Var}}(\mu)\); a linear weighted mean has sandwich variance \(\sum c_i^2\sigma_i^2\) with \(c_i=w_i/\sum w\), which is not used here.

**Normalizer mismatch.** Weighted mean uses \(\sum w\); \(\mathrm{MSE}\) uses \(M\). Normalized weighted MSE \(\sum w(y-\mu)^2/\sum w\) is larger by \(M/\sum w\ge 1\) and should be **reported separately**, not silently swapped in. Dividing by \(M\) **shrinks** the matching target and **raises** \(\kappa\) relative to standard weighted MOM. Saturating \(w\) already nearly zeros \(n=1\) points (\(w\approx 0.01\)); effective sample size is \(\sum w\), often \(\ll M\). Do not replace saturating weights, \(\mathrm{MSE}/M\), or add noise subtraction in the primary fit.

**Other issues (do not change the baseline).** \(\kappa<0\) when \(\mathrm{MSE}\ge\mu(1-\mu)\) (Bernoulli bound); the map \(\mathrm{MSE}\mapsto\kappa\) is non-Lipschitz near that edge—floor/handle it. \(\kappa=2\) is \(\mathrm{Beta}(1,1)\) **only if** \(\mu=1/2\); otherwise \(\mathrm{Beta}(2\mu,2(1-\mu))\) is not uniform.

### 2. Variant-only LOO
Historical rule: drop only `var==target`; other substitutions at the same residue remain. Average **posterior means with equal variant weights**. That is the baseline.

Shared full-data \((\mu,\kappa)\) is **not** guaranteed \(O(1/M)\) leakage. Influence on \(\mu\) is \(w_i/\sum w\), which is \(O(1)\) if one well-counted variant dominates. Influence on \(\mathrm{MSE}\) is quadratic and worse for boundary \(y\in\{0,1\}\) with large \(w\). When \(\mathrm{MSE}\) is near \(\mu(1-\mu)\), one point can swing \(\kappa\) from large to \(\approx 0\). Indirect path: the target shifts \((\alpha,\beta)\), which shifts **every** donor posterior mean that later enters the spatial average. Production prior on all labeled classes is coherent; for **evaluative** LOO it is leakage. Refit-without-target is an evaluation correction, not a silent production change.

### 3. Kernel
Compact sine: \(K=1\) on \(d\le a-\pi\), \(\tfrac12-\tfrac12\sin((d-a)/2)\) on \((a-\pi,a+\pi)\), else \(0\). \(C^1\) at the knots. With \(a=3\): \(K(3)=1/2\), \(K(0)\approx 0.998747\), cutoff \(a+\pi\approx 6.14159\)Å. Because \(a-\pi<0\), there is **no physical plateau**; decay starts at contact. Historical \(3.14\) vs \(\pi\) is a declared tiny correction. Sigmoid sensitivity, unit at \(0\) and half at \(h\): \(2/(1+e^{\log 3\cdot d/h})\). Do not replace the compact sine as primary.

### 4. Geometry, copies, polymer
Documented “centroid” has **unrecovered atom selection**. New **explicit primary**: side-chain heavy-atom mass-weighted center, Gly C\(\alpha\) fallback; C\(\alpha\)–C\(\alpha\) a named sensitivity. This is not claimed as the historical producer. Do not make min-heavy-atom primary. Side-chain vs C\(\alpha\) distances have **no universal order**; at a 6Å cutoff the choice is first-order.

Variant ID survives chain expansion; strip the target from **all** copies. Collapse donor copies by **max kernel** (nearest eligible copy), then **average densities** across equivalent target contexts. Reuse **one** donor posterior draw across copies (no duplicated clinical counts). Strict same-chain / same-IDR polymer; no mixed ordered–IDR or cross-chain distances. Missing experimental residue is not automatically disorder. Monomer-first: explicit isoform/canonical map, then numbered oligomers.

Polymer \(3.8\sqrt{N}\)Å with the \(a=3\) compact kernel admits only \(N\le 2\) (\(N=2\Rightarrow 5.37\)Å; \(N=3\Rightarrow 6.58\)Å \(>\) cutoff). Preserve this; **flag no-donor** cases rather than broadening \(a\), \(K\), or the polymer.

### Conditional variance of the spatial estimate
Given \((\alpha,\beta)\) and independent donor Betas, the linear combination \(\hat Y=\sum c_i\theta_i\) (\(\sum c_i=1\)) has
\[
\mathbb{E}[\hat Y]=\sum c_i m_i,\qquad \mathrm{Var}(\hat Y\mid\alpha,\beta,\mathrm{data})=\sum_i c_i^2 v_i,
\]
\(m_i=\mathbb{E}[\theta_i]\), \(v_i=m_i(1-m_i)/(\alpha+\beta+n_i+1)\). Equal-weight averaging of **means** is the plug-in \(\sum c_i m_i\); that number is fixed given data and is **not** \(\mathrm{Var}(\hat Y)\). A mixture of densities (oligomer copies) has the same mean but extra between-copy variance. \(\sum c_i^2 v_i\) omits: (i) dependence through shared \((\hat\mu,\hat\kappa)\), (ii) hyperparameter uncertainty, (iii) cohort/ascertainment, (iv) geometry/atom-selection/IDR calls. Empirical-Bayes plug-in therefore understates predictive uncertainty.

### Required vs optional
**Required for a faithful, well-posed run:** \(\kappa\) floor; report \(\mathrm{MSE}/M\) and normalized \(\sum w(y-\mu)^2/\sum w\); exact \(0/1\) \(y\); saturating \(w\) and \(\mathrm{MSE}/M\) as primary; variant-ID LOO including all copies; max-kernel donor collapse + one draw; equal-weight (polymer) / kernel-weight (structure) as specified; flag zero donors; same-chain/same-IDR only; isoform map before oligomers; do not label \(\kappa=2\) uniform unless \(\mu=1/2\); declare SC-centroid as new primary, not historical identity.

**Optional sensitivities:** \(.0005/.9995\); sigmoid \(K\); C\(\alpha\) distances; \(\pi\) vs \(3.14\); noise-subtracted latent \(V\); unweighted moments; prior refit excluding the target (evaluation); broader polymer/kernel.

### Execution checklist
1. Collect count-bearing classes; \(y=A/n\) exact; \(w=1-1/(n+0.01)\).
2. \(\mu=\sum wy/\sum w\); \(\mathrm{MSE}=\sum w(y-\mu)^2/M\); store normalized weighted MSE.
3. \(\kappa=\mu(1-\mu)/\mathrm{MSE}-1\) with a declared floor; \(\alpha,\beta=\mu\kappa,(1-\mu)\kappa\).
4. Posteriors \(\mathrm{Beta}(\alpha+A,\beta+U)\) for all labeled variants.
5. LOO: drop target ID on every chain copy; keep other residue substitutions.
6. Distances: SC heavy mass center (Gly C\(\alpha\)); C\(\alpha\) run as sensitivity.
7. Structure: compact sine \(a=3\) (optional sigmoid); collapse donors by max \(K\).
8. Disorder: \(3.8\sqrt{N}\) only on the same IDR/chain segment; \(N\le 2\) can contribute; else flag no donor.
9. Average equally (polymer) or by \(K\) (structure); average copy-densities; reuse one donor draw.
10. Quote \(\sum c_i m_i\) as the point estimate; if reporting variance, use \(\sum c_i^2 v_i\) **and** state it is conditional on plug-in \((\alpha,\beta)\), not full uncertainty.
