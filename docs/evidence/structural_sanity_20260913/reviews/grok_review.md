# Grok independent methodological review

Requested model: grok-4.6, high reasoning. Actual model and usage are retained in grok_raw.json. Reviewer saw the saved aggregate brief only and did not inspect source files.

**Scientific-method review (aggregate math and design only; not independent verification of source data).**

**(1) Verified-from-provided-math conclusions and uncertainties**

The GCK missense prior is internally consistent: \(\mu=\alpha/(\alpha+\beta)=1.081159/2.918755=.370418\), \(S=\alpha+\beta=2.918755\). Posterior means match the stated Beta update: \(A=0,U=1\) gives \(1.081159/3.918755=.275894\); \(A=1,U=0\) gives \(2.081159/3.918755\approx.5311\). \(w(n=1)=1-1/1.01=.009901\) is correct. Kernel \(k(d)=2/(1+e^{\ln 3\cdot d/h})\) is 1 at \(d=0\) and 0.5 at \(d=h=3\); no cutoff is a fixed choice, not an error.

The GCK median-density jump \(.182394\to.374072\) with **byte-identical spatial weights** is a prior-floor effect: the new median sits on \(\mu\approx.370\). That is expected when many donors are weakly counted (385/634 population-only; \(A=0,U=1\) floor \(.276\)) and the spatial feature is a normalized weighted mean of those posteriors. It is **not** evidence that nearly all residues cause MODY, nor a claim about substitution completeness. Singleton weight \(\approx1.53\%\) despite a large unit share is the intended \(w(n)=1-1/(n+.01)\) behavior; I cannot confirm “52.4% of units” from the 246 population-only \(n=1\) figure alone (\(246/634\approx38.8\%\)), so treat 52.4% as an all-singleton claim that still needs a unit-level \(n\) audit.

Hold-out MAE deltas are negligible on the stated labels: GCK AM-only \(.124636\) vs AM+density \(.123400\) (\(\Delta\approx.0012\)); clinical 242 almost unchanged (\(.111226\to.111203\)). BRCA2 has only 44 AM-common targets, all literature-archive fallbacks, tightly clustered posteriors, and no density gain. These metrics reuse posterior means from the **same** count collection; they do not show independent disease-outcome calibration. Historical \(v=\sum w(y-\mu)^2/M\) (not \(\sum w\)) inflates \(S\) relative to a weighted-variance estimator; keep it as the named primary, with a diagnostic \(v_w=\sum w(y-\mu)^2/\sum w\).

BRCA2 “no polymer” is a **manifest gap**, not a disorder-negative finding: 127/6656 missense structurally supported, 72/3418 residues with COM, two disjoint peptide frames, zero polymer by construction. Missing coordinates \(\neq\) IDR.

**(2) Safe BRCA2 hybrid geometry**

Keep the two experimental frames as **disjoint bound-3D contexts** (no cross-frame Euclidean distances; peptide copies along a filament are assembly copies, not extra variants). Add a separate **canonical-IDR polymer layer** only for explicitly validated contiguous IDR segments with provenance (source, chain, start–end on the canonical sequence). Intra-segment distances: \(3.8\sqrt{|i-j|}\) only inside one validated segment on one molecular chain. Mixed ordered/IDR and cross-chain IDR pairs remain undefined (not zero, not Euclidean, not polymer).

Prefer **explicitly labeled alternatives** (bound-3D vs free-polymer) over a fused coordinate object: each supported context normalizes its own donors, then contexts average equally—the existing rule. Do not invent a full-length assembly. Do not treat gap residues as IDR. Do not emit density 0 for unsupported sequence.

Tests: (a) no pair with coordinates from different frames; (b) each canonical IDR interval counted once, not once per PDB copy; (c) missing-coord residues stay “unsupported,” never auto-IDR; (d) polymer pairs fail if segment IDs differ or chain IDs differ; (e) bound-peptide fragments do not inherit partner-filament geometry; (f) provenance fields survive into the density audit.

**(3) GCK diagnosis**

- **Prior floor/weighting:** Median density \(\approx\mu\) after prior swap, with unchanged kernels, is global shrinkage, not a residue-wise MODY map. High-n clinical units dominate \(\mu\); singletons barely move the prior.
- **Sparse counts:** \(A=670\) vs literature \(U=72\) vs gnomAD \(U=4369\). Most spatial averages are of near-prior, \(A=0\) posteriors plus a few high-\(A\) donors. Neighborhood maps will look “hot” wherever any donors exist.
- **Endpoint mix:** Pooling MODY/hyperglycemia with activating/hypoglycemia makes \(A\) a mixed-mechanism count. Local means are not MODY risk.
- **Local signal:** Tiny MAE deltas, especially on the 242 clinical targets, argue against a strong extra spatial component on these labels. Variant-only LOO must stay; residue-collapsing would further inflate apparent spatial support.

**(4) Ranked immediate checks**

1. BRCA2 provenance audit: frame IDs, missing≠IDR, no cross-frame distances, IDR intervals unique on canonical sequence, unsupported \(\neq 0\).
2. GCK diagnostic split of \(A\) by endpoint (MODY vs activating) **without** changing the primary pooled prior; report density vs \(\mu\) residuals, not raw histograms.
3. Confirm exact-variant LOO across aliases/copies for every training target (variant-only, not residue-only).
4. Named prior diagnostic: \(v\) vs \(v_w\), singleton share of \(M\) vs of \(\sum w\), and density median minus \(\mu\).
5. Hold out labels that are **not** the same-collection posterior means (or report that calibration remains circular).
