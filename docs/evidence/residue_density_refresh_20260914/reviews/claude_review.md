## Review: generic residue‑neighborhood plotting framework

**Standing caveat.** I have inspected no data, code, or fitted prior. Everything below is derived from the brief and its synthetic arithmetic. Nothing here validates any implementation; the checks are written so they can *fail*.

### Overall

The design is defensible **as a descriptive feature** and is explicit about that. The main risks are not conceptual but (a) arithmetic edge cases in the historical moment‑matching formula, (b) shrinkage‑induced flatness being misread as a biological or geometric result, and (c) label circularity if catalogue classifications rather than patient counts feed A. The framework's own disclaimers ("no independent clinical calibration by construction") should be treated as load‑bearing, not boilerplate.

### 1. Fitting formula, orientation, exclusion, aggregation — errors vs conventions

Separate three categories: **arithmetic bugs**, **stated conventions** (weighted‑mean/MSE, prior held fixed, target‑excluded, population‑as‑unaffected), and **untested assumptions**.

Bugs to trap (P0):
- **w = 1 − 1/(n+0.01) is negative for n = 0** (≈ −99) and ≈ 0.5 at n = 1. Any eligible unit with n = 0 corrupts mu and v silently, possibly flipping signs. *Test:* assert min(w) > 0 and n ≥ 1 for every eligible unit; inject a synthetic n = 0 row and confirm it is rejected, not absorbed.
- **kappa = mu(1−mu)/v − 1 can be ≤ 0** whenever v ≥ mu(1−mu)/2, producing non‑positive alpha/beta. Because v is an unnormalized MSE divided by M while mu is a *weighted* mean, the two are on inconsistent scales — that is the declared convention, but it makes the degenerate branch reachable. *Test:* assert alpha > 0, beta > 0; construct a synthetic high‑overdispersion set (many n = 1 units, mixed 0/1) and record whether kappa goes negative. Report the guard, do not quietly clip.
- **Orientation.** *Test:* build synthetic data with 90% affected; assert alpha/(alpha+beta) ≈ 0.9 and that a variant with A=5,U=0 has higher p than A=0,U=5. Independently assert p is monotone increasing in A and decreasing in U. A swap survives most reproduce‑the‑code tests but fails this.
- **Scale invariance sanity:** duplicating every observation should leave mu roughly unchanged but *increase* kappa; if kappa is insensitive to doubling n, the weighting is not doing what is claimed.

Target exclusion and same‑residue alternatives (P0/P1): the stated rule leaves a sibling variant at the *same residue* at d = 0, W = 1. *Test:* report, per plotted point, the fraction of normalized weight coming from d = 0 donors; if that median exceeds ~0.5 the "neighborhood" feature is largely a same‑residue feature and should be labeled as such. Separately, *test for record identity collisions*: two HGVS strings mapping to the same underlying patient/report must not appear as two donors. Falsifier: shuffle residue coordinates while keeping counts; if D_i ranks barely move, geometry contributes nothing and the panel is a count plot.

Coordinate/isoform mapping (P1): *Test* round‑trip residue numbering — take 20 residues, map sequence→structure→sequence, assert identity and matching amino‑acid identity at each position; assert no negative or off‑by‑one indices at construct termini; assert multi‑chain/multi‑model structures are not double‑counted. Falsifier: a systematic ±1 offset usually shows as a one‑residue shift in the D_i profile versus the raw‑count profile.

### 2. Flatness: shrinkage vs kernel range

Normalized W sums to 1, so if all p_j are near the prior mean, D_i ≈ mu everywhere regardless of bandwidth. This is shrinkage flatness, not geometric smoothing, and the synthetic Beta(1,9) example makes it concrete: every A=0/U=1 donor is exactly 9.09%, so any weighting returns 9.09%.

Report these per plotted point (P0) — they are cheap and jointly diagnostic:
1. **raw kernel mass** Σ_j K(d_ij) (unnormalized);
2. **nearest‑donor distance** and **number of donors within 3 Å / 6 Å**;
3. **normalized‑weight radii** r50, r90;
4. **effective donors** 1/ΣW²; **effective exposure** Σ W_j·n_j;
5. **prior‑mass fraction** Σ_j W_j·(alpha+beta)/(alpha+beta+n_j) — the share of D_i that is literally the prior;
6. **top‑3 contribution decomposition**;
7. **context disagreement** (spread of D_i across contexts when averaging equally);
8. **bandwidth sensitivity**: recompute at half and double the 3 Å half‑distance.

Discriminator: high prior‑mass fraction **and** bandwidth‑invariant ranks ⇒ shrinkage dominance (more kernel tuning cannot help). Rank changes under bandwidth ⇒ range dominance. Note the brief's warning is correct: raw K decays like ≈exp(−0.366 d), yet after normalization an isolated residue's weight can be carried entirely by a 30 Å donor with raw mass ~1e‑4.

Also (P1): the polymer surrogate d = 3.8·√|i−j| is an *expected* random‑coil separation, not a contact distance; K(E[d]) ≠ E[K(d)] by Jensen, and 38 Å polymer is not physically comparable to 38 Å measured. *Test:* plot D_i computed from 3D‑only donors against polymer‑inclusive D_i; if they disagree beyond the declared tolerance, polymer donors must be a separate, separately colored series.

### 3. Biological reasonableness

Neighborhood sharing is plausible only for *mechanistically coherent* neighborhoods. Tolerated and disruptive substitutions coexist at a sensitive residue; loss‑ and gain‑of‑function variants can cluster in the same pocket with **opposite** phenotypes, and averaging p over them erases sign. *Test:* stratify donors by mechanism/phenotype direction and check whether stratified D_i separate; if they do, the pooled feature is uninterpretable. Germline vs somatic rows, age‑dependent penetrance (an "unaffected" carrier is unaffected *at a censoring age*), sex‑specific risk, and case‑series A vs population‑reference U mean A/n is an ascertainment ratio, not a risk.

Most important: **a catalogue classification is not a patient count**. If pathogenic labels were assigned partly via PM1/PS1/PP3, the label is already a function of the neighborhood, and the feature predicts it tautologically. *Falsifier (P0):* recompute after dropping rows whose evidence includes PM1/PS1/PP3; if the signal collapses, report circularity.

### 4. Plot semantics — pass/fail

P0: own posterior and neighborhood feature never share an unlabeled axis; residues with multiple variant units either show all units or declare aggregation, weighting (count‑weighted vs equal‑variant) and multiplicity in the panel; unsupported positions render as gaps with broken lines, never 0; 3D vs polymer vs unsupported use distinct encodings; no interval is drawn for D_i (the model supplies none) — between‑donor spread, if shown, is labeled "between‑variant spread," not uncertainty; diagnostics (raw‑fraction, raw‑K pool) are labeled diagnostic. P1: hand‑recompute one plotted D_i from its listed donors and assert agreement to 1e‑6.

### 5. Ownership, endpoints, evaluation — before any risk tuning

Order: (1) type and audit source rows; (2) deduplicate repeated people, cohorts and re‑published families — genuine relatives are real observations but dependent, so report an effective sample size; (3) only then evaluate. Held‑out splits must be by **family/source/cohort**, never by variant. Internal agreement with the same fitted posterior labels establishes nothing.

External checks and their limits: deep mutational scanning and biophysical assays validate *molecular disruption* and can refute a hot‑spot claim, but not penetrance; electrophysiology/trafficking establishes direction (LoF vs GoF), not risk; ΔΔG/burial is orthogonal structural support only; only an independent cohort with observed endpoints, analyzed with calibration curves and a Brier decomposition, can address disease‑risk calibration. Until that exists, the primary output is a descriptive feature and all alternatives stay labeled diagnostics.
