This is a conceptual statistical review of the provided modeling framework. It is based entirely on the mathematical properties described in the prompt. No real dataset, source records, or implementations have been inspected, and this review does not constitute verification of any actual data.

**Arithmetic Check**
The arithmetic provided in the synthetic examples is exactly correct. With \(\alpha=1\) and \(\beta=9\), a variant with \(A=0\) and \(U=1\) yields a posterior \(p_j = 1/11 \approx 9.09\%\). A normalized weighted average of such identical values strictly equals 9.09%. Pooling 100 such unaffected observations under a single prior yields \((1+0)/(1+9+100) = 1/110 \approx 0.909\%\). Finally, achieving a posterior of \(\le 0.1\%\) (or 0.001) with \(A=0\) requires \(1/(10+U) \le 0.001\), which simplifies correctly to \(U \ge 990\).

**1. Components and Convex Constraints**
The model separates four mechanisms: the *kernel range* (decay of raw influence over distance), *normalization* (forcing weights to sum to 1), *shared per-variant shrinkage* (the beta-binomial prior pulling individual estimates toward 10%), and *variant-versus-carrier aggregation* (averaging posterior probabilities versus pooling underlying counts).

Because the feature \(D_i = \sum W_{ij} p_j\) uses normalized geometric weights (\(\sum W_{ij} = 1\), \(W_{ij} \ge 0\)), it is mathematically constrained to be a convex combination of the local donor posteriors. Consequently, \(D_i\) cannot be lower than the minimum \(p_j\) in the neighborhood. This does not create a *universal* positive floor for all possible counts—if a donor had \(U \ge 990\), the local estimate could drop to 0.1%. However, it does create a strict positive floor dictated by the *actual* sparsity of local counts. An all-zero-affected eligible segment cleanly isolates this phenomenon: it demonstrates that the apparent "nonzero risk" is purely an artifact of the shrinkage prior and small individual \(U_j\) counts, entirely independent of the kernel's range.

**2. Contrasting Estimators**
The three estimators target different estimands and rely on different assumptions:
*   **\(\sum(W \cdot p)\)**: Averages individually regularized risk estimates. It assumes variants are independent entities whose risks must be shrunk individually. Because \(W\) normalizes to 1, this estimator is invariant to arbitrary raw-K scaling, but it fails to accumulate regional evidence (averaging 9.09% yields 9.09%).
*   **\(\sum(W \cdot A/n)\)**: Averages raw empirical fractions. Its estimand is unregularized spatial risk. It avoids the prior's floor but is highly volatile, treating \(A=0/U=1\) as exactly 0%, which is statistically unsafe.
*   **Kernel-pooled \(\frac{\alpha+\sum(K \cdot A)}{\alpha+\beta+\sum(K \cdot n)}\)**: Pools underlying *evidence* rather than estimates. Its estimand is the regional risk, assuming local observations are drawn from a shared spatial pool. This allows unaffected counts to accumulate (e.g., dropping the risk to 0.909%). However, it is highly sensitive to arbitrary raw-K rescaling; multiplying \(K\) by 10 artificially inflates the effective sample size.

Pooling observations across different coordinate contexts introduces the danger of counting the same individuals repeatedly, which artificially inflates statistical confidence. Regional evidence is only conceptually transferable to an individual variant if the underlying pathogenic mechanism (e.g., regional structural destabilization) is shared spatially, rather than strictly variant-specific.

**3. Support-Aware Diagnostics**
A positive-tail kernel with normalized weights cannot be evaluated solely by its raw far-tail decay. If a local neighborhood is empty, normalization will scale up the microscopic raw weight of a distant donor to 100% influence. To diagnose this, the model requires support-aware metrics:
*   **Absolute (unnormalized) kernel mass**: Detects when the model is operating in a regional vacuum.
*   **Nearest eligible donor distance**: Flags estimates driven entirely by the far tail.
*   **Effective donor number**: The inverse sum of squared normalized weights.
*   **Raw pooled A/U evidence**: Contextualizes the actual local sample size.

If the absolute kernel mass is negligible, the model should explicitly abstain or fall back to a labeled baseline estimator. An alternative estimator must be explicitly justified; normalization must not be silently allowed to inflate distant noise just to force an estimate.

**4. Source-Count Audit**
A statistical estimator is only as reliable as its inputs. Variant database catalogue entries, clinical pathogenicity classifications, and raw carrier counts are distinctly different data types. Count-based estimators require the latter. A source-count audit must confirm valid observation ownership (ensuring cases and controls meet strict, distinct clinical endpoints) and rigorous deduplication. If the same individual is sequenced multiple times across different datasets, failing to deduplicate artificially inflates \(n\), directly corrupting the posterior.

**5. Recommended Order of Operations**
1.  **Audit Data:** Verify deduplication and valid clinical endpoints for all \(A\) and \(U\) counts before any modeling occurs.
2.  **Calculate Baselines:** Compute the individual \(p_j\) for all variants to identify the mathematical floor imposed by local \(U\) counts.
3.  **Implement Diagnostics:** Track unnormalized kernel mass and effective donor numbers to identify regions where normalization inappropriately amplifies the far tail.
4.  **Evaluate Estimators:** Compare the current \(\sum(W \cdot p)\) behavior against a kernel-pooled evidence model, explicitly separating mathematical deductions (e.g., the convex combination bound) from biological hypotheses (e.g., optimal kernel width).

A narrower kernel will not fix a flat spatial curve if local counts are sparse, and zero observed cases do not establish a 0.1% risk without massive unaffected exposure.
