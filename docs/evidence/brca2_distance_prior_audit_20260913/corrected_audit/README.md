# Independent corrected BRCA2 audit

The corrected union, empirical priors, posterior orientation, selected structural
calculations and complete residue aggregation pass independent checks. The
[executable audit](audit_corrected.py) and [receipt](audit.json) preserve the
checks and final source hashes. The corrected runner was read but not modified.

All 24,851 eligible population alleles occur exactly once in the rebuilt union,
preserving all 15,877,471 gnomAD carrier observations across consequences.
Removing 3,644 clinical keys with no remaining carrier evidence released 782
population alleles to separate genomic units; two alleles were instead assigned
to restored real patient keys. Each unit's population count was independently
recomputed from its member alleles. No A=0/U=0 units remain.

The source quarantine and five patient restorations were applied independently
to the original clinical-key count ledger and match the authoritative corrected
ledger exactly. The separate class fits were recomputed directly from the
historical weighted mean/MSE formula:

| Corrected class | Units | Affected | Literature unaffected | gnomAD unaffected | α | β |
|---|---:|---:|---:|---:|---:|---:|
| Missense | 4,584 | 2,085 | 1,482 | 1,516,774 | 0.202612130 | 4.219485782 |
| Nonsense | 453 | 1,174 | 108 | 13,375 | 3.014199618 | 1.772592652 |

Every saved posterior satisfies **α posterior = α empirical + affected** and
**β posterior = β empirical + literature unaffected + gnomAD carriers**. The
missense prior mean is 4.581810% and its strength is 4.422098. Alpha and beta
are not reversed.

Four released alleles explain the 511-carrier reduction within the strict
nonsense class. Their population annotations are canonical-transcript
frameshifts, with nucleotide length changes not divisible by three. The old
clinical protein keys ending in X had assigned them to nonsense. All four
remain in the complete union and correctly stay outside both strict nonsense
and missense fits. The [allele trace](released_frameshift_alleles.csv) includes
the before/after unit, transcript, HGVS and sequence-length evidence.

| Old clinical key | Genomic allele | Length change | gnomAD carriers |
|---|---|---:|---:|
| L760X | 13-32336632-TTTATA-T | −5 | 1 |
| C2473X | 13-32355268-AGT-A | −2 | 6 |
| S3366X | 13-32398608-C-CT | +1 | 501 |
| R3370X | 13-32398619-C-CT | +1 | 3 |

Eight targets cover experimental, AlphaFold, polymer and zero-affected-donor
cases. Fresh enumeration with PPA's reference backend yielded 27 supported
contexts. The audit directly computes K(d)=2/[1+exp(log(3)·d/3)], verifies
target-identity exclusion and one vote per distinct donor within each context,
and reproduces the posterior average, observed-fraction average and prior/count
decomposition. It separately uses unnormalized K to pool A and U, adds one prior
within each context, and averages supported contexts equally. The maximum
absolute discrepancy across checked quantities is **1.14×10⁻¹³**. See the
[context calculations](selected_context_formula_checks.csv) and
[target calculations](selected_target_formula_checks.csv).

All five corrected distribution summaries were recomputed, including exact-zero
and ≤0.1% counts. All 3,418 rows in the saved residue table match independently
aggregated variant scores, counts, support and the previous residue curve,
including missing values. The primary supports 4,371 targets at 2,525 residues:
3,041 polymer and 1,330 structured targets. The 33 supported targets at 22
residues with no affected donor have exactly zero observed-fraction scores;
their posterior-neighborhood scores range from 0.199771% to 3.736785%.
Residues 865–866 contain six variants with A=0/U=20 and a corrected mean
neighborhood of 3.061727%, versus 10.785819% previously. A zero-affected singleton
has posterior mean 3.736785%; the fixed prior needs at least 199 unaffected
observations for an individual zero-affected posterior mean ≤0.1%.

This verifies calculations and the specified source correction. It does not
certify every remaining clinical observation, person-level independence or
disease calibration, and it does not rerun prediction-model outer LOO. The
[parent report](../README.md) identifies the remaining source/endpoint queues.

Reproduction from the repository root:

```bash
OPENBLAS_NUM_THREADS=1 ../BayesianPenetranceEstimator/.venv/bin/python \
  docs/evidence/brca2_distance_prior_audit_20260913/corrected_audit/audit_corrected.py
```

All 18 hashes in the finalized corrected input manifest matched on this run.
Audit script SHA-256:
`ef8f06f632eeebd54c4492533c35da81ad36ab2ddbdd0540c493292ebd06032f`.
Corrected runner SHA-256:
`d6fccca9ae904a3be2357d637d6924cadfdd7a1601e9f80b661d95e586793ff4`.
