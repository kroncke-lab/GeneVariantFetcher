# HNF1A endpoint sanity corrections, 14 September 2026

The bounded review found that the HNF1A input mixed type-2-diabetes studies,
somatic liver-tumor variants, liver adenoma status and MODY counts. Ten exact
observations are corrected before the residue map is rebuilt. Across the
complete literature input, A changes from **762 to 633**, and U from **115 to
83**. The new table has 234 clinical keys, of which 229 retain clinical evidence.
All observed population alleles must remain in the rebuilt union.

| Source and variant | Frozen A/U | Current MODY contribution | Reason |
| --- | ---: | ---: | --- |
| PMID 24915262, E508K | 52/12 | Excluded | Explicit type-2-diabetes case/control study |
| PMID 28116330, I27L | 38/23 | Excluded | Explicit type-2-diabetes case/control study |
| PMID 20172480, I27L | 31/0 | Excluded | Abstract-only mixed MODY3/MODY2/type-2 cohort; variant ownership and MODY endpoint assignment unresolved |
| PMID 30121369, L214Q and E32X | 1/0 each | Excluded | Explicit somatic tumor variants |
| PMID 29101032, Q511L | 1/0 | Excluded | Tumor-only variant, absent from adjacent non-tumorous liver |
| PMID 38133737, K205E | 1/0 | Excluded | Healthy liver is wild type; no pathogenic germline HNF1A variant identified |
| PMID 26631547, S247T | 1/0 | Excluded | Somatic hepatocellular-carcinoma allele |
| PMID 31483937, R272H | 5/0 | **4/1** | Five carriers, four with diabetes and one with no diabetes onset at follow-up |
| PMID 14598263, R229X | 2/1 | **0/3** | Three germline carriers with normal glycemia; two had liver adenomas |

“Excluded” means the source does not establish a count for the requested
endpoint. It is not an observed A=0/U=0 and does not create an unaffected person.
The four restored unaffected observations have explicit source evidence of no
diabetes at the reported follow-up; they are not claimed to remain unaffected
throughout life. Missense and nonsense remain separate downstream classes.

## Source interpretation

[PMID 24915262](https://pubmed.ncbi.nlm.nih.gov/24915262/) reports E508K in a
type-2-diabetes case/control study. The 52/12 counts are real for that endpoint;
their removal from a MODY fit is an endpoint restriction, not a claim that the
paper or its carrier counts are wrong. Similarly,
[PMID 28116330](https://pubmed.ncbi.nlm.nih.gov/28116330/) Table 1 provides the
38 and 23 I27L carriers in diabetic and control groups. The broader-diabetes
sensitivity input restores these two established case/control observations
while preserving all other repairs.

[PMID 20172480](https://pubmed.ncbi.nlm.nih.gov/20172480/) is available in the
frozen source only as an abstract. It partitions 31 subjects into 10 with MODY3,
15 with MODY2 and six with type-2 diabetes/I27L. That is insufficient to assign
all 31 to HNF1A I27L and the requested MODY endpoint. We quarantine the original
observation rather than asserting that its correct MODY count is six. This
source remains unresolved for a future full-text review.

The four missense somatic observations and the E32X nonsense observation are
traced to the exact primary source text. Somatic tumor alterations are not
germline MODY carrier observations. The source phrases, original DB identities,
source files, line locators and decisions are in `observation_decisions.csv`.

[PMID 31483937](https://pubmed.ncbi.nlm.nih.gov/31483937/) describes the R272H
family in its Results: the index patient developed diabetes eight years after
his liver presentation; his mother, one sister and brother also had diabetes.
The other genetically confirmed sister had no diabetes onset at follow-up.
All five had liver lesions, explaining the misleading frozen 5/0 partition.
The separate [PMID 39408812](https://pubmed.ncbi.nlm.nih.gov/39408812/) R272H
observation is retained: that patient is explicitly described as having MODY
with a confirmed heterozygous germline variant. An adenoma keyword alone is not
grounds to exclude a patient.

In [PMID 14598263](https://pubmed.ncbi.nlm.nih.gov/14598263/), Family A comprises
three genetically confirmed R229X carriers. A1 and A2 had liver adenomatosis;
A3 did not. The clinical descriptions establish normal fasting glucose in all
three, with normal glycated hemoglobin also reported for A2. The frozen 2/1
partition describes adenoma status, so the diabetes partition is corrected to
0/3. The family's unrelated ancestor with late-onset diabetes is not assigned
this genotype. The same paper's G55fs family has a more complicated history
including gestational diabetes and elevated glycated hemoglobin; that
frameshift observation is outside the requested missense/nonsense fits and is
not reclassified in this bounded review.

## Coverage, artifacts and limits

The initial source screen selected the three largest original affected-count
contributors and already-flagged somatic/adenoma source rows. It covered 18
original observations with A155/U49; nine were on canonical missense keys with
A131/U35. Remaining flagged sources were checked for context: coexisting
adenoma or a somatic second hit does not itself invalidate a clinically
documented germline MODY carrier. No claim is made that every remaining source
is MODY-specific, independent, or free of diagnosis-based count assumptions.

- `clinical_input.csv.gz`: complete authoritative clinical-key ledger, using
  the original identity-flag schema with refresh accounting fields.
- `observation_decisions.csv`: ten exact `(key, PMID, variant ID)` decisions,
  original counts, removals, restorations and source evidence.
- `corrected_observations.csv.gz`: retained observations plus explicit source
  replacements, independently summed to check every clinical key.
- `clinical_key_deltas.csv`: affected-key view of the complete ledger.
- `broader_diabetes_sensitivity_clinical_input.csv.gz`: optional alternate
  endpoint input that restores the two established type-2 case/control rows;
  this is not the primary MODY input and has not been fitted here.
- `adjudicate.py` and `checks.json`: reproducible arithmetic, source text
  checks, hashes and output receipt.

Drop zero-clinical-evidence keys before the population union is rebuilt, then
rejoin retained observations and all observed gnomAD alleles by identity.
Do not duplicate population alleles or keep empty former protein aggregates.
The primary ledger SHA-256 is
`7bb88e32ab4a71287d6d92a410824600b4ff74f4ca8487f0e8b791d046f1b3e9`.
This source task does not change the empirical formula or structural method,
and does not establish a clinically calibrated MODY probability map.
