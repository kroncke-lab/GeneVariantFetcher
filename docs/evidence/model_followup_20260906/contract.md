# Frozen bounded-roster experiment

This is an opened-source diagnostic, not a new gold benchmark or a production
count-recovery promotion. Eight packets represent seven previously opened papers.
Five positive packets contain 139 row appearances (91 distinct study people;
the 48-person DOC is tested twice). Three negative packets have no reconstructable
identified-person roster. Aggregate evidence in a negative packet is not absent
scalar-count evidence; refusing to fabricate people is the task here.

Astra low and medium receive identical source, prompt, schema, output cap and
180-second direct HTTP deadline within each packet. Cap is 4096 for positive
packets and 1024 for negative packets. One draw per cell, random interleaving
(seed 20260906), maximum two concurrent Astra requests. No SDK or content retries.
Failures and length stops remain in the operational denominator. This measures
performance under these budgets, not intrinsic reasoning capacity. The DOC's
original-grid versus shifted-Markdown pair is a separate representation ablation.
The structured input expands only actual original HTML rowspans, retaining blank
cells. Neither arm receives reference records or benchmark gold values.

The source reference and source packets are hashed before dispatch. Reference
records were prepared from original DOC geometry, article tables and prose by
this investigator; they are not independent blinded human adjudication. No new
gold comparison is allowed until all planned model outputs are locked. Programmatic
comparison normalizes whitespace, prefix notation, No Dx spacing and identifier
colon/hyphen typography; it does not normalize a phenotype into another endpoint.

Each person record must carry source row locator(s), target-gene variant(s),
family and person IDs, mutation membership, literal phenotype, index-case status,
and endpoint classification. The packet identifies its study/cohort and endpoint
at the relevant baseline or recorded diagnosis time. Validation reports source
row identity, missing/extra people, source-location errors, variant binding,
genotype, literal phenotype, index status and endpoint differences separately.
Reported totals cannot repair incorrect person records. Reject unknown source
locators, conflicting duplicate person records, wrong gene binding, invented
people and schema violations. Count each unique carrier once per variant;
co-carriers may contribute once to each of their variants. Noncarriers do not
contribute. Unknown endpoints remain explicit. A subset is never silently called
a complete paper count: complete reference roster coverage and correct genotype/
variant binding are required for a complete packet-derived count.

Endpoint contracts:

- MYBPC3 20433692: table-defined diagnosis (numeric age positive; No Dx negative;
  unknown mark unknown). This does not certify ECG-normality. Suggestive ECG/prose
  family discrepancies remain limitations. Blank original family cells require
  explicit prose joins; index cases are already table people.
- MYBPC3 21302287: clinical HCM at enrollment, Table 3(a), all listed related and
  unrelated patients. Unrelated-proband cohort totals are a different denominator.
  Copy reported protein notation; no invented correction of article discrepancies.
- RYR2 25435091: initial diagnostic EST/adrenaline CPVT. Earlier symptoms and
  asymptomatic follow-up are different timepoints/endpoints.
- RYR2 30403697: documented VT and/or SCA at presentation. This is explicitly a
  different endpoint from all clinical CPVT or any symptom. Do not label its
  negative values as generally unaffected. Only numbered Table 1 subjects enter;
  family-history relatives remain outside this defined cohort.

Primary outcomes: complete valid JSON response rate, exact person roster,
reference-field agreement, source-bound person/variant coverage, and deterministic
packet-derived C/endpoint-positive/endpoint-negative/unknown totals. Report raw
and source-validated results separately. Existing gold agreement, if calculated,
is a secondary diagnostic with endpoint/cohort discrepancies visible and is never
the acceptance criterion. No percentage-of-human or population gain claim.

The old campaign ledger remains immutable. New total reserve ceiling is $4.90,
within its remaining original $150 envelope; input byte bounds include overhead
and Astra cache-write premium. Text inputs stay below long-context thresholds.
Unknown charges retain reservations, including full documented 128k Grok output
uncertainty for early legacy max_tokens probes. Explicit documented cap spelling
is used for subsequent probes. CLI consultations are excluded as authorized.
