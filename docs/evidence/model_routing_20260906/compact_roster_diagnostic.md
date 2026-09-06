# Why the compact roster is useful, and what it does not prove

This is a source-only diagnostic on MYBPC3 20433692, selected after an observed
Astra output-budget failure and before any new gold scores. The three CLI
readings suggested separating patient transcription from clinical aggregation.
No diagnostic values are merged into benchmark predictions.

The frozen input repeats one 4,488-character clinical table twice, byte for byte.
The compact prompt uses one copy with source line IDs and returns table cells
in arrays. It keeps family/person IDs, genotype status, diagnosis cell text,
index-case markers and whether grouping was inherited from a merged cell.

Astra low returned all 48 table people, with 48 unique family/person pairs,
in 49.8 seconds. Its 3,031 input and 1,463 output tokens correspond to a
$0.10346 proxy at the campaign rates. All line IDs, grouping and emitted cells
match a mechanical audit of the converted text; the three wrapped fragments
were excluded correctly. The 20 starred index cases are already included.
The roster has 44 genotype-positive and four genotype-negative people.

For comparison, the two medium-effort whole-paper calls consumed 64,000
reasoning tokens and emitted no visible extraction. Their returned-usage
proxy is $3.54630 before the small repair call. These are different tasks:
source selection, output schema, effort and scope all changed together. This
is an operational diagnostic, not a causal estimate of the benefit of low
reasoning or a per-correct-count benchmark.

The transcription preserves 34 carrier rows with numeric diagnosis age,
nine with No Dx, and one unknown. Those groups are **not new accepted A/U
counts**. The prose separately names six healthy carriers and four carriers
with suggestive ECG findings, including inconsistent family identifiers.
The subsequent original-DOC layout audit inspected all three rendered pages
and expanded the HTML rowspans. The six compared fields match on all 48 people,
with a qualification: two noncarrier family cells are explicitly blank; their
H49 membership is supported by the main-text prose, not by a merged table cell.
See `original_layout_audit.json`. The table footnote explicitly defines No Dx
as "Unaffected or healthy"; those nine carriers include three whom the prose
calls suggestive on ECG. Thus a table-defined unaffected endpoint, no diagnostic
HCM, and ECG-normal/healthy are different predicates. Preserve that definition
and the prose discrepancy instead of silently choosing one. These checks do not
establish reference-count agreement. Agreeing CLIs are not reference adjudicators.

Next implementation to test: produce bounded, source-linked person records;
retain explicit cohort, variant and genotype ownership; validate source rows
and deduplicate family/person IDs; classify phenotype with the paper's endpoint;
aggregate only validated records in code. Keep contradictions and unknowns
visible. This complements literal count evidence instead of bypassing its
validator with unverified model sums.
