# Recommendations after source and tranche checks

These are source-supported next changes. Except for the verified one-paper
source upgrade, they were not implemented during the paired runs and their
future gains are not measured results. Keep the same acceptance thresholds
and use Azure for routine experiments; no Anthropic call is needed for the
work below.

| Priority | Change | Evidence and success check |
| --- | --- | --- |
| 1 | Discover and validate article-specific components, then compare their quality when updating the cache. | PMID 25163546's actual article-page Supplementary Data link yielded a 53-page roster PDF; two unrelated cached files had falsely represented supplement availability. A file-count tie initially prevented the corpus builder from copying the two real files. Require article/DOI linkage, successful format parsing, table/figure identity and file hashes. A valid table must outrank an unrelated larger file or an equal number of files. |
| 2 | Continue clinical-table reading after identity extraction succeeds. | PMID 30059973 already has Tables 11/14 in the manuscript. Stop treating successful identity extraction as completion of carrier/A/U extraction. Route only the still-missing fields to the clinical component, with its headers, footnotes and cohort context. This remains an untested implementation. |
| 3 | Preserve table grids and join explicit patient IDs. | PMID 20433692's original DOC has merged mutation/family cells; its folded text shifts columns and duplicates the supplement. Reconstruct the grid, then deduplicate family + person + genotype before reading diagnosis. Extend to 21302287 and 30403697 only when IDs and variant ownership are explicit. Affected non-carriers must not become variant carriers. |
| 4 | Abstain independently by field. | PMID 18929323 explicitly supplies carrier counts 13 and 6. The candidate quotes them but clears its structured carrier fields because the A/U split is unavailable. Preserve source-supported carrier integers while leaving A/U unknown; do not allocate aggregate symptoms across variants. |
| 5 | Resolve endpoint conflicts before claiming reference-matching gains. | PMID 25814417 already produces 97 affected / 62 unaffected / 26 unknown from its patient rows. Reference 73/106 uses a different endpoint/denominator; this is not additional parser upside. PMID 27566755's COUNT column likewise does not by itself specify the clinical affected endpoint. Keep evidence and uncertainty visible. |

The downloader should keep a small component inventory for each paper: article
identity, expected component label, discovered link or embedded PDF page range,
retrieval result, format validation, source hash, converted representation and
whether the component contains identities, counts, phenotypes or only aggregate
statistics. Supplement files may legitimately lack per-variant counts, as the
new 25163546 roster does. The correct outcome is a recorded evidence gap and
field-specific abstention, not an invented count of one.

Follow body references such as “Online Table” to embedded PDF pages as well as
external links. Include the full table header hierarchy, caption, footnotes and
endpoint definitions in extraction input. Retain original DOC/PDF assets so
conversion can be repaired without downloading or paying for a second model
read. Repeatedly downloading the main body cannot repair a missing roster.
Store retained archives with members in their standard extraction directory,
and use the normal fold sentinels. The local 25163546 reuse check preserves
all 20 source-listed cDNA notations, folds exactly two components, and is
byte-identical on a second fold. This avoids making the raw archive produce
a second copy of tables that were already extracted.

Spend API budget only after local checks establish useful evidence. Try
deterministic format and table reconstruction first; then a bounded Azure call
on the relevant component and unresolved fields. Use vision only when text/grid
recovery demonstrably loses the evidence. Cache by source hash, component and
extraction protocol. Repeated full-paper calls are a poor default for all-NA
rows whose required component is absent.

Measure three separate outcomes: acquisition of the required component,
correct field capture given that component, and end-to-end performance across
all registered papers. Report identity TP/FP/FN, A/U supply and exact supplied
values, and omission-inclusive MAE together. Preserve unknown phenotype and
distinct cohort denominators; never infer unaffected as total minus affected
without an explicit complete partition. Keep incomplete-reference extras
visible in official scoring and separately inspect their source support.

Operational follow-up: normalize or reject malformed
`extraction_metadata.total_variants_found` before arithmetic, and always emit a
failed RUN_STATUS for a process-level exception. One tranche-03 baseline job
needed an unchanged-runtime operational retry for this type mismatch; no
production fix was inserted between the measured arms.

The second tranche adds a necessary constraint on priorities 2–5: bind each
count to its **unit and cohort ownership** before retaining it as a person
count. SCN5A 20129283 H558R's Number = 408 matches gold but appears in a table
captioned as reference alleles, using previously reported healthy controls.
A table-cell provenance stamp or an outlier exemption cannot establish people,
phenotype endpoint or new-study attribution. Preserve the raw cell and its
caption, mark the ambiguity, and adjudicate it before claiming that the large
unaffected improvement reflects correct clinical counts. This single row drives
most of the new pooled A/U gain; remaining-case improvement is much smaller.
