# BRCA2 source adjudication for the 14 September residue refresh

The two known BRCA2 missense queues are now resolved for the current
breast/ovarian hereditary-disease fit: their unsupported affected assignments
are quarantined, with the original observations and source information retained.
Unknown germline status or carrier phenotype is not converted to an unaffected
count. This is an additional source correction, not a change to the empirical
prior formula, gnomAD policy, variant unit or structural kernel.

The starting clinical-key ledger is the correction published at `c8e34c96`,
SHA-256 `8202b4021bb0466f00c31fbcdc1a50bbd6474c1831eb1f06d4ee80fb13722cf6`.
All five genuine PMID 40664060 breast-cancer observations restored in that
correction remain unchanged: R118H, P2276T, V2076I, S28N and W2619C.

## Decisions and exact observation scope

| Source | Type | Original observations | A quarantined | U changed | A restored |
| --- | --- | ---: | ---: | ---: | ---: |
| PMID 33054725 | Missense | 205 | 205 | 0 | 0 |
| PMID 33054725 | Nonsense | 25 | 25 | 0 | 0 |
| PMID 36385461 | Missense | 23 | 23 | 0 | 0 |
| PMID 36385461 | Nonsense | 255 | 255 | 0 | 0 |
| **Total** | | **508** | **508** | **0** | **0** |

The previous canonical missense queue contained 227 of these observations:
205 from PMID 33054725 and 22 from PMID 36385461. The remaining missense row
already failed the frozen canonical eligibility checks. The exact same source
proof also resolves the listed nonsense records; the two variant classes
remain separate in downstream fitting. Other variant types and papers are
unchanged by this bounded task.

Each decision identifies the original `(key, PMID, observation variant ID)`,
the frozen DB variant ID, affected-fact source row, exact full-text line(s),
original workbook and sheet, and one-based Excel row. All **508 affected fact
quotes match the complete archived primary source and original supplement
workbook rows**. This is not a PMID-only subtraction.

## PMID 33054725: tumor samples with unavailable germline status

The [primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7556962/) is a
pan-tumor study. Every selected missense/nonsense affected fact comes from
Supplementary Table 1, the TCGA cohort, in
`12885_2020_7481_MOESM3_ESM.xlsx`, Sheet1. These are real tumor sample rows;
the paper's limitations explicitly state that germline mutation data were
unavailable for the TCGA cases. The selected rows also describe cancers other
than breast or ovarian cancer. They therefore do not establish affected
germline carriers for the current endpoint. It would be inaccurate to label
all of them proven somatic or invented patients.

For example, A1439T is the record for `TCGA-BR-4184-01`, stomach adenocarcinoma.
The parser inferred A=1 from that row. Its actual tumor observation is retained
in the evidence, while its contribution to this hereditary-endpoint fit is
removed. The 205 missense records span 18 other tumor categories, including
bladder, colorectal, melanoma, lung and stomach cancers.

The separate SYSUCC workbook was checked for possible valid restorations.
Its seven explicitly germline BRCA2 entries have stomach, pancreatic or
colorectal endpoints: three frameshift entries and four nonsense entries.
None is missense or an established breast/ovarian case. The workbook has two
columns called “Cancer Type”; the earlier one contains repeated TCGA labels
that conflict with the SYSUCC cohort context. We preserve both fields in
`SYSUCC_germline_review.csv` and use the SYSUCC-specific endpoint field, not
the repeated ovarian labels. No affected or unaffected restoration qualifies
from this separate cohort.

## PMID 36385461: classification is not a carrier phenotype

The [primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10098510/) compiles
Asian variants from publications and genomic databases. Exact affected-fact
tracing corrects the preliminary September 13 queue description: the selected
affected facts originate in **Supplementary Table S7B**, labeled `T24` by the
extractor, in `IJC-152-1159-s011.xlsx`, BRCA2 sheet. They are not sourced from
the shorter Table S3C inventory used as an initial warning example.

S7B does contain pooled carrier totals. Its “Cancer”, “Non-cancer” and “Total”
columns are numbers tested across contributing cohorts, not counts of affected
and unaffected carriers of the row's variant. The parser used the clinical
classification column as a phenotype and inferred one affected carrier. For
five nonsense records, the variant-paper location points to a main-text table,
but the actual affected=1 fact still traces to the S7B row; the fact-level
identity controls the decision.

| Variant | Source carrier total | Cancer subjects tested | Non-cancer subjects tested | Source Excel row | Frozen A |
| --- | ---: | ---: | ---: | ---: | ---: |
| D2723G | 3 | 63,828 | 37,086 | 179 | 1 |
| I2675V | 692 | 80,248 | 271,840 | 189 | 1 |
| R2336H | 73 | 20,085 | 600 | 270 | 1 |
| V211L | 2 | 133 | 0 | 887 | 1 |

The companion complete inventory, Table S2, also states that 1 is substituted
when the carrier number is unavailable. We do not assume every displayed 1 is
an imputed value. We preserve the reported total, all tested denominators and
reference codes for every selected row. Those fields alone do not prove the
individual carriers' germline status and breast/ovarian endpoint, nor resolve
overlap with original papers already represented in the fit. Even where no
non-cancer subjects were tested, “cancer” is not a confirmed endpoint-specific
carrier count. Consequently the unsupported frozen A=1 is quarantined; the
pooled totals are not reassigned to A or U. Reinstatement would require the
underlying clinical evidence and overlap adjudication, not a classification.

## Files and rebuilding contract

- `BRCA2_literature_identity_corrected.csv.gz` is the **authoritative complete
  clinical-key table** for the refreshed union. It preserves the previous
  schema and adds `refresh_*` accounting fields. It has 6,110 rows, of which
  2,138 retain clinical evidence.
- `BRCA2_clinical_key_deltas.csv.gz` is the affected-key view of the same table.
- `BRCA2_observation_decisions.csv.gz` records all 508 exact decisions, source
  fields, exclusions and zero restorations. Unknown source A/U remain missing;
  the zero contribution to the fit is an exclusion, not an observed zero.
- `selected_source_records.csv.gz` preserves the original count provenance and
  affected fact quotes. `preserved_40664060_restorations.csv` records the five
  previously established clinical observations.
- `SYSUCC_germline_review.csv` records the separate-cohort restoration check.
  `decision_summary.csv` and `checks.json` provide compact receipts and hashes.

Before rebuilding the clinical/population union, drop rows marked
`drop_clinical_key_zero_clinical_evidence` (equivalently `n_literature == 0`).
Their observed gnomAD alleles must remain and be reassigned to genomic
population units where appropriate. Do not merely subtract counts from the
old protein aggregate while retaining its old grouping. Do not duplicate
the five genuine clinical restorations or remove another paper's observation
because it shares a DB variant ID. The source audit does not alter gnomAD
counts or perform the union/model/structural rebuild.

`adjudicate.py` reproduces the decisions from the frozen observations, full
primary text, original XLSX files and a read-only SQLite connection. It asserts
the baseline clinical-ledger hash, original source-DB hash, exact source matches,
all observation identities, arithmetic, unchanged previous restorations and
unchanged unaffected counts. It uses pandas, numpy and openpyxl. The generated
clinical ledger has SHA-256
`573c091ac1c88f8c5733f91b27681887fe9f1b856f4395b62eef476d2d8b4929`.

The known 227-row missense queue is resolved **for inclusion in this fit**.
Its excluded carrier endpoints are still unknown, and a future source review
may recover valid observations. This does not validate every remaining BRCA2
paper, every other variant type, or the other four genes. The refreshed plot
should be described as using adjudicated source exclusions, not as a fully
ascertained or clinically calibrated penetrance map.
