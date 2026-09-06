# Observed failure mechanisms

These observations come from source files, request/response traces and code.
They are distinct from the final gold-agreement metrics.

- **A requested setting was not an effective setting.** LiteLLM silently
  discarded reasoning effort for the new deployment names. The corrected path
  explicitly allows the parameter, with actual SDK HTTP-body tests. Earlier
  pilots remain abandoned unscored; their labels cannot rank reasoning levels.
- **The model can consume its output budget without emitting an answer.**
  Corrected Astra-medium calls on MYBPC3 20433692 and RYR2 25814417 exhausted
  32k as reasoning and emitted no visible extraction. The same happened on
  other large-source calls during the run. This does not establish an intrinsic
  model limit; the current prompt, effort and cap are part of the outcome.
- **Retries can repeat the same failure and block downstream work.** The
  production empty-response retry, JSON repair and failed-paper retry pass
  are retained in this test. Hidden SDK retries add another request layer.
  Repeated full-paper calls are costly even when they return no new evidence.
  Budget waits affect run elapsed time and must not be attributed wholly to
  inference latency.
- **A failed reader can prevent table hints from being merged.** In MYBPC3
  20433692, eleven source-filtered table hints were available before the first
  Astra attempt. Its JSON repair returned an empty object and extraction
  raised before `_merge_table_variants_with_overflow_qc`. The candidates are
  not thereby proven correct, but their fate deserves a separate, validated
  fallback test; a reader failure need not mean the source contained nothing.
- **A complete-looking body can lack the relevant table.** RYR2 18929323's
  supplied text contains explicit variant carrier totals, which both current
  readers capture, but cites a patient-characteristics Table 1 whose body is
  absent. Its symptomatic/asymptomatic prose totals do not identify variants.
  Controls, repeat measurements and variant carriers are different groups.
- **Correct table transcription is not yet a clinical count.** The MYBPC3
  compact probe transcribes 48 people accurately in the audited fields. Its
  No Dx footnote and prose support different phenotype distinctions, and two
  blank family cells require prose for their join. Repeated copies of the same
  table are one cohort. Index cases already occur among the people.
- **The acceptance contract intentionally limits gains.** The independent
  reader may propose patient-row-derived values, but the existing literal
  validator rejects those as primary fills. It cannot add missing identities
  or overwrite existing nonnull values. For example, the independent reader
  proposes the three current-study intronic variants IVS6+5G>A, IVS11-9G>A and
  IVS29+5G>A from MYBPC3 20433692, but the fresh baseline retains only the 13
  protein-labelled variants, so those proposals cannot enter this overlay.
  This observation is source-based, not a gold-match judgment.
  Raw proposals and accepted values
  measure different stages; neither model agreement nor a good raw score can
  silently waive the evidence contract.

The follow-up should isolate bounded reading, source retrieval/representation,
validated person aggregation and error handling. It should not assume that
selecting a stronger model will fix all four at once. See
`roster_followup_design.md` for the unimplemented experiment design.
