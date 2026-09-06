# Source-only CLI diagnostic disposition

All three CLIs received the same frozen MYBPC3 20433692 source without gold,
baseline predictions or tools. These are fallible diagnostic readings, not
reference adjudications or additional scored benchmark arms.

Their useful shared findings are supported by the supplied source: distinguish
the literature-pooled Table 3 from current-study patients; retain family plus
person IDs; distinguish genotype-positive carriers from affected non-carriers;
and preserve disagreement between prose and table family IDs. A short person
roster followed by deterministic aggregation is a concrete alternative to
repeating long scalar-count JSON and whole-paper reasoning.

Qualifications to their individual claims:

- Claude called the R502Q table people "five relative rows plus" index cases.
  Stars already identify index cases among those five people; adding them again
  would double count. The source-roster prompt explicitly prevents this.
- Agy discussed a 32k context window. The observed failure was a 32k **output
  and reasoning budget**, not context overflow. The successful compact probe
  changes source selection, task, schema and effort together; it does not
  isolate a low-versus-medium model effect.
- Neither "No Dx" nor the lack of non-carrier rows proves asymptomatic carrier
  counts. Source prose distinguishes six healthy carriers from four with
  suggestive ECG findings. The diagnosis-cell transcription retains that
  uncertainty and does not assign an accepted phenotype scalar.
- The table text appears twice in the input with the same hash. This is one
  table with duplicate representations, not two cohorts. The compact probe uses
  one copy; all original scientific comparison arms retain the frozen source.
- The compact source-transcription audit checks all 48 person rows, line IDs,
  genotype/status cells and group carry-forward against the converted text.
  It cannot independently validate the original DOC's geometry or resolve
  clinical definitions. There is no scalar gold score for this diagnostic.

These findings support testing an evidence roster plus aggregation next. They
do not justify accepting raw model sums through the unchanged literal validator.

Follow-up source audit before scoring: the original DOC was rendered on all
three pages and its HTML rowspans checked. All six compared fields agree for
48 people after two blank family cells are bound using main-text prose. The
footnote explicitly defines No Dx as "Unaffected or healthy," while the prose
distinguishes ECG-suggestive carriers from healthy carriers. This strengthens
the layout check without resolving which phenotype endpoint the benchmark
should represent. See `../original_layout_audit.json`.
