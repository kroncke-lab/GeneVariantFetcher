# Follow-up design supported by the source diagnostic

This is an unimplemented next experiment. TASKS.md owns its priority and gate;
this document does not enable a recovery stage or approve a new count rule.

Use a capable model for bounded evidence reading: one count-bearing table or
short narrative section with its caption, column headers, footnotes, study-scope
sentences and a stable source locator. Deduplicate identical representations
before dispatch while keeping every original artifact in the source manifest.
Prefer structured table parsing and rowspan expansion when the original
DOC/DOCX/HTML exposes that information. Use a model for unresolved evidence
reading; routine cell expansion and arithmetic can remain deterministic.
Escalate ambiguous page geometry to the original page image. A paper with a
missing supplement needs retrieval first; another reader cannot recover unseen
patient rows.

The intermediate record should contain study/cohort, variant assertion, person
and family identifiers, genotype status, literal phenotype evidence, endpoint
and timepoint, index-case status, source coordinates, and an explicit statement
of which fields came from merged cells or another passage. An empty family cell
is different from a rowspan. Preserve both candidate bindings until evidence
supports the join. Do not attach a whole-cohort count to every variant or treat
an index person as an additional relative.

Separate transcription, endpoint interpretation, and arithmetic:

1. Transcribe a compact person roster without emitting scalar counts. Start
   with Astra low and a small explicit output budget. Validate schema, row
   coverage, source locators, unique person identities and genotype membership.
   A malformed or truncated response gets a smaller task, with a bounded retry;
   repeating the same exhausted reasoning call is not a recovery strategy.
2. Read the paper's phenotype definition and map literal evidence to it. Keep
   diagnostic disease, suggestive/subclinical findings, asymptomatic status,
   healthy/normal status and unknown assessment distinct. Require the definition
   used for an affected/unaffected field to travel with that assertion. Use a
   separate short Astra call only where the definition or join is unresolved.
   A second model's agreement is a review signal, not evidence.
3. Aggregate validated records in code. Count genotype-positive people within
   the intended study and variant group. Check affected/unaffected subsets
   against carriers and preserve unclassified people. Keep shared/co-carrier
   identities available to prevent double counting across representations.
   Derived values need their own reviewed evidence contract and trust path;
   they must not be disguised as literal published integers.
4. Keep raw output, source-validated roster, endpoint decisions, accepted sums,
   rejections and contradictions separately inspectable. Evaluate each boundary
   to distinguish a reading error from a deliberately conservative gate.

MYBPC3 20433692 illustrates why this matters. Its table footnote calls No Dx
"Unaffected or healthy," while prose separately identifies six healthy carriers
and four with suggestive ECG findings but no diagnostic HCM. A correct roster
can support different endpoint-defined totals. Neither a generic "No Dx means
normal" rule nor forcing all clinical evidence into a binary label is adequate.
Two noncarrier family cells are blank in the original table; main-text prose
supplies their H49 membership. Source joins need that evidence explicitly.

For the next opened-paper test, freeze the roster and endpoint contract before
scoring. Include a table with merged cells, a narrative family, an ambiguous
endpoint, a shared cohort, a count-free source and an unavailable supplement.
Compare accepted and rejected derived values against source adjudication as
well as the existing reference. Measure additional exact values, wrong supplied
values, positive-reference coverage, count-bearing identity extras, missing
identities, API cost per additional exact field, timeout rate and latency. Keep
previously disputed reference rows visible and report their influence.

Only after that contract passes opened calibration should a fixed implementation
enter the repository's next unopened discovery/confirmation sequence. No
percentage of human performance follows from this single-table diagnostic.

For the effort ablation suggested by the final Claude review, hold the bounded
source, roster schema, prompt and output cap constant and compare Astra low
with medium. Freeze both outputs before evaluating them. Measure row coverage,
source bindings, completeness and cost separately from accepted phenotype sums.
This isolates effort on that task; it does not itself explain the whole-paper
failures. Separate source-scope and fallback tests should change only their own
factor, with the same scorer and contract.
