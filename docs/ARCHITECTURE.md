# GeneVariantFetcher — Technical Architecture

A deep dive into GVF's pipeline architecture, module responsibilities, and extension points.

## Pipeline Overview

GVF is a multi-stage pipeline that transforms a gene symbol into a structured database of variants and patient data extracted from the biomedical literature.

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                         GeneVariantFetcher Pipeline                              │
└─────────────────────────────────────────────────────────────────────────────────┘

INPUT: Gene Symbol (e.g., "KCNH2")
  │
  ├──────────────────────────────────────────────────────────────────────────────┐
  │ STEP 0: Synonym Discovery (optional)                                         │
  │   • NCBI Gene Database → gene aliases                                         │
  │   • Example: KCNH2 → ["HERG", "LQT2", "Kv11.1"]                               │
  └──────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STEP 1: PMID Discovery                                                           │
│   • PubMind API → {gene}_pmids_pubmind.txt                                       │
│   • PubMed E-Utilities → {gene}_pmids_pubmed.txt                                 │
│   OUTPUT: {gene}_pmids.txt (merged, deduplicated)                                │
└──────────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STEP 1.5: Fetch Abstracts                                                        │
│   • PubMed E-Utilities (efetch)                                                  │
│   OUTPUT: abstract_json/{PMID}.json                                              │
└──────────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STEP 1.6: Filter Papers (Two-Tier)                                               │
│   Scope gate: deterministic non-human target-ortholog rejection (always on)      │
│   Tier 1: KeywordFilter (regex, fast) — ~65 clinical/variant keywords            │
│   Tier 2: InternFilter (LLM, provider-aware model) — relevance classification     │
│   OUTPUT: pmid_status/filtered_out.csv                                           │
└──────────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STEP 2: Download Full-Text (PMCHarvester)                                        │
│   Sources (priority order):                                                      │
│     1. PMC OA (Open Access) → BioC XML                                           │
│     2. Publisher APIs (Elsevier, Springer, Wiley) → XML/PDF                      │
│     3. Unpaywall → OA PDF links                                                  │
│   OUTPUT: pmc_fulltext/{PMID}_FULL_CONTEXT.md                                    │
└──────────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STEP 2.5: Data Scout (optional)                                                  │
│   • Identifies high-value data zones (tables, methods sections)                  │
│   • Creates condensed context for LLM extraction                                 │
│   OUTPUT: pmc_fulltext/{PMID}_DATA_ZONES.md                                      │
└──────────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STEP 3: Variant Extraction (ExpertExtractor; Azure staging route shown)           │
│   • Input: DATA_ZONES.md > FULL_CONTEXT.md > abstract                            │
│   • Cheap paper census estimates variant/count ranges for escalation only         │
│   • Kimi routes candidate tables; deterministic parser extracts table rows         │
│   • Grok 4.3 runs primary full-text extraction when table parsing is insufficient  │
│   • GPT-5.4 / DeepSeek / Kimi verify compact claim cards for high-risk outputs    │
│   • Pre-scan: Regex scanner on FULL_CONTEXT.md (not condensed text)              │
│     - Detects concatenated gene+variant (HERGG604S, KCNH2A561V)                  │
│     - Unicode arrow normalization (→, ➔, ⟶)                                      │
│   OUTPUT: extractions/{gene}_PMID_{pmid}.json                                    │
└──────────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STEP 4: Aggregation                                                              │
│   • HGVS variant name normalization                                              │
│   • Fuzzy matching for variant deduplication                                     │
│   • Cross-paper evidence aggregation; never derives penetrance from counts       │
│   OUTPUT: {gene}_penetrance_summary.json                                         │
└──────────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STEP 5: SQLite Migration                                                         │
│   • Normalized relational schema                                                 │
│   • Indexed for efficient querying                                               │
│   OUTPUT: {gene}.db                                                              │
└──────────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
FINAL OUTPUTS:
  • {gene}.db                       (SQLite database)
  • {gene}_workflow_summary.json    (run statistics)
  • {gene}_penetrance_summary.json  (aggregated data)
  • run_manifest.json               (execution metadata)
```

`gvf-run` is the current turnkey orchestrator. After migration it runs additive
recovery layers, Step 3.45 figure-count adoption, VariantFeatures enrichment and
false-positive quarantine, the default-on per-fact trust gate, source
QC/recovery and corpus sync, then reporting and optional review publication.
VariantFeatures remains warning-only for ordinary exploratory runs; collaborator
runs can make it mandatory with `--require-vf-enrich`, and publication always
requires both enrichment and quarantine. Publication also requires the exact
pinned PMID manifest and is refused after any failed stage. The older
`python -m cli.automated_workflow` entry point remains a lower-level
compatibility path; it is not registered as a `gvf` subcommand.

Explicit PMID manifests bypass recall-oriented Tier 1/Tier 2 filtering, but not
the deterministic paper-scope gate. A title that explicitly names a non-human
ortholog of the target gene is rejected before harvest/model calls. If full text
reveals the same condition later, extraction writes
`paper_scope_exclusion_reason=nonhuman_target_gene_ortholog`; SQLite migration
persists it, source replay refuses it even when forced, and ClinVar, PubTator,
and figure recovery exclude the PMID and purge legacy evidence links while
retaining the paper/extraction metadata as an audit record.

Penetrance and phenotype partitions are source claims, not arithmetic outputs.
Extraction retains a penetrance percentage only when a verbatim source quote
contains that exact percentage and names the variant. Aggregation never computes
`affected / carriers`, and claim verification never completes an
affected/unaffected partition from diagnosis, enrollment, or subtraction.

Variant-paper provenance separates origin from corroboration:
`variant_papers.source_layer` is one primary enum, while
`observed_source_layers` is the ordered, de-duplicated set of every lane that
observed the link. Consumers use the primary for diagnostic stratification and
the observed set for membership questions such as figure-count adoption.
Legacy comma-joined rows are read origin-first and normalized on replay.

Run-manifest provenance uses a versioned combined extractor digest. Its code
component hashes executable Python tokens (comments and cosmetic whitespace do
not churn it); its configuration component hashes the allowlisted resolved
model/routing/scientific knobs. The manifest records both components, the
algorithm and Python minor, plus the pre-v2 raw-byte digest for historical
auditability. Runtime alias and reference assets are loaded from the owned
`gvf_data` package through `importlib.resources`.

## Module Responsibilities

### Core Pipeline Modules

| Module | Location | Responsibility |
|--------|----------|----------------|
| **CLI Entry** | `cli/__init__.py` | Typer app, command registration |
| **Turnkey Orchestrator** | `cli/gvf_run.py` | Current end-to-end orchestration, recovery, trust, corpus sync, and reporting |
| **Legacy Workflow** | `cli/automated_workflow.py` | Lower-level compatibility orchestration |
| **Pipeline Steps** | `pipeline/steps.py` | Individual step implementations |

### Discovery & Filtering

| Module | Location | Responsibility |
|--------|----------|----------------|
| **Discovery** | `gene_literature/discovery.py` | PMID collection from multiple sources |
| **PubMind Fetcher** | `gene_literature/pubmind_fetcher.py` | PubMind variant-specific search |
| **PubMed Client** | `gene_literature/pubmed_client.py` | NCBI E-Utilities interface |
| **Synonym Finder** | `gene_literature/synonym_finder.py` | NCBI Gene synonym lookup |
| **Filters** | `pipeline/filters.py` | Keyword + LLM paper filtering |

### Harvesting (Full-Text Acquisition)

| Module | Location | Responsibility |
|--------|----------|----------------|
| **Orchestrator** | `harvesting/orchestrator.py` | Multi-source download coordination |
| **PMC API** | `harvesting/pmc_api.py` | PubMed Central BioC access |
| **Elsevier API** | `harvesting/elsevier_api.py` | ScienceDirect downloads |
| **Springer API** | `harvesting/springer_api.py` | SpringerLink + Nature downloads |
| **Wiley API** | `harvesting/wiley_api.py` | Wiley Online Library access |
| **Format Converters** | `harvesting/format_converters.py` | XML/PDF/DOCX → Markdown; dense PDF pages use a bounded PyMuPDF fast path before pdfminer layout sorting |

### Extraction & Analysis

| Module | Location | Responsibility |
|--------|----------|----------------|
| **Extraction** | `pipeline/extraction.py` | LLM-based variant extraction |
| **Prompts** | `pipeline/prompts.py` | LLM prompt templates |
| **Data Scout** | `pipeline/data_scout.py` | High-value zone identification |
| **Aggregation** | `pipeline/aggregation.py` | Cross-paper data combination |

### Utilities

| Module | Location | Responsibility |
|--------|----------|----------------|
| **Variant Normalizer** | `utils/variant_normalizer.py` | HGVS standardization |
| **Variant Scanner** | `utils/variant_scanner.py` | Regex-based variant detection (runs on full text, not condensed) |
| **LLM Utils** | `utils/llm_utils.py` | OpenAI/LiteLLM interface |
| **SQLite Migration** | `harvesting/migrate_to_sqlite.py` | JSON → SQLite conversion |

## Data Flow Detail

### Stage 1: PMID → Text

```
Gene Symbol
    │
    ├─→ PubMind API ──────────→ PMIDs (gene-variant focused)
    │
    ├─→ PubMed E-Utils ───────→ PMIDs (broader search)
    │
    └─→ Merge & Deduplicate ──→ {gene}_pmids.txt
                                    │
                                    ▼
                            Abstract Fetch (E-Utils)
                                    │
                                    ▼
                            abstract_json/{PMID}.json
                                    │
                                    ▼
                            Relevance Filtering
                            (Keyword + LLM tiers)
                                    │
                                    ▼
                            Filtered PMID list
```

### Stage 2: PMID → Full Text

```
Filtered PMIDs
    │
    ├─→ PMC OA Check ─────────→ BioC XML (if available)
    │
    ├─→ DOI Resolution ───────→ Publisher identification
    │       │
    │       ├─→ Elsevier API ─→ ScienceDirect XML
    │       ├─→ Springer API ─→ SpringerLink HTML/PDF
    │       └─→ Wiley API ────→ Wiley XML/PDF
    │
    └─→ Unpaywall ────────────→ OA PDF location
            │
            ▼
    Format Conversion (XML/PDF/HTML → Markdown)
            │
            ▼
    Supplement Processing (Excel, Word, PDFs)
            │
            ▼
    {PMID}_FULL_CONTEXT.md (unified document)
```

### Stage 3: Text → Variants

```
{PMID}_FULL_CONTEXT.md
    │
    ├─→ Variant Scanner ──────→ Pre-detected variant patterns
    │   (runs on FULL text)      - Concatenated: HERGG604S, KCNH2A561V
    │                            - Unicode arrows normalized (→, ➔, ⟶)
    │
    ├─→ Data Scout (optional)─→ {PMID}_DATA_ZONES.md
    │                           (condensed high-value sections)
    │
    └─→ ExpertExtractor ──────→ LLM extraction
            │                   (provider-aware model cascade)
            │
            ▼
    Structured JSON extraction
    {
      "paper_metadata": {...},
      "variants": [
        {
          "protein_notation": "p.Gly628Ser",
          "cdna_notation": "c.1883G>A",
          "clinical_significance": "pathogenic",
          "penetrance_data": {...},
          "functional_data": {...}
        }
      ],
      "extraction_metadata": {...}
    }
```

### Stage 4: Variants → Database

```
extractions/{gene}_PMID_*.json
    │
    ├─→ Variant Normalization ─→ Standardized names
    │   (frameshift, nonsense, protein notation)
    │
    ├─→ Fuzzy Matching ────────→ Variant deduplication
    │   (handles notation variants)
    │
    └─→ Aggregation ───────────→ {gene}_penetrance_summary.json
            │
            ▼
    SQLite Migration
            │
            ▼
    {gene}.db (normalized relational schema)
```

## External API Integration

### API Rate Limits & Authentication

| API | Auth Required | Rate Limit | Notes |
|-----|---------------|------------|-------|
| **PubMed E-Utils** | Email (key optional) | 3/sec (10/sec with key) | NCBI_EMAIL required |
| **PubMind** | None | Courteous use | Variant-focused search |
| **PMC OA** | None | 3/sec | Free full-text |
| **Elsevier** | API key | 5/sec | Requires registration |
| **Springer** | API key | Variable | Free for researchers |
| **Wiley** | API key | Variable | TDM agreement needed |
| **Unpaywall** | Email | 100k/day | OA link resolution |
| **LLM provider** | API key | Per plan | One provider required for extraction — Azure AI, Anthropic, or OpenAI |

### Configuration Surface

GVF has two authoritative configuration surfaces:

- `config/settings.py` contains typed, validated fields and provider/model
  resolution.
- `.env.example` documents direct environment reads in runtime modules and
  operator scripts, including defaults and behavioral consequences.

`tests/unit/test_env_documentation.py` scans shipped Python with the AST and
fails when a direct environment read appears in neither surface. This prevents
an undocumented stale shell value from silently changing paper selection,
source policy, recovery, or scoring acceptance.

### Model Provider And Reasoning Effort

`config/settings.py` resolves the effective model for each stage. The shipped
provider default is Anthropic. Set `MODEL_PROVIDER=azure` to use the measured
Azure staging route below; explicit per-tier model environment variables always
win. The lower-level `gvf extract` command also exposes `--model-provider`.
A separate final per-paper sniff test with GPT-5.6 Sol at `xhigh` (Step 3.8) is
retained but parked/default-off. When enabled, it records exact fact/field findings and does not
replace routine Tier 2 or mutate extracted counts. Step 3.9 deterministically
composes source-verified, objective contradictions into a field-level trusted
projection. Weak unsupported-count findings remain advisory. Sonnet/Opus are
reserved for optional exception-adjudication and
hard-case escalation over compact claim cards.

Measured Azure staging routing:

| Stage | Model |
|-------|-------|
| Tier 2 triage | `azure_ai/gpt-5.4` (`azure_ai/gpt-5.4-nano` only if deployed on the same endpoint) |
| Cheap paper census | deterministic regex/table/count pass; stored as `extraction_metadata.paper_census` |
| Table routing | `azure_ai/Kimi-K2.6-1`; falls back on empty/bad routes |
| Tier 3 extraction | `azure_ai/grok-4.3` |
| Internal claim QA/debate | `azure_ai/gpt-5.4`, `azure_ai/DeepSeek-V4-Pro`, `azure_ai/Kimi-K2.6-1` |
| Final per-paper sniff test (Step 3.8) | `azure_ai/gpt-5.6-sol` at `xhigh`; **PARKED — default off since 2026-07-26** (cost/latency not justified by measured effect). Code and tests retained; revive with `PAPER_FINAL_CHECK_ENABLED=1` *and* `PAPER_FINAL_CHECK_GATE_ENABLED=1` |
| Final-check trust composer (Step 3.9) | Deterministic; **PARKED with Step 3.8** (alone it can only refuse stale findings and fail acceptance). When enabled: source-verified objective count/phenotype contradictions only; weak unsupported-count findings stay advisory and raw counts stay unchanged |
| Optional exception-adjudication queue | `FINAL_ADJUDICATOR_MODELS` (`anthropic/claude-sonnet-5` by default) |
| Optional hard-case escalation | `FINAL_ARBITER_MODEL` (`anthropic/claude-opus-4-8` by default) |

The paper census is deliberately approximate. It produces ranges for variant
rows, unique variants, carriers, affected, and unaffected counts plus risk flags
such as denominator-like columns or missing table bodies. These values are never
used as extracted facts. They only raise review priority, trigger targeted
fallback, and help the adjudication dashboard explain why a paper or claim was
escalated.

OpenAI-style reasoning models expose a reasoning-effort knob. GVF can set it per
pipeline stage through environment variables. Routine-stage efforts default to
unset; the canonical Step 3.8 per-paper check defaults to `xhigh`:

| Variable | Stage |
|----------|-------|
| `TIER2_REASONING_EFFORT` | Tier 2 relevance filtering in `pipeline/filters.py` |
| `TIER3_REASONING_EFFORT` | Tier 3 extraction and compact claim adjudication in `pipeline/extraction.py` |
| `TABLE_ROUTER_REASONING_EFFORT` | Clinical table classification in `pipeline/table_router.py` |
| `VISION_REASONING_EFFORT` | Figure and pedigree extraction in `harvesting/figure_text_extractor.py`, `harvesting/figure_variant_reader.py`, and `pipeline/pedigree_extractor.py` |
| `PAPER_FINAL_CHECK_REASONING_EFFORT` | Parked final per-paper sniff test in `pipeline/paper_final_check.py` (`xhigh` when enabled) |
| `PAPER_FINAL_CHECK_GATE_ENABLED` | Enable exact, source-quoted final-check fact/field trust composition (default `false` — parked with `PAPER_FINAL_CHECK_ENABLED`); weak/unsupported-only reasons remain advisory |
| `PAPER_FINAL_CHECK_GATE_MIN_SEVERITY` | Minimum applied severity (`high` by default) |
| `PAPER_FINAL_CHECK_GATE_REQUIRE_SOURCE_GROUNDED` | Keep DB-only flags advisory (default `true`) |
| `FINAL_ADJUDICATOR_REASONING_EFFORT` | Optional exception-adjudication queue when overridden to an OpenAI-style model |
| `FINAL_ARBITER_REASONING_EFFORT` | Optional hard-case queue when overridden to an OpenAI-style model |

Valid values are `none`, `minimal`, `low`, `medium`, `high`, and `xhigh`; `max`
is accepted by GVF as an alias for `xhigh`. This compatibility mapping mattered
in the Luna shadow: the tested Azure deployment rejected literal `max` and
reported `xhigh` as its maximum accepted value. Validation lives in
`config/settings.py`.
Chat-completions calls use
`utils/llm_utils.build_reasoning_effort_kwargs`; Responses API calls use
`utils/llm_utils.build_responses_reasoning_param`. Both helpers no-op for models
without an OpenAI-style effort knob.

Treat reasoning effort as a secondary lever. Source acquisition, supplement
folding, extraction logic, and matcher behavior usually move recall more than a
non-default effort value. Change one stage at a time and re-score before keeping
the setting.

Before a long run, verify that the active Azure endpoint has the configured
deployment names:

```bash
.venv/bin/python scripts/smoke_azure_models.py --include-final
```

### LLM execution and decision traces

`gvf-run` configures `utils/llm_trace.py` for the **whole run lifetime** — not
just extraction — so post-extraction stages (count recovery, per-paper final
check, claim verification) are traced on every path, including `--skip extract`.
The recorder preserves exact textual requests, secret-safe parameters, provider
response envelopes, timing/usage, errors, and provider-exposed reasoning
*content* (token counts are reported separately and never imply exposed
reasoning). Scope is thread-local, which is what keeps the ~20 concurrent Tier-2
filter workers that share one filter instance from cross-attributing one paper's
evidence to another. Decision events carry `accepted_response_trace_id`,
`discarded_trace_ids` and `attempt_trace_links`, so retries, JSON repairs, parse
failures and fallbacks are distinguishable from the accepted response.

Run isolation matters because a `--skip extract` run reuses an older run's
directory: it writes to `llm_traces_<gvf run id>/` and its own report file, and a
manifest rebuild refuses to relabel another run's records. Integrity is reported
at one of three honest levels (`generated_now` / `write_time_verified` /
`locked`) by cross-checking each record against the write-time digest in
`trace_index.jsonl`.

The extraction-blinded paper evaluation applies a stronger integrity contract:
route, route-decision, extraction, and final-decision traces are required for
every paper, hashed into `llm_traces/trace_manifest.json`, and locked (together
with `trace_index.jsonl`) before gold scoring. `utils/llm_trace_html.py` renders
those same records into the self-contained `llm_trace_report.html` per-paper
browser timeline — with a route-coverage panel, bounded embedded bodies, and an
explicit omissions list; `scripts/build_llm_trace_html.py` is a thin CLI over it.
The benchmark hashes and locks that report before gold scoring as well. See
`docs/LLM_TRACING.md`.

### Compact count-claim verification

`pipeline/claim_verifier.py` verifies one count-bearing variant against a small
source evidence card and fails closed field by field. GPT-5.6 at `xhigh` gets a
64k output allowance because hidden reasoning consumes that budget before the
small JSON response; ordinary efforts/models retain the 2.5k budget. Evidence
routing expands one-letter/three-letter protein aliases and centers excerpts on
the variant mention when a converted paper collapses a long paragraph onto one
line. Generic `total - affected = unaffected` completion is disabled for
card-aware verification because nested cohorts can have genotyped people who
never entered phenotype follow-up.

The verifier currently runs inside `ExpertExtractor`, before downstream figure
and recovery layers settle. The next production change is to apply the same
compact, risk-ranked check after all layers merge; until that is measured, Luna
use here remains a shadow/evaluation route.

Adjudicated count gold is resolved centrally by `utils/gold_standard.py`.
`gold_v2_status` is an exact closed vocabulary: any unknown value raises rather
than silently dropping a row. A populated status makes all three `gold_v2_*`
fields authoritative, including explicit nulls, and exclusion states remove the
row consistently from the paper evaluator, general comparison scorer, and
claim-verification pilot.

### Additive count recovery (Step 3.55, default OFF)

Every other count stage is *subtractive*: `carrier_guard` NULLs implausible
counts, `count_outlier_guard` clears statistical outliers, `count_classifier`
refuses counts whose provenance is not per-variant. `pipeline/count_recovery.py`
is the only stage that **fills** a count the extractor never emitted — measured
on the locked 48-paper set, production left `total_carriers_observed` NULL on
58.4% of the gold variants it had correctly identified.

Contract:

- **Fill only.** It asks only about variants already stored for the paper, and
  writes only into NULL columns. An existing count is never overwritten.
- **Reasoning-safe budget.** GPT-5.6 `xhigh` calls receive a 64k output allowance;
  other efforts retain the 8192 recovery budget. A 2-paper Luna probe showed
  broad gap filling was low-yield on the fixed set (0/162 completed gaps), so
  this stage remains default OFF.
- **Paper-derived only.** ClinVar/PubTator linkage-only variants are excluded by
  default (`--include-linkage-rows` is an inspection mode), because asking for a
  count the paper never states is meaningless.
- **Quote-grounded, role-proven.** Every accepted value needs a verbatim quote
  from the source text the model actually saw, which locally binds the number to
  the requested variant *and* proves a per-variant role. Denominator grammars
  (`X of Y`, `X/Y`, `among Y`, `N cases of <phenotype>`), allele counts,
  population frequencies, family counts and measurements are rejected in both
  the prose and table branches; prose with two or more candidate integers fails
  closed unless the requested one carries an explicit carrier label locally. The
  model also declares a structured `count_role`, and a non-per-variant
  declaration is a veto.
- **Role into the trust path.** Accepted values write `count_role` +
  `evidence_locator` into `variant_papers.count_provenance` (merged, so a stale
  `cohort_total` cannot spuriously fire `count_is_total`) and land as
  `trust_tier='quarantine'`; trust-gate rule `recovered_count_unverified` (tg4)
  keeps any recovery-sourced count without a proven role quarantined. Step 3.55
  **refuses** to run when `--skip trust-gate` is in effect.
- **Durable and reversible.** The DB is copied to `*.before_count_recovery.db`,
  each paper commits under its own `SAVEPOINT`, and every write is logged to
  `count_recovery_log` with quote, role, locator, model and effort.
- **Auditable.** Each batch runs in a `count_recovery` trace scope and emits a
  `count_recovery_decision` event.

`COUNT_RECOVERY_ENABLED` stays **off** until a clean v2 measurement from
untouched baseline DB copies shows it improves count recall without an
unacceptable MAE or attribution regression (see `TASKS.md`).

### Gold-free trust gate (Step 3.7) and traceability

Before the LLM final check, a deterministic **gold-free** gate
(`pipeline/trust_gate.py`, rule generation **tg4**) soft-quarantines
structurally-implausible counts into a two-tier (`trusted` / `quarantine`)
projection. It **never** deletes a row or NULLs a raw value — quarantined fields
are only masked from the trusted projection the scorer/report read by default,
and can be re-tiered after a rule change. Rules are gene-class-agnostic so they
transfer from cardiac missense to BRCA case-control:

| Reason | Fires when |
|--------|-----------|
| `arith_inconsistent` | affected+unaffected+uncertain exceed the total, a fully-specified partition ≠ total, or any count is negative |
| `negative_count` | any carrier/affected/unaffected/uncertain < 0 |
| `count_is_total` | the "carrier" count is really a cohort denominator / screened-N |
| `population_count` | population allele magnitude (gnomAD/MAF label, or above the population ceiling) |
| `paper_outlier` | carrier count is a wild within-paper outlier (median × k) |
| `study_type_mismatch` | clinical carrier counts from a review/functional/GWAS paper |
| `implied_unaffected_zero` | a *derived* unaffected=0 (affected=total, zero not sourced) in an affirmatively cohort/biobank/case-control/family-cascade study — masks **only** the unaffected field; dormant on proband/unknown-design papers |
| `recovered_count_unverified` | a count written by Step 3.55 whose provenance does not name the per-variant role for the field it filled — recovery lands rows as `quarantine` on purpose, so this gate is what *promotes* them |

Each reason maps to the count fields it masks (`REASON_FIELDS`); most mask the
whole fact, `implied_unaffected_zero` masks only `unaffected`. The scorer reads
`penetrance_data.field_trust` per field (`cli/compare_variants.py`), falling back
to the row-level `trust_tier` when field masks are absent.

**Traceability.** `variant_papers.source_notation` stores the verbatim variant
string as written in each paper (before normalization), captured from the LLM
schema and the deterministic regex-table parser, so a curator can trace a
normalized ID back to the source and catch normalization errors.

**Design-aware discovery.** PubMed discovery includes a dedicated
penetrance/segregation/cascade/unaffected-carrier query lane
(`gene_literature/discovery.py`), and extraction priority scores that signal
(`pipeline/extraction_priority.py`), so carrier-first penetrance cohorts are not
out-prioritized by high-volume patient-first case series.

### Retry & Circuit Breaker

GVF implements exponential backoff with jitter:

```python
@dataclass
class RetryConfig:
    max_retries: int = 3
    base_delay: float = 1.0
    max_delay: float = 60.0
    exponential_base: float = 2.0
    jitter: bool = True
```

Circuit breakers prevent cascading failures:
- Skip extractions on files < 500 chars
- Skip files with < 30% alphanumeric content
- Timeout extraction after 1200 seconds

## Checkpoint/Resume System

GVF supports resumable jobs via the checkpoint system:

```python
# Pipeline steps (in order)
class PipelineStep(Enum):
    PENDING = "pending"
    DISCOVERING_SYNONYMS = "discovering_synonyms"
    FETCHING_PMIDS = "fetching_pmids"
    FETCHING_ABSTRACTS = "fetching_abstracts"
    FILTERING_PAPERS = "filtering_papers"
    DOWNLOADING_FULLTEXT = "downloading_fulltext"
    SCOUTING_DATA = "scouting_data"
    EXTRACTING_VARIANTS = "extracting_variants"
    AGGREGATING_DATA = "aggregating_data"
    MIGRATING_DATABASE = "migrating_database"
    COMPLETED = "completed"
    FAILED = "failed"
```

Checkpoints are stored in `~/.gvf_jobs/{job_id}/checkpoint.json` with atomic writes and file locking.

## Extension Points

### Adding a New Publisher API

1. Create `harvesting/{publisher}_api.py`:
   ```python
   class PublisherAPI:
       def can_handle(self, doi: str) -> bool:
           """Return True if this publisher can handle the DOI."""
           pass

       def download(self, doi: str, output_dir: Path) -> Optional[Path]:
           """Download and return path to content, or None on failure."""
           pass
   ```

2. Register in `harvesting/orchestrator.py`:
   ```python
   self.publishers = [
       ElsevierAPI(),
       SpringerAPI(),
       WileyAPI(),
       YourNewPublisherAPI(),  # Add here
   ]
   ```

### Adding a New Gene

There is no central gene-config file or `GENE_CONFIGS` dict — per-gene wiring is
intentionally minimal:

- Add the canonical protein length to `PROTEIN_LENGTHS` in
  `utils/variant_normalizer.py`.
- Add runtime aliases/query aliases to `BUILTIN_GENE_METADATA` in
  `utils/gene_metadata.py`, and optionally add a variant alias map at
  `gvf_data/{gene_lower}_variant_aliases.json`. The cardiac-synonyms JSON is a
  historical cold-start benchmark input, not the runtime registry.

See `docs/NEW_GENE_RUNBOOK.md` for the full add-a-gene flow.

### Customizing Extraction Prompts

Modify `pipeline/prompts.py`:
```python
EXTRACTION_PROMPT = """
Your custom extraction instructions here.
Focus on specific data types relevant to your use case.
"""
```

### Adding New Output Formats

Extend `harvesting/migrate_to_sqlite.py` or create a new exporter:
```python
def export_to_csv(db_path: Path, output_dir: Path):
    """Export SQLite database to CSV files."""
    pass
```

## Performance Considerations

### Typical Resource Usage

| Gene Size | Papers | RAM | Disk | Time |
|-----------|--------|-----|------|------|
| Small (rare) | 10-50 | 2 GB | 500 MB | 10-20 min |
| Medium | 50-200 | 4 GB | 2 GB | 30-60 min |
| Large (BRCA1) | 200-500 | 8 GB | 5 GB | 1-3 hours |

### Optimization Tips

1. **Parallel downloads**: Orchestrator fetches 3 papers concurrently by default
2. **SSD storage**: Significantly improves SQLite migration speed
3. **API keys**: Publisher keys unlock 2-3x more papers
4. **Data Scout**: `--scout-first` reduces LLM token usage for long papers

## Known Limitations

1. **Gene-specific tuning**: Best performance on cardiac genes (KCNH2, SCN5A, etc.)
2. **Variant aliases**: Only KCNH2 has a comprehensive alias dictionary
3. **PDF extraction**: Tables in scanned PDFs may not extract cleanly
4. **LLM hallucination**: Occasional false positives in extraction (validate critical findings)

## Further Reading

- [OUTPUT_FORMAT.md](OUTPUT_FORMAT.md) — Database schema and file formats
- [RECALL_STATUS.md](RECALL_STATUS.md) — Current recall baseline and blockers
- [API_KEYS.md](API_KEYS.md) — Obtaining API credentials
