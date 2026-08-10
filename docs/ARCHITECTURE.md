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
│ STEP 3: Variant Extraction (ExpertExtractor)                                     │
│   • Input: DATA_ZONES.md > FULL_CONTEXT.md > abstract                            │
│   • Cheap deterministic census estimates ranges/risk for escalation only          │
│   • Fixed-width, vertical, and large Markdown tables have deterministic fast paths │
│   • An LLM routes other candidate tables; deterministic code reads their cells     │
│   • Tier 3 LLM full-text extraction fills gaps when table paths are insufficient   │
│   • High-risk results receive compact evidence-card verification/adjudication      │
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
│   • Cross-paper penetrance aggregation                                           │
│   OUTPUT: {gene}_penetrance_summary.json                                         │
└──────────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STEP 5: SQLite Migration                                                         │
│   • Normalized relational schema (8 core tables)                                 │
│   • Indexed for efficient querying                                               │
│   OUTPUT: {gene}.db                                                              │
└──────────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────────────────────────┐
│ STEPS 3.5–4.6 in gvf-run: recovery, trust, source QA, publication                │
│   • ClinVar / PubTator / figure recovery layers add attributable rows             │
│   • Optional count recovery fills only NULL fields (default off)                  │
│   • Deterministic trust gate creates trusted/quarantine field projections          │
│   • Per-paper reasoning review + composer are parked together (default off)        │
│   • Source QA/recovery, corpus sync, report and publish finish the run             │
└──────────────────────────────────────────────────────────────────────────────────┘
  │
  ▼
FINAL OUTPUTS:
  • {gene}.db                       (SQLite database)
  • {gene}_workflow_summary.json    (run statistics)
  • {gene}_penetrance_summary.json  (aggregated data)
  • run_manifest.json               (execution metadata)
```

## Module Responsibilities

### Core Pipeline Modules

| Module | Location | Responsibility |
|--------|----------|----------------|
| **CLI Entry** | `cli/__init__.py` | Typer app, command registration |
| **Automated Workflow** | `cli/automated_workflow.py` | End-to-end orchestration |
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
| **Format Converters** | `harvesting/format_converters.py` | XML/PDF/DOCX → Markdown |

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

### Extraction decision hierarchy

The word "tier" is used in two related places. Discovery tiers decide which
papers deserve extraction; extraction layers decide how facts are read from a
selected paper.

1. **Discovery and relevance.** Tier 1 is a fail-open keyword/regex screen.
   Tier 2 is an LLM relevance classifier. Supplying `--pmid-file`, as the fixed
   benchmarks do, intentionally skips discovery and both filters.
2. **Source resolution.** Cached consolidated context is preferred, then the
   harvester tries PMC/publisher/OA routes and folds convertible supplements.
   Data Scout may make a smaller high-value context, but the deterministic
   variant scanner always sees the original full context.
3. **Deterministic extraction first where structure is strong.** Regex scanning,
   PDF-linearized table reconstruction, large-Markdown-table parsing,
   fixed-width tables, and vertical tables run without an LLM. These paths
   normalize cells and counts using explicit rules.
4. **LLM-assisted routing, deterministic reading.** The table router classifies
   candidate tables and maps their columns. Code—not the routing model—then
   walks every selected row/cell. An empty, malformed, or unusable route falls
   through instead of suppressing Tier 3.
5. **Tier 3 full-text extraction.** A provider-routed model reads source context
   when deterministic/table routes are insufficient. Regex and router results
   are merged with its structured output. A deterministic paper census supplies
   only ranges and risk flags; it never becomes an extracted fact.
6. **Risk-gated judgment.** High-risk outputs are converted to compact evidence
   cards for a second-model verifier/adjudicator. This is not a second full-paper
   extraction on every paper. Production currently opens this lane for any risk
   score above threshold. The experimental, default-off
   `ENABLE_TIER3_REASON_CLASS_ROUTING` policy narrows that decision to
   count-semantic/precision reasons; it classifies completeness-only,
   missing-count, and source blockers for their appropriate recovery lanes.
   Those route labels are provenance today, not newly implemented recovery
   stages.
7. **Post-extraction layers.** JSON normalization/migration is followed by
   ClinVar, PubTator and figure recovery. Optional count recovery fills only
   NULL count fields with quote/role evidence and is currently off. The
   deterministic trust gate is on and masks questionable fields only in the
   trusted projection. The per-paper final review and its deterministic composer
   are parked together/off; source QA and recovery then finish the run.

### Model provider and reasoning effort

`config/settings.py` resolves each stage independently. Environment variables
can override any stage, so code defaults are not proof of what a particular run
used. The authoritative record is that run's `run_manifest.json`, which stores
the resolved models, reasoning efforts, enablement/policy switches, git state,
dependency hash, and a hash over behavior-defining extraction files. Exact
redacted requests, responses, repairs, fallbacks, and accepted decisions live in
the run's `llm_traces/` tree.

The active workstation routing recorded when the 2026-08-08 fixed-48 rerun
started is:

| Stage | Resolved route / state |
|-------|------------------------|
| Tier 1 relevance | deterministic keyword/regex; skipped for explicit benchmark PMIDs |
| Tier 2 relevance | `azure_ai/gpt-5.6-sol`; skipped for explicit benchmark PMIDs |
| Paper census | deterministic only |
| Table routing | `azure_ai/Kimi-K2.6-1` |
| Tier 3 full-text extraction | `azure_ai/grok-4.3` |
| High-risk Tier 3 adjudication | `azure_ai/gpt-5.6-sol` over compact evidence cards |
| Reason-class verifier routing | **off**; production retains the above-threshold rule |
| Early debate candidates | `azure_ai/gpt-5.6-sol`, `azure_ai/DeepSeek-V4-Pro`, `azure_ai/Kimi-K2.6-1` |
| Count recovery (Step 3.55) | `azure_ai/gpt-5.6-sol` at `high`, **off** |
| Trust gate (Step 3.7) | deterministic, **on** |
| Final paper review (Step 3.8) | `azure_ai/gpt-5.6-sol` at `xhigh`, **parked/off** |
| Final-review composer (Step 3.9) | deterministic, **parked/off** |

GPT-5.6 Luna and Terra were both connectivity-smoked successfully against the
configured Azure endpoint on 2026-08-08. They are registered as optional
deployment aliases, not production defaults. Luna is the first candidate for
high-volume, low-judgment work such as Tier 2 relevance and extraction-priority
triage; it should not be assigned count attribution, high-risk adjudication, or
final review until a fixed-set A/B shows no quality regression. Terra is the
intermediate candidate. Change one route at a time and score both the curated
fixture and the locked fixed-48 set before retaining it.

The fixed-48 run itself cannot evaluate Luna for Tier 2: an explicit PMID
manifest bypasses discovery and both relevance tiers. Its 279 calls instead
landed in Kimi table routing (13), Grok extraction (43), Sol source-grounded
claim verification (146), Sol figure reading (76), and Sol paper adjudication
(1). The two large Sol lanes require attribution or visual judgment, so the
first Luna experiment belongs on a separately labeled relevance/priority set,
not as an unmeasured swap inside this fixed-paper replay.

The paper census is deliberately approximate. It produces ranges for variant
rows, unique variants, carriers, affected, and unaffected counts plus risk flags
such as denominator-like columns or missing table bodies. These values are never
used as extracted facts. They only raise review priority, trigger targeted
fallback, and help the adjudication dashboard explain why a paper or claim was
escalated.

OpenAI-style reasoning models expose a reasoning-effort knob. GVF can set it per
pipeline stage through environment variables. Routine-stage efforts default to
unset; Step 3.8 is parked, but would use `xhigh` when explicitly enabled:

| Variable | Stage |
|----------|-------|
| `TIER2_REASONING_EFFORT` | Tier 2 relevance filtering in `pipeline/filters.py` |
| `TIER3_REASONING_EFFORT` | Tier 3 extraction and compact claim adjudication in `pipeline/extraction.py` |
| `ENABLE_TIER3_REASON_CLASS_ROUTING` | Experimental reason-class verifier gate (default `false`); completeness/missing-count/source risks are classified but do not open claim verification |
| `TABLE_ROUTER_REASONING_EFFORT` | Clinical table classification in `pipeline/table_router.py` |
| `VISION_REASONING_EFFORT` | Figure and pedigree extraction in `harvesting/figure_text_extractor.py`, `harvesting/figure_variant_reader.py`, and `pipeline/pedigree_extractor.py` |
| `PAPER_FINAL_CHECK_REASONING_EFFORT` | Parked final per-paper review in `pipeline/paper_final_check.py` (`xhigh` when enabled) |
| `PAPER_FINAL_CHECK_GATE_ENABLED` | Enable exact, source-quoted final-check fact/field trust composition (default `false` — parked with `PAPER_FINAL_CHECK_ENABLED`); weak/unsupported-only reasons remain advisory |
| `PAPER_FINAL_CHECK_GATE_MIN_SEVERITY` | Minimum applied severity (`high` by default) |
| `PAPER_FINAL_CHECK_GATE_REQUIRE_SOURCE_GROUNDED` | Keep DB-only flags advisory (default `true`) |
| `FINAL_ADJUDICATOR_REASONING_EFFORT` | Optional exception-adjudication queue when overridden to an OpenAI-style model |
| `FINAL_ARBITER_REASONING_EFFORT` | Optional hard-case queue when overridden to an OpenAI-style model |

Valid values are `none`, `minimal`, `low`, `medium`, `high`, and `xhigh`; `max`
is accepted as an alias for `xhigh`. Validation lives in `config/settings.py`.
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
.venv/bin/python scripts/smoke_azure_models.py \
  azure_ai/gpt-5.6-luna azure_ai/gpt-5.6-terra
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

Known current limitation: some extraction fallbacks still copy a generic
`patients.count` or carrier-only table total into `affected` when a variant is
disease-associated. That is not equivalent to an explicit phenotype split and
was visible in the 2026-08-08 fixed-48 replay (KCNQ1 PMID 18713323). Until the
tracked deterministic fix and legacy-row trust backstop are measured and
promoted, consumers should require an explicit affected label/provenance before
treating such an affected count as high-confidence.

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
- Optionally add synonyms to `config/cardiac_gene_synonyms.json` and a variant
  alias map at `data/{gene_lower}_variant_aliases.json`.

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
