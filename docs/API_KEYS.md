# GeneVariantFetcher — API Keys Guide

How to obtain API keys for each service GVF uses, and what works without them.

## Quick Reference

| API | Required? | Free? | Coverage Impact |
|-----|-----------|-------|-----------------|
| **LLM provider** | ✅ Yes | Varies | Anthropic, OpenAI, or Azure AI is required for extraction |
| **Anthropic** | One LLM option | No (pay-per-use) | Claude extraction and vision-capable tiers |
| **OpenAI** | One LLM option | No (pay-per-use) | OpenAI extraction tiers |
| **Azure AI** | One LLM option | Varies | Azure-hosted model deployments |
| **NCBI (Email)** | ✅ Yes | Yes | Required for PubMed |
| **NCBI API Key** | Optional | Yes | 3x faster rate limits |
| **Elsevier API key** | Optional | Yes | Metadata/API access for ScienceDirect |
| **Elsevier Insttoken** | Recommended | Institutional | Highest-leverage unlock for subscription full text |
| **Springer** | Optional | Depends on access plan | Springer/Nature routes; article access varies |
| **Wiley** | Optional | Depends on TDM access | Wiley routes; subscription/access may be required |
| **PMC** | No publisher key | Yes | Available archive articles; not all PubMed records |

## Required: One LLM Provider Key

**Purpose:** Powers LLM-based variant extraction (the core of GVF).
GVF accepts Anthropic, OpenAI, or Azure AI credentials through LiteLLM. A key
does not select its provider: `MODEL_PROVIDER` must match the credential. The
shipped default is `anthropic`.

### Common Options

- Anthropic: create an API key in the Anthropic console and set `ANTHROPIC_API_KEY`.
- OpenAI: create an API key at [platform.openai.com](https://platform.openai.com/) and set `OPENAI_API_KEY`.
- Azure AI: configure your Azure AI Foundry endpoint and key with `AZURE_AI_API_KEY` and `AZURE_AI_API_BASE`.

### Configuration

```bash
# In your .env file
MODEL_PROVIDER=anthropic
ANTHROPIC_API_KEY=your-anthropic-key
# or:
# MODEL_PROVIDER=openai
# OPENAI_API_KEY=sk-your-key-here
# or:
# MODEL_PROVIDER=azure
# AZURE_AI_API_KEY=your-azure-ai-key
# AZURE_AI_API_BASE=https://your-resource.services.ai.azure.com

# Or as environment variable
export ANTHROPIC_API_KEY=your-anthropic-key
```

To use deployments on multiple Azure resources, export `AZURE_AI_MODEL_ROUTES`
as a JSON object in the process environment. For example:

```bash
export AZURE_AI_MODEL_ROUTES='{"gpt-6-astra":{"api_base":"https://second-resource.services.ai.azure.com/openai/v1","api_key_env":"SECOND_AZURE_API_KEY"}}'
```

Supply `SECOND_AZURE_API_KEY` through your existing secret-management environment.
The route contains the variable's name, never its credential value. These
arbitrary key variables must be exported; Settings does not load them from
`.env`. Unlisted models use the ordinary `AZURE_AI_API_BASE` and
`AZURE_AI_API_KEY`. A selected route with a missing key or malformed endpoint
fails before dispatch. Both chat and Azure Responses vision paths honor the
per-deployment route. Selecting a deployment does not change stage defaults.

### Illustrative Cost Range

| Gene Size | Papers Extracted | Approx Cost |
|-----------|-----------------|-------------|
| Small (20 papers) | ~15 | $0.50-2.00 |
| Medium (100 papers) | ~60 | $2.00-8.00 |
| Large (300 papers) | ~150 | $5.00-20.00 |

*These are rough historical planning ranges, not current provider pricing.
Verify deployed-model prices and measure actual token usage before budgeting a
large run.*

---

## Required: NCBI Email

**Purpose:** Compliance with NCBI E-Utilities terms of service

This is just your email address, not an API key. NCBI uses it to contact you if your usage causes problems.

### Configuration

```bash
# Provided via command line
gvf extract KCNH2 --email you@institution.edu --output ./results
```

---

## Recommended: NCBI API Key

**Purpose:** Increases rate limits from 3/sec to 10/sec

### How to Get It

1. Go to [ncbi.nlm.nih.gov](https://www.ncbi.nlm.nih.gov/)
2. Click **Log in** → **Register** if needed
3. Go to your account settings
4. Find **API Key Management**
5. Click **Create an API Key**

### Configuration

```bash
# In your .env file
NCBI_API_KEY=your-ncbi-api-key
```

---

## Recommended: Elsevier API Key + Insttoken

**Purpose:** Access ScienceDirect content (Cell, Lancet, many journals)

### How to Get It

1. Go to [dev.elsevier.com](https://dev.elsevier.com/)
2. Click **Get Started** → Create an account
3. Verify your email
4. Create a new application:
   - Name: "GeneVariantFetcher Research"
   - Description: "Academic text mining for genetic variant extraction"
5. Copy the **API Key**

### Requirements

- Must be affiliated with an institution that has Elsevier access
- Use limited to non-commercial research
- Must comply with their [Text and Data Mining policy](https://www.elsevier.com/about/policies/text-and-data-mining)

### Insttoken

The API key alone does not unlock many Vanderbilt-subscribed ScienceDirect
articles. Institutional access can require an `X-ELS-Insttoken` header,
configured as `ELSEVIER_INSTTOKEN`. Current measured blockers live in
[RECALL_STATUS.md](RECALL_STATUS.md), not this credential guide.

Obtain it through Vanderbilt library/e-resources support or an approved
ScienceDirect TDM workflow. Treat it like a credential.

### Configuration

```bash
# In your .env file
ELSEVIER_API_KEY=your-elsevier-key
ELSEVIER_INSTTOKEN=your-x-els-insttoken
```

### Coverage

The API key and institutional token enable configured Elsevier routes when
article access is available. Measure newly recovered bodies and supplements
on the actual PMID set; a credential does not imply a fixed coverage gain.

---

## Recommended: Springer Nature API Key

**Purpose:** Access Springer and Nature content

### How to Get It

1. Go to [dev.springernature.com](https://dev.springernature.com/)
2. Click **Sign Up** → Create account
3. Verify your institutional email
4. Go to **Applications** → **Create Application**
5. Fill out the form:
   - Application Name: "GVF Research"
   - Use Case: "Text mining for genetic variant research"
6. Copy the **API Key**

### Requirements

- Institutional affiliation preferred
- Non-commercial research use
- Rate limits apply (generous for academic use)

### Configuration

```bash
# In your .env file
SPRINGER_API_KEY=your-springer-key
```

### Coverage

Access depends on the article and enabled API product. Check the download
ledger and missing-supplement worklist for the measured effect on your cohort.

---

## Optional: Wiley API Key

**Purpose:** Access Wiley Online Library content

### How to Get It

1. Go to [onlinelibrary.wiley.com/library-info/resources/text-and-datamining](https://onlinelibrary.wiley.com/library-info/resources/text-and-datamining)
2. Review their TDM (Text and Data Mining) policy
3. Click **Request Access** or contact tdm@wiley.com
4. Provide:
   - Your institutional affiliation
   - Research purpose
   - Expected usage volume
5. Wiley will provide API credentials

### Requirements

- Institutional subscription typically required
- TDM agreement must be signed
- Non-commercial research only

### Configuration

```bash
# In your .env file
WILEY_API_KEY=your-wiley-key
```

### Coverage

Access depends on the article and TDM entitlement. Measure usable bodies and
supplements separately; publisher presence alone does not prove count coverage.

---

## Free: PubMed Central (PMC)

**No API key needed!**

PMC is a free full-text archive distinct from PubMed's citation database.
It includes participating-journal content, selected deposits and author
manuscripts; it does not contain every PubMed record or every open-access
article. See the [PMC FAQ](https://pmc.ncbi.nlm.nih.gov/about/faq/).

For automated retrieval, PMC provides article datasets through designated
services. Public readability is distinct from inclusion in the
[Open Access Subset](https://pmc.ncbi.nlm.nih.gov/tools/openftlist/).
Body availability also does not guarantee that the needed supplement is present.

---

## What Works Without Publisher Keys

GVF can discover candidate PMIDs, reuse cached sources, and try configured
public full-text routes. Extraction still requires an LLM provider. Missing
bodies or count-bearing supplements remain source gaps; neither discovery
completeness nor a fixed download percentage is guaranteed.

`--no-source-recovery` skips the later recovery pass. It does not restrict
ordinary harvesting to PMC. Start with the selected provider and NCBI email,
then use the actual blocked-paper worklist to decide which publisher access
would help. Keep provider-token usage and source availability as separate
measurements.

---

## Configuration Summary

Create a `.env` file in your GeneVariantFetcher directory:

```bash
# === REQUIRED ===
NCBI_EMAIL=your-email@example.org

# === REQUIRED: one LLM provider ===
ANTHROPIC_API_KEY=your-anthropic-key
# OPENAI_API_KEY=sk-your-openai-key
# AZURE_AI_API_KEY=your-azure-ai-key

# === RECOMMENDED ===
NCBI_API_KEY=your-ncbi-key
ELSEVIER_API_KEY=your-elsevier-key
ELSEVIER_INSTTOKEN=your-elsevier-insttoken
SPRINGER_API_KEY=your-springer-key

# === OPTIONAL ===
WILEY_API_KEY=your-wiley-key
CORE_API_KEY=your-core-key

# Private review-gold pipeline (project-owner provisioned)
GVF_REVIEW_GOLD_TOKEN=your-machine-token
GVF_REVIEW_GOLD_URL=https://variantbrowser.org/review/api/gold-standard/
```

`GVF_REVIEW_GOLD_TOKEN` is not an LLM/provider key. It is the shared
machine-to-machine credential for reading lead-approved gold from the private
Variant_Browser Azure review database. The same value is stored as the Azure App
Service setting and the GVF GitHub Actions secret; rotate both together.

---

## Troubleshooting

### "Invalid API key" errors

1. Check for extra whitespace in your `.env` file
2. Verify the key hasn't expired
3. Ensure you're using the correct key type (API key, not secret key where applicable)

### "Rate limit exceeded"

1. Add NCBI_API_KEY for higher limits
2. GVF has built-in rate limiting; this usually indicates a configuration issue
3. Wait a few minutes and retry

### "Unauthorized" from publisher APIs

1. Verify your institution has a subscription
2. Check if you need to be on VPN/campus network
3. Some APIs require IP registration — contact the publisher

### Papers still missing

1. Check `pmc_fulltext/paywalled_missing.csv` for blocked papers
2. Some papers require institutional access regardless of API keys
3. Very recent papers may not be available yet

---

## Measuring Cost and Coverage

Use the chosen run's source ledger and LLM traces. Report the number of attempted
papers, usable bodies, available count-bearing sections, model tokens and the
price basis/date. Credential combinations do not have a universal paper-coverage
percentage or dollar cost per gene. Dated observations remain in
[PROTOCOL_COST_EVAL.md](PROTOCOL_COST_EVAL.md); current validation and budget
constraints remain in [TASKS.md](../TASKS.md).

---

## Next Steps

- [QUICKSTART.md](QUICKSTART.md) — Get running with your keys
- [ARCHITECTURE.md](ARCHITECTURE.md) — Understand how GVF uses these APIs
- [OUTPUT_FORMAT.md](OUTPUT_FORMAT.md) — What GVF produces
