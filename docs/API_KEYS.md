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
| **Springer** | Optional | Yes | +10-15% paper coverage |
| **Wiley** | Optional | Yes (TDM) | +5-10% paper coverage |
| **PMC** | No key needed | Yes | ~30% of papers |

## Required: One LLM Provider Key

**Purpose:** Powers LLM-based variant extraction (the core of GVF).
GVF accepts Anthropic, OpenAI, or Azure AI credentials through LiteLLM.

### Common Options

- Anthropic: create an API key in the Anthropic console and set `ANTHROPIC_API_KEY`.
- OpenAI: create an API key at [platform.openai.com](https://platform.openai.com/) and set `OPENAI_API_KEY`.
- Azure AI: configure your Azure AI Foundry endpoint and key with `AZURE_AI_API_KEY` and `AZURE_AI_API_BASE`.

### Configuration

```bash
# In your .env file
ANTHROPIC_API_KEY=your-anthropic-key
# or OPENAI_API_KEY=sk-your-key-here
# or AZURE_AI_API_KEY=your-azure-ai-key

# Or as environment variable
export ANTHROPIC_API_KEY=your-anthropic-key
```

### Cost Estimate

| Gene Size | Papers Extracted | Approx Cost |
|-----------|-----------------|-------------|
| Small (20 papers) | ~15 | $0.50-2.00 |
| Medium (100 papers) | ~60 | $2.00-8.00 |
| Large (300 papers) | ~150 | $5.00-20.00 |

*Costs vary based on paper length, provider, and extraction complexity.*

---

## Required: NCBI Email

**Purpose:** Compliance with NCBI E-Utilities terms of service

This is just your email address, not an API key. NCBI uses it to contact you if your usage causes problems.

### Configuration

```bash
# Provided via command line
gvf extract KCNH2 --email you@institution.edu
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
articles. The current recall blocker is the institutional `X-ELS-Insttoken`
header value, configured as `ELSEVIER_INSTTOKEN`.

Obtain it through Vanderbilt library/e-resources support or an approved
ScienceDirect TDM workflow. Treat it like a credential.

### Configuration

```bash
# In your .env file
ELSEVIER_API_KEY=your-elsevier-key
ELSEVIER_INSTTOKEN=your-x-els-insttoken
```

### Coverage

Elsevier publishes ~20% of biomedical literature. With the API key plus
institutional token, expect:
- +15-20% more papers downloaded
- Access to Cell, Lancet, American Journal of Human Genetics, many more

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

Springer Nature publishes Nature, Scientific Reports, European Journal of Human Genetics, and more:
- +10-15% more papers downloaded
- Critical for Nature family journals

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

Wiley publishes Human Mutation, Clinical Genetics, Genetic Epidemiology:
- +5-10% more papers
- Important for genetics-focused journals

---

## Free: PubMed Central (PMC)

**No API key needed!**

PMC provides free access to ~30% of PubMed articles. GVF uses this by default.

### What's in PMC?

- All NIH-funded research (after 12-month embargo)
- All articles from open access journals
- Author-deposited manuscripts

### Limitations

- Not all journals deposit in PMC
- Embargo periods mean recent papers may be missing
- Supplemental materials not always available

---

## What Works Without Keys

### PMC-Only Mode

Without any publisher keys, GVF still:
- ✅ Discovers all relevant PMIDs
- ✅ Downloads ~30% of papers (PMC open access)
- ✅ Extracts variants from available papers
- ✅ Creates SQLite database

### What You Miss

- ~70% of papers are behind paywalls
- Some high-impact journals (Cell, Nature) limited without keys
- Fewer variants extracted overall

### Recommendation

For best results, obtain at least:
1. **One LLM provider key** (required) — Anthropic (default), OpenAI, or Azure AI
2. **Elsevier API key + `ELSEVIER_INSTTOKEN`** (biggest current impact)
3. **Springer** (free, good coverage)

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
```

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

## Cost-Benefit Summary

| Configuration | Paper Coverage | Est. Cost/Gene |
|--------------|----------------|----------------|
| LLM key only | ~30% | $2-10 |
| + Elsevier | ~50% | $2-10 |
| + Springer | ~60% | $2-10 |
| + Wiley | ~65% | $2-10 |
| Full setup | ~70% | $2-10 |

*Paper coverage varies by gene. Cost is primarily LLM usage.*

---

## Next Steps

- [QUICKSTART.md](QUICKSTART.md) — Get running with your keys
- [ARCHITECTURE.md](ARCHITECTURE.md) — Understand how GVF uses these APIs
- [OUTPUT_FORMAT.md](OUTPUT_FORMAT.md) — What GVF produces
