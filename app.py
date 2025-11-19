import streamlit as st
import requests
import pandas as pd
import json
from io import StringIO
import plotly.express as px
import plotly.graph_objects as go
import os
from openai import OpenAI
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception
from concurrent.futures import ThreadPoolExecutor, as_completed

st.set_page_config(page_title="LitVar2 Gene Variant Search", layout="wide")

# OpenAI client setup using Replit AI Integrations
# the newest OpenAI model is "gpt-5" which was released August 7, 2025.
# do not change this unless explicitly requested by the user
AI_INTEGRATIONS_OPENAI_API_KEY = os.environ.get("AI_INTEGRATIONS_OPENAI_API_KEY")
AI_INTEGRATIONS_OPENAI_BASE_URL = os.environ.get("AI_INTEGRATIONS_OPENAI_BASE_URL")

openai_client = OpenAI(
    api_key=AI_INTEGRATIONS_OPENAI_API_KEY,
    base_url=AI_INTEGRATIONS_OPENAI_BASE_URL
)

st.title("NCBI LitVar2 Gene Variant Search")
st.markdown("""
Search for genetic variants and associated PubMed publications by gene name using the NCBI LitVar2 API.
""")

@st.cache_data(ttl=3600)
def search_gene_variants(gene_name):
    """
    Search for variants associated with a gene using NCBI LitVar2 API.
    Returns variant data from the search endpoint.
    Cached for 1 hour to improve performance.
    """
    base_url = "https://www.ncbi.nlm.nih.gov/research/bionlp/litvar/api/v1"
    search_url = f"{base_url}/entity/search/{gene_name}"
    
    try:
        response = requests.get(search_url, timeout=30)
        response.raise_for_status()
        data = response.json()
        return data, None
    except requests.exceptions.RequestException as e:
        return None, f"API Error: {str(e)}"
    except json.JSONDecodeError:
        return None, "Error: Unable to parse API response"

@st.cache_data(ttl=3600)
def get_pmids_for_rsids(rsids):
    """
    Get PMIDs for a list of rsIDs using the rsids2pmids endpoint.
    Batches requests to avoid URL length limits (max 100 rsIDs per request).
    Cached for 1 hour to improve performance.
    """
    if not rsids:
        return {}
    
    rsids = tuple(rsids)
    
    base_url = "https://www.ncbi.nlm.nih.gov/research/bionlp/litvar/api/v1"
    rsid_to_pmids = {}
    batch_size = 100
    
    for i in range(0, len(rsids), batch_size):
        batch = rsids[i:i + batch_size]
        rsids_param = ','.join(batch)
        pmids_url = f"{base_url}/public/rsids2pmids?rsids={rsids_param}"
        
        try:
            response = requests.get(pmids_url, timeout=60)
            response.raise_for_status()
            data = response.json()
            
            if isinstance(data, list):
                for item in data:
                    if isinstance(item, dict) and 'pmids' in item and 'rsid' in item:
                        rsid_to_pmids[item['rsid']] = item['pmids']
        except Exception as e:
            continue
    
    return rsid_to_pmids

def extract_variant_data(api_response, gene_name):
    """
    Extract variant information from the API response and fetch PMIDs.
    Returns a list of dictionaries containing variant details with PMIDs.
    """
    variants_data = []
    
    if not api_response or not isinstance(api_response, list):
        return variants_data
    
    rsids = []
    for item in api_response:
        if isinstance(item, dict) and 'rsid' in item:
            rsids.append(item['rsid'])
    
    rsid_to_pmids = get_pmids_for_rsids(rsids)
    
    for item in api_response:
        variant_info = extract_single_variant(item, rsid_to_pmids, gene_name)
        if variant_info:
            variants_data.append(variant_info)
    
    return variants_data

def extract_single_variant(variant_item, rsid_to_pmids, gene_name):
    """
    Extract information from a single variant item.
    """
    if not isinstance(variant_item, dict):
        return None
    
    variant_info = {}
    
    rsid = variant_item.get('rsid', 'N/A')
    variant_info['rsID'] = rsid
    
    variant_info['Variant_ID'] = variant_item.get('id', 'N/A')
    
    hgvs = variant_item.get('hgvs', '')
    if hgvs:
        variant_info['HGVS'] = hgvs
    else:
        variant_info['HGVS'] = 'N/A'
    
    hgvs_prot = variant_item.get('hgvs_prot', '')
    if hgvs_prot:
        variant_info['HGVS_Protein'] = hgvs_prot
    else:
        variant_info['HGVS_Protein'] = 'N/A'
    
    gene_data = variant_item.get('gene', {})
    if isinstance(gene_data, dict):
        variant_info['Gene'] = gene_data.get('name', gene_name)
    else:
        variant_info['Gene'] = gene_name
    
    pmids = rsid_to_pmids.get(rsid, [])
    variant_info['PMID_Count'] = len(pmids)
    variant_info['PMIDs'] = ', '.join(map(str, pmids)) if pmids else ''
    
    variant_info['Publications_Count'] = variant_item.get('pmids_count', 0)
    
    return variant_info

def get_all_pmids(variants_data):
    """
    Extract all unique PMIDs from the variants data.
    """
    all_pmids = set()
    for variant in variants_data:
        if variant.get('PMIDs'):
            pmid_list = variant['PMIDs'].split(', ')
            all_pmids.update(pmid_list)
    
    return sorted(list(all_pmids), key=lambda x: int(x) if x.isdigit() else 0)

def create_visualizations(variants_data):
    """
    Create visualization charts for variant distribution and publication trends.
    """
    df = pd.DataFrame(variants_data)
    
    fig_pmid_dist = px.histogram(
        df, 
        x='PMID_Count',
        title='Distribution of PMID Counts per Variant',
        labels={'PMID_Count': 'Number of PMIDs', 'count': 'Number of Variants'},
        nbins=30
    )
    fig_pmid_dist.update_layout(showlegend=False)
    
    top_variants = df.nlargest(15, 'PMID_Count')[['rsID', 'HGVS', 'PMID_Count']]
    top_variants['Label'] = top_variants.apply(
        lambda x: f"{x['rsID']} ({x['HGVS']})" if x['HGVS'] != 'N/A' else x['rsID'], 
        axis=1
    )
    
    fig_top_variants = px.bar(
        top_variants,
        x='PMID_Count',
        y='Label',
        orientation='h',
        title='Top 15 Variants by Publication Count',
        labels={'PMID_Count': 'Number of Publications', 'Label': 'Variant'}
    )
    fig_top_variants.update_layout(yaxis={'categoryorder': 'total ascending'})
    
    return fig_pmid_dist, fig_top_variants

@st.cache_data(ttl=3600)
def fetch_pubmed_abstracts(pmids):
    """
    Fetch PubMed article abstracts using NCBI E-utilities API.
    Returns a dictionary mapping PMID to article data (title + abstract).
    Cached for 1 hour to improve performance.
    """
    if not pmids:
        return {}
    
    pmids = tuple(pmids)
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    pmid_to_text = {}
    
    batch_size = 200
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i + batch_size]
        params = {
            'db': 'pubmed',
            'id': ','.join(batch),
            'retmode': 'xml',
            'rettype': 'abstract'
        }
        
        try:
            response = requests.get(base_url, params=params, timeout=60)
            response.raise_for_status()
            
            import xml.etree.ElementTree as ET
            root = ET.fromstring(response.content)
            
            for article in root.findall('.//PubmedArticle'):
                pmid_elem = article.find('.//PMID')
                if pmid_elem is not None:
                    pmid = pmid_elem.text
                    
                    title_elem = article.find('.//ArticleTitle')
                    title = title_elem.text if title_elem is not None and title_elem.text else ""
                    
                    abstract_parts = article.findall('.//AbstractText')
                    abstract = ' '.join([
                        (part.text if part.text else '') for part in abstract_parts
                    ])
                    
                    pmid_to_text[pmid] = {
                        'title': title,
                        'abstract': abstract,
                        'full_text': f"{title}\n\n{abstract}"
                    }
        except Exception as e:
            continue
    
    return pmid_to_text

def is_rate_limit_error(exception):
    """Check if the exception is a rate limit error."""
    error_msg = str(exception)
    return (
        "429" in error_msg
        or "RATELIMIT_EXCEEDED" in error_msg
        or "quota" in error_msg.lower()
        or "rate limit" in error_msg.lower()
        or (hasattr(exception, "status_code") and exception.status_code == 429)
    )

def extract_individuals_from_article(pmid, article_text, gene_name):
    """
    Extract individual-level variant and phenotype data from a PubMed article
    using biomedical text extraction with LLM.
    """
    extraction_prompt = f"""You are a biomedical text-extraction engine that reads scientific literature and outputs highly structured JSON records at the level of individual persons.

Your goal:
Given a PubMed article (PMID: {pmid}) about the gene {gene_name}, extract every variant-specific individual mentioned in the article, with classification of affected vs unaffected based on age-appropriate penetrance.

Article Text:
{article_text[:15000]}

Follow these extraction rules:

1. Identify all individuals mentioned (proband, case, patient, subject, family member, etc.)
2. Extract genetic variant(s) for each individual using HGVS notation when present
3. Extract phenotypes for each individual (convert to HPO codes when possible)
4. Determine affected vs unaffected status using clinical logic
5. Extract age, sex, ancestry, and timing information
6. Include sentence-level evidence for each extracted fact

Output Format:
Return a JSON array of individuals. Each individual should be a JSON object with these fields:
- individual_id: Unique ID like "PMID_{pmid}_case1"
- pmid: "{pmid}"
- gene: "{gene_name}"
- variants: Array of variant objects with hgvs_c, hgvs_p, rsid, genomic fields
- age: Number or null
- sex: "male", "female", or null
- phenotypes_hpo: Array of HPO codes
- phenotypes_raw: Array of raw text phenotypes
- affected_status: "affected", "unaffected", or "uncertain"
- evidence: Object with "sentence" and "section" fields
- data_from: "narrative", "table", or "supplement"

IMPORTANT: 
- Only output valid JSON
- If you cannot link a variant to a specific individual, do not create the entry
- Return empty array [] if no individuals can be extracted

Output only the JSON array, no other text."""

    @retry(
        stop=stop_after_attempt(3),
        wait=wait_exponential(multiplier=1, min=2, max=30),
        retry=retry_if_exception(is_rate_limit_error),
        reraise=True
    )
    def call_openai():
        # the newest OpenAI model is "gpt-5" which was released August 7, 2025.
        # do not change this unless explicitly requested by the user
        response = openai_client.chat.completions.create(
            model="gpt-5",
            messages=[{"role": "user", "content": extraction_prompt}],
            max_completion_tokens=8192,
            response_format={"type": "json_object"}
        )
        return response.choices[0].message.content or ""
    
    try:
        json_response = call_openai()
        individuals = json.loads(json_response)
        
        if isinstance(individuals, dict) and 'individuals' in individuals:
            return individuals['individuals']
        elif isinstance(individuals, list):
            return individuals
        else:
            return []
    except Exception as e:
        st.warning(f"Extraction error for PMID {pmid}: {str(e)}")
        return []

def batch_extract_individuals(pmids, gene_name, article_data, max_workers=2):
    """
    Process multiple PMIDs concurrently with rate limiting and automatic retries.
    """
    def process_pmid(pmid):
        if pmid not in article_data or not article_data[pmid].get('full_text'):
            return []
        
        article_text = article_data[pmid]['full_text']
        return extract_individuals_from_article(pmid, article_text, gene_name)
    
    all_individuals = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_pmid, pmid): pmid for pmid in pmids}
        for future in as_completed(futures):
            try:
                individuals = future.result()
                all_individuals.extend(individuals)
            except Exception as e:
                continue
    
    return all_individuals

st.sidebar.header("Search Parameters")

search_mode = st.sidebar.radio(
    "Search Mode",
    ["Single Gene", "Multiple Genes"],
    help="Choose to search one gene or multiple genes at once"
)

if search_mode == "Single Gene":
    gene_input = st.sidebar.text_input(
        "Enter Gene Symbol", 
        placeholder="e.g., BRCA1, TP53, CFTR",
        help="Enter the official gene symbol (HGNC nomenclature recommended)"
    )
else:
    gene_input = st.sidebar.text_area(
        "Enter Gene Symbols (one per line)", 
        placeholder="BRCA1\nTP53\nCFTR",
        help="Enter one gene symbol per line",
        height=100
    )

st.sidebar.markdown("---")
st.sidebar.subheader("Filter Options")

min_pmids = st.sidebar.number_input(
    "Min PMIDs per variant",
    min_value=0,
    value=0,
    help="Filter variants with at least this many PMIDs"
)

search_button = st.sidebar.button("Search", type="primary", use_container_width=True)

if search_button and gene_input:
    if search_mode == "Single Gene":
        genes_to_search = [gene_input.strip()]
    else:
        genes_to_search = [g.strip() for g in gene_input.strip().split('\n') if g.strip()]
    
    if not genes_to_search:
        st.warning("Please enter at least one gene symbol.")
    else:
        all_variants_data = []
        failed_genes = []
        
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        for idx, gene_name in enumerate(genes_to_search):
            status_text.text(f"Processing {gene_name} ({idx + 1}/{len(genes_to_search)})...")
            
            api_response, error = search_gene_variants(gene_name)
            
            if error:
                failed_genes.append((gene_name, error))
            elif not api_response:
                failed_genes.append((gene_name, "No variant data found"))
            else:
                variants_data = extract_variant_data(api_response, gene_name)
                if variants_data:
                    all_variants_data.extend(variants_data)
                else:
                    failed_genes.append((gene_name, "No variants found"))
            
            progress_bar.progress((idx + 1) / len(genes_to_search))
        
        progress_bar.empty()
        status_text.empty()
        
        if failed_genes:
            with st.expander(f"⚠️ {len(failed_genes)} gene(s) failed", expanded=False):
                for gene, reason in failed_genes:
                    st.warning(f"{gene}: {reason}")
        
        if all_variants_data:
            if min_pmids > 0:
                filtered_variants = [v for v in all_variants_data if v.get('PMID_Count', 0) >= min_pmids]
                if filtered_variants != all_variants_data:
                    st.info(f"Filtered from {len(all_variants_data)} to {len(filtered_variants)} variants (min {min_pmids} PMIDs)")
                variants_data = filtered_variants
            else:
                variants_data = all_variants_data
            
            if not variants_data:
                st.warning(f"No variants found matching the filter criteria (min {min_pmids} PMIDs)")
            else:
                gene_names_str = ", ".join(genes_to_search[:3]) + ("..." if len(genes_to_search) > 3 else "")
                if len(genes_to_search) == 1:
                    gene_display = genes_to_search[0]
                else:
                    gene_display = f"{len(genes_to_search)} genes ({gene_names_str})"
                
                st.success(f"Found {len(variants_data)} variant(s) across {gene_display}")
                gene_name = "_".join(genes_to_search) if len(genes_to_search) <= 3 else f"{len(genes_to_search)}_genes"
                
                all_pmids = get_all_pmids(variants_data)
                total_pmids = len(all_pmids)
                
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Total Variants", len(variants_data))
                with col2:
                    st.metric("Unique PMIDs", total_pmids)
                with col3:
                    total_associations = sum(v['PMID_Count'] for v in variants_data)
                    st.metric("Total Associations", total_associations)
                
                tab1, tab2, tab3, tab4 = st.tabs(["📊 Visualizations", "📋 Variant Details", "📚 PMID List", "🧬 Individual Extractions"])
                
                with tab1:
                    st.subheader("Data Visualizations")
                    fig_pmid_dist, fig_top_variants = create_visualizations(variants_data)
                    
                    col_viz1, col_viz2 = st.columns(2)
                    with col_viz1:
                        st.plotly_chart(fig_pmid_dist, use_container_width=True)
                    with col_viz2:
                        st.plotly_chart(fig_top_variants, use_container_width=True)
                    
                    st.markdown("""
                    **Insights:**
                    - The left chart shows how PMIDs are distributed across variants
                    - The right chart highlights the most-studied variants based on publication count
                    """)
                
                with tab2:
                    st.subheader("Variant Details")
                    df = pd.DataFrame(variants_data)
                    
                    df_display = df.copy()
                    df_display['dbSNP_Link'] = df_display['rsID'].apply(
                        lambda x: f"https://www.ncbi.nlm.nih.gov/snp/{x}" if x != 'N/A' else ''
                    )
                    
                    df_display['First_5_PMIDs'] = df_display['PMIDs'].apply(
                        lambda x: ', '.join(x.split(', ')[:5]) if x else ''
                    )
                    
                    column_order = ['Gene', 'dbSNP_Link', 'Variant_ID', 'HGVS', 'HGVS_Protein', 'PMID_Count', 'Publications_Count', 'First_5_PMIDs']
                    available_columns = [col for col in column_order if col in df_display.columns]
                    df_display_final = df_display[available_columns]
                    
                    st.dataframe(
                        df_display_final,
                        use_container_width=True,
                        height=400,
                        column_config={
                            "dbSNP_Link": st.column_config.LinkColumn(
                                "rsID",
                                help="Click to view variant in dbSNP database"
                            ),
                            "First_5_PMIDs": st.column_config.TextColumn(
                                "First 5 PMIDs",
                                help="First 5 PMIDs (see PMID List tab for clickable links to all PMIDs)"
                            ),
                        }
                    )
                    
                    st.markdown("**Clickable PMID Links (first 20 variants):**")
                    for idx, row in df.head(20).iterrows():
                        pmid_str = row.get('PMIDs')
                        if pmid_str and isinstance(pmid_str, str):
                            pmids = pmid_str.split(', ')[:5]
                            pmid_links = ' • '.join([f"[{p}](https://pubmed.ncbi.nlm.nih.gov/{p}/)" for p in pmids])
                            st.markdown(f"**{row['rsID']}**: {pmid_links}")
                    
                    csv = df[['Gene', 'rsID', 'Variant_ID', 'HGVS', 'HGVS_Protein', 'PMID_Count', 'Publications_Count', 'PMIDs']].to_csv(index=False)
                    st.download_button(
                        label="Download Variant Data (CSV)",
                        data=csv,
                        file_name=f"{gene_name}_variants.csv",
                        mime="text/csv"
                    )
                
                with tab3:
                    st.subheader("All Unique PMIDs")
                    st.info(f"Total unique PMIDs: {total_pmids}")
                    
                    if all_pmids:
                        pmid_links = [f"[{pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)" for pmid in all_pmids[:50]]
                        pmid_display = ', '.join(all_pmids)
                        
                        col_a, col_b = st.columns([2, 1])
                        with col_a:
                            st.text_area("PMID List (all)", pmid_display, height=150)
                        with col_b:
                            st.markdown("**Quick Links (first 50):**")
                            st.markdown(" • ".join(pmid_links), unsafe_allow_html=True)
                            
                            st.download_button(
                                label="Download PMIDs (TXT)",
                                data='\n'.join(all_pmids),
                                file_name=f"{gene_name}_pmids.txt",
                                mime="text/plain"
                            )
                
                with tab4:
                    st.subheader("Individual-Level Data Extraction")
                    st.markdown("""
                    Extract structured individual-level variant and phenotype data from PubMed articles using AI-powered biomedical text extraction.
                    """)
                    
                    st.info("This feature uses Replit AI Integrations (charged to your credits) to analyze PubMed articles and extract individual patient data including variants, phenotypes, and clinical status.")
                    
                    max_pmids_extract = st.number_input(
                        "Number of articles to extract (max 20)",
                        min_value=1,
                        max_value=20,
                        value=min(5, len(all_pmids)),
                        help="Processing uses AI credits - start small for testing"
                    )
                    
                    if st.button("🚀 Extract Individual Data", type="primary"):
                        pmids_to_process = all_pmids[:max_pmids_extract]
                        
                        with st.spinner(f"Fetching {len(pmids_to_process)} PubMed abstracts..."):
                            article_data = fetch_pubmed_abstracts(pmids_to_process)
                        
                        st.success(f"Fetched {len(article_data)} articles successfully")
                        
                        if article_data:
                            with st.spinner(f"Extracting individual-level data from {len(article_data)} articles... This may take a few minutes."):
                                progress_extraction = st.progress(0)
                                extracted_individuals = []
                                
                                for idx, pmid in enumerate(pmids_to_process):
                                    if pmid in article_data:
                                        individuals = extract_individuals_from_article(
                                            pmid, 
                                            article_data[pmid]['full_text'],
                                            genes_to_search[0] if len(genes_to_search) == 1 else gene_name
                                        )
                                        extracted_individuals.extend(individuals)
                                    
                                    progress_extraction.progress((idx + 1) / len(pmids_to_process))
                                
                                progress_extraction.empty()
                            
                            if extracted_individuals:
                                st.success(f"✅ Extracted {len(extracted_individuals)} individual records from {len(article_data)} articles")
                                
                                df_individuals = pd.DataFrame(extracted_individuals)
                                
                                if 'variants' in df_individuals.columns:
                                    df_individuals['variants_str'] = df_individuals['variants'].apply(
                                        lambda x: json.dumps(x) if isinstance(x, (list, dict)) else str(x)
                                    )
                                
                                if 'evidence' in df_individuals.columns:
                                    df_individuals['evidence_sentence'] = df_individuals['evidence'].apply(
                                        lambda x: x.get('sentence', '') if isinstance(x, dict) else ''
                                    )
                                
                                if 'phenotypes_raw' in df_individuals.columns:
                                    df_individuals['phenotypes_str'] = df_individuals['phenotypes_raw'].apply(
                                        lambda x: ', '.join(x) if isinstance(x, list) else str(x)
                                    )
                                
                                display_columns = []
                                for col in ['individual_id', 'pmid', 'gene', 'age', 'sex', 'affected_status', 
                                           'phenotypes_str', 'variants_str', 'evidence_sentence']:
                                    if col in df_individuals.columns:
                                        display_columns.append(col)
                                
                                st.dataframe(
                                    df_individuals[display_columns] if display_columns else df_individuals,
                                    use_container_width=True,
                                    height=400
                                )
                                
                                json_output = json.dumps(extracted_individuals, indent=2)
                                st.download_button(
                                    label="📥 Download Extraction Results (JSON)",
                                    data=json_output,
                                    file_name=f"{gene_name}_individual_extractions.json",
                                    mime="application/json"
                                )
                                
                                csv_output = df_individuals.to_csv(index=False)
                                st.download_button(
                                    label="📥 Download Extraction Results (CSV)",
                                    data=csv_output,
                                    file_name=f"{gene_name}_individual_extractions.csv",
                                    mime="text/csv"
                                )
                                
                                with st.expander("📊 Extraction Statistics"):
                                    col_stat1, col_stat2, col_stat3 = st.columns(3)
                                    with col_stat1:
                                        affected_count = len([i for i in extracted_individuals if i.get('affected_status') == 'affected'])
                                        st.metric("Affected Individuals", affected_count)
                                    with col_stat2:
                                        unique_variants = set()
                                        for i in extracted_individuals:
                                            if isinstance(i.get('variants'), list):
                                                for v in i['variants']:
                                                    if isinstance(v, dict) and v.get('hgvs_c'):
                                                        unique_variants.add(v['hgvs_c'])
                                        st.metric("Unique Variants", len(unique_variants))
                                    with col_stat3:
                                        articles_with_data = len(set([i.get('pmid') for i in extracted_individuals if i.get('pmid')]))
                                        st.metric("Articles with Data", articles_with_data)
                            else:
                                st.warning("No individual-level data could be extracted from the selected articles. This may occur if articles don't contain specific patient-level variant information.")
                        else:
                            st.error("Failed to fetch article abstracts from PubMed")
                    
                    with st.expander("ℹ️ About Individual Extraction"):
                        st.markdown("""
                        **What this extracts:**
                        - Individual patients/subjects with variant data
                        - Genetic variants in HGVS notation (cDNA and protein)
                        - Phenotypes with HPO codes when available
                        - Clinical status (affected/unaffected)
                        - Age, sex, and other demographics
                        - Evidence sentences supporting each extraction
                        
                        **Important notes:**
                        - Uses AI (GPT-5) to analyze article text
                        - Works best with case reports and detailed studies
                        - Limited to abstracts (full text not available via API)
                        - May not extract data from purely aggregate studies
                        - Costs are charged to your Replit credits
                        """)

elif search_button and not gene_input:
    st.warning("Please enter a gene symbol to search.")

st.sidebar.markdown("---")
st.sidebar.markdown("""
### About
This application uses the [NCBI LitVar2 API](https://www.ncbi.nlm.nih.gov/research/litvar2/) 
to retrieve genetic variant information and associated literature.

**Data Source:** NCBI/NLM/NIH  
**API:** LitVar2 REST API
""")

if not search_button:
    st.info("👈 Enter a gene symbol in the sidebar and click Search to begin.")
    
    with st.expander("ℹ️ How to use this tool"):
        st.markdown("""
        1. Enter a gene symbol (e.g., BRCA1, TP53, CFTR) in the sidebar
        2. Click the **Search** button
        3. View the results including:
           - Number of variants found
           - Associated PMIDs (PubMed IDs)
           - Variant details with rsIDs and HGVS notation
        4. Download the results as CSV or TXT files
        
        **Note:** The API searches for variants associated with the gene symbol 
        and returns literature references from PubMed.
        """)
