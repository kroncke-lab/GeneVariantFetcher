import streamlit as st
import requests
import pandas as pd
import json
from io import StringIO
import plotly.express as px
import plotly.graph_objects as go

st.set_page_config(page_title="LitVar2 Gene Variant Search", layout="wide")

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
                
                tab1, tab2, tab3 = st.tabs(["📊 Visualizations", "📋 Variant Details", "📚 PMID List"])
                
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
                        if row.get('PMIDs'):
                            pmids = row['PMIDs'].split(', ')[:5]
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
