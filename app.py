import streamlit as st
import requests
import pandas as pd
import json
from io import StringIO

st.set_page_config(page_title="LitVar2 Gene Variant Search", layout="wide")

st.title("NCBI LitVar2 Gene Variant Search")
st.markdown("""
Search for genetic variants and associated PubMed publications by gene name using the NCBI LitVar2 API.
""")

def search_gene_variants(gene_name):
    """
    Search for variants associated with a gene using NCBI LitVar2 API.
    Returns variant data from the search endpoint.
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

def get_pmids_for_rsids(rsids):
    """
    Get PMIDs for a list of rsIDs using the rsids2pmids endpoint.
    Batches requests to avoid URL length limits (max 100 rsIDs per request).
    """
    if not rsids:
        return {}
    
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

st.sidebar.header("Search Parameters")
gene_input = st.sidebar.text_input(
    "Enter Gene Symbol", 
    placeholder="e.g., BRCA1, TP53, CFTR",
    help="Enter the official gene symbol (HGNC nomenclature recommended)"
)

search_button = st.sidebar.button("Search", type="primary", use_container_width=True)

if search_button and gene_input:
    with st.spinner(f"Searching for variants in {gene_input}..."):
        gene_name = gene_input.strip()
        api_response, error = search_gene_variants(gene_name)
        
        if error:
            st.error(error)
        elif not api_response:
            st.warning(f"No variant data found for gene: {gene_name}")
        else:
            with st.spinner("Fetching PMID data for variants..."):
                variants_data = extract_variant_data(api_response, gene_name)
            
            if not variants_data:
                st.warning(f"No variants found for gene: {gene_name}")
            else:
                st.success(f"Found {len(variants_data)} variant(s) for gene: {gene_name}")
                
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
                
                st.subheader("Variant Details")
                df = pd.DataFrame(variants_data)
                column_order = ['Gene', 'rsID', 'Variant_ID', 'HGVS', 'HGVS_Protein', 'PMID_Count', 'Publications_Count', 'PMIDs']
                available_columns = [col for col in column_order if col in df.columns]
                df = df[available_columns]
                
                st.dataframe(df, use_container_width=True, height=400)
                
                st.subheader("All Unique PMIDs")
                st.info(f"Total unique PMIDs: {total_pmids}")
                
                if all_pmids:
                    pmid_display = ', '.join(all_pmids)
                    st.text_area("PMID List", pmid_display, height=150)
                    
                    st.download_button(
                        label="Download PMIDs (TXT)",
                        data='\n'.join(all_pmids),
                        file_name=f"{gene_name}_pmids.txt",
                        mime="text/plain"
                    )
                
                csv = df.to_csv(index=False)
                st.download_button(
                    label="Download Variant Data (CSV)",
                    data=csv,
                    file_name=f"{gene_name}_variants.csv",
                    mime="text/csv"
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
