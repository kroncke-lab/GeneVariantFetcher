import requests

EMAIL = "brett.kroncke@gmail.com"

def get_oa_url(doi):
    url = f"https://api.unpaywall.org/v2/{doi}?email={EMAIL}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        if data and data['is_oa']:
            if data['best_oa_location'] and data['best_oa_location']['url_for_pdf']:
                return data['best_oa_location']['url_for_pdf']
            elif data['best_oa_location']:
                return data['best_oa_location']['url']
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data for DOI {doi}: {e}")
    return None