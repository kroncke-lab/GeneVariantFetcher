#!/usr/bin/env python3
"""Quick test of AHA DOI resolution."""

import requests

session = requests.Session()
session.headers.update({
    'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36'
})

# Test AHA DOIs
test_dois = [
    '10.1161/01.cir.102.10.1178',
    '10.1161/hc0702.105124',
    '10.1161/CIRCEP.109.882159',
]

for doi in test_dois:
    url = f'https://doi.org/{doi}'
    print(f'\n{"="*60}')
    print(f'Testing: {url}')
    
    try:
        resp = session.get(url, allow_redirects=True, timeout=30)
        print(f'Status: {resp.status_code}')
        print(f'Final URL: {resp.url}')
        print(f'Content length: {len(resp.text)} chars')
        
        # Check if it's a full article page
        if 'full text' in resp.text.lower() or 'abstract' in resp.text.lower():
            print('✓ Looks like an article page')
        if 'sign in' in resp.text.lower() or 'login' in resp.text.lower():
            print('⚠ May require login')
        
        # Check for Cloudflare or bot detection
        if 'cloudflare' in resp.text.lower():
            print('⚠ Cloudflare protection detected')
        if 'captcha' in resp.text.lower():
            print('⚠ CAPTCHA required')
        if 'robot' in resp.text.lower() or 'bot' in resp.text.lower():
            print('⚠ Bot detection triggered')
            
        # Check for abstract in the 403 response
        if '<div class="abstractSection' in resp.text or 'Abstract' in resp.text:
            print('✓ Abstract may be available in response')
            # Try to extract abstract
            from bs4 import BeautifulSoup
            soup = BeautifulSoup(resp.text, 'html.parser')
            abstract = soup.find('div', class_='abstractSection')
            if abstract:
                print(f'Abstract found: {abstract.get_text()[:200]}...')
            
    except Exception as e:
        print(f'Error: {type(e).__name__}: {e}')
