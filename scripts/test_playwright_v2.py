#!/usr/bin/env python3
"""Test Playwright with wait for Cloudflare challenge."""

from playwright.sync_api import sync_playwright
import time

print("Testing Playwright with AHA DOI (with wait)...")

try:
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        context = browser.new_context(
            user_agent='Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36'
        )
        page = context.new_page()
        
        # Navigate to AHA DOI
        page.goto('https://doi.org/10.1161/01.cir.102.10.1178', timeout=30000)
        
        print(f'Initial URL: {page.url}')
        print(f'Initial Title: {page.title()}')
        
        # Wait for Cloudflare challenge to resolve (if present)
        if 'Just a moment' in page.title():
            print('Cloudflare challenge detected, waiting up to 15 seconds...')
            for i in range(15):
                time.sleep(1)
                title = page.title()
                print(f'  [{i+1}s] Title: {title}')
                if 'Just a moment' not in title:
                    print('  Challenge resolved!')
                    break
        
        print(f'\nFinal URL: {page.url}')
        print(f'Final Title: {page.title()}')
        content = page.content()
        print(f'Content length: {len(content)} chars')
        
        # Check what we got
        if 'cloudflare' in content.lower():
            print('⚠ Still on Cloudflare challenge page')
        elif 'abstract' in content.lower():
            print('✓ Found abstract in content')
            # Try to find the abstract text
            abstract_section = page.query_selector('.abstractSection')
            if abstract_section:
                print(f'Abstract: {abstract_section.inner_text()[:500]}...')
        
        browser.close()
        
except Exception as e:
    print(f'\n❌ Failed: {e}')
