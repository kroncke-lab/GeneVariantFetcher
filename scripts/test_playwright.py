#!/usr/bin/env python3
"""Test if Playwright works on this system."""

from playwright.sync_api import sync_playwright

print("Testing Playwright with AHA DOI...")

try:
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        
        # Test with an AHA DOI that was getting 403 with requests
        page.goto('https://doi.org/10.1161/01.cir.102.10.1178', timeout=30000)
        
        print(f'Final URL: {page.url}')
        print(f'Title: {page.title()}')
        content = page.content()
        print(f'Content length: {len(content)} chars')
        
        # Check if we got the actual article page
        if 'Circulation' in page.title() or 'ahajournals' in page.url:
            print('✓ Successfully reached AHA article page!')
        
        # Check for abstract
        if 'abstract' in content.lower():
            print('✓ Page contains abstract')
        
        browser.close()
        print('\n✅ Playwright is working!')
        
except Exception as e:
    print(f'\n❌ Playwright test failed: {e}')
