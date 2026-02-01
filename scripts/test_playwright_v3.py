#!/usr/bin/env python3
"""Verify we can extract actual content from AHA via Playwright."""

from playwright.sync_api import sync_playwright

print("Testing Playwright AHA content extraction...")

with sync_playwright() as p:
    browser = p.chromium.launch(headless=True)
    context = browser.new_context(
        user_agent='Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36'
    )
    page = context.new_page()
    
    page.goto('https://doi.org/10.1161/01.cir.102.10.1178', timeout=30000)
    
    print(f'Title: {page.title()}')
    print(f'URL: {page.url}')
    
    # Try to find abstract
    abstract = page.query_selector('.abstractSection')
    if abstract:
        text = abstract.inner_text()
        print(f'\n✅ ABSTRACT FOUND ({len(text)} chars):')
        print(text[:800])
        print('...')
    else:
        print('\n⚠ No .abstractSection found')
        # Try alternate selectors
        for selector in ['#abstract', '.abstract', '[data-widget-name="abstract"]']:
            elem = page.query_selector(selector)
            if elem:
                print(f'Found with {selector}: {elem.inner_text()[:200]}...')
                break
    
    # Check if full text is available
    fulltext = page.query_selector('.article__body')
    if fulltext:
        print(f'\n✅ FULL TEXT available ({len(fulltext.inner_text())} chars)')
    
    browser.close()
