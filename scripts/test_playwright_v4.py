#!/usr/bin/env python3
"""Full content check for AHA via Playwright."""

from playwright.sync_api import sync_playwright

print("Testing full AHA content extraction...")

with sync_playwright() as p:
    browser = p.chromium.launch(headless=True)
    page = browser.new_page()
    
    page.goto('https://doi.org/10.1161/01.cir.102.10.1178', timeout=30000)
    
    print(f'Title: {page.title()}')
    
    # Get all text content
    body = page.query_selector('body')
    full_text = body.inner_text() if body else ""
    
    print(f'\nTotal page text: {len(full_text)} chars')
    
    # Check for key content indicators
    checks = [
        ('Abstract', 'abstract' in full_text.lower()),
        ('Methods', 'method' in full_text.lower()),
        ('Results', 'result' in full_text.lower()),
        ('Discussion', 'discussion' in full_text.lower()),
        ('References', 'reference' in full_text.lower()),
        ('Mutations', 'mutation' in full_text.lower()),
        ('LQTS', 'lqts' in full_text.lower() or 'long-qt' in full_text.lower()),
    ]
    
    print('\nContent checks:')
    for name, found in checks:
        print(f'  {"✓" if found else "✗"} {name}')
    
    # Print a sample
    print(f'\nFirst 1500 chars of page text:')
    print('-' * 50)
    print(full_text[:1500])
    
    browser.close()
