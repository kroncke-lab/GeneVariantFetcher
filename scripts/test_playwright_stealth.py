#!/usr/bin/env python3
"""Test Playwright with stealth settings to bypass bot detection."""

from playwright.sync_api import sync_playwright
import time

print("Testing Playwright with stealth settings...")

with sync_playwright() as p:
    browser = p.chromium.launch(
        headless=True,
        args=[
            '--disable-blink-features=AutomationControlled',
            '--disable-dev-shm-usage',
            '--no-first-run',
            '--disable-infobars',
        ]
    )
    
    context = browser.new_context(
        user_agent='Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36',
        viewport={'width': 1920, 'height': 1080},
        locale='en-US',
        timezone_id='America/New_York',
    )
    
    page = context.new_page()
    
    # Mask webdriver property
    page.add_init_script("""
        Object.defineProperty(navigator, 'webdriver', {
            get: () => undefined
        });
    """)
    
    page.goto('https://doi.org/10.1161/01.cir.102.10.1178', timeout=60000, wait_until='networkidle')
    
    # Wait for challenge if present
    for i in range(10):
        title = page.title()
        if 'Just a moment' not in title and 'Verifying' not in title:
            break
        print(f'  Waiting for challenge... ({i+1}s)')
        time.sleep(1)
    
    print(f'\nTitle: {page.title()}')
    print(f'URL: {page.url}')
    
    body = page.query_selector('body')
    text = body.inner_text() if body else ""
    print(f'Content: {len(text)} chars')
    
    if 'mutation' in text.lower() or 'lqts' in text.lower():
        print('✅ SUCCESS - Got actual article content!')
        print(f'\nFirst 500 chars:\n{text[:500]}')
    elif 'cloudflare' in text.lower() or 'verifying' in text.lower():
        print('❌ Still blocked by Cloudflare')
    else:
        print(f'⚠ Unknown result. First 300 chars:\n{text[:300]}')
    
    browser.close()
