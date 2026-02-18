#!/usr/bin/env python3
"""
Semi-Manual Paper Download Assistant

For the ~10% of papers that are truly CAPTCHA'd or paywalled,
this tool assists human-in-the-loop downloading via institutional proxy.

Requirements:
    pip install pandas playwright nest_asyncio
    playwright install chromium

Usage:
    python semi_manual_download.py --input papers_to_download.csv --output downloads/
"""

import os
import pandas as pd
import shutil
import time
from pathlib import Path
from playwright.async_api import async_playwright
import nest_asyncio
import asyncio
import argparse

nest_asyncio.apply()


class SemiManualDownloadAssistant:
    """
    Opens a browser for each paper, navigates to proxy URL,
    waits for manual download, then organizes the file.
    """
    
    def __init__(self, download_folder, output_folder='downloads'):
        self.download_folder = Path(download_folder).expanduser()
        self.output_folder = Path(output_folder)
        self.output_folder.mkdir(parents=True, exist_ok=True)
        self.results = []
        
        # Vanderbilt proxy URLs - customize for your institution
        self.proxy_urls = {
            'Wiley': 'https://onlinelibrary-wiley-com.proxy.library.vanderbilt.edu',
            'Elsevier': 'https://www-sciencedirect-com.proxy.library.vanderbilt.edu',
            'Springer': 'https://link-springer-com.proxy.library.vanderbilt.edu',
            'Oxford': 'https://academic-oup-com.proxy.library.vanderbilt.edu',
            'Nature': 'https://www-nature-com.proxy.library.vanderbilt.edu',
            'Cell': 'https://www-cell-com.proxy.library.vanderbilt.edu',
            'NEJM': 'https://www-nejm-org.proxy.library.vanderbilt.edu',
            'BMJ': 'https://www-bmj-com.proxy.library.vanderbilt.edu',
            'AHA': 'https://www-ahajournals-org.proxy.library.vanderbilt.edu',
        }
    
    def get_proxy_url(self, publisher, doi):
        """Build institutional proxy URL for a paper."""
        if publisher in self.proxy_urls and doi:
            return f"{self.proxy_urls[publisher]}/doi/{doi}"
        return None
    
    def wait_for_download(self, pmid, timeout=120):
        """Watch download folder for new PDF, rename and move it."""
        print(f"\n    Waiting for download...")
        start_time = time.time()
        initial_pdfs = set(self.download_folder.glob("*.pdf"))
        
        while time.time() - start_time < timeout:
            current_pdfs = set(self.download_folder.glob("*.pdf"))
            new_pdfs = current_pdfs - initial_pdfs
            
            if new_pdfs:
                new_pdf = list(new_pdfs)[0]
                time.sleep(3)  # Wait for download to complete
                dest = self.output_folder / f"PMID_{pmid}.pdf"
                shutil.move(str(new_pdf), str(dest))
                print(f"    ✓ Downloaded: PMID_{pmid}.pdf")
                return dest
            time.sleep(2)
        
        print(f"    ✗ Timeout")
        return None
    
    async def assist_download(self, df_papers):
        """Main loop: open browser, navigate, wait for download."""
        print(f"Semi-Manual Download Assistant")
        print(f"==============================")
        print(f"Papers to process: {len(df_papers)}")
        print(f"Output folder: {self.output_folder}")
        print(f"\nInstructions:")
        print(f"  1. Browser will open to paper URL")
        print(f"  2. Log in if prompted (VUNet for Vanderbilt)")
        print(f"  3. Click PDF download button")
        print(f"  4. Press ENTER when ready for next paper")
        print(f"  5. Type 'skip' to skip a paper")
        print(f"  6. Type 'quit' to exit\n")
        
        playwright = await async_playwright().start()
        browser = await playwright.chromium.launch(headless=False)
        context = await browser.new_context()
        page = await context.new_page()
        
        for idx, (i, paper) in enumerate(df_papers.iterrows()):
            print(f"\n[{idx+1}/{len(df_papers)}] PMID: {paper['pmid']}")
            print(f"    Title: {paper.get('title', 'N/A')[:70]}...")
            print(f"    Publisher: {paper.get('publisher', 'Unknown')}")
            
            result = {
                'pmid': paper['pmid'],
                'title': paper.get('title', ''),
                'publisher': paper.get('publisher', ''),
                'success': False,
                'skipped': False,
                'error': None
            }
            
            try:
                # Try proxy URL first, fallback to PubMed
                proxy_url = self.get_proxy_url(
                    paper.get('publisher', ''), 
                    paper.get('doi', '')
                )
                
                if proxy_url:
                    print(f"    URL: {proxy_url}")
                    await page.goto(proxy_url, timeout=60000)
                else:
                    pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{paper['pmid']}/"
                    print(f"    URL: {pubmed_url} (no proxy available)")
                    await page.goto(pubmed_url, timeout=60000)
                
                # Wait for user action
                downloaded = self.wait_for_download(paper['pmid'])
                result['success'] = downloaded is not None
                
            except Exception as e:
                result['error'] = str(e)
                print(f"    Error: {e}")
            
            self.results.append(result)
            
            # Save progress after each paper
            pd.DataFrame(self.results).to_csv(
                self.output_folder / 'download_progress.csv', 
                index=False
            )
            
            if idx < len(df_papers) - 1:
                user_input = input("    Press ENTER for next (or 'skip'/'quit'): ").strip().lower()
                if user_input == 'quit':
                    print("\nExiting early...")
                    break
                elif user_input == 'skip':
                    result['skipped'] = True
                    continue
        
        await browser.close()
        await playwright.stop()
        
        # Summary
        success_count = sum(r['success'] for r in self.results)
        skip_count = sum(r.get('skipped', False) for r in self.results)
        print(f"\n{'='*50}")
        print(f"COMPLETE!")
        print(f"  Downloaded: {success_count}/{len(self.results)}")
        print(f"  Skipped: {skip_count}")
        print(f"  Failed: {len(self.results) - success_count - skip_count}")
        print(f"  Progress saved to: {self.output_folder / 'download_progress.csv'}")
        
        return self.results


async def main():
    parser = argparse.ArgumentParser(description='Semi-manual paper download assistant')
    parser.add_argument('--input', '-i', required=True, help='CSV file with papers (needs pmid column)')
    parser.add_argument('--output', '-o', default='downloads', help='Output folder for PDFs')
    parser.add_argument('--downloads', '-d', default='~/Downloads', help='Browser downloads folder')
    args = parser.parse_args()
    
    # Load papers
    df = pd.read_csv(args.input)
    
    if 'pmid' not in df.columns:
        print("Error: Input CSV must have 'pmid' column")
        return
    
    # Setup and run
    assistant = SemiManualDownloadAssistant(
        download_folder=args.downloads,
        output_folder=args.output
    )
    await assistant.assist_download(df)


if __name__ == '__main__':
    asyncio.run(main())
