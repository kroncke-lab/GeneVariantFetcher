#!/usr/bin/env python3
"""
Reprocess papers that had Excel conversion failures.
Now that openpyxl is installed, we can convert Excel supplements properly.
"""

import sys
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from harvesting.format_converters import FormatConverter

OUTPUT_DIR = Path("/mnt/temp2/kronckbm/gvf_output/KCNH2/20260128_210249/pmc_fulltext")
converter = FormatConverter()


def reprocess_paper(pmid: str) -> int:
    """Reprocess a single paper's Excel supplements."""
    fulltext_path = OUTPUT_DIR / f"{pmid}_FULL_CONTEXT.md"
    supp_dir = OUTPUT_DIR / f"{pmid}_supplements"

    if not fulltext_path.exists() or not supp_dir.exists():
        return 0

    # Read current content
    content = fulltext_path.read_text()

    # Check if has Excel errors
    if "Error converting Excel file" not in content:
        return 0

    print(f"\nReprocessing PMID {pmid}...")

    # Find Excel files
    excel_files = list(supp_dir.glob("*.xlsx")) + list(supp_dir.glob("*.xls"))

    updated = 0
    for excel_file in excel_files:
        try:
            # Convert Excel to markdown
            markdown = converter.excel_to_markdown(excel_file)

            if markdown and len(markdown) > 50:
                # Find the error placeholder in content and replace
                error_marker = "# SUPPLEMENTAL FILE"
                filename = excel_file.name

                # Build new section
                new_section = f"# SUPPLEMENTAL FILE: {filename}\n\n{markdown}\n"

                # Replace error with actual content
                old_pattern = f"# SUPPLEMENTAL FILE [^#]*{filename}[^#]*\\[Error converting Excel file[^\\]]*\\]"
                import re

                content, count = re.subn(
                    old_pattern, new_section, content, flags=re.DOTALL
                )

                if count == 0:
                    # Try simpler replacement
                    if filename in content:
                        # Find line with filename and next error line
                        lines = content.split("\n")
                        new_lines = []
                        skip_next_error = False
                        for i, line in enumerate(lines):
                            if filename in line and "SUPPLEMENTAL" in line:
                                new_lines.append(f"# SUPPLEMENTAL FILE: {filename}")
                                new_lines.append("")
                                new_lines.append(markdown)
                                skip_next_error = True
                            elif skip_next_error and "Error converting" in line:
                                skip_next_error = False
                                continue
                            else:
                                new_lines.append(line)
                        content = "\n".join(new_lines)

                updated += 1
                print(f"  ✓ Converted {filename} ({len(markdown)} chars)")
        except Exception as e:
            print(f"  ✗ Failed {excel_file.name}: {e}")

    if updated > 0:
        # Write updated content
        fulltext_path.write_text(content)
        print(f"  Updated FULL_CONTEXT.md with {updated} Excel conversions")

    return updated


def main():
    # Find all papers with Excel errors
    failed_pmids = []
    for f in OUTPUT_DIR.glob("*_FULL_CONTEXT.md"):
        content = f.read_text()
        if "Error converting Excel file" in content:
            pmid = f.name.replace("_FULL_CONTEXT.md", "")
            failed_pmids.append(pmid)

    print(f"Found {len(failed_pmids)} papers with Excel conversion failures")

    total_updated = 0
    for pmid in failed_pmids:
        total_updated += reprocess_paper(pmid)

    print(f"\n{'=' * 50}")
    print(f"Total Excel files converted: {total_updated}")
    print(f"Papers updated: {len([p for p in failed_pmids if reprocess_paper(p) > 0])}")


if __name__ == "__main__":
    main()
