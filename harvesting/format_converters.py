"""
Format Converters Module

Converts various file formats (XML, Excel, Word, PDF) to markdown
for unified LLM processing.
"""

from pathlib import Path
from typing import Optional
import xml.etree.ElementTree as ET
import pandas as pd
from bs4 import BeautifulSoup
import re

try:
    from markitdown import MarkItDown
    MARKITDOWN_AVAILABLE = True
except ImportError:
    MARKITDOWN_AVAILABLE = False
    print("WARNING: markitdown not available. Will use basic conversion for docx/xlsx.")


class FormatConverter:
    """Converts various file formats to markdown."""

    def __init__(self):
        """Initialize format converter."""
        self.markitdown = MarkItDown() if MARKITDOWN_AVAILABLE else None

    def xml_to_markdown(self, xml_content: str) -> str:
        """
        Convert PubMed Central XML to markdown.

        Args:
            xml_content: XML string content from PMC

        Returns:
            Markdown formatted text
        """
        try:
            root = ET.fromstring(xml_content)

            markdown = "# MAIN TEXT\n\n"

            title_elem = root.find(".//article-title")
            if title_elem is not None:
                title = ''.join(title_elem.itertext()).strip()
                markdown += f"## {title}\n\n"

            abstract_elem = root.find(".//abstract")
            if abstract_elem is not None:
                markdown += "### Abstract\n\n"
                abstract_text = ''.join(abstract_elem.itertext()).strip()
                markdown += f"{abstract_text}\n\n"

            def process_section(sec, level: int = 3):
                """Recursively process <sec> elements preserving hierarchy."""
                nonlocal markdown

                title_elem = sec.find("title")
                if title_elem is not None:
                    sec_title = ''.join(title_elem.itertext()).strip()
                    markdown += f"{'#' * level} {sec_title}\n\n"

                for child in sec:
                    tag = child.tag.split('}')[-1]  # Handle optional namespaces
                    if tag == "p":
                        para_text = ''.join(child.itertext()).strip()
                        if para_text:
                            markdown += f"{para_text}\n\n"
                    elif tag == "sec":
                        process_section(child, min(level + 1, 6))

            body_elem = root.find(".//body")
            if body_elem is not None:
                sections = body_elem.findall("./sec")

                if sections:
                    for sec in sections:
                        process_section(sec)
                else:
                    # Some PMC XML files use <p> directly under <body> without <sec> wrappers
                    for p in body_elem.findall(".//p"):
                        para_text = ''.join(p.itertext()).strip()
                        if para_text:
                            markdown += f"{para_text}\n\n"

            return markdown
        except Exception as e:
            print(f"  Error parsing XML: {e}")
            return "# MAIN TEXT\n\n[Error parsing XML content]\n\n"

    def pmc_html_to_markdown(self, html_content: str) -> str:
        """
        Convert PMC HTML (when XML lacks body text) to markdown.

        Args:
            html_content: HTML string from a PMC article page

        Returns:
            Markdown formatted text (best effort)
        """
        soup = BeautifulSoup(html_content, 'html.parser')
        markdown = "# MAIN TEXT\n\n"

        # Title
        title_elem = soup.find('h1')
        if title_elem:
            title = title_elem.get_text(strip=True)
            markdown += f"## {title}\n\n"

        # Abstract
        abstract_section = soup.find('section', class_=lambda c: c and 'abstract' in c)
        if abstract_section:
            markdown += "### Abstract\n\n"
            for p in abstract_section.find_all('p'):
                text = p.get_text(" ", strip=True)
                if text:
                    markdown += f"{text}\n\n"

        def process_section(section, level: int = 3):
            """Recursively walk PMC <section> blocks."""
            nonlocal markdown
            heading = section.find(['h2', 'h3', 'h4'], recursive=False)
            if heading:
                sec_title = heading.get_text(" ", strip=True)
                if sec_title:
                    markdown += f"{'#' * level} {sec_title}\n\n"

            for child in section.children:
                if getattr(child, "name", None) is None:
                    continue
                tag = child.name.lower()
                if tag == "p":
                    text = child.get_text(" ", strip=True)
                    if text:
                        markdown += f"{text}\n\n"
                elif tag in {"ul", "ol"}:
                    for li in child.find_all("li", recursive=False):
                        text = li.get_text(" ", strip=True)
                        if text:
                            markdown += f"- {text}\n"
                    markdown += "\n"
                elif tag == "section":
                    process_section(child, min(level + 1, 6))

        # Main body: PMC pages typically wrap text in section.body.main-article-body
        body_section = soup.find('section', class_=lambda c: c and 'main-article-body' in c)
        if body_section:
            for sec in body_section.find_all('section', recursive=False):
                # Skip duplicate abstract sections we already handled
                sec_classes = sec.get('class', [])
                if any('abstract' == c for c in sec_classes):
                    continue
                process_section(sec)
        else:
            # Fallback: just pull all reasonably long paragraphs
            for p in soup.find_all('p'):
                text = p.get_text(" ", strip=True)
                if len(text) > 50:
                    markdown += f"{text}\n\n"

        return markdown

    def excel_to_markdown(self, file_path: Path) -> str:
        """
        Convert Excel file to markdown tables.

        Args:
            file_path: Path to Excel file (.xlsx or .xls)

        Returns:
            Markdown formatted tables
        """
        try:
            xl_file = pd.ExcelFile(file_path)
            markdown = ""

            for sheet_name in xl_file.sheet_names:
                df = pd.read_excel(file_path, sheet_name=sheet_name)

                if df.empty:
                    continue

                markdown += f"#### Sheet: {sheet_name}\n\n"

                if len(df) > 100:
                    markdown += f"*Note: Showing first 100 rows of {len(df)} total rows*\n\n"
                    df_display = df.head(100)
                else:
                    df_display = df

                markdown += df_display.to_markdown(index=False)
                markdown += "\n\n"

            return markdown
        except Exception as e:
            print(f"    Error converting Excel {file_path}: {e}")
            return f"[Error converting Excel file: {e}]\n\n"

    def docx_to_markdown(self, file_path: Path) -> str:
        """
        Convert Word document to markdown.

        Args:
            file_path: Path to Word document (.docx)

        Returns:
            Markdown formatted text
        """
        if self.markitdown:
            try:
                result = self.markitdown.convert(str(file_path))
                return result.text_content
            except Exception as e:
                print(f"    Error converting DOCX with markitdown {file_path}: {e}")
                return f"[Error converting DOCX file: {e}]\n\n"
        else:
            try:
                from docx import Document
                doc = Document(file_path)
                text = "\n\n".join([para.text for para in doc.paragraphs if para.text.strip()])
                return text + "\n\n"
            except Exception as e:
                print(f"    Error converting DOCX {file_path}: {e}")
                return f"[Error converting DOCX file: {e}]\n\n"

    def doc_to_markdown(self, file_path: Path) -> str:
        """
        Convert legacy Word document (.doc) to markdown.

        Args:
            file_path: Path to Word document (.doc)

        Returns:
            Markdown formatted text
        """
        # MarkItDown can handle .doc files directly
        if self.markitdown:
            try:
                result = self.markitdown.convert(str(file_path))
                if result and result.text_content:
                    return result.text_content
            except Exception as e:
                print(f"    Warning: markitdown failed for .doc file {file_path.name}: {e}")

        # Try antiword as a fallback (if installed)
        try:
            import subprocess
            result = subprocess.run(
                ['antiword', str(file_path)],
                capture_output=True,
                text=True,
                timeout=60
            )
            if result.returncode == 0 and result.stdout.strip():
                return result.stdout + "\n\n"
        except FileNotFoundError:
            pass  # antiword not installed
        except Exception as e:
            print(f"    Warning: antiword fallback failed for {file_path.name}: {e}")

        # Try catdoc as another fallback (if installed)
        try:
            import subprocess
            result = subprocess.run(
                ['catdoc', str(file_path)],
                capture_output=True,
                text=True,
                timeout=60
            )
            if result.returncode == 0 and result.stdout.strip():
                return result.stdout + "\n\n"
        except FileNotFoundError:
            pass  # catdoc not installed
        except Exception as e:
            print(f"    Warning: catdoc fallback failed for {file_path.name}: {e}")

        # Try LibreOffice conversion as final fallback
        try:
            import subprocess
            import tempfile
            with tempfile.TemporaryDirectory() as tmpdir:
                result = subprocess.run(
                    ['soffice', '--headless', '--convert-to', 'txt:Text', '--outdir', tmpdir, str(file_path)],
                    capture_output=True,
                    text=True,
                    timeout=120
                )
                if result.returncode == 0:
                    txt_file = Path(tmpdir) / (file_path.stem + '.txt')
                    if txt_file.exists():
                        text = txt_file.read_text(encoding='utf-8', errors='ignore')
                        if text.strip():
                            return text + "\n\n"
        except FileNotFoundError:
            pass  # LibreOffice not installed
        except Exception as e:
            print(f"    Warning: LibreOffice fallback failed for {file_path.name}: {e}")

        # Lightweight heuristic extraction directly from the binary when no
        # conversion tools are available. Many legacy .doc files store
        # human-readable text alongside formatting bytes; decoding and
        # filtering printable ranges often recovers table contents well enough
        # for downstream parsing.
        try:
            raw_bytes = file_path.read_bytes()
            # Decode with a permissive codec, then drop non-printable noise
            decoded = raw_bytes.decode("latin-1", errors="ignore")
            cleaned = re.sub(r"[^\x09\x0A\x0D\x20-\x7E]", "\n", decoded)
            # Collapse excessive whitespace/newlines while preserving row-ish
            # structure so tables remain legible to the extractor.
            normalized_lines = [line.strip() for line in cleaned.splitlines() if line.strip()]
            merged = "\n".join(normalized_lines)
            if len(merged) > 200:
                return merged + "\n\n"
        except Exception as e:
            print(f"    Warning: heuristic .doc extraction failed for {file_path.name}: {e}")

        # Final fallback - indicate manual review needed
        return f"[Legacy .doc file available at: {file_path.name} - text extraction failed, manual review required]\n\n"

    def pdf_to_markdown(self, file_path: Path) -> str:
        """
        Convert PDF to markdown.

        Args:
            file_path: Path to PDF file

        Returns:
            Markdown formatted text or placeholder
        """
        # First verify the file is a valid PDF
        try:
            with open(file_path, 'rb') as f:
                header = f.read(8)
                if not header.startswith(b'%PDF'):
                    print(f"    Warning: {file_path.name} is not a valid PDF file")
                    return f"[Invalid PDF file: {file_path.name}]\n\n"
        except Exception as e:
            print(f"    Error reading PDF header {file_path}: {e}")
            return f"[Error reading PDF file: {file_path.name}]\n\n"

        # Try markitdown first
        if self.markitdown:
            try:
                result = self.markitdown.convert(str(file_path))
                if result and result.text_content:
                    return result.text_content
            except NameError as e:
                # Handle internal markitdown errors like 'excel_to_markdown' not defined
                print(f"    Warning: markitdown internal error for {file_path.name}: {e}")
            except Exception as e:
                print(f"    Warning: markitdown failed for {file_path.name}: {e}")

        # Try PyMuPDF (fitz) as fallback if available
        try:
            import fitz  # PyMuPDF
            doc = fitz.open(str(file_path))
            text_content = []
            for page_num, page in enumerate(doc, 1):
                text = page.get_text()
                if text.strip():
                    text_content.append(f"### Page {page_num}\n\n{text.strip()}")
            doc.close()
            if text_content:
                return "\n\n".join(text_content) + "\n\n"
        except ImportError:
            pass  # PyMuPDF not installed
        except Exception as e:
            print(f"    Warning: PyMuPDF fallback failed for {file_path.name}: {e}")

        # Try pdfplumber as another fallback
        try:
            import pdfplumber
            text_content = []
            with pdfplumber.open(str(file_path)) as pdf:
                for page_num, page in enumerate(pdf.pages, 1):
                    text = page.extract_text()
                    if text and text.strip():
                        text_content.append(f"### Page {page_num}\n\n{text.strip()}")
            if text_content:
                return "\n\n".join(text_content) + "\n\n"
        except ImportError:
            pass  # pdfplumber not installed
        except Exception as e:
            print(f"    Warning: pdfplumber fallback failed for {file_path.name}: {e}")

        # Final fallback - just indicate the file exists
        return f"[PDF file available at: {file_path.name} - text extraction failed, manual review required]\n\n"
