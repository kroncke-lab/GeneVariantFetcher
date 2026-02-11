"""
Format Converters Module

Converts various file formats (XML, Excel, Word, PDF) to markdown
for unified LLM processing.

Supports optional image extraction from PDFs for downstream vision analysis.
"""

import base64
import logging
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from bs4 import BeautifulSoup

logger = logging.getLogger(__name__)

try:
    from markitdown import MarkItDown

    MARKITDOWN_AVAILABLE = True
except ImportError:
    MARKITDOWN_AVAILABLE = False
    print("WARNING: markitdown not available. Will use basic conversion for docx/xlsx.")

try:
    import olefile

    OLEFILE_AVAILABLE = True
except ImportError:
    OLEFILE_AVAILABLE = False


class FormatConverter:
    """Converts various file formats to markdown."""

    def __init__(self):
        """Initialize format converter."""
        self.markitdown = MarkItDown() if MARKITDOWN_AVAILABLE else None

    def _convert_tsv_to_markdown_tables(self, text: str) -> str:
        """
        Convert tab-separated content to markdown tables where possible.

        This helps LLMs better understand tabular data extracted from documents.
        Only converts sections that look like actual tables (multiple columns, consistent structure).
        """
        lines = text.split("\n")
        result_lines = []
        in_table = False
        table_lines = []

        for line in lines:
            # Check if line has tabs (potential table row)
            if "\t" in line:
                cols = line.split("\t")
                # Consider it a table row if it has 2+ non-empty columns
                non_empty_cols = [c for c in cols if c.strip()]
                if len(non_empty_cols) >= 2:
                    if not in_table:
                        in_table = True
                        table_lines = []
                    table_lines.append(cols)
                    continue

            # Line without tabs or single column - end current table if any
            if in_table and table_lines:
                # Convert accumulated table to markdown
                markdown_table = self._table_lines_to_markdown(table_lines)
                result_lines.append(markdown_table)
                table_lines = []
                in_table = False

            result_lines.append(line)

        # Handle table at end of text
        if in_table and table_lines:
            markdown_table = self._table_lines_to_markdown(table_lines)
            result_lines.append(markdown_table)

        return "\n".join(result_lines)

    def _table_lines_to_markdown(self, table_lines: list) -> str:
        """Convert list of column lists to markdown table format."""
        if not table_lines:
            return ""

        # Find max columns and normalize all rows
        max_cols = max(len(row) for row in table_lines)

        # Build markdown table
        md_lines = []
        for i, row in enumerate(table_lines):
            # Pad row to max columns
            padded_row = row + [""] * (max_cols - len(row))
            # Clean up cells (remove excessive whitespace, escape pipes)
            cleaned_row = [cell.strip().replace("|", "\\|") for cell in padded_row]
            md_lines.append("| " + " | ".join(cleaned_row) + " |")

            # Add separator after first row (header)
            if i == 0:
                md_lines.append("|" + "|".join(["---"] * max_cols) + "|")

        return "\n".join(md_lines)

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
                title = "".join(title_elem.itertext()).strip()
                markdown += f"## {title}\n\n"

            abstract_elem = root.find(".//abstract")
            if abstract_elem is not None:
                markdown += "### Abstract\n\n"
                abstract_text = "".join(abstract_elem.itertext()).strip()
                markdown += f"{abstract_text}\n\n"

            def process_section(sec, level: int = 3):
                """Recursively process <sec> elements preserving hierarchy."""
                nonlocal markdown

                title_elem = sec.find("title")
                if title_elem is not None:
                    sec_title = "".join(title_elem.itertext()).strip()
                    markdown += f"{'#' * level} {sec_title}\n\n"

                for child in sec:
                    tag = child.tag.split("}")[-1]  # Handle optional namespaces
                    if tag == "p":
                        para_text = "".join(child.itertext()).strip()
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
                        para_text = "".join(p.itertext()).strip()
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
        soup = BeautifulSoup(html_content, "html.parser")
        markdown = "# MAIN TEXT\n\n"

        # Title
        title_elem = soup.find("h1")
        if title_elem:
            title = title_elem.get_text(strip=True)
            markdown += f"## {title}\n\n"

        # Abstract
        abstract_section = soup.find("section", class_=lambda c: c and "abstract" in c)
        if abstract_section:
            markdown += "### Abstract\n\n"
            for p in abstract_section.find_all("p"):
                text = p.get_text(" ", strip=True)
                if text:
                    markdown += f"{text}\n\n"

        def process_section(section, level: int = 3):
            """Recursively walk PMC <section> blocks."""
            nonlocal markdown
            heading = section.find(["h2", "h3", "h4"], recursive=False)
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
        body_section = soup.find(
            "section", class_=lambda c: c and "main-article-body" in c
        )
        if body_section:
            for sec in body_section.find_all("section", recursive=False):
                # Skip duplicate abstract sections we already handled
                sec_classes = sec.get("class", [])
                if any("abstract" == c for c in sec_classes):
                    continue
                process_section(sec)
        else:
            # Fallback: just pull all reasonably long paragraphs
            for p in soup.find_all("p"):
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

                # NOTE: Previously truncated to 100 rows which silently dropped variant data
                # from Excel supplements. Increased to 10000 to preserve full tables.
                # Fix applied 2026-02-10 for GVF pipeline accuracy.
                if len(df) > 10000:
                    markdown += (
                        f"*Note: Showing first 10000 rows of {len(df)} total rows*\n\n"
                    )
                    df_display = df.head(10000)
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
                text = "\n\n".join(
                    [para.text for para in doc.paragraphs if para.text.strip()]
                )
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
                print(
                    f"    Warning: markitdown failed for .doc file {file_path.name}: {e}"
                )

        # Try LibreOffice HTML conversion first - best for preserving table structure
        try:
            import subprocess
            import tempfile

            from bs4 import BeautifulSoup

            with tempfile.TemporaryDirectory() as tmpdir:
                result = subprocess.run(
                    [
                        "soffice",
                        "--headless",
                        "--convert-to",
                        "html",
                        "--outdir",
                        tmpdir,
                        str(file_path),
                    ],
                    capture_output=True,
                    text=True,
                    timeout=120,
                )
                if result.returncode == 0:
                    html_file = Path(tmpdir) / (file_path.stem + ".html")
                    if html_file.exists():
                        html = html_file.read_text(encoding="utf-8", errors="ignore")
                        if html.strip():
                            # Parse HTML and convert tables to markdown
                            soup = BeautifulSoup(html, "html.parser")
                            markdown_parts = []

                            # Extract all tables as markdown
                            for table in soup.find_all("table"):
                                rows = table.find_all("tr")
                                md_rows = []
                                for row in rows:
                                    cells = row.find_all(["td", "th"])
                                    cell_texts = [
                                        c.get_text(strip=True).replace("|", "/")
                                        for c in cells
                                    ]
                                    if any(cell_texts):
                                        md_rows.append(
                                            "| " + " | ".join(cell_texts) + " |"
                                        )

                                if md_rows:
                                    # Add separator after header
                                    if len(md_rows) > 1:
                                        num_cols = md_rows[0].count("|") - 1
                                        separator = "|" + "---|" * num_cols
                                        md_rows.insert(1, separator)
                                    markdown_parts.append("\n".join(md_rows))

                            # Also get non-table text
                            for table in soup.find_all("table"):
                                table.decompose()
                            body_text = soup.get_text(separator="\n", strip=True)

                            if markdown_parts:
                                result_text = (
                                    body_text + "\n\n" + "\n\n".join(markdown_parts)
                                )
                            else:
                                result_text = body_text

                            if result_text.strip():
                                print(
                                    f"    ✓ Extracted via LibreOffice HTML ({len(result_text)} chars, {len(markdown_parts)} tables)"
                                )
                                return result_text + "\n\n"
        except ImportError:
            # BeautifulSoup not available, fall back to text conversion
            try:
                import subprocess
                import tempfile

                with tempfile.TemporaryDirectory() as tmpdir:
                    result = subprocess.run(
                        [
                            "soffice",
                            "--headless",
                            "--convert-to",
                            "txt:Text",
                            "--outdir",
                            tmpdir,
                            str(file_path),
                        ],
                        capture_output=True,
                        text=True,
                        timeout=120,
                    )
                    if result.returncode == 0:
                        txt_file = Path(tmpdir) / (file_path.stem + ".txt")
                        if txt_file.exists():
                            text = txt_file.read_text(encoding="utf-8", errors="ignore")
                            if text.strip():
                                print(
                                    f"    ✓ Extracted text via LibreOffice ({len(text)} chars)"
                                )
                                return text + "\n\n"
            except Exception:
                pass
        except FileNotFoundError:
            pass  # LibreOffice not installed
        except Exception as e:
            print(
                f"    Warning: LibreOffice conversion failed for {file_path.name}: {e}"
            )

        # Try antiword as fallback (if installed)
        # First try with -t flag for tab-delimited output (better for tables)
        try:
            import subprocess

            # Try tab-delimited output first (better for tables)
            result = subprocess.run(
                ["antiword", "-t", str(file_path)],
                capture_output=True,
                text=True,
                timeout=60,
            )
            if result.returncode == 0 and result.stdout.strip():
                text = result.stdout
                # Convert tab-delimited tables to markdown format for better LLM extraction
                text = self._convert_tsv_to_markdown_tables(text)
                print(
                    f"    ✓ Extracted text via antiword (tab-delimited, {len(text)} chars)"
                )
                return text + "\n\n"

            # If -t flag fails or produces empty output, try standard output
            result = subprocess.run(
                ["antiword", str(file_path)], capture_output=True, text=True, timeout=60
            )
            if result.returncode == 0 and result.stdout.strip():
                print(f"    ✓ Extracted text via antiword ({len(result.stdout)} chars)")
                return result.stdout + "\n\n"
        except FileNotFoundError:
            pass  # antiword not installed
        except Exception as e:
            print(f"    Warning: antiword fallback failed for {file_path.name}: {e}")

        # Try catdoc as another fallback (if installed)
        try:
            import subprocess

            result = subprocess.run(
                ["catdoc", str(file_path)], capture_output=True, text=True, timeout=60
            )
            if result.returncode == 0 and result.stdout.strip():
                return result.stdout + "\n\n"
        except FileNotFoundError:
            pass  # catdoc not installed
        except Exception as e:
            print(f"    Warning: catdoc fallback failed for {file_path.name}: {e}")

        # Try OLE compound document parsing as final fallback

        # Try OLE compound document parsing (more reliable than raw bytes)
        if OLEFILE_AVAILABLE:
            try:
                ole = olefile.OleFileIO(str(file_path))
                text_parts = []

                # Try to extract text from WordDocument stream (main text storage)
                for stream_name in ole.listdir():
                    try:
                        stream_data = ole.openstream(stream_name).read()
                        # Word stores text as UTF-16LE in some streams
                        # Try to decode as UTF-16 first, then fallback to latin-1
                        decoded_text = None
                        for encoding in ["utf-16-le", "utf-16", "latin-1", "cp1252"]:
                            try:
                                decoded_text = stream_data.decode(
                                    encoding, errors="ignore"
                                )
                                # Filter to keep only printable text
                                cleaned = re.sub(
                                    r"[^\x09\x0A\x0D\x20-\x7E\u00A0-\u00FF]",
                                    " ",
                                    decoded_text,
                                )
                                cleaned = re.sub(r"\s+", " ", cleaned).strip()
                                # Only keep if it looks like meaningful text (has words)
                                if len(cleaned) > 50 and re.search(
                                    r"[a-zA-Z]{3,}", cleaned
                                ):
                                    text_parts.append(cleaned)
                                    break
                            except:
                                continue
                    except:
                        continue

                ole.close()

                if text_parts:
                    # Join all extracted text and clean up
                    full_text = "\n".join(text_parts)
                    # Remove duplicate content (OLE files often have redundant streams)
                    lines = full_text.split("\n")
                    seen = set()
                    unique_lines = []
                    for line in lines:
                        line_clean = line.strip()
                        if line_clean and line_clean not in seen:
                            seen.add(line_clean)
                            unique_lines.append(line_clean)

                    result = "\n".join(unique_lines)
                    if len(result) > 200:
                        print(f"    ✓ Extracted {len(result)} chars via OLE parsing")
                        return result + "\n\n"

            except Exception as e:
                print(f"    Warning: OLE parsing failed for {file_path.name}: {e}")

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
            normalized_lines = [
                line.strip() for line in cleaned.splitlines() if line.strip()
            ]
            merged = "\n".join(normalized_lines)
            if len(merged) > 200:
                return merged + "\n\n"
        except Exception as e:
            print(
                f"    Warning: heuristic .doc extraction failed for {file_path.name}: {e}"
            )

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
            with open(file_path, "rb") as f:
                header = f.read(8)
                if not header.startswith(b"%PDF"):
                    print(f"    Warning: {file_path.name} is not a valid PDF file")
                    return f"[Invalid PDF file: {file_path.name}]\n\n"
        except Exception as e:
            print(f"    Error reading PDF header {file_path}: {e}")
            return f"[Error reading PDF file: {file_path.name}]\n\n"

        # Try markitdown first (requires markitdown[pdf] extra for pdfminer.six)
        if self.markitdown:
            try:
                result = self.markitdown.convert(str(file_path))
                if (
                    result
                    and result.text_content
                    and len(result.text_content.strip()) > 100
                ):
                    print(
                        f"    ✓ Extracted PDF via markitdown ({len(result.text_content)} chars)"
                    )
                    return result.text_content
                elif result and result.text_content:
                    print(
                        f"    Warning: markitdown returned minimal content ({len(result.text_content)} chars)"
                    )
            except NameError as e:
                # Handle internal markitdown errors like 'excel_to_markdown' not defined
                print(
                    f"    Warning: markitdown internal error for {file_path.name}: {e}"
                )
            except ImportError as e:
                # Missing pdfminer.six - markitdown[pdf] not installed
                print(
                    f"    Warning: markitdown PDF support not installed (need markitdown[pdf]): {e}"
                )
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
        logger.debug(f"PDF conversion for {file_path} complete")
        return f"[PDF file available at: {file_path.name} - text extraction failed, manual review required]\n\n"

    def pdf_to_markdown_with_images(
        self,
        file_path: Path,
        output_dir: Optional[Path] = None,
        extract_images: bool = True,
        min_image_size: int = 100,
    ) -> Tuple[str, List[Dict]]:
        """
        Convert PDF to markdown AND extract images for vision analysis.

        Args:
            file_path: Path to PDF file
            output_dir: Directory to save extracted images. If None, images are
                       returned as base64 in the metadata.
            extract_images: If True, extract images from the PDF
            min_image_size: Minimum width/height in pixels to extract (filters icons)

        Returns:
            Tuple of (markdown_text, list of image metadata dicts)

            Image metadata dict contains:
            - page: Page number (1-indexed)
            - index: Image index on that page
            - width: Image width in pixels
            - height: Image height in pixels
            - size_bytes: Size of image data
            - ext: File extension (png, jpeg, etc.)
            - path: Path to saved image (if output_dir provided)
            - base64: Base64-encoded image data (if no output_dir)
        """
        extracted_images: List[Dict] = []

        # First verify the file is a valid PDF
        try:
            with open(file_path, "rb") as f:
                header = f.read(8)
                if not header.startswith(b"%PDF"):
                    logger.warning(f"{file_path.name} is not a valid PDF file")
                    return f"[Invalid PDF file: {file_path.name}]\n\n", []
        except Exception as e:
            logger.error(f"Error reading PDF header {file_path}: {e}")
            return f"[Error reading PDF file: {file_path.name}]\n\n", []

        # Try PyMuPDF (fitz) - best for image extraction
        try:
            import fitz  # PyMuPDF

            doc = fitz.open(str(file_path))
            text_content = []

            # Create figures directory if extracting images
            figures_dir = None
            if extract_images and output_dir:
                figures_dir = output_dir / "figures"
                figures_dir.mkdir(exist_ok=True)

            for page_num, page in enumerate(doc, 1):
                # Extract text
                text = page.get_text()
                if text.strip():
                    text_content.append(f"### Page {page_num}\n\n{text.strip()}")

                # Extract images
                if extract_images:
                    for img_index, img in enumerate(page.get_images(full=True)):
                        try:
                            xref = img[0]
                            base_image = doc.extract_image(xref)
                            image_bytes = base_image["image"]
                            image_ext = base_image.get("ext", "png")

                            # Get dimensions
                            width = base_image.get("width", 0)
                            height = base_image.get("height", 0)

                            # Skip tiny images (icons, bullets, decorations)
                            if width < min_image_size or height < min_image_size:
                                continue

                            # Skip very small file sizes (likely low-value)
                            if len(image_bytes) < 1000:
                                continue

                            img_filename = f"fig_p{page_num}_{img_index}.{image_ext}"
                            img_meta: Dict = {
                                "page": page_num,
                                "index": img_index,
                                "width": width,
                                "height": height,
                                "size_bytes": len(image_bytes),
                                "ext": image_ext,
                                "filename": img_filename,
                            }

                            if figures_dir:
                                img_path = figures_dir / img_filename
                                img_path.write_bytes(image_bytes)
                                img_meta["path"] = img_path
                                logger.debug(
                                    f"Extracted image: {img_filename} ({width}x{height})"
                                )
                            else:
                                # Store as base64 if no output dir
                                img_meta["base64"] = base64.b64encode(
                                    image_bytes
                                ).decode()

                            extracted_images.append(img_meta)

                            # Add placeholder in markdown for reference
                            text_content.append(
                                f"\n[Figure {len(extracted_images)} from page {page_num}: "
                                f"{img_filename} ({width}x{height}px)]\n"
                            )

                        except Exception as e:
                            logger.warning(
                                f"Failed to extract image {img_index} from page {page_num}: {e}"
                            )

            doc.close()

            if text_content:
                markdown = "\n\n".join(text_content) + "\n\n"
                if extracted_images:
                    logger.info(
                        f"Extracted {len(extracted_images)} images from {file_path.name}"
                    )
                return markdown, extracted_images

        except ImportError:
            logger.warning(
                "PyMuPDF not installed - falling back to text-only extraction"
            )
        except Exception as e:
            logger.warning(f"PyMuPDF extraction failed for {file_path.name}: {e}")

        # Fall back to text-only extraction
        text_only = self.pdf_to_markdown(file_path)
        return text_only, []

    def get_image_as_base64_url(self, image_path: Path) -> str:
        """
        Convert an image file to a base64 data URL for vision API calls.

        Args:
            image_path: Path to image file

        Returns:
            Data URL string (e.g., "data:image/png;base64,...")
        """
        image_bytes = image_path.read_bytes()
        b64 = base64.b64encode(image_bytes).decode()

        ext = image_path.suffix.lower()
        mime_type = {
            ".png": "image/png",
            ".jpg": "image/jpeg",
            ".jpeg": "image/jpeg",
            ".gif": "image/gif",
            ".webp": "image/webp",
            ".tiff": "image/tiff",
            ".bmp": "image/bmp",
        }.get(ext, "image/png")

        return f"data:{mime_type};base64,{b64}"
