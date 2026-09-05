"""Bounded, metadata-bound recovery of public article copies.

Used by the initial harvester and paywall recovery, including papers PubMed
marks free whose publisher URL fails. No PMID, gene, DOI or PDF URL exceptions.
Repository identifiers need not contain a DOI: exact DOI index metadata plus
the downloaded article's title are the two independent identity witnesses.
This recovers bodies, not a claim that all supplements/tables are complete.
"""

from __future__ import annotations

import hashlib
import ipaddress
import json
import re
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote, urljoin, urlparse

import fitz
import requests
from bs4 import BeautifulSoup

from harvesting.browser_html.content_quality import validate_article_content
from harvesting.content_validation import validate_content_quality
from harvesting.scholar_pdf_fallback import _compact_identity_text, _normalized_doi
from harvesting.unpaywall_api import UnpaywallClient
from pipeline.source_quality import SUPPLEMENT_SURFACE_STATUS_MARKER
from utils.doi import DOI_RE, clean_doi

MAX_BYTES = 30 * 1024 * 1024
MAX_ASSET_REQUESTS = 14


def _title_keys(title: str) -> set[str]:
    """Some indexes omit a subtitle; allow only an exact, specific main title."""
    keys = {_compact_identity_text(title)}
    if ":" in title:
        main = _compact_identity_text(title.split(":", 1)[0])
        if len(main) >= 40:
            keys.add(main)
    return keys


@dataclass
class RepositoryPDFResult:
    success: bool = False
    pmid: str = ""
    doi: str = ""
    title: str = ""
    markdown: str = ""
    source_url: str = ""
    pdf_path: str = ""
    pdf_sha256: str = ""
    page_count: int = 0
    identity_reason: str = ""
    candidate: dict = field(default_factory=dict)
    attempts: list[dict] = field(default_factory=list)
    error: str = ""


def _public_url(url: str) -> bool:
    parsed = urlparse(url)
    host = parsed.hostname or ""
    if parsed.scheme not in {"http", "https"} or not host or parsed.username:
        return False
    if host == "localhost" or host.endswith((".local", ".localhost")):
        return False
    try:
        return ipaddress.ip_address(host).is_global
    except ValueError:
        return "." in host


def article_identity(pages: list[str], title: str, doi: str) -> tuple[bool, str]:
    """Match the article header, not its bibliography or a repository cover."""
    expected = _title_keys(title)
    if not any(len(key) >= 20 for key in expected):
        return False, "no sufficiently specific title from identifier metadata"
    for page in pages[:3]:
        # Repository citation/licence covers identify a record, not its attached
        # article. Verify the article page as the independent witness.
        low = page.lower()
        if (
            ("to cite this version" in low and "hal id:" in low)
            or "this is the peer reviewed version of" in low
            or ("document version" in low and "link to publication" in low)
            or ("general rights" in low and "download date" in low)
        ):
            continue
        # Do not split on the word 'Background' inside an article title.
        heading = re.search(
            r"(?im)^\s*(?:abstract|introduction|references|bibliography)(?:\s*$|\s*[:—–-])",
            page,
        )
        header = page[: min(heading.start() if heading else 3500, 3500)]
        if not any(key in _compact_identity_text(header) for key in expected):
            continue
        # Explicit different article DOI defeats a same-title/version match.
        found = [clean_doi(match.group(0)).lower() for match in DOI_RE.finditer(header)]
        if found and doi not in found:
            return False, "article header has a conflicting DOI"
        return True, "exact DOI metadata linkage and article-header title match"
    return (
        False,
        "requested title absent from article header (cover/references excluded)",
    )


class RepositoryPDFRecovery:
    def __init__(self, *, email: str, session=None, unpaywall=None, max_seconds=120):
        self.session = session or requests.Session()
        self.email = email
        self.unpaywall = unpaywall or UnpaywallClient(email, session=self.session)
        self.max_seconds = max_seconds

    def recover(self, *, pmid: str, output_dir: Path, doi=None, title=None):
        result = RepositoryPDFResult(
            pmid=str(pmid), doi=_normalized_doi(doi), title=title or ""
        )
        self.deadline = time.monotonic() + self.max_seconds
        self.asset_requests = 0
        self.visited = set()
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        try:
            self._metadata(result)
            if not result.doi:
                result.error = "no DOI resolved from bibliographic metadata"
                return result
            for provider in (
                self._unpaywall_candidates,
                self._openalex_candidates,
                self._hal_candidates,
            ):
                if time.monotonic() >= self.deadline:
                    break
                for candidate in provider(result):
                    recovered = self._candidate(candidate, result, output_dir)
                    if recovered:
                        result.success = True
                        return result
            result.error = (
                "no indexed copy passed download, identity and article-body gates"
            )
        except (requests.RequestException, ValueError, OSError) as exc:
            result.error = str(exc)
        finally:
            audit = asdict(result)
            audit.pop("markdown")
            audit["recorded_at"] = datetime.now(timezone.utc).isoformat()
            (output_dir / f"{pmid}_repository_recovery.json").write_text(
                json.dumps(audit, indent=2) + "\n"
            )
        return result

    def _json(self, url, params, result, provider):
        if time.monotonic() >= self.deadline:
            return {}
        try:
            response = self.session.get(
                url,
                params=params,
                timeout=min(20, max(1, self.deadline - time.monotonic())),
                headers={
                    "Accept": "application/json",
                    "User-Agent": f"GeneVariantFetcher/1.0 (mailto:{self.email})",
                },
            )
            response.raise_for_status()
            data = response.json()
            if not isinstance(data, dict):
                raise ValueError("metadata response is not an object")
            result.attempts.append(
                {"provider": provider, "url": response.url, "status": "metadata_ok"}
            )
            return data
        except (requests.RequestException, ValueError) as exc:
            result.attempts.append(
                {"provider": provider, "status": "metadata_failed", "reason": str(exc)}
            )
            return {}

    def _metadata(self, result):
        if result.doi and result.title:
            return
        data = self._json(
            "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
            {
                "query": f"EXT_ID:{result.pmid} AND SRC:MED",
                "format": "json",
                "pageSize": 1,
            },
            result,
            "europe_pmc",
        )
        for record in data.get("resultList", {}).get("result", []):
            if str(record.get("id")) != result.pmid or record.get("source") != "MED":
                continue
            resolved = _normalized_doi(record.get("doi"))
            if result.doi and resolved and result.doi != resolved:
                raise ValueError("PMID metadata conflicts with supplied DOI")
            result.doi = result.doi or resolved
            result.title = result.title or record.get("title", "")

    def _bound(self, data, result):
        title = data.get("title") or data.get("display_name") or ""
        doi_matches = _normalized_doi(data.get("doi")) == result.doi
        if not result.title and doi_matches:
            result.title = title
        bound = bool(
            doi_matches
            and title
            and _title_keys(title).intersection(_title_keys(result.title))
        )
        result.attempts.append(
            {
                "status": "metadata_identity_match"
                if bound
                else "metadata_identity_mismatch",
                "observed_doi": data.get("doi"),
                "observed_title": title,
            }
        )
        return bound

    def _unpaywall_candidates(self, result):
        data, error = self.unpaywall.find_open_access(result.doi)
        result.attempts.append(
            {
                "provider": "unpaywall",
                "status": "metadata_ok" if data else "metadata_failed",
                "reason": error or "",
            }
        )
        if not data or not self._bound(data, result):
            return []
        locations = data.get("oa_locations") or [data.get("best_oa_location") or {}]
        candidates = self._locations(
            locations, "unpaywall", f"https://api.unpaywall.org/v2/{result.doi}"
        )
        result.attempts.append(
            {
                "provider": "unpaywall",
                "status": "indexed_locations",
                "candidates": candidates,
            }
        )
        return candidates

    def _openalex_candidates(self, result):
        url = "https://api.openalex.org/works/https://doi.org/" + quote(
            result.doi, safe="/"
        )
        data = self._json(url, {}, result, "openalex")
        if not self._bound(data, result):
            return []
        candidates = self._locations(data.get("locations", []), "openalex", url)
        result.attempts.append(
            {
                "provider": "openalex",
                "status": "indexed_locations",
                "candidates": candidates,
            }
        )
        return candidates

    def _hal_candidates(self, result):
        url = "https://api.archives-ouvertes.fr/search/"
        value = result.doi.replace("\\", "\\\\").replace('"', '\\"')
        data = self._json(
            url,
            {
                "q": f'doiId_s:"{value}"',
                "wt": "json",
                "rows": 5,
                "fl": "doiId_s,title_s,fileMain_s,uri_s,halId_s",
            },
            result,
            "hal",
        )
        candidates = []
        for doc in data.get("response", {}).get("docs", []):
            titles = doc.get("title_s") or []
            if self._bound(
                {"doi": doc.get("doiId_s"), "title": titles[0] if titles else ""},
                result,
            ) and doc.get("fileMain_s"):
                candidates.append(
                    {
                        "url": doc["fileMain_s"],
                        "provider": "hal",
                        "metadata_url": url,
                        "repository_record": doc.get("uri_s"),
                        "kind": "pdf",
                    }
                )
        result.attempts.append(
            {"provider": "hal", "status": "indexed_locations", "candidates": candidates}
        )
        return candidates

    @staticmethod
    def _locations(locations, provider, metadata_url):
        candidates = []
        for loc in locations:
            repository = (
                loc.get("host_type") == "repository"
                or (loc.get("source") or {}).get("type") == "repository"
            )
            pdf = loc.get("pdf_url") or loc.get("url_for_pdf")
            landing = loc.get("landing_page_url") or loc.get("url_for_landing_page")
            common = {
                "provider": provider,
                "metadata_url": metadata_url,
                "repository_record": landing,
                "version": loc.get("version"),
                "repository": repository,
            }
            if pdf:
                candidates.append({**common, "url": pdf, "kind": "pdf"})
            # Metadata-only repository records can link publicly readable PDFs
            # even when an index's OA flag or direct PDF field is stale.
            if landing and repository and "pubmed.ncbi.nlm.nih.gov" not in landing:
                candidates.append({**common, "url": landing, "kind": "landing"})
        return sorted(
            candidates, key=lambda c: (c["kind"] != "pdf", not c["repository"])
        )

    def _fetch(self, url, *, redirect_chain=()):
        if (
            not _public_url(url)
            or self.asset_requests >= MAX_ASSET_REQUESTS
            or time.monotonic() >= self.deadline
        ):
            raise ValueError("invalid URL or recovery request/time budget exhausted")
        if url in redirect_chain or len(redirect_chain) >= 5:
            raise ValueError("redirect loop or hop limit exceeded")
        self.asset_requests += 1
        # Inspect each redirect before requesting it; do not follow local URLs.
        response = self.session.get(
            url,
            timeout=min(20, max(1, self.deadline - time.monotonic())),
            headers={
                "User-Agent": f"GeneVariantFetcher/1.0 (mailto:{self.email})",
                "Accept": "application/pdf,text/html;q=0.9,*/*;q=0.8",
            },
            stream=True,
            allow_redirects=False,
        )
        try:
            if response.is_redirect:
                location = response.headers.get("Location")
                if not location:
                    raise ValueError("redirect is missing Location")
                return self._fetch(
                    urljoin(url, location), redirect_chain=(*redirect_chain, url)
                )
            response.raise_for_status()
            if int(response.headers.get("Content-Length") or 0) > MAX_BYTES:
                raise ValueError("payload exceeds recovery size limit")
            chunks = []
            size = 0
            for chunk in response.iter_content(64 * 1024):
                size += len(chunk)
                if size > MAX_BYTES or time.monotonic() > self.deadline:
                    raise ValueError("payload exceeds recovery size/time limit")
                chunks.append(chunk)
            return b"".join(chunks), response.url
        finally:
            response.close()

    def _candidate(self, candidate, result, output_dir):
        url = candidate["url"]
        if url in self.visited:
            return False
        self.visited.add(url)
        try:
            body, final_url = self._fetch(url)
            if not body.lstrip().startswith(b"%PDF-"):
                if candidate["kind"] != "landing":
                    raise ValueError("response is not a PDF (magic bytes)")
                soup = BeautifulSoup(body, "html.parser")
                links = []
                for meta in soup.select(
                    'meta[name="citation_pdf_url"], meta[name="eprints.document_url"]'
                ):
                    if meta.get("content"):
                        links.append(urljoin(final_url, meta["content"]))
                for anchor in soup.select("a[href]"):
                    href = anchor["href"]
                    label = anchor.get_text(" ", strip=True).lower()
                    if ".pdf" in urlparse(href).path.lower() or label in {
                        "download",
                        "download pdf",
                        "full text",
                        "view/open",
                    }:
                        links.append(urljoin(final_url, href))
                result.attempts.append(
                    {
                        "provider": candidate["provider"],
                        "url": url,
                        "final_url": final_url,
                        "status": "landing",
                        "pdf_links": len(set(links)),
                    }
                )
                for link in dict.fromkeys(links):
                    if self._candidate(
                        {
                            **candidate,
                            "url": link,
                            "kind": "pdf",
                            "discovered_on": final_url,
                        },
                        result,
                        output_dir,
                    ):
                        return True
                return False
            with fitz.open(stream=body, filetype="pdf") as doc:
                if doc.is_encrypted or len(doc) > 500:
                    raise ValueError("encrypted or oversized PDF")
                pages = [p.get_text(sort=True) for p in doc]
            ok, reason = article_identity(pages, result.title, result.doi)
            if not ok:
                raise ValueError(reason)
            text = "\n\n".join(pages)
            # Reflow only the validation view; preserve page/line boundaries in
            # the source consumed by extraction and retain the original PDF.
            validation_text = re.sub(r"(?<!\n)\n(?!\n)", " ", text)
            for gate in (validate_content_quality, validate_article_content):
                valid, why = gate(validation_text)
                if not valid:
                    raise ValueError(why)
            digest = hashlib.sha256(body).hexdigest()
            pdf_path = output_dir / f"{result.pmid}_repository_{digest[:12]}.pdf"
            pdf_path.write_bytes(body)
            result.pdf_path = str(pdf_path.resolve())
            result.pdf_sha256 = digest
            result.page_count = len(pages)
            result.source_url = final_url
            result.candidate = candidate
            result.identity_reason = reason
            result.markdown = "\n\n".join(
                f"<!-- PDF page {i} -->\n\n{page}" for i, page in enumerate(pages, 1)
            )
            result.attempts.append(
                {
                    "provider": candidate["provider"],
                    "url": url,
                    "final_url": final_url,
                    "status": "accepted",
                    "pdf_sha256": digest,
                    "pages": len(pages),
                    "identity_reason": reason,
                    "discovered_on": candidate.get("discovered_on"),
                }
            )
            return True
        except (requests.RequestException, ValueError, OSError, RuntimeError) as exc:
            result.attempts.append(
                {
                    "provider": candidate["provider"],
                    "url": url,
                    "status": "rejected",
                    "reason": str(exc),
                }
            )
            return False


def write_repository_source(result, output_dir, *, gene=None, supplement_markdown=""):
    """Persist body + evidence in the ordinary corpus-builder input layout."""
    if not result.success:
        raise ValueError("cannot persist an unsuccessful recovery")
    output_dir = Path(output_dir)
    text = (
        f"<!-- {SUPPLEMENT_SURFACE_STATUS_MARKER} unavailable -->\n\n"
        f"# Repository article recovery\n\nPMID: {result.pmid} | DOI: {result.doi}\n"
        f"Source: {result.source_url}\n\n{result.markdown}{supplement_markdown}"
    )
    path = output_dir / f"{result.pmid}_FULL_CONTEXT.md"
    path.write_text(text)
    (output_dir / f"{result.pmid}_CLEANED.md").write_text(text)
    artifact_path = output_dir / f"{result.pmid}_artifacts.json"
    artifact = json.loads(artifact_path.read_text()) if artifact_path.exists() else {}
    artifact.update(
        {
            "pmid": result.pmid,
            "doi": result.doi,
            "title": result.title,
            "gene_symbol": gene or artifact.get("gene_symbol"),
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "main_text": {
                **(artifact.get("main_text") or {}),
                "source": "repository_pdf",
                "chars": len(result.markdown),
                "source_url": result.source_url,
                "pdf_path": result.pdf_path,
                "pdf_sha256": result.pdf_sha256,
                "pages": result.page_count,
                "converter": "PyMuPDF sorted page text; original PDF retained",
                "source_identity_verified": True,
                "identity_reason": result.identity_reason,
                "discovery": result.candidate,
            },
            "supplement_surface_status": "unavailable",
        }
    )
    for key in ("figures", "supplements", "supplement_links"):
        artifact.setdefault(key, [])
    artifact.setdefault("notes", []).append(
        "Repository body recovered; separate supplement surface remains unverified."
    )
    artifact_path.write_text(json.dumps(artifact, indent=2) + "\n")
    return path, text
