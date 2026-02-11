"""
Pydantic models for GeneVariantFetcher GUI API.

Request and response models for the FastAPI endpoints.
"""

from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field

# =============================================================================
# Job Models
# =============================================================================


class JobCreateRequest(BaseModel):
    """Request to create a new pipeline job."""

    gene_symbol: str = Field(..., description="Gene symbol (e.g., BRCA1, SCN5A)")
    email: str = Field(..., description="Email for NCBI E-utilities")
    output_dir: str = Field(..., description="Output directory for results")

    # Pipeline limits
    max_pmids: int = Field(100, ge=1, le=20000, description="Maximum PMIDs to fetch")
    max_papers_to_download: int = Field(
        50, ge=1, le=1000, description="Maximum papers to download"
    )

    # Synonym settings
    auto_synonyms: bool = Field(False, description="Auto-discover gene synonyms")
    synonyms: List[str] = Field(default_factory=list, description="Manual synonyms")

    # Specific PMIDs for testing (skip discovery if provided)
    specific_pmids: List[str] = Field(
        default_factory=list, description="Specific PMIDs for testing (skips discovery)"
    )

    # Advanced settings
    tier_threshold: int = Field(1, ge=0, le=10, description="Model cascade threshold")
    use_clinical_triage: bool = Field(False, description="Use clinical triage filter")
    enable_tier1: bool = Field(True, description="Enable keyword filtering")
    enable_tier2: bool = Field(True, description="Enable LLM filtering")
    use_pubmind: bool = Field(True, description="Use PubMind source")
    use_pubmed: bool = Field(True, description="Use PubMed source")
    use_europepmc: bool = Field(False, description="Use Europe PMC source")
    tier2_confidence_threshold: float = Field(0.5, ge=0.0, le=1.0)
    scout_enabled: bool = Field(True, description="Enable DATA_ZONES condensation")


class JobResponse(BaseModel):
    """Response with job information."""

    job_id: str
    gene_symbol: str
    current_step: str
    current_step_display: str
    created_at: str
    updated_at: str
    started_at: Optional[str]
    completed_at: Optional[str]
    error_message: Optional[str]
    is_resumable: bool
    output_dir: str
    statistics: Dict[str, Any] = Field(default_factory=dict)


class JobListResponse(BaseModel):
    """Response with list of jobs."""

    jobs: List[JobResponse]
    incomplete_count: int


# =============================================================================
# Settings Models
# =============================================================================


class EnvSettingsResponse(BaseModel):
    """Response with environment settings."""

    settings: Dict[str, Optional[str]]
    env_file_exists: bool
    env_file_path: str


class EnvSettingsUpdateRequest(BaseModel):
    """Request to update environment settings."""

    settings: Dict[str, str]


# =============================================================================
# Directory Browser Models
# =============================================================================


class DirectoryEntry(BaseModel):
    """A directory entry for the file browser."""

    name: str
    path: str
    is_dir: bool
    is_readable: bool = True


class DirectoryListResponse(BaseModel):
    """Response with directory contents."""

    current_path: str
    parent_path: Optional[str]
    entries: List[DirectoryEntry]
    can_create: bool = True


# =============================================================================
# Pipeline Analysis Models
# =============================================================================


class PmidFilteredEntry(BaseModel):
    """A PMID that was deliberately filtered out with its reason."""

    pmid: str
    reason: str


class PmidPaywalledEntry(BaseModel):
    """A PMID that was paywalled or could not be downloaded."""

    pmid: str
    reason: str
    url: Optional[str] = None


class PipelineStageBreakdown(BaseModel):
    """Breakdown of PMIDs at each pipeline stage for debugging."""

    # Discovery stage
    discovered_count: int = 0
    discovered_pmids: List[str] = Field(default_factory=list)

    # Harvesting stage (abstract fetch)
    harvested_count: int = 0
    harvested_pmids: List[str] = Field(default_factory=list)
    harvest_gap_count: int = 0  # discovered - harvested - filtered
    harvest_gap_pmids: List[str] = Field(
        default_factory=list
    )  # PMIDs lost between discovery and harvest

    # Filtering stage
    filtered_out_count: int = 0
    filtered_out_entries: List[PmidFilteredEntry] = Field(default_factory=list)

    # Download stage
    downloaded_count: int = 0
    downloaded_pmids: List[str] = Field(default_factory=list)
    paywalled_count: int = 0
    paywalled_entries: List[PmidPaywalledEntry] = Field(default_factory=list)
    download_gap_count: int = 0  # harvested - downloaded - paywalled - filtered
    download_gap_pmids: List[str] = Field(
        default_factory=list
    )  # PMIDs lost during download

    # Status flags
    has_discovery_data: bool = False
    has_harvest_data: bool = False
    has_filter_data: bool = False
    has_download_data: bool = False

    # Interruption detection
    likely_interrupted: bool = False
    interruption_stage: Optional[str] = None
    interruption_details: Optional[str] = None


class FolderAnalysis(BaseModel):
    """Analysis of an existing folder with downloaded papers."""

    path: str
    is_valid: bool
    gene_symbol: Optional[str] = None
    detected_gene_symbols: List[str] = Field(default_factory=list)

    # File counts
    full_context_count: int = 0
    data_zones_count: int = 0
    extraction_count: int = 0

    # Lists of PMIDs found
    full_context_pmids: List[str] = Field(default_factory=list)
    data_zones_pmids: List[str] = Field(default_factory=list)
    extraction_pmids: List[str] = Field(default_factory=list)

    # Status flags
    has_fulltext: bool = False
    has_scouted: bool = False
    scouting_complete: bool = False
    has_extractions: bool = False

    # Messages
    status_message: str = ""
    recommendations: List[str] = Field(default_factory=list)

    # Pipeline stage breakdown for debugging
    pipeline_breakdown: Optional[PipelineStageBreakdown] = None


# =============================================================================
# Folder Job Models
# =============================================================================


class FolderJobRequest(BaseModel):
    """Request to create a job from an existing folder."""

    folder_path: str = Field(..., description="Path to folder with downloaded papers")
    gene_symbol: str = Field(..., description="Gene symbol for extraction")
    email: str = Field(..., description="Email for NCBI (if needed)")

    # Options
    run_scout: bool = Field(False, description="Run Data Scout on unscouted files")
    skip_already_extracted: bool = Field(
        True, description="Skip PMIDs that already have extractions"
    )

    # Advanced settings
    tier_threshold: int = Field(1, ge=0, le=10)
    scout_enabled: bool = Field(True)
    scout_min_relevance: float = Field(0.3, ge=0.0, le=1.0)


class ResumeJobRequest(BaseModel):
    """Request to resume an interrupted pipeline from a folder."""

    folder_path: str = Field(
        ..., description="Path to folder with partial pipeline results"
    )
    gene_symbol: str = Field(..., description="Gene symbol for the pipeline")
    email: str = Field(..., description="Email for NCBI E-utilities")

    # Resume stage - where to resume from
    resume_stage: str = Field(
        ..., description="Stage to resume from: 'downloading' or 'extraction'"
    )

    # PMIDs to process (for downloading stage)
    pmids_to_download: List[str] = Field(
        default_factory=list,
        description="PMIDs to download (missing from interrupted download)",
    )

    # Pipeline limits for downloading
    max_papers_to_download: int = Field(
        500, ge=1, le=1000, description="Maximum papers to download"
    )

    # Options for extraction
    run_scout: bool = Field(True, description="Run Data Scout on downloaded files")
    skip_already_extracted: bool = Field(
        True, description="Skip PMIDs that already have extractions"
    )

    # Advanced settings
    tier_threshold: int = Field(1, ge=0, le=10)
    scout_enabled: bool = Field(True)
    scout_min_relevance: float = Field(0.3, ge=0.0, le=1.0)


# =============================================================================
# Browser Fetch Models
# =============================================================================


class BrowserFetchRequest(BaseModel):
    """Request to start browser-based fetching of paywalled papers."""

    csv_path: str = Field(..., description="Path to paywalled_missing.csv")
    gene_symbol: str = Field(..., description="Gene symbol for file naming")
    headless: bool = Field(True, description="Run browser in headless mode")
    max_papers: Optional[int] = Field(None, description="Maximum papers to process")
    use_claude: bool = Field(False, description="Use Claude to find download links")
    wait_for_captcha: bool = Field(
        False, description="Wait up to 5 minutes for manual CAPTCHA completion"
    )
    retry_failures_only: bool = Field(
        False, description="Only process papers with browser_failed status"
    )
    auto_process: bool = Field(
        True,
        description="Automatically convert to markdown and run scout after download",
    )
    fallback_to_manual: bool = Field(
        True,
        description="Auto-switch to manual download mode when Cloudflare loops detected",
    )


class BrowserFetchStatus(BaseModel):
    """Status of a browser fetch task."""

    task_id: str
    status: str  # "running", "completed", "failed"
    progress: int = 0
    total: int = 0
    successful: int = 0
    failed: int = 0
    current_pmid: Optional[str] = None
    current_step: Optional[str] = None  # Current step being executed
    error: Optional[str] = None
    output_dir: Optional[str] = None
    recent_logs: List[str] = Field(default_factory=list)  # Last N log lines
