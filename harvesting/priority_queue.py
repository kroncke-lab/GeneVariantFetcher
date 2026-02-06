"""
Priority Paper Queue - Manual acquisition system for high-value papers.

Papers that fail automated download are queued for manual acquisition,
prioritized by importance (citation count, relevance score, etc.)
"""

import json
import logging
from dataclasses import dataclass, asdict, field
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import List, Optional, Dict, Any

logger = logging.getLogger(__name__)


class Priority(Enum):
    """Paper acquisition priority levels."""
    CRITICAL = 1    # Must have - key paper for analysis
    HIGH = 2        # Important - likely contains valuable data
    MEDIUM = 3      # Useful - would improve coverage
    LOW = 4         # Nice to have - not essential


class Status(Enum):
    """Paper acquisition status."""
    PENDING = "pending"           # Waiting for manual acquisition
    IN_PROGRESS = "in_progress"   # Someone is working on it
    ACQUIRED = "acquired"         # Successfully obtained
    UNAVAILABLE = "unavailable"   # Confirmed not available anywhere
    SKIPPED = "skipped"           # Decided not to acquire


@dataclass
class QueuedPaper:
    """A paper in the manual acquisition queue."""
    pmid: str
    doi: Optional[str] = None
    title: Optional[str] = None
    journal: Optional[str] = None
    year: Optional[int] = None
    
    # Priority and status
    priority: Priority = Priority.MEDIUM
    status: Status = Status.PENDING
    
    # Failure info
    failure_reason: str = ""
    failed_sources: List[str] = field(default_factory=list)
    
    # Manual acquisition
    acquisition_notes: str = ""
    acquired_via: Optional[str] = None  # "author", "interlibrary", "purchase", etc.
    acquired_by: Optional[str] = None   # Who acquired it
    acquired_date: Optional[str] = None
    
    # Metadata
    citation_count: Optional[int] = None
    relevance_score: Optional[float] = None  # 0-1 from LLM filter
    has_supplements: Optional[bool] = None
    supplement_count: Optional[int] = None
    
    # Timestamps
    added_date: str = field(default_factory=lambda: datetime.now().isoformat())
    updated_date: str = field(default_factory=lambda: datetime.now().isoformat())
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        d = asdict(self)
        d['priority'] = self.priority.name
        d['status'] = self.status.value
        return d
    
    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> 'QueuedPaper':
        """Create from dictionary."""
        d = d.copy()
        d['priority'] = Priority[d.get('priority', 'MEDIUM')]
        d['status'] = Status(d.get('status', 'pending'))
        return cls(**d)


class PriorityQueue:
    """
    Manages the manual acquisition queue.
    
    Usage:
        queue = PriorityQueue(Path("manual_acquisition_queue.json"))
        queue.add_paper(pmid="12345", failure_reason="Paywall", priority=Priority.HIGH)
        
        for paper in queue.get_pending(priority=Priority.HIGH):
            print(f"Need to acquire: {paper.pmid}")
    """
    
    def __init__(self, queue_file: Path):
        """
        Initialize queue.
        
        Args:
            queue_file: Path to JSON file for persistence
        """
        self.queue_file = Path(queue_file)
        self.papers: Dict[str, QueuedPaper] = {}
        self._load()
    
    def _load(self):
        """Load queue from file."""
        if self.queue_file.exists():
            try:
                with open(self.queue_file) as f:
                    data = json.load(f)
                self.papers = {
                    pmid: QueuedPaper.from_dict(paper_data)
                    for pmid, paper_data in data.get('papers', {}).items()
                }
                logger.info(f"Loaded {len(self.papers)} papers from queue")
            except Exception as e:
                logger.error(f"Failed to load queue: {e}")
                self.papers = {}
    
    def _save(self):
        """Save queue to file."""
        data = {
            'updated': datetime.now().isoformat(),
            'papers': {pmid: paper.to_dict() for pmid, paper in self.papers.items()}
        }
        
        # Ensure parent directory exists
        self.queue_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(self.queue_file, 'w') as f:
            json.dump(data, f, indent=2)
    
    def add_paper(
        self,
        pmid: str,
        failure_reason: str,
        priority: Priority = Priority.MEDIUM,
        doi: Optional[str] = None,
        title: Optional[str] = None,
        failed_sources: Optional[List[str]] = None,
        **kwargs
    ) -> QueuedPaper:
        """
        Add a paper to the queue.
        
        Args:
            pmid: PubMed ID
            failure_reason: Why automated download failed
            priority: Acquisition priority
            doi: Digital Object Identifier
            title: Paper title
            failed_sources: List of sources that were tried
            **kwargs: Additional QueuedPaper fields
            
        Returns:
            The queued paper
        """
        if pmid in self.papers:
            # Update existing entry
            paper = self.papers[pmid]
            paper.failure_reason = failure_reason
            if failed_sources:
                paper.failed_sources = list(set(paper.failed_sources + failed_sources))
            paper.updated_date = datetime.now().isoformat()
            logger.info(f"Updated paper {pmid} in queue")
        else:
            # Add new entry
            paper = QueuedPaper(
                pmid=pmid,
                doi=doi,
                title=title,
                priority=priority,
                failure_reason=failure_reason,
                failed_sources=failed_sources or [],
                **kwargs
            )
            self.papers[pmid] = paper
            logger.info(f"Added paper {pmid} to queue (priority: {priority.name})")
        
        self._save()
        return paper
    
    def update_status(
        self,
        pmid: str,
        status: Status,
        acquired_via: Optional[str] = None,
        acquired_by: Optional[str] = None,
        notes: Optional[str] = None
    ):
        """Update paper status."""
        if pmid not in self.papers:
            raise ValueError(f"Paper {pmid} not in queue")
        
        paper = self.papers[pmid]
        paper.status = status
        paper.updated_date = datetime.now().isoformat()
        
        if status == Status.ACQUIRED:
            paper.acquired_date = datetime.now().isoformat()
            paper.acquired_via = acquired_via
            paper.acquired_by = acquired_by
        
        if notes:
            paper.acquisition_notes = notes
        
        self._save()
        logger.info(f"Updated {pmid} status to {status.value}")
    
    def get_pending(
        self,
        priority: Optional[Priority] = None,
        limit: Optional[int] = None
    ) -> List[QueuedPaper]:
        """
        Get pending papers, optionally filtered by priority.
        
        Returns papers sorted by priority (CRITICAL first), then by date added.
        """
        papers = [
            p for p in self.papers.values()
            if p.status == Status.PENDING
            and (priority is None or p.priority == priority)
        ]
        
        # Sort by priority value (lower = higher priority), then by date
        papers.sort(key=lambda p: (p.priority.value, p.added_date))
        
        if limit:
            papers = papers[:limit]
        
        return papers
    
    def get_stats(self) -> Dict[str, Any]:
        """Get queue statistics."""
        total = len(self.papers)
        by_status = {}
        by_priority = {}
        
        for paper in self.papers.values():
            status = paper.status.value
            priority = paper.priority.name
            
            by_status[status] = by_status.get(status, 0) + 1
            by_priority[priority] = by_priority.get(priority, 0) + 1
        
        return {
            'total': total,
            'by_status': by_status,
            'by_priority': by_priority,
            'pending_critical': sum(1 for p in self.papers.values() 
                                   if p.status == Status.PENDING and p.priority == Priority.CRITICAL),
            'pending_high': sum(1 for p in self.papers.values()
                               if p.status == Status.PENDING and p.priority == Priority.HIGH)
        }
    
    def export_for_manual_work(self, output_file: Path, status: Status = Status.PENDING) -> int:
        """
        Export papers to a CSV for manual work.
        
        Returns number of papers exported.
        """
        import csv
        
        papers = [p for p in self.papers.values() if p.status == status]
        papers.sort(key=lambda p: (p.priority.value, p.added_date))
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                'PMID', 'DOI', 'Title', 'Priority', 'Failure Reason',
                'Failed Sources', 'Notes', 'PubMed URL'
            ])
            
            for p in papers:
                writer.writerow([
                    p.pmid,
                    p.doi or '',
                    p.title or '',
                    p.priority.name,
                    p.failure_reason,
                    ', '.join(p.failed_sources),
                    p.acquisition_notes,
                    f"https://pubmed.ncbi.nlm.nih.gov/{p.pmid}/"
                ])
        
        logger.info(f"Exported {len(papers)} papers to {output_file}")
        return len(papers)
