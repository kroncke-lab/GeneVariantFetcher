"""
Pedigree Extractor Module

Uses vision-capable LLMs to extract family/carrier information from
pedigree diagrams in scientific papers.

Pedigree figures are common in case reports and family studies but often
contain information not described in the text (e.g., asymptomatic carriers
shown as half-filled symbols).
"""

import base64
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

from litellm import completion

from config.settings import get_settings
from utils.llm_utils import parse_llm_json_response

logger = logging.getLogger(__name__)


# Prompt for detecting if an image is a pedigree diagram
PEDIGREE_DETECTION_PROMPT = """Analyze this figure from a medical genetics paper.

Is this a family pedigree diagram? Pedigree diagrams show family trees with:
- Squares (males) and circles (females)
- Filled shapes (affected individuals) and empty shapes (unaffected)
- Half-filled shapes (carriers)
- Lines connecting family members
- Generation labels (I, II, III, etc.)

Answer with JSON only:
{
    "is_pedigree": true or false,
    "confidence": 0.0 to 1.0,
    "reason": "brief explanation of why this is or isn't a pedigree"
}"""


# Prompt for extracting family structure from a pedigree
PEDIGREE_EXTRACTION_PROMPT = """Extract the family structure from this pedigree diagram.

Standard pedigree symbols:
- Square = male, Circle = female
- Filled/solid = affected by the condition
- Empty/unfilled = unaffected
- Half-filled (left or right) = carrier/heterozygote
- Diagonal line through symbol = deceased
- Arrow or "P" = proband (index case)
- Diamond = sex unknown
- Horizontal line between symbols = mating/partnership
- Vertical lines = descent to offspring

For each individual visible in the pedigree, extract:
1. Their position (generation Roman numeral + number, e.g., "II-3")
2. Sex (male/female/unknown)
3. Affection status (affected/unaffected/carrier/unknown)
4. Whether they are the proband
5. Any genotype labels shown near the symbol (e.g., "+/−", "WT/mut", "N/M")
6. Any other annotations (age, deceased, etc.)

Return as JSON:
{
    "individuals": [
        {
            "id": "II-1",
            "sex": "male",
            "affection_status": "affected",
            "is_proband": false,
            "is_carrier": false,
            "genotype_label": "+/−",
            "is_deceased": false,
            "age_if_shown": null,
            "notes": "any other annotations"
        }
    ],
    "total_generations": 3,
    "inheritance_pattern": "autosomal dominant" or "autosomal recessive" or "X-linked" or "unknown",
    "family_notes": "any additional observations about the pedigree"
}

Be thorough - count every individual shown, even those without detailed annotations.
If you cannot determine a value, use null rather than guessing."""


class PedigreeExtractor:
    """
    Extract carrier/family information from pedigree figures using vision models.

    This class provides methods to:
    1. Detect whether an image is a pedigree diagram
    2. Extract structured family data from pedigree images
    3. Process all figures from a paper to find and analyze pedigrees
    """

    def __init__(
        self,
        model: Optional[str] = None,
        detection_confidence_threshold: float = 0.7,
    ):
        """
        Initialize the PedigreeExtractor.

        Args:
            model: Vision-capable model to use. If None, uses settings.vision_model.
                   Must support image input (e.g., gpt-4o, gemini-1.5-pro, claude-3-5-sonnet).
            detection_confidence_threshold: Minimum confidence to consider an image
                                           a pedigree (0.0-1.0).
        """
        settings = get_settings()
        self.model = model or settings.vision_model
        self.detection_threshold = detection_confidence_threshold

        logger.info(f"PedigreeExtractor initialized with model={self.model}")

    def _image_to_base64_url(self, image_path: Path) -> str:
        """
        Convert an image file to a base64 data URL for vision API calls.

        Args:
            image_path: Path to the image file

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
        }.get(ext, "image/png")

        return f"data:{mime_type};base64,{b64}"

    def _call_vision_model(
        self,
        prompt: str,
        image_path: Path,
        max_tokens: int = 2000,
    ) -> Optional[Dict]:
        """
        Make a vision API call with an image.

        Args:
            prompt: Text prompt for the model
            image_path: Path to the image file
            max_tokens: Maximum tokens in response

        Returns:
            Parsed JSON response, or None on failure
        """
        try:
            image_url = self._image_to_base64_url(image_path)

            response = completion(
                model=self.model,
                messages=[
                    {
                        "role": "user",
                        "content": [
                            {"type": "text", "text": prompt},
                            {"type": "image_url", "image_url": {"url": image_url}},
                        ],
                    }
                ],
                temperature=0,
                max_tokens=max_tokens,
                response_format={"type": "json_object"},
            )

            result_text = response.choices[0].message.content
            return parse_llm_json_response(result_text)

        except Exception as e:
            logger.error(f"Vision API call failed for {image_path.name}: {e}")
            return None

    def is_pedigree(self, image_path: Path) -> tuple[bool, float, str]:
        """
        Check if an image is a pedigree diagram.

        Args:
            image_path: Path to the image file

        Returns:
            Tuple of (is_pedigree, confidence, reason)
        """
        result = self._call_vision_model(
            PEDIGREE_DETECTION_PROMPT,
            image_path,
            max_tokens=200,
        )

        if result is None:
            return False, 0.0, "Detection failed"

        is_ped = result.get("is_pedigree", False)
        confidence = result.get("confidence", 0.0)
        reason = result.get("reason", "")

        return is_ped, confidence, reason

    def extract_pedigree(self, image_path: Path) -> Optional[Dict]:
        """
        Extract family structure from a pedigree image.

        Args:
            image_path: Path to the pedigree image

        Returns:
            Dict containing:
            - individuals: List of family members with their attributes
            - total_generations: Number of generations shown
            - inheritance_pattern: Detected inheritance pattern
            - family_notes: Additional observations

            Returns None on failure.
        """
        result = self._call_vision_model(
            PEDIGREE_EXTRACTION_PROMPT,
            image_path,
            max_tokens=4000,
        )

        if result is None:
            return None

        # Validate required fields
        if "individuals" not in result:
            logger.warning(f"Extraction result missing 'individuals' for {image_path}")
            return None

        return result

    def process_figures_directory(
        self,
        figures_dir: Path,
        detect_only: bool = False,
    ) -> List[Dict]:
        """
        Process all figures in a directory, finding and extracting pedigrees.

        Args:
            figures_dir: Directory containing extracted figure images
            detect_only: If True, only detect pedigrees without full extraction

        Returns:
            List of results for each identified pedigree, containing:
            - image: Filename of the image
            - confidence: Detection confidence
            - (if not detect_only) Full extraction data
        """
        if not figures_dir.exists():
            logger.warning(f"Figures directory does not exist: {figures_dir}")
            return []

        # Find all image files
        image_extensions = {".png", ".jpg", ".jpeg", ".gif", ".webp"}
        image_files = [
            f
            for f in sorted(figures_dir.iterdir())
            if f.suffix.lower() in image_extensions
        ]

        if not image_files:
            logger.info(f"No images found in {figures_dir}")
            return []

        logger.info(f"Processing {len(image_files)} images from {figures_dir}")
        results = []

        for img_path in image_files:
            # First, detect if this is a pedigree
            is_ped, confidence, reason = self.is_pedigree(img_path)

            if is_ped and confidence >= self.detection_threshold:
                logger.info(
                    f"Pedigree detected: {img_path.name} (confidence: {confidence:.2f})"
                )
                print(
                    f"  Found pedigree: {img_path.name} (confidence: {confidence:.0%})"
                )

                if detect_only:
                    results.append(
                        {
                            "image": img_path.name,
                            "image_path": str(img_path),
                            "is_pedigree": True,
                            "confidence": confidence,
                            "detection_reason": reason,
                        }
                    )
                else:
                    # Full extraction
                    extraction = self.extract_pedigree(img_path)
                    if extraction:
                        individuals = extraction.get("individuals", [])
                        affected = sum(
                            1
                            for ind in individuals
                            if ind.get("affection_status") == "affected"
                        )
                        carriers = sum(
                            1 for ind in individuals if ind.get("is_carrier")
                        )
                        print(
                            f"    Extracted: {len(individuals)} individuals, "
                            f"{affected} affected, {carriers} carriers"
                        )

                        results.append(
                            {
                                "image": img_path.name,
                                "image_path": str(img_path),
                                "confidence": confidence,
                                "detection_reason": reason,
                                **extraction,
                            }
                        )
                    else:
                        logger.warning(
                            f"Extraction failed for pedigree {img_path.name}"
                        )
                        results.append(
                            {
                                "image": img_path.name,
                                "image_path": str(img_path),
                                "is_pedigree": True,
                                "confidence": confidence,
                                "detection_reason": reason,
                                "extraction_error": "Full extraction failed",
                            }
                        )
            else:
                logger.debug(
                    f"Not a pedigree: {img_path.name} "
                    f"(is_ped={is_ped}, confidence={confidence:.2f})"
                )

        return results

    def summarize_for_extraction(self, pedigree_results: List[Dict]) -> str:
        """
        Create a text summary of pedigree findings for inclusion in variant extraction.

        This summary can be appended to the paper text before sending to the
        ExpertExtractor to incorporate pedigree-derived carrier information.

        Args:
            pedigree_results: Results from process_figures_directory()

        Returns:
            Markdown-formatted summary of pedigree findings
        """
        if not pedigree_results:
            return ""

        lines = [
            "\n\n# PEDIGREE ANALYSIS (extracted from figures)\n",
            "*The following carrier information was extracted from pedigree figures using vision analysis.*\n",
        ]

        for idx, result in enumerate(pedigree_results, 1):
            image_name = result.get("image", "Unknown")
            individuals = result.get("individuals", [])
            inheritance = result.get("inheritance_pattern", "unknown")
            family_notes = result.get("family_notes", "")

            lines.append(f"\n## Pedigree {idx} ({image_name})\n")

            if inheritance and inheritance != "unknown":
                lines.append(f"**Inheritance pattern:** {inheritance}\n")

            if family_notes:
                lines.append(f"**Notes:** {family_notes}\n")

            if individuals:
                lines.append(
                    "\n| ID | Sex | Status | Proband | Carrier | Genotype | Notes |"
                )
                lines.append(
                    "|-----|-----|--------|---------|---------|----------|-------|"
                )

                for ind in individuals:
                    ind_id = ind.get("id", "?")
                    sex = ind.get("sex", "?")
                    status = ind.get("affection_status", "?")
                    is_proband = "Yes" if ind.get("is_proband") else ""
                    is_carrier = "Yes" if ind.get("is_carrier") else ""
                    genotype = ind.get("genotype_label") or ""
                    notes_parts = []
                    if ind.get("is_deceased"):
                        notes_parts.append("deceased")
                    if ind.get("age_if_shown"):
                        notes_parts.append(f"age {ind['age_if_shown']}")
                    if ind.get("notes"):
                        notes_parts.append(ind["notes"])
                    notes = "; ".join(notes_parts)

                    lines.append(
                        f"| {ind_id} | {sex} | {status} | {is_proband} | "
                        f"{is_carrier} | {genotype} | {notes} |"
                    )

                # Summary counts
                affected_count = sum(
                    1
                    for ind in individuals
                    if ind.get("affection_status") == "affected"
                )
                unaffected_count = sum(
                    1
                    for ind in individuals
                    if ind.get("affection_status") == "unaffected"
                )
                carrier_count = sum(1 for ind in individuals if ind.get("is_carrier"))

                lines.append(
                    f"\n**Summary:** {len(individuals)} individuals total, "
                    f"{affected_count} affected, {unaffected_count} unaffected, "
                    f"{carrier_count} carriers\n"
                )

        return "\n".join(lines)


def extract_pedigrees_from_paper(
    figures_dir: Path,
    model: Optional[str] = None,
    detect_only: bool = False,
) -> List[Dict]:
    """
    Convenience function to extract pedigree data from a paper's figures.

    Args:
        figures_dir: Directory containing figure images (e.g., {pmid}_figures/)
        model: Vision model to use (defaults to settings.vision_model)
        detect_only: If True, only detect pedigrees without full extraction

    Returns:
        List of pedigree extraction results
    """
    extractor = PedigreeExtractor(model=model)
    return extractor.process_figures_directory(figures_dir, detect_only=detect_only)
