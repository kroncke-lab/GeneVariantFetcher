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
import os
from pathlib import Path
from typing import Dict, List, Optional

import requests

from config.settings import get_settings
from utils.llm_trace import (
    OUTCOME_DISCARDED,
    OUTCOME_ERROR,
    OUTCOME_PARSE_FAILED,
    OUTCOME_PARSED,
    attempt_link_summary,
    capture_llm_call,
    last_llm_trace,
    llm_attempt_ledger,
    llm_trace_scope,
    note_llm_attempt,
    note_llm_outcome,
    record_trace_event,
)
from utils.llm_utils import (
    azure_responses_api_url,
    build_reasoning_effort_kwargs,
    build_responses_reasoning_param,
    litellm_completion,
    normalize_azure_ai_api_base,
    parse_llm_json_response,
)

logger = logging.getLogger(__name__)


# Models that use the Azure AI Foundry Responses API rather than chat
# completions. gpt-5 family deployments (gpt-5, gpt-5-codex, gpt-5.3-codex-1)
# return "operation unsupported" on /chat/completions and must be called via
# /openai/v1/responses.
_RESPONSES_API_PREFIXES = (
    "gpt-5",
    "azure_ai/gpt-5",
)


def _uses_responses_api(model: str) -> bool:
    name = model.lower()
    return any(name.startswith(p) for p in _RESPONSES_API_PREFIXES)


def _strip_provider_prefix(model: str) -> str:
    """Drop a 'azure_ai/' prefix to get the bare deployment name."""
    if model.startswith("azure_ai/"):
        return model[len("azure_ai/") :]
    return model


def _call_azure_responses_api_vision(
    *,
    deployment: str,
    prompt: str,
    image_data_url: str,
    max_output_tokens: int,
) -> Optional[Dict]:
    """Call Azure AI Foundry Responses API with a single image input.

    Uses the api-version=v1 route since Azure-hosted gpt-5 deployments only
    accept that version for the Responses API. Returns the parsed JSON the
    model emitted, or None on transport / non-200 errors. Caller is
    responsible for retry semantics.
    """
    base = normalize_azure_ai_api_base()
    key = os.environ.get("AZURE_AI_API_KEY", "")
    if not base or not key:
        logger.error(
            "Responses API call requires AZURE_AI_API_BASE and AZURE_AI_API_KEY"
        )
        return None

    url = azure_responses_api_url(base)
    body = {
        "model": deployment,
        "input": [
            {
                "role": "user",
                "content": [
                    {"type": "input_text", "text": prompt},
                    {"type": "input_image", "image_url": image_data_url},
                ],
            }
        ],
        "max_output_tokens": max_output_tokens,
        **build_responses_reasoning_param(
            deployment, get_settings().vision_reasoning_effort
        ),
    }

    def send_request() -> dict:
        response = requests.post(
            url,
            headers={"api-key": key, "Content-Type": "application/json"},
            json=body,
            timeout=120,
        )
        if response.status_code != 200:
            raise RuntimeError(
                f"Responses API returned {response.status_code}: {response.text[:300]}"
            )
        try:
            return response.json()
        except json.JSONDecodeError as exc:
            raise RuntimeError(
                f"Responses API returned invalid JSON: {response.text[:300]}"
            ) from exc

    try:
        data, _trace = capture_llm_call(
            provider="azure_openai_responses_http",
            requested_model=deployment,
            resolved_model=deployment,
            request={"endpoint": url, "body": body},
            call=send_request,
        )
    except (requests.RequestException, RuntimeError) as e:
        logger.error("Responses API request failed: %s", e)
        return None

    # The Responses API returns an "output" array. Find the first output_text
    # block inside the first message-shaped item.
    text_payload: Optional[str] = None
    for item in data.get("output", []) or []:
        if item.get("type") == "message":
            for content in item.get("content", []) or []:
                if content.get("type") == "output_text":
                    text_payload = content.get("text")
                    break
        if text_payload:
            break

    if text_payload is None:
        logger.warning(
            f"Responses API returned no output_text; keys={list(data.keys())}"
        )
        return None

    return parse_llm_json_response(text_payload)


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
        *,
        gene: Optional[str] = None,
        pmid: Optional[str] = None,
    ):
        """
        Initialize the PedigreeExtractor.

        Args:
            model: Vision-capable model to use. If None, uses settings.vision_model.
                   Must support image input (e.g., gpt-4o, gemini-1.5-pro, claude-3-5-sonnet).
            detection_confidence_threshold: Minimum confidence to consider an image
                                           a pedigree (0.0-1.0).
            gene: Target gene symbol, used only to scope the vision trace.
            pmid: PubMed ID, used only to scope the vision trace. A pedigree
                  prompt contains neither, so without these the trace lands in
                  ``_unscoped/`` and a curator cannot tell which paper a
                  pedigree read belongs to.
        """
        settings = get_settings()
        self.model = model or settings.get_vision_model()
        self.detection_threshold = detection_confidence_threshold
        self.gene = gene
        self.pmid = pmid
        # Why the last vision call produced no usable result, so the decision
        # event can distinguish empty output / unparseable JSON / a failed call.
        self._last_vision_status: Optional[str] = None
        self._last_vision_error: Optional[str] = None

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
        *,
        stage: str = "pedigree_vision",
        role: str = "pedigree_vision",
    ) -> Optional[Dict]:
        """
        Make a vision API call with an image.

        Routes Azure gpt-5 family deployments through the Responses API
        (`/openai/v1/responses?api-version=v1`) since chat completions returns
        "operation unsupported" for those models. Other vision-capable models
        continue to use the standard chat-completions path via LiteLLM.

        ``stage``/``role`` separate detection from extraction in the trace: they
        ask different questions of the same model about the same image, and
        collapsing them left a curator unable to tell which call produced a
        detection confidence and which produced a family structure.

        Args:
            prompt: Text prompt for the model
            image_path: Path to the image file
            max_tokens: Maximum tokens in response

        Returns:
            Parsed JSON response, or None on failure
        """
        responses_path = _uses_responses_api(self.model)
        with llm_trace_scope(
            gene=self.gene,
            pmid=self.pmid,
            stage=stage,
            component=self.__class__.__name__,
            operation="responses_api" if responses_path else "chat_completions",
            image_path=str(image_path),
        ):
            trace_id: Optional[str] = None
            try:
                image_url = self._image_to_base64_url(image_path)
                if responses_path:
                    parsed = _call_azure_responses_api_vision(
                        deployment=_strip_provider_prefix(self.model),
                        prompt=prompt,
                        image_data_url=image_url,
                        # gpt-5 family uses reasoning tokens — give it a generous
                        # budget so the visible output isn't pre-empted.
                        max_output_tokens=max(max_tokens * 4, 4096),
                    )
                    trace_id = note_llm_attempt(last_llm_trace(), role=role)
                    note_llm_outcome(
                        trace_id,
                        OUTCOME_PARSED if parsed else OUTCOME_PARSE_FAILED,
                    )
                    self._last_vision_status = (
                        "parsed" if parsed else "unparseable_response"
                    )
                    return parsed

                response = litellm_completion(
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
                    **build_reasoning_effort_kwargs(
                        self.model, get_settings().vision_reasoning_effort
                    ),
                )
                trace_id = note_llm_attempt(last_llm_trace(), role=role)
                result_text = response.choices[0].message.content
                if not (result_text or "").strip():
                    note_llm_outcome(trace_id, OUTCOME_DISCARDED)
                    self._last_vision_status = "empty_output"
                    logger.error(
                        f"Vision API returned no content for {image_path.name}"
                    )
                    return None
                parsed = parse_llm_json_response(result_text)
                note_llm_outcome(trace_id, OUTCOME_PARSED)
                self._last_vision_status = "parsed"
                return parsed

            except Exception as e:
                if trace_id is None:
                    trace_id = note_llm_attempt(last_llm_trace(), role=role)
                # A JSON failure on a successful provider call is a parse failure;
                # anything else means the call itself did not come back.
                if isinstance(e, (json.JSONDecodeError, ValueError)):
                    note_llm_outcome(trace_id, OUTCOME_PARSE_FAILED)
                    self._last_vision_status = "unparseable_response"
                else:
                    note_llm_outcome(trace_id, OUTCOME_ERROR)
                    self._last_vision_status = "call_failed"
                self._last_vision_error = f"{type(e).__name__}: {e}"
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
        with llm_attempt_ledger():
            self._reset_vision_status()
            result = self._call_vision_model(
                PEDIGREE_DETECTION_PROMPT,
                image_path,
                max_tokens=200,
                stage="pedigree_detection",
                role="pedigree_detect",
            )

            if result is None:
                self._record_pedigree_decision(
                    "pedigree_detection_decision",
                    image_path,
                    {"is_pedigree": False, "confidence": 0.0},
                    decision_source=self._last_vision_status or "call_failed",
                )
                return False, 0.0, "Detection failed"

            is_ped = result.get("is_pedigree", False)
            confidence = result.get("confidence", 0.0)
            reason = result.get("reason", "")
            self._record_pedigree_decision(
                "pedigree_detection_decision",
                image_path,
                {
                    "is_pedigree": bool(is_ped),
                    "confidence": confidence,
                    "detection_threshold": self.detection_threshold,
                    "passes_threshold": bool(
                        is_ped
                        and isinstance(confidence, (int, float))
                        and confidence >= self.detection_threshold
                    ),
                    "rationale": reason,
                },
                decision_source="model_detection",
            )
            return is_ped, confidence, reason

    def _reset_vision_status(self) -> None:
        self._last_vision_status: Optional[str] = None
        self._last_vision_error: Optional[str] = None

    def _record_pedigree_decision(
        self,
        event_type: str,
        image_path: Path,
        data: Dict,
        *,
        decision_source: str,
    ) -> None:
        record_trace_event(
            event_type,
            {
                "image": image_path.name,
                "image_path": str(image_path),
                "model": self.model,
                "decision_source": decision_source,
                "error": getattr(self, "_last_vision_error", None),
                **data,
                **attempt_link_summary(),
            },
        )

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
        with llm_attempt_ledger():
            self._reset_vision_status()
            result = self._call_vision_model(
                PEDIGREE_EXTRACTION_PROMPT,
                image_path,
                max_tokens=4000,
                stage="pedigree_extraction",
                role="pedigree_extract",
            )

            if result is None:
                self._record_pedigree_decision(
                    "pedigree_extraction_decision",
                    image_path,
                    {"individual_count": 0},
                    decision_source=self._last_vision_status or "call_failed",
                )
                return None

            # Validate required fields
            if "individuals" not in result:
                logger.warning(
                    f"Extraction result missing 'individuals' for {image_path}"
                )
                self._record_pedigree_decision(
                    "pedigree_extraction_decision",
                    image_path,
                    {"individual_count": 0, "returned_keys": sorted(result)[:20]},
                    decision_source="missing_individuals_field",
                )
                return None

            individuals = result.get("individuals")
            self._record_pedigree_decision(
                "pedigree_extraction_decision",
                image_path,
                {
                    "individual_count": (
                        len(individuals) if isinstance(individuals, list) else 0
                    ),
                    "total_generations": result.get("total_generations"),
                    "inheritance_pattern": result.get("inheritance_pattern"),
                    "rationale": result.get("family_notes"),
                },
                decision_source="model_extraction",
            )
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
