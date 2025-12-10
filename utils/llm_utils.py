"""
LLM utilities for making structured calls to language models.

This module provides a base class and utilities for making LLM calls with
consistent error handling, retry logic, and JSON response parsing.
"""

import json
import logging
from typing import Dict, Any, Optional, List
from litellm import completion

from .retry_utils import llm_retry

logger = logging.getLogger(__name__)


class BaseLLMCaller:
    """
    Base class for components that make LLM API calls.

    This class encapsulates common LLM calling patterns including:
    - Consistent parameter configuration (model, temperature, max_tokens)
    - Retry logic for transient failures
    - JSON response parsing and validation
    - Error handling and logging

    Attributes:
        model: The LLM model identifier (e.g., "gpt-4o", "claude-3-opus")
        temperature: Sampling temperature for response generation (0.0-1.0)
        max_tokens: Maximum number of tokens in the response

    Example:
        >>> class MyFilter(BaseLLMCaller):
        >>>     def __init__(self):
        >>>         super().__init__(model="gpt-4o", temperature=0.3)
        >>>
        >>>     def analyze(self, text: str) -> dict:
        >>>         prompt = f"Analyze this text: {text}"
        >>>         return self.call_llm_json(prompt)
    """

    def __init__(
        self,
        model: str = "gpt-4o",
        temperature: float = 0.3,
        max_tokens: float = 4000
    ):
        """
        Initialize the LLM caller with model parameters.

        Args:
            model: The LLM model to use (default: "gpt-4o")
            temperature: Sampling temperature (default: 0.3)
            max_tokens: Maximum response length (default: 4000)
        """
        self.model = model
        self.temperature = temperature
        self.max_tokens = max_tokens
        logger.info(
            f"Initialized {self.__class__.__name__} with model={model}, "
            f"temperature={temperature}, max_tokens={max_tokens}"
        )

    def _attempt_json_repair(self, raw_text: str) -> Optional[Dict[str, Any]]:
        """
        Best-effort repair for malformed/truncated JSON emitted by the LLM.

        Uses the same model with a constrained prompt to fix the structure while
        trimming incomplete trailing items. Returns None on failure.

        For large variant extractions, uses higher token limits to preserve data.
        """
        # Use a larger char limit for repair - variant tables can be very large
        # Each variant is ~300-500 chars of JSON, so 150 variants needs ~75K chars
        max_repair_chars = 80000
        repair_chunk = raw_text[:max_repair_chars]

        repair_prompt = (
            "The following text is intended to be a JSON object but is malformed or truncated. "
            "Return a valid JSON object that keeps ALL intact content and discards only incomplete tail fragments. "
            "CRITICAL: Preserve all special characters in protein notation, especially asterisks (*) for stop codons "
            "(e.g., 'p.Arg412*' or 'p.Gly24fs*58'). "
            "Do not add new information. Respond with JSON only.\n\n"
            f"{repair_chunk}"
        )

        # Use higher token limit for repair to handle large variant tables
        # Default repair needs more tokens for 100+ variant extractions
        repair_max_tokens = min(self.max_tokens, 16000)

        try:
            response = completion(
                model=self.model,
                messages=[
                    {"role": "system", "content": "You fix malformed JSON. Preserve all data including special characters like asterisks (*) in protein notation. Respond with JSON only."},
                    {"role": "user", "content": repair_prompt},
                ],
                temperature=0,
                max_tokens=repair_max_tokens,
                response_format={"type": "json_object"},
            )
            repaired_text = response.choices[0].message.content
            return parse_llm_json_response(repaired_text)
        except Exception as repair_exc:
            logger.error(f"JSON repair attempt failed: {repair_exc}")
            return None

    @llm_retry
    def call_llm_json_with_status(
        self,
        prompt: str,
        system_message: Optional[str] = None,
        response_format: Optional[Dict[str, Any]] = None
    ) -> tuple:
        """
        Call the LLM and parse JSON, returning (data, was_truncated, raw_text).

        This variant returns additional metadata needed for continuation extraction.
        """
        messages = []
        if system_message:
            messages.append({"role": "system", "content": system_message})
        messages.append({"role": "user", "content": prompt})

        if response_format is None:
            response_format = {"type": "json_object"}

        response = completion(
            model=self.model,
            messages=messages,
            temperature=self.temperature,
            max_tokens=self.max_tokens,
            response_format=response_format
        )

        result_text = response.choices[0].message.content
        finish_reason = getattr(response.choices[0], "finish_reason", None)
        was_truncated = finish_reason == "length"

        try:
            result_data = parse_llm_json_response(result_text)
            return result_data, was_truncated, result_text
        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse LLM JSON response: {e}")
            if was_truncated:
                logger.warning("LLM response was cut off due to max_tokens; attempting repair.")
            repaired = self._attempt_json_repair(result_text)
            if repaired is not None:
                logger.info("JSON repair succeeded after initial parse failure.")
                return repaired, was_truncated, result_text
            raise

    @llm_retry
    def call_llm_json(
        self,
        prompt: str,
        system_message: Optional[str] = None,
        response_format: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        Call the LLM with a prompt and parse the JSON response.

        This method handles the full lifecycle of an LLM call:
        1. Constructs the message array with optional system message
        2. Calls the LLM API with retry logic
        3. Extracts and parses the JSON response
        4. Validates the response and handles errors

        Args:
            prompt: The user prompt to send to the LLM
            system_message: Optional system message to set context
            response_format: Optional response format specification
                           (default: {"type": "json_object"})

        Returns:
            Parsed JSON response as a dictionary

        Raises:
            json.JSONDecodeError: If the response is not valid JSON
            Exception: For other LLM API errors

        Example:
            >>> caller = BaseLLMCaller()
            >>> result = caller.call_llm_json("What is 2+2?")
            >>> print(result)
        """
        # Build messages array
        messages = []
        if system_message:
            messages.append({"role": "system", "content": system_message})
        messages.append({"role": "user", "content": prompt})

        # Use default JSON response format if not specified
        if response_format is None:
            response_format = {"type": "json_object"}

        logger.debug(f"Calling LLM with model={self.model}, prompt length={len(prompt)}")

        try:
            # Make the LLM API call
            response = completion(
                model=self.model,
                messages=messages,
                temperature=self.temperature,
                max_tokens=self.max_tokens,
                response_format=response_format
            )

            # Extract response text
            result_text = response.choices[0].message.content
            finish_reason = getattr(response.choices[0], "finish_reason", None)

            # Parse JSON response
            result_data = parse_llm_json_response(result_text)

            logger.debug(f"LLM call successful, response keys: {list(result_data.keys())}")
            return result_data

        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse LLM JSON response: {e}")
            logger.error(f"Response text: {result_text[:500]}")
            if finish_reason == "length":
                logger.warning("LLM response was cut off due to max_tokens; attempting repair.")
            repaired = self._attempt_json_repair(result_text)
            if repaired is not None:
                logger.info("JSON repair succeeded after initial parse failure.")
                return repaired
            raise
        except Exception as e:
            logger.error(f"LLM API call failed: {e}")
            raise

    def call_llm_text(
        self,
        prompt: str,
        system_message: Optional[str] = None
    ) -> str:
        """
        Call the LLM and return the raw text response.

        Use this method when you need free-form text rather than structured JSON.

        Args:
            prompt: The user prompt to send to the LLM
            system_message: Optional system message to set context

        Returns:
            Raw text response from the LLM

        Example:
            >>> caller = BaseLLMCaller()
            >>> text = caller.call_llm_text("Write a haiku about code")
            >>> print(text)
        """
        # Build messages array
        messages = []
        if system_message:
            messages.append({"role": "system", "content": system_message})
        messages.append({"role": "user", "content": prompt})

        logger.debug(f"Calling LLM for text response, prompt length={len(prompt)}")

        try:
            response = completion(
                model=self.model,
                messages=messages,
                temperature=self.temperature,
                max_tokens=self.max_tokens
            )

            result_text = response.choices[0].message.content
            logger.debug(f"LLM text response received, length={len(result_text)}")
            return result_text

        except Exception as e:
            logger.error(f"LLM API call failed: {e}")
            raise


def parse_llm_json_response(response_text: str) -> Dict[str, Any]:
    """
    Parse a JSON response from an LLM.

    This function handles common edge cases in LLM JSON responses:
    - Strips markdown code blocks (```json ... ```)
    - Handles extra whitespace
    - Provides clear error messages on parse failures

    Args:
        response_text: Raw text response from the LLM

    Returns:
        Parsed JSON as a dictionary

    Raises:
        json.JSONDecodeError: If the response cannot be parsed as JSON

    Example:
        >>> response = '```json\\n{"result": "success"}\\n```'
        >>> data = parse_llm_json_response(response)
        >>> print(data)
        {'result': 'success'}
    """
    # Strip markdown code blocks if present
    cleaned_text = response_text.strip()
    if cleaned_text.startswith("```json"):
        cleaned_text = cleaned_text[7:]  # Remove ```json
    if cleaned_text.startswith("```"):
        cleaned_text = cleaned_text[3:]  # Remove ```
    if cleaned_text.endswith("```"):
        cleaned_text = cleaned_text[:-3]  # Remove trailing ```

    cleaned_text = cleaned_text.strip()

    try:
        return json.loads(cleaned_text)
    except json.JSONDecodeError as e:
        logger.error(f"Failed to parse JSON: {e}")
        logger.error(f"Response text (first 500 chars): {cleaned_text[:500]}")
        raise


def create_structured_prompt(
    task_description: str,
    input_data: str,
    output_format: Dict[parse_llm_json_response, parse_llm_json_response],
    examples: Optional[List[Dict[str, Any]]] = None
) -> str:
    """
    Create a well-structured prompt for LLM calls.

    This utility helps build consistent prompts with clear structure:
    - Task description
    - Input data
    - Expected output format
    - Optional examples

    Args:
        task_description: Clear description of what the LLM should do
        input_data: The actual data to process
        output_format: Dictionary describing expected output fields
        examples: Optional list of example inputs/outputs

    Returns:
        Formatted prompt string

    Example:
        >>> prompt = create_structured_prompt(
        >>>     task_description="Extract gene names from text",
        >>>     input_data="The BRCA1 gene is important",
        >>>     output_format={"genes": "list of gene names found"},
        >>>     examples=[{"input": "TP53 mutation", "output": {"genes": ["TP53"]}}]
        >>> )
    """
    prompt_parts = [
        "# Task",
        task_description,
        "",
        "# Input",
        input_data,
        "",
        "# Expected Output Format",
        "Return a JSON object with the following structure:",
    ]

    for field, description in output_format.items():
        prompt_parts.append(f"- {field}: {description}")

    if examples:
        prompt_parts.extend([
            "",
            "# Examples",
        ])
        for i, example in enumerate(examples, 1):
            prompt_parts.append(f"\nExample {i}:")
            prompt_parts.append(f"Input: {example.get('input', 'N/A')}")
            prompt_parts.append(f"Output: {json.dumps(example.get('output', {}), indent=2)}")

    return "\n".join(prompt_parts)
