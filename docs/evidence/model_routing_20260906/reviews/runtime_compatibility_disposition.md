# Runtime compatibility review disposition

Claude, Grok and Agy reviewed the patch through their CLIs with a supplied
source diff and tools disabled. Their receipts are retained. Reviews are
adversarial suggestions, not independent validation of model accuracy.

- Added a specific Grok 4.6 clamp regression: 32k remains distinct from older
  Grok's 15k policy. The exact-name match already precedes substring fallback.
- Added request-capture coverage for Astra figure text, figure variants and
  pedigree Responses calls. It checks the routed host/key, bare deployment,
  `max_output_tokens`, `reasoning.effort`, and absence of chat/sampling fields.
  Existing body builders already used this shape, so no runtime change was
  needed. Figure variants delegate to the shared figure text client.
- Actual SDK HTTP tests already verify the Chat Completions effort parameter
  after LiteLLM transformation. The live Astra probe also passed; Grok's later
  availability failure is not cured or disproved by offline transport tests.
- `json` was already imported. Smoke's `grok-4` hint includes Grok 4.6. Missing
  selected-route credentials raise the exact configured variable name before
  the older generic missing-default message. Those reported defects were not
  reproduced in the full source.
- Route keys are documented bare deployment names. A provider-prefixed key is
  an unselected entry, not a malformed selected route. No alias guessing was
  added during the frozen experiment. Malformed JSON intentionally fails the
  opted-in process, and missing selected keys never fall back to another key.
- Arbitrary HTTPS hosts remain allowed for an explicitly operator-controlled
  endpoint map; the map and key environment are equally trusted. No guarantee
  against hostile process-environment access is claimed. In particular Agy's
  assertion that this avoids all process-dump leakage is too strong: the key
  itself remains in the process environment. Nothing writes it to this report.
- Astra omits unsupported temperature settings; temperature zero cannot promise
  deterministic inference. `none`/`minimal` map to `low` as documented, so these
  two explicitly supported models do not promise disabled reasoning.

No acceptance rules, sources, primary prompts or production defaults changed
in response to this review. New-model vision accuracy remains untested live.
