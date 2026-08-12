# Audit erratum — 2026-08-12

The locked predictions and reported extraction metrics in this historical run
remain unchanged. A repository audit found three provenance/reporting defects in
the generated presentation:

- `predictions.json` explicitly records `telemetry_available: false`; the
  `report.md` statements claiming exactly zero tokens and zero elapsed time are
  placeholders, not measured cost.
- The lock has no LLM trace-manifest digest. The report's generic trace-artifact
  list must not be read as evidence that per-call traces were attached to this
  evaluation lock.
- The current `selection.json` SHA-256 is
  `4db3a283d673bfd37c2f7d5e662c99e0b9a821540c775d14deb180ea6910b596`,
  matching `LOCK.json`; the older `report.json` embeds the stale value
  `af3f53af4db31b69036c80ce8361e3272e15e7e5828c57e524023d838d7ed744`.

The report generator now renders missing telemetry and missing locked traces as
unavailable, and external schema-1 projections disclose that their source
workflow may have produced read-only gold scorecards before the projection lock.
Historical locked files were not silently rewritten.
