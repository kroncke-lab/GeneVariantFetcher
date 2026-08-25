# Runtime Data

This owned Python package contains small tracked inputs used by GVF at runtime.
It avoids publishing the collision-prone top-level package name ``data`` and is
not generated pipeline output.

## Contents

- `kcnh2_variant_aliases.json` - gene-specific variant alias mappings loaded by
  `utils/variant_normalizer.py`. **Benchmark contamination — being removed.**
  Its own metadata declares `"source": "Gold standard Excel + generated forms"`
  (856 variants / 4,766 aliases), and `_lookup_alias` consults it at runtime for
  every gene, so KCNH2 — one of the four gold-120 genes — normalizes with access
  to the evaluation answer key. Re-matching all 554 pairs in the
  `20260824_aha_table_sourcebound_gold118` lock with the lookup stubbed out
  gives 554/554, so no *score-time* result depends on it; extraction-time
  influence cannot be recovered from already-normalized predictions. An alias
  map is legitimate infrastructure — the requirement is that its contents come
  from public notation resources rather than from a gold standard. Do not cite
  KCNH2 metrics as evidence of gene-general behaviour until this is rebuilt.
  Tracked in `TASKS.md`; detail in
  `docs/evidence/generalization_consult_20260825.md` §1.
- `reference_sequences/exon_maps.json` - NCBI coding-exon protein spans for
  structural matching (`utils/structural_alleles.py`).
- `reference_sequences/` - canonical RefSeq protein FASTA files,
  `manifest.json`, and `exon_maps.json`, used by
  `pipeline/reference_validation.py` and `utils/structural_alleles.py`.

## Maintenance

- Add new `<gene>_variant_aliases.json` files here when a gene needs stable
  alias normalization. Keep the filename lower-case. **Never build one from a
  gold standard, an evaluation fixture, or any answer key** — these files are
  read at runtime, so a gold-derived map silently makes that gene's benchmark
  numbers non-general. Record the real provenance in the file's `metadata.source`
  and prefer public notation resources (HGVS grammar, RefSeq, ClinVar exports).
- Update reference FASTA files through `scripts/fetch_reference_sequences.py` so
  fetched sequences are checked against the manifest and expected protein
  lengths.
- Do not put run outputs, downloaded papers, databases, or metrics here. Those
  belong under ignored directories such as `results/`, `validation_runs/`,
  `recall_metrics/`, or `corpus/`.
