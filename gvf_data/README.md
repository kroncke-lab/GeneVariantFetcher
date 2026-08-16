# Runtime Data

This owned Python package contains small tracked inputs used by GVF at runtime.
It avoids publishing the collision-prone top-level package name ``data`` and is
not generated pipeline output.

## Contents

- `kcnh2_variant_aliases.json` - gene-specific variant alias mappings loaded by
  `utils/variant_normalizer.py`.
- `reference_sequences/exon_maps.json` - NCBI coding-exon protein spans for
  structural matching (`utils/structural_alleles.py`).
- `reference_sequences/` - canonical RefSeq protein FASTA files,
  `manifest.json`, and `exon_maps.json`, used by
  `pipeline/reference_validation.py` and `utils/structural_alleles.py`.

## Maintenance

- Add new `<gene>_variant_aliases.json` files here when a gene needs stable
  alias normalization. Keep the filename lower-case.
- Update reference FASTA files through `scripts/fetch_reference_sequences.py` so
  fetched sequences are checked against the manifest and expected protein
  lengths.
- Do not put run outputs, downloaded papers, databases, or metrics here. Those
  belong under ignored directories such as `results/`, `validation_runs/`,
  `recall_metrics/`, or `corpus/`.
