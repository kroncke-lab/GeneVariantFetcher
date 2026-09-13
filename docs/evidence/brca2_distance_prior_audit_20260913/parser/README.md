# Exact catalog-row reproduction before the parser correction

At clean HEAD `1e366f45d0624983df82a1e97b959157079ba941`, the deterministic
markdown parser still converted the ClinVar catalog entry for BRCA2 A1043V
into one carrier and one affected person, with `clinical_significance` set to
`pathogenic`. The source row has no patient identifier, patient count, or
phenotype observation. Its germline classification says conflicting
classifications of pathogenicity.

The frozen source is PMID 40664060, catalog header line 10679 and A1043V row
16755 in the existing cleaned paper. The archived observation is variant_id
12496 in the BRCA2_0 run database. `source_excerpt.md` preserves the exact
header, separator, and row. `reproduction.json` records the source hash,
unmodified parser/router hashes, archived fact provenance, and results.

The three `Somatic clinical impact...` headers match the markdown parser's bare
`clinical` heuristic, although none identifies a person. An unescaped pipe in
Condition(s) also misaligns later cells, but repairing that delimiter alone
still yields affected=1: literal `nan` annotation cells are not recognized as
unknown phenotype. Replacing those cells with `unknown` removes affected=1
while still inventing one carrier. Removing the three annotation columns
turns off the false patient-row gate. The table router's deterministic mapping
and its row-level test already reject this same snippet.

`reproduce_catalog_row.py` is an offline, read-only demonstration of this
pre-correction behavior. It uses no model calls and reads the database through
SQLite `mode=ro`. Its assertion intentionally records the faulty baseline; it
will fail once that specific parser behavior is corrected. The parser-fix
regression tests should verify the post-correction behavior instead. This
demonstrates a current parser route, not a rerun of the entire live extraction
pipeline or proof about every other table in the paper.
