# Manual acquisition worklist — 2026-09-08

Brett's constraint: EZproxy is no longer usable, so remaining source gaps are
closed by a person downloading papers by hand, and a person can do roughly one
to two hundred, not thousands. This directory ranks the candidates so that
effort goes where the gold says the rows are.

Built by [`scripts/recall_audit/rank_manual_acquisition.py`](../../../scripts/recall_audit/rank_manual_acquisition.py)
from the 2026-09-03 source-presence sweep, PubMed metadata (ESummary) and
Unpaywall open-access status. Nothing under `corpus/` was written.

## What the list is out of

The pipeline's own discovery step (PubMed search, Tier 1/2 relevance filtering,
then acquisition) has touched **5,939** PMIDs for the four cardiac genes across
all runs, as recorded in `corpus/INDEX.csv`: KCNH2 1,310, KCNQ1 2,383, SCN5A
1,588, RYR2 658. Of those, 4,942 (83%) have a usable full text on disk and
**996 (17%) are stubs**, an abstract or landing page only. The stubs split into
251 gold papers and 746 papers outside the gold standard; ten more gold PMIDs
were never fetched at all.

Two lists therefore exist here, and they answer different questions:

- **Gold-backed list** (`manual_acquisition_worklist.csv`, 363 papers ranked,
  250 listed): every gold paper with rows behind the acquisition ceiling,
  ranked by *measured* yield. This is what moves the benchmark. It includes the
  251 gold stubs and the gold papers whose body is on disk but whose variants
  sit in a supplement we never fetched.
- **Gold-free list** (`nogold/manual_acquisition_worklist_nogold_top200.csv`,
  the 746 non-gold stubs ranked by *predicted* yield): what a new gene would
  get, since no gold exists to rank by. The predictor is the abstract-only
  acquisition expected value (`scripts/acquisition_ev/predict_yield.py`); on
  the cardiac gold it was validated to recover roughly half of all gold
  carriers when the top 20% of papers by score are acquired, about twice an
  abstract-length baseline. Its top 200 is 88 paywalled PDFs, 66 open-access
  PDFs that a bot never fetched, and 46 papers without a DOI (mostly older or
  non-English journals).

## What "yield" means here

The sweep classified every gold row of the mixed-gold inventory by whether its
variant string exists in anything we hold on disk. Four classes are unreachable
by any reading protocol until the source is acquired: `source_absent`,
`text_absent_stub_body` (landing page / abstract only), `text_absent_garbled_body`
(PDF rendered as glyph codes) and `text_absent_substitution` (a body is on disk
but the variant is in none of it, which in practice means a missing supplement).
A paper's **hard yield** is the number of such gold rows summed over every gene
that cites it, because one download serves all of them. The two undecidable
classes (figure images on disk, non-searchable notation) are reported as
**possible yield** and never ranked on.

Genes: KCNH2, KCNQ1, SCN5A, RYR2 and the BRCA2 collaborator set. Yield is gold
identity presence; a downloaded paper can still lack a count-bearing table.

## Result

| | papers | hard-ceiling gold rows | share of all 1,314 |
| --- | ---: | ---: | ---: |
| Papers with any hard-ceiling row | 363 | 1,314 | 100% |
| Top 25 | 25 | 578 | 44.0% |
| Top 50 | 50 | 778 | 59.2% |
| Top 100 | 100 | 986 | 75.0% |
| **Top 150 (recommended)** | 150 | 1,097 | 83.5% |
| Top 200 | 200 | 1,151 | 87.6% |
| Top 250 (listed) | 250 | 1,201 | 91.4% |

The list is steep: the first paper alone (PMID 23631430, a clinical-laboratory
LQTS referral series cited by three genes) holds 122 rows, and the top 25 hold
44% of everything. Beyond about 150 papers each additional download buys one or
two rows.

Access classes in the recommended top 150, from what the corpus holds plus
Unpaywall:

| class | papers | hard rows | what a person does |
| --- | ---: | ---: | --- |
| `free_pdf_not_fetched` | 42 | 226 | An open-access copy exists but nothing usable is on disk: a bot wall or fetch failure. Download the PDF at the listed URL, plus supplements. |
| `free_supplements` | 15 | 160 | Open-access body is on disk, variants are not in it: fetch the supplementary tables from the OA page or PMC. |
| `paywalled_pdf` | 49 | 369 | No OA copy, no usable body: institutional browser access for PDF and supplements. |
| `paywalled_supplements` | 31 | 276 | Body on disk, variants absent, no OA copy: supplements via institutional access. |
| `no_doi` | 11 | 32 | PubMed has no DOI (older papers): resolve by title first. |
| `garbled_pdf` | 2 | 34 | The body is glyph codes: re-download a clean PDF. |

So 57 of the 150 papers (386 rows) need no institutional access at all, and
those are the ones to do first. The 80 paywalled papers (645 rows) need a
browser session with library access; Wolters Kluwer (Circulation and the AHA
journals), Elsevier and Oxford dominate that half. Two-thirds of the list is
from the 2000s.

## How to use the list

1. Open `manual_acquisition_worklist.csv` (250 rows, one paper each) or the
   readable `manual_acquisition_worklist.md`. Columns give rank, PMID, genes,
   hard and possible yield, cumulative share, access class, the specific
   manual action, corpus state, DOI/PubMed/OA URLs, journal, year and title.
2. Work down the ranks, doing the `free_*` classes first within any batch.
3. Place each downloaded PDF and supplement under the corpus layout the
   corpus builder expects (`corpus/<GENE>/<PMID>/`) and run
   `scripts/build_source_corpus.py` so the fold and index are updated; then
   `scripts/refresh_run_db.py` replays affected runs. Nothing in the worklist
   writes to the corpus by itself.

Regenerate (network: PubMed ESummary in two batches, one Unpaywall lookup per
DOI, cached in `metadata_cache.json`):

```bash
.venv/bin/python scripts/recall_audit/rank_manual_acquisition.py \
  --out-dir docs/evidence/manual_acquisition_20260908 --max-papers 250
```

For a gene without gold the same report ranks a PMID list by the abstract-only
acquisition expected value (`scripts/acquisition_ev/predict_yield.py`):

```bash
.venv/bin/python scripts/recall_audit/rank_manual_acquisition.py \
  --pmid-file <pmids.txt> --gene <GENE> --out-dir <dir>
```

## Limits

Unpaywall status is as of the run date and "bot wall" is inferred from an OA
copy existing while nothing usable is on disk; the harvester's own failure
reason is not recorded per paper. Hard yield counts gold rows, so it favours
papers with many variants over papers with rich per-variant counts; the
`positive_carrier_rows` column is the closest proxy for count value. Papers
outside the gold inventory are not ranked here; use the gold-free mode.

Files: `manual_acquisition_worklist.csv`, `manual_acquisition_worklist.md`,
`summary.json` (capture curve and class counts), `metadata_cache.json`;
`nogold/<GENE>_nongold_stubs.txt` (the 746 non-gold stub PMIDs by gene),
`nogold/<GENE>/` (per-gene gold-free rankings), and
`nogold/manual_acquisition_worklist_nogold_top200.{csv,md}` (merged by predicted
expected value; `nogold/run_all.sh` regenerates them).
