# General repository PDF fallback — 2026-09-05

The manual recovery is now implemented in the normal downloader and the
standalone paywall-recovery tool. Both requested papers succeed from PMID alone,
without supplied repository URLs. A fixed panel of **14 other gold-standard
papers** yielded **one additional successful PDF download**; the other 13 remain
explicit failures. This is acquisition evidence, not a new extraction score.

## How the manual downloads worked

The previous audit searched the article title/DOI, followed an author-repository
record, downloaded its PDF, verified the article identity and inspected the
count-bearing page. It used:

- PMID **20031634**, DOI `10.1161/CIRCGENETICS.109.853374`:
  [HAL author deposit](https://univ-brest.hal.science/hal-00750425/document).
  The publisher returned 403, CiteSeer timed out and CORE returned 403.
- PMID **25163546**, DOI `10.1093/eurheartj/ehu301`:
  [UCL manuscript](https://discovery.ucl.ac.uk/id/eprint/1447311/1/Syrris_Haas%20Eur%20Heart%20J.pdf).

The PDFs were converted with `pdftotext -layout` into the standard flat
`<PMID>_FULL_CONTEXT.md`, `<PMID>_CLEANED.md` and artifact-sidecar layout. The
existing `scripts/build_source_corpus.py` dry run was inspected, then `--apply`
copied the source upgrades into the external `corpus/` and updated its index.
No SQLite count rows were patched. The [prior audit](../phenotype_failure_panel_20260905/README.md)
retains its original PDF hashes and the 13-row Table 1 transcription.

## Why normal execution missed the copies

The [historical failure rows](historical_failures.json) come from the immutable
`20260903_protocol_cont120_01_candidate` production artifacts.

1. **The free-access flag bypassed repository discovery.**
   `initialize_free_text_access()` returned immediately when PubMed marked a
   paper free. Unpaywall was only consulted for papers *not* marked free. Both
   examples took the publisher path; a failed publisher did not send them back
   to the OA-index lookup.
2. **Publisher/browser failure became an abstract.** For 20031634 DOI
   resolution failed and PubMed supplied no usable free URL. The error then
   incorrectly said “No DOI or URL” even though its DOI was present. For
   25163546 the OUP article URL returned 403. Both browser strategies returned
   empty bodies. The normal harvester then wrote abstract-only context.
3. **Even the old Unpaywall success path lost its converted text.** It fetched
   one PDF, converted it, discarded the markdown and returned the PDF URL as
   `free_url`. The next stage treated that URL as HTML. It also did not keep
   trying alternate copies after the selected PDF failed, and did not follow
   repository landing pages without direct PDF URLs.
4. **The separate Scholar path did not cover these repository identifiers.**
   It was not the normal harvester's last resort, and its source gate requires
   the PDF URL itself to contain the requested DOI/PMID. HAL and UCL record
   URLs do not. That gate is unchanged; the new path instead verifies the
   index's exact DOI relationship and the attached article's header.

## What changed

`harvesting/repository_pdf_fallback.py` is shared by:

- `PMCHarvester`: after publisher and any enabled browser fallback fail; also when
  PMC has no body/DOI or the returned PMC content fails validation.
- `scripts/fetch_paywalled.py`: before the final Scholar attempt, including
  no-matching-strategy and supplement-only outcomes. Already obtained
  supplement text and artifact metadata are preserved when adding a body.

The fallback resolves missing bibliographic fields by exact PMID through
Europe PMC, tries **all Unpaywall locations**, then **OpenAlex locations**, then
HAL's DOI search. A repository landing page may supply a public PDF even when
the index's direct-PDF field or OA flag is stale. Links come from indexed
records and their HTML; no paper-specific URL guesses or gold identities are
used. See [OpenAlex locations](https://help.openalex.org/data/locations/) and
[HAL search documentation](https://api.archives-ouvertes.fr/docs/search).

Acceptance requires actual PDF magic bytes, readable substantive article text,
an exact DOI metadata relationship and a matching article title near the start
of the attached article. Repository citation/licence covers and bibliographies
cannot supply that title witness; an explicit conflicting header DOI rejects
the attachment. Exact, sufficiently specific main titles may match when one
index omits a colon-delimited subtitle. Fuzzy topic matching is not used.

Requests are bounded by a two-minute recovery budget, 14 asset requests
(redirects included), 20-second individual request timeouts and a 30 MiB
payload limit. Unpaywall retains its existing 30-second request timeout and
single empty-response retry. No browser challenge bypass or new model call is
part of this fallback.

Successful output retains the PDF, SHA-256, source URL, discovery chain,
identity decision and page boundaries. Conversion uses PyMuPDF's sorted page
text directly; it never routes the PDF back through the HTML scraper. The
artifact and FULL_CONTEXT explicitly retain **body-only / supplement surface
unavailable** status, so corpus reuse cannot mistake a recovered body for a
complete supplement set. The original PDF remains available for table/figure
inspection; searchable text alone does not prove table-cell fidelity.

## Live validation

The later [implementation campaign](../recall_campaign_20260905/README.md)
preserves nested caption metadata and caption text, handles missing/looping
redirects, cleans header DOI render glue, and preserves enabled browser
completion before this fallback. The 3/16 result below remains the original
fixed source-only measurement; fresh extraction results are reported separately
in that campaign. Some acquired supplement text never proves that the whole
supplement surface is complete.

The panel was frozen before the network test: the two requested PMIDs plus the
first 14 other distinct PMIDs from the existing two-cohort ranked problem-paper
report. [panel.json](panel.json) records the selection. No gold variant or
phenotype value was used as a discovery input. This intentionally difficult
panel is calibration evidence, not a representative estimate of all papers.

| PMID | Gene | Final result |
|---|---|---|
| 20031634 | SCN5A | HAL PDF, 7 pages; 46,499 characters of page text |
| 25163546 | SCN5A | Amsterdam UMC PDF, 16 pages; 158,001 characters of page text |
| 29650123 | KCNH2 | Repisalud author PDF, 13 pages; 43,551 characters of page text |
| 15898185 | SCN5A | Indexed handle returned 404 |
| 27114410 | RYR2 | No usable indexed PDF/record candidate |
| 31520628 | KCNQ1 | Linked PDF returned 403; other records had no usable PDF |
| 14678125 | KCNQ1 | No usable indexed PDF/record candidate |
| 19926015 | RYR2 | Publisher/PMC PDFs blocked; repository records had no PDF |
| 22677073 | RYR2 | PMC blocked; other repository record had no PDF |
| 29121719 | KCNH2 | No usable indexed PDF/record candidate |
| 24667783 | KCNQ1 | Publisher URL returned HTML; repository/PMC copies blocked |
| 29350269 | RYR2 | Indexed ZORA record returned 404 |
| 19687231 | KCNQ1 | Indexed PMC copy blocked |
| 22966897 | SCN5A | No usable indexed PDF/record candidate |
| 10528853 | BRCA1 | BMJ and PMC copies blocked |
| 17905336 | KCNH2 | Repository record had no PDF link |

Thus **3/16 PDFs downloaded**, comprising **2/2 requested reproductions and
1/14 other papers**. These are downloads, not 3 new corpus recoveries. The
corpus already had a good body for 29650123; its missing mutation roster
remains unresolved. Its existing corpus source was retained. The larger
25163546 published copy upgrades the prior UCL manuscript rendering; the
already-recovered 20031634 rendering is also retained by the corpus builder.

The two requested cases additionally ran through the real
`PMCHarvester.download_pmid()` entry point, from fresh output directories,
with ordinary publisher failures and no injected download URLs. Both returned
repository body text: [normal-entrypoint results](normal_entrypoint_results.json).
This verifies more than calling the new helper directly.

PMID 20031634's new source retains Table 1's headers, 13 family rows and the
printed total **115 carriers / 54 BrS-ECG+ / 55 not BrS-ECG+ / 6 undetermined**;
eight ECG-positive noncarriers are a separate column. The downloaded PDF hash
matches the visually verified HAL paper from the previous audit. The Amsterdam
UMC copy's article body/methods were also checked visually. Neither check
claims that all clinical supplements have now been acquired.

## Reproduction and evidence

```bash
.venv/bin/python scripts/recall_audit/evaluate_repository_fallback.py \
  --panel docs/evidence/repository_fallback_20260905/panel.json \
  --out-dir validation_runs/repository_fallback_20260905/new_replay \
  --email brett.kroncke@gmail.com --workers 3
```

Use a fresh directory; the driver refuses to overwrite evidence. Then run the
ordinary corpus builder, inspect its dry run, and apply eligible upgrades:

```bash
.venv/bin/python scripts/build_source_corpus.py \
  --roots validation_runs/repository_fallback_20260905/new_replay/sources \
  --gene SCN5A --assume-gene SCN5A
# Repeat the same command with --apply after reviewing the dry run.
```

Durable outputs: [results.csv](results.csv), [manifest.json](manifest.json),
[discovery_traces.json](discovery_traces.json), and the historical/entry-point
records linked above. PDFs, converted source, rendered checks, full test logs
and corpus-sync logs remain under
`validation_runs/repository_fallback_20260905/`. Earlier development replays
are retained; `panel_final` is the reported runtime and cohort.

The acceptance surface includes alternate-location continuation, real PDF
checks, exact metadata binding, title/subtitle handling, misleading covers,
wrong attachments, citing-paper rejection, public landing links, missing DOI
recovery, private redirect rejection, both entry points, invalid PMC bodies,
and preservation of existing supplements. The complete offline suite passes
**2,865 tests**; Ruff and whitespace checks pass.
[verification.json](verification.json) binds the final runtime, verifies the
historical lock/gold inputs, and checks the corpus upgrade; the
[corpus-sync record](corpus_sync.txt) preserves dry-run/apply outcomes.
No scored extraction arm, gold standard, locked prediction, count guard or
accepted recall/MAE changes in this work. The
[current phenotype figure](../../figures/evaluated_phenotype_counts/phenotype_count_recovery_stratified.png)
and [run manifest](../../figures/evaluated_phenotype_counts/phenotype_count_recovery_stratified.json)
remain the previously evaluated protocol evidence, not a score for this source
fallback.
