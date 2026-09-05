**Do not promote the bundled candidate on this $100, and do not treat paired 120 discovery+confirmation as the registered second tranche.** Merge only acquisition/QC that cannot change accepted outputs on already-complete frozen papers. Keep (a)(b)(c) off the headline until a mechanism check, a single paired frozen-120 contrast, and an actual unopened set exist. Acquisition of 20031634/25163546 is a different claim and must stay out of the extraction delta.

## What is actually being tested

The registered gate (recall Δ ≥ 1pp, 95% cluster-bootstrap recall LB ≥ −1pp, precision LB ≥ −2pp, then the identical candidate on the next unopened tranche) is a two-stage promotion rule. **Re-running the same 120 is a noise replicate, not a tranche.** The rules were written from named misses (R594Q, a lost paired coordinate, excluded-from-figure prose, Q2958R). Confirmation on that same gold is almost powerless against overfitting.

Cluster-bootstrap on one run captures **which papers are in the 120**, not LLM jitter. You already know jitter is ±1–2 TP. If 1pp is a few TPs, a single A120 can pass the written LB and fail the next seed. Deterministic source-grounded replay is required for the point estimate; replication is required for residual model noise. Those are not interchangeable.

Frozen source bytes and no cached DB/predictions are the right scored-arm constraints. The snapshotted current-prechange runtime is the baseline **code**, not a cached score. Both arms must re-extract from the same bytes. If that freeze predates 20031634’s 13 rows, keep those papers out of the extraction contrast or freeze a post-acquisition corpus and run **both** arms on it so acquisition is not credited to (a)(b)(c).

## Ownership and precision bugs in the three rules

**(a) Nameless rows with complete one-letter `source_notation` (R594Q) gain protein identity.** This is matcher/identity fill, not new extraction. Precision rests entirely on gene and range guards. Complete `R594Q` is still not unique across proteins; without a gene already bound on the row (not inferred from the paper), this can attach the right token to the wrong protein. If “model rows” are constructs/in vitro/computational rather than patient alleles, filling protein identity imports non-clinical rows into the clinical identity space. Equivalent notations (`p.R594Q`, `Arg594Gln`) are a recall hole, not a precision control. **Score it as an identity-layer change.** If it only fires when gene is already bound and residue is in-range, it is the least dangerous of the three; it is still not an extraction win.

**(b) One protein column paired with one labelled nucleotide-change/cDNA column; bare `14803 G>A` → `c.14803G>A`.** Ownership is correct (table emission; the generic cleaner already accepted well-formed cDNA). The precision bugs are interpretation, not pairing:

- Forcing `c.` is a transcript claim. A “nucleotide change” header can be genomic, mitochondrial, or non-cDNA numbering. Rejecting columns labelled genomic is necessary and not sufficient; **unlabelled or weakly labelled position columns should not coerce `c.`**.
- “Exactly one protein” / “exactly one cDNA column” is table-schema, not row alignment. Merged cells, `14803 G>A / 14804 C>T`, and allele pairs in one cell still pass a one-and-one column check.
- Rejecting multi-protein, multi-allele, cell, genomic, and BRCA-legacy tables is a precision hedge that can also drop true paired coordinates. BRCA-legacy as a blanket reject is a gene-specific exception: either it is a documented FP factory, or it is an unmeasured recall leak.

This rule can be made merge-adjacent only for the **reject** side, and only after calibration shows no lost TPs on the motivating paired-coordinate paper. The `c.` coercion is the part that must not ship on schema pairing alone.

**(c) Identity-only recovery of `VARIANT-GENE [n cases]` in prose that says excluded from Figure/Table.** This is the precision hazard and a count-recovery rule in identity clothing. `[n cases]` is a count; if `n` never becomes affected/unaffected/carrier, say so in the emitter and test it. Remaining bugs:

- “Excluded from Figure/Table” often means **studied but not plotted**. If gold includes those alleles, recovery is recall. If gold is figure/table-only, this is systematic over-emission. That is a gold-definition dependency, not a small regex edge.
- Exact adjacent gene is not binding. Two genes in the same clause will attach the integer to the wrong gene.
- Token+paragraph evidence is auditable and not a uniqueness key. The same allele can be table-extracted and prose-repeated; without dedupe against table identities you double-count.
- You already know the rule does **not** get uncounted Q2958R and does **not** take “excluded from study.” So the motivating FN that would justify a headline move is still out of scope. Do not promote a partial discourse heuristic because it is adjacent to a named miss.

**(c) should not be in the same candidate as (a) and (b).** It can spend the entire −2pp precision allowance on a few papers and still look acceptable in a global LB.

## Acquisition is a separate, still-unextracted claim

General-repository fallback plus DOI+header-title match, cover excluded, body-only/retry status is the only clearly merge-safe packet, **if** it cannot overwrite a paper that already has usable body/supplement bytes.

Facts you have: 20031634 now has 13 source variant rows (old run 0 TP / 13 FN); 25163546 body landed, supplement roster did not; 29650123 already had body, no new roster; **no model has extracted the new sources**. Source-on-disk is not TP. 25163546’s gold may still live in the missing supplement; scoring body-only as a full paper will look like extraction failure. Keep body-only/retry papers out of frozen scored arms or mark them ineligible.

## Discovery design does not have the power the gate needs

A120 at $10–12 and a $50–60 “paired 120 discovery+confirmation” is four to six full arms on the **same** papers. That can estimate seed noise. It cannot satisfy “next unopened tranche.” It also cannot separate (a), (b), (c), and acquisition.

Calibration `<$5` on opened motivating sources is the only high-power spend in the proposal — **as mechanism tests**, not as a recall point. If those PMIDs are inside gold-120, you have already peeked for rule-writing; spending them again does not create independence.

Power against a 1pp recall move with ±1–2 TP model noise is inadequate without deterministic replay. Power against a precision failure localized to (c) is inadequate with a global −2pp LB. Power against paper-specific overfit is approximately zero if discovery and confirmation share gold-120.

## Best practical use of $100

Treat $100 as **one paired frozen contrast, cheap mechanism kills, noise bound, and a wall for a true holdout** — not as two headline 120s.

1. **Freeze now (cost $0).** Snapshot candidate code vs current-prechange runtime. Freeze source bytes **excluding** newly acquired text from the extraction contrast, or freeze post-acquisition bytes and give them to both arms. Record body-only/retry. No cached predictions.

2. **Opened calibration, cap ~$5, sequential kill switches, not a bundle.** On the motivating frozen papers only, deterministic replay if the path is regex/matcher; paid LLM only if (c) is a model call. Pass/fail per rule:
   - (a) fires iff complete one-letter notation **and** gene already bound **and** residue in range; does not fire on nameless in-vitro/model constructs.
   - (b) emits the lost paired coordinate **without** `c.` unless the column label is explicitly cDNA/`c.`/transcript; no emit on multi-cell/multi-allele; BRCA-legacy reject does not drop a true pair on that paper.
   - (c) emits identity with original token+paragraph, gene string-equal adjacent, positive integer literal present, **no** phenotype/carrier fields, **no** “excluded from study,” **no** Q2958R-if-uncounted, **no** identity that already exists from the table.
   Kill or rewrite the failing rule before any A120. If (c) cannot be shown non-duplicative, drop it from the candidate.

3. **One paired A120, ~$22–24.** Baseline = snapshotted prechange. Candidate = **surviving** rules only. Same frozen bytes. Deterministic replay for regex/identity; if any LLM remains, lock prompt hash and decoding. Report extraction Δ only. Do not add 20031634’s new rows to this Δ. Apply the registered cluster-bootstrap **and** treat ±1–2 TP as a separate model-noise band: if Δ ≥ 1pp but the lower noise band is < 1pp, you do not have a discovery pass.

4. **Ablation on delta PMIDs only (~$0–10).** Replay papers that changed. Attribute TPs/FPs to (a), (b), or (c). If (c) owns new FPs or (a) owns TPs that are matcher-only, unbundle before any promotion language.

5. **Replicate / reserve ~$25–40, in this order.**
   - If discovery fails the gate or sits inside noise: **stop**. Spend at most one candidate replicate ($10–12) to confirm it was noise, then keep the rest unspent.
   - If discovery clears the gate **and** ablation says the gain is (a) and/or a non-coercing (b): spend one candidate replicate on the same frozen 120 to bound jitter. That still is not tranche 2.
   - Remaining money: either extract 20031634 (and 25163546 body) as a **labelled acquisition arm** with no headline claim, or buy the smallest **unopened** gold you have (other gene, unused PMIDs). If no unopened gold exists, **you cannot complete the registered rule with this $100.** Do not spend the reserve to reconfirm gold-120.

6. **Merge path that preserves the accepted headline.**
   - Merge now: repository fallback, DOI+title PDF match, cover exclusion, body-only/retry status, write guard that does not replace complete papers.
   - Behind a flag, off default: (a) with gene-bound+range; (b) pairing **without** `c.` coercion; (c) not on by default.
   - Headline / default `gvf-run` path unchanged until an unopened tranche passes the written rule on the identical candidate.

## Uncertainties (do not fill with a score)

- How many TPs 1pp is on this denominator, so whether ±1–2 TP swallows the gate.
- Whether cluster-bootstrap is by PMID (correct) or by variant (too optimistic).
- Whether the prechange snapshot predates 20031634/25163546 bytes.
- Whether (a)/(c) are deterministic or LLM; only (b)’s regex side is clearly replayable.
- Whether gold-120 includes “excluded from figure/table” alleles.
- Whether 25163546’s roster is in the missing supplement.
- Whether any true unopened tranche exists; if not, the promotion rule is not completable this campaign.
- What “nameless model rows” are (missing fields vs biochemical models).

**Implementable line:** spend ≤$5 to kill or unbundle (c) and the `c.` coercion; spend ~$24 on one frozen paired 120 of the remainder; use ~$12 only if that Δ is above both the 1pp line and the jitter band; put acquisition in a separate labelled arm; leave $25+ unspent unless a real holdout exists. Anything that calls two gold-120 runs “discovery and confirmation” spends the budget on a replicate and still leaves the headline unprotected.
