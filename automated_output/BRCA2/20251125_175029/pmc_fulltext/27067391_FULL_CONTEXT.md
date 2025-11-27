# MAIN TEXT

## A unified analytic framework for prioritization of non-coding variants of uncertain significance in heritable breast and ovarian cancer

### Abstract

BackgroundSequencing of both healthy and disease singletons yields many novel and low frequency variants of uncertain significance (VUS). Complete gene and genome sequencing by next generation sequencing (NGS) significantly increases the number of VUS detected. While prior studies have emphasized protein coding variants, non-coding sequence variants have also been proven to significantly contribute to high penetrance disorders, such as hereditary breast and ovarian cancer (HBOC). We present a strategy for analyzing different functional classes of non-coding variants based on information theory (IT) and prioritizing patients with large intragenic deletions.MethodsWe captured and enriched for coding and non-coding variants in genes known to harbor mutations that increase HBOC risk. Custom oligonucleotide baits spanning the complete coding, non-coding, and intergenic regions 10 kb up- and downstream of ATM, BRCA1, BRCA2, CDH1, CHEK2, PALB2, and TP53 were synthesized for solution hybridization enrichment. Unique and divergent repetitive sequences were sequenced in 102 high-risk, anonymized patients without identified mutations in BRCA1/2. Aside from protein coding and copy number changes, IT-based sequence analysis was used to identify and prioritize pathogenic non-coding variants that occurred within sequence elements predicted to be recognized by proteins or protein complexes involved in mRNA splicing, transcription, and untranslated region (UTR) binding and structure. This approach was supplemented by in silico and laboratory analysis of UTR structure.Results15,311 unique variants were identified, of which 245 occurred in coding regions. With the unified IT-framework, 132 variants were identified and 87 functionally significant VUS were further prioritized. An intragenic 32.1 kb interval in BRCA2 that was likely hemizygous was detected in one patient. We also identified 4 stop-gain variants and 3 reading-frame altering exonic insertions/deletions (indels).ConclusionsWe have presented a strategy for complete gene sequence analysis followed by a unified framework for interpreting non-coding variants that may affect gene expression. This approach distills large numbers of variants detected by NGS to a limited set of variants prioritized as potential deleterious changes.Electronic supplementary materialThe online version of this article (doi:10.1186/s12920-016-0178-5) contains supplementary material, which is available to authorized users.

### Background

Advances in NGS have enabled panels of genes, whole exomes, and even whole genomes to be sequenced for multiple individuals in parallel. These platforms have become so cost-effective and accurate that they are beginning to be adopted in clinical settings, as evidenced by recent FDA approvals [1, 2]. However, the overwhelming number of gene variants revealed in each individual has challenged interpretation of clinically significant genetic variation [3–5].

After common variants, which are rarely pathogenic, are eliminated, the number of VUS in the residual set remains substantial. Assessment of pathogenicity is not trivial, considering that nearly half of the unique variants are novel, and cannot be resolved using published literature and variant databases [6]. Furthermore, loss-of-function variants (those resulting in protein truncation are most likely to be deleterious) represent a very small proportion of identified variants. The remaining variants are missense and synonymous variants in the exon, single nucleotide changes, or in frame insertions or deletions in intervening and intergenic regions. Functional analysis of large numbers of these variants often cannot be performed, due to lack of relevant tissues, and the cost, time, and labor required for each variant. Another problem is that in silico protein coding prediction tools exhibit inconsistent accuracy and are thus problematic for clinical risk evaluation [7–9]. Consequently, many HBOC patients undergoing genetic susceptibility testing will receive either an inconclusive (no BRCA variant identified) or an uncertain (BRCA VUS) result. The former has been reported in up to 80 % of cases and depends on the number of genes tested [10]. The occurrence of uncertain BRCA mutations varies greatly (as high as 46 % in African American populations and as low as 2.1 %) among tested individuals depending on the laboratory and the patient’s ethnicity [11–13]. The inconsistency in diagnostic yield is significant, considering that HBOC accounts for 5–10 % of all breast/ovarian cancer [14, 15].

One strategy to improve variant interpretation in patients is to reduce the full set of variants to a manageable list of potentially pathogenic variants. Evidence for pathogenicity of VUS in genetic disease is often limited to amino acid coding changes [16, 17], and mutations affecting splicing, transcriptional activation, and mRNA stability tend to be underreported [18–24]. Splicing errors are estimated to represent 15 % of disease-causing mutations [25], but may be much higher [26, 27]. The impact of a single nucleotide change in a recognition sequence can range from insignificant to complete abolition of a protein binding site. Aberrant splicing events causing frameshifts often disrupt protein function; in-frame changes are dependent on gene context. The complexity of interpretation of non-coding sequence variants benefits from computational approaches [28] and direct functional analyses [29–33] that may each support evidence of pathogenicity.

Ex vivo transfection assays developed to determine the pathogenicity of VUS predicted to lead to splicing aberrations (using in silico tools) have been successful in identifying pathogenic sequence variants [34, 35]. IT-based analysis of splicing variants has proven to be robust and accurate (as determined by functional assays for mRNA expression or binding assays) at analyzing splice site (SS) variants, including splicing regulatory factor binding sites (SRFBSs), and in distinguishing them from polymorphisms in both rare and common diseases [36–39]. However, IT can be applied to any sequence recognized and bound by another factor [40], such as with transcription factor binding sites (TFBSs) and RNA-binding protein binding sites (RBBSs). IT is used as a measure of sequence conservation and is more accurate than consensus sequences [41]. The individual information (Ri) of a base is related to thermodynamic entropy, and therefore free energy of binding, and is measured on a logarithmic scale (in bits). By comparing the change in information (ΔRi) for a nucleotide variation of a bound sequence, the resulting change in binding affinity is ≥ 2ΔRi, such that a 1 bit change in information will result in at least a 2-fold change in binding affinity [42].

IT measures nucleotide sequence conservation and does not provide information on effects of variants on mRNA secondary (2°) structure, nor can it accurately predict effects of amino acid sequence changes. Associations of structural changes in untranslated regions (UTR) of mRNA with disease justifies including predicted effects of these changes on 2° structure in the comprehensive analysis of sequence variants [43]. Other in silico methods have attempted to address these deficiencies. For example, Halvorsen et al. (2010) introduced an algorithm called SNPfold, which computes the potential effect of a single nucleotide variant (SNV) on mRNA 2° structure [20]. Predictions made by SNPfold can be tested by the SHAPE assay (Selective 2’-Hydroxyl Acylation analyzed by Primer Extension) [44], which provides evidence for sequence variants that lead to structural changes in mRNA by detection of covalent adducts in mRNA.

The implications of improved VUS interpretation are particularly relevant for HBOC due to its incidence and the adoption of panel testing for these individuals [45, 46]. It has been suggested that patients with a high risk profile receiving uninformative results would imply that deleterious variants lie in untested regions of BRCA1/2, untested genes, or are unrecognized [47, 48]. This is also supported by studies where families with linkage to BRCA1/2 had no detectable pathogenic mutation (however it is noteworthy that detection rates of BRCA mutations in families with documented linkage to these loci appears to vary by ascertainment, inclusion criteria, and technology used to identify the mutations) [49, 50]. The concept of non-BRCA gene association has been demonstrated by the identification of low-to-moderate risk HBOC genes, and variants within coding and non-coding regions affecting splicing and regulatory factor binding [51, 52]. Consequently, VUS, which include rare missense changes, other coding and non-coding changes in all of these genes, greatly outnumber the catalog of known deleterious mutations [53].

Here, we develop and evaluate IT-based models to predict potential non-coding sequence mutations in SSs, TFBSs, and RBBSs in 7 genes sequenced in their entirety. These models were used to analyze 102 anonymous HBOC patients who did not exhibit known BRCA1/2 coding mutations at the time of initial testing, despite meeting the criteria for BRCA genetic testing. The genes are: ATM, BRCA1, BRCA2, CDH1, CHEK2, PALB2, and TP53, and have been reported to harbor mutations that increase HBOC risk [54–76]. We apply these IT-based methods to analyze variants in the complete sequences of coding, non-coding, and up- and downstream regions of the 7 genes. In this study, we established and applied a unified IT-based framework, first filtering out common variants, then to “flag” potentially deleterious ones. Then, using context-specific criteria and information from the published literature, we prioritized likely candidates.

### Methods

Nucleic acid hybridization capture reagents designed from genomic sequences generally avoid repetitive sequence content to avoid cross hybridization [77]. Complete gene sequences harbor numerous repetitive sequences, and an excess of denatured C0t-1 DNA is usually added to hybridization to prevent inclusion of these sequences [78]. RepeatMasker software completely masks all repetitive and low-complexity sequences [79]. We increased sequence coverage in complete genes with capture probes by enriching for both single-copy and divergent repeat (>30 % divergence) regions, such that, under the correct hybridization and wash conditions, all probes hybridize only to their correct genomic locations [77]. This step was incorporated into a modified version of Gnirke and colleagues’ (2009) in-solution hybridization enrichment protocol, in which the majority of library preparation, pull-down, and wash steps were automated using a BioMek® FXP Automation Workstation (Beckman Coulter, Mississauga, Canada) [80].

Genes ATM (RefSeq: NM_000051.3, NP_000042.3), BRCA1 (RefSeq: NM_007294.3, NP_009225.1), BRCA2 (RefSeq: NM_000059.3, NP_000050.2), CDH1 (RefSeq: NM_004360.3, NP_004351.1), CHEK2 (RefSeq: NM_145862.2, NP_665861.1), PALB2 (RefSeq: NM_024675.3, NP_078951.2), and TP53 (RefSeq: NM_000546.5, NP_000537.3) were selected for capture probe design by targeting single copy or highly divergent repeat regions (spanning 10 kb up- and downstream of each gene relative to the most upstream first exon and most downstream final exon in RefSeq) using an ab initio approach [77]. If a region was excluded by ab initio but lacked a conserved repeat element (i.e. divergence > 30 %) [79], the region was added back into the probe-design sequence file. Probe sequences were selected using PICKY 2.2 software [81]. These probes were used in solution hybridization to capture our target sequences, followed by NGS on an Illumina Genome Analyzer IIx (Additional file 1: Methods).

Genomic sequences from both strands were captured using overlapping oligonucleotide sequence designs covering 342,075 nt among the 7 genes (Fig. 1). In total, 11,841 oligonucleotides were synthesized from the transcribed strand consisting of the complete, single copy coding, and flanking regions of ATM (3513), BRCA1 (1587), BRCA2 (2386), CDH1 (1867), CHEK2 (889), PALB2 (811), and TP53 (788). Additionally, 11,828 antisense strand oligos were synthesized (3497 ATM, 1591 BRCA1, 2395 BRCA2, 1860 CDH1, 883 CHEK2, 826 PALB2, and 776 TP53). Any intronic or intergenic regions without probe coverage are most likely due to the presence of conserved repetitive elements or other paralogous sequences.Fig. 1Capture Probe Coverage over Sequenced Genes. The genomic structure of the 7 genes chosen are displayed with the UCSC Genome Browser. Top row for each gene is a custom track with the “dense” visualization modality selected with black regions indicating the intervals covered by the oligonucleotide capture reagent. Regions without probe coverage contain conserved repetitive sequences or correspond to paralogous sequences that are unsuitable for probe design

Capture Probe Coverage over Sequenced Genes. The genomic structure of the 7 genes chosen are displayed with the UCSC Genome Browser. Top row for each gene is a custom track with the “dense” visualization modality selected with black regions indicating the intervals covered by the oligonucleotide capture reagent. Regions without probe coverage contain conserved repetitive sequences or correspond to paralogous sequences that are unsuitable for probe design

For regions lacking probe coverage (of ≥ 10 nt, N = 141; 8 in ATM, 26 in BRCA1, 10 in BRCA2, 29 in CDH1, 36 in CHEK2, 15 in PALB2, and 17 in TP53), probes were selected based on predicted Tms similar to other probes, limited alignment to other sequences in the transcriptome (<10 times), and avoidance of stable, base-paired 2° structures (with unaFOLD) [82, 83]. The average coverage of these sequenced regions was 14.1–24.9 % lower than other probe sets, indicating that capture was less efficient, though still successful.

Genomic DNA from 102 patients previously tested for inherited breast/ovarian cancer without evidence of a predisposing genetic mutation, was obtained from the Molecular Genetics Laboratory (MGL) at the London Health Sciences Centre in London, Ontario, Canada. Patients qualified for genetic susceptibility testing as determined by the Ontario Ministry of Health and Long-Term Care BRCA1 and BRCA2 genetic testing criteria [84] (see Additional file 2). The University of Western Ontario research ethics board (REB) approved this anonymized study of these individuals to evaluate the analytical methods presented here. BRCA1 and BRCA2 were previously analyzed by Protein Truncation Test (PTT) and Multiplex Ligation-dependent Probe Amplification (MLPA). The exons of several patients (N = 14) had also been Sanger sequenced. No pathogenic sequence change was found in any of these individuals. In addition, one patient with a known pathogenic BRCA variant was re-sequenced by NGS as a positive control.

Variant analysis involved the steps of detection, filtering, IT-based and coding sequence analysis, and prioritization (Fig. 2). Sequencing data were demultiplexed and aligned to the specific chromosomes of our sequenced genes (hg19) using both CASAVA (Consensus Assessment of Sequencing and Variation; v1.8.2) [85] and CRAC (Complex Reads Analysis and Classification; v1.3.0) [86] software. Alignments were prepared for variant calling using Picard [87] and variant calling was performed on both versions of the aligned sequences using the UnifiedGenotyper tool in the Genome Analysis Toolkit (GATK) [88]. We used the recommended minimum phred base quality score of 30, and results were exported in variant call format (VCF; v4.1). A software program was developed to exclude variants called outside of targeted capture regions and those with quality scores < 50. Variants flagged by bioinformatic analysis (described below) were also assessed by manually inspecting the reads in the region using the Integrative Genomics Viewer (IGV; version 2.3) [89, 90] to note and eliminate obvious false positives (i.e. variant called due to polyhomonucleotide run dephasing, or PCR duplicates that were not eliminated by Picard). Finally, common variants (≥1 % allele frequency based on dbSNP 142 or > 5 individuals in our study cohort) were not prioritized.Fig. 2Framework for the Identification of Potentially Pathogenic Variants. Integrated laboratory processing and bioinformatic analysis procedures for comprehensive complete gene variant determination and analysis. Intermediate datasets resulting from filtering are represented in yellow and final datasets in green. Non-bioinformatic steps, such as sample preparation are represented in blue and prediction programs in purple. Sequencing analysis yields base calls for all samples. CASAVA [85] and CRAC [86] were used to align these sequencing results to hg19. GATK [88] was used to call variants from this data against GRCh37 release of the reference human genome. Variants with a quality score < 50 and/or call confidence score < 30 were eliminated along with variants falling outside of our target regions. SNPnexus [112–114] was used to identify the genomic location of the variants. Nonsense and indels were noted and prediction tools were used to assess the potential pathogenicity of missense variants. The Shannon Pipeline [91] evaluated the effect of a variant on natural and cryptic SSs, as well as SRFBSs. ASSEDA [38] was used to predict the potential isoforms as a result of these variants. PWMs for 83 TFs were built using an information weight matrix generator based on Bipad [106]. Mutation Analyzer evaluated the effect of variants found 10 kb upstream up to the first intron on protein binding. Bit thresholds (R
i values) for filtering variants on software program outputs are indicated. Variants falling within the UTR sequences were assessed using SNPfold [20], and the most probable variants that alter mRNA structure (p < 0.1) were then processed using mFold to predict the effect on stability [83]. All UTR variants were scanned with a modified version of the Shannon Pipeline, which uses PWMs computed from nucleotide frequencies for 28 RBPs in RBPDB [109] and 76 RBPs in CISBP-RNA [110]. All variants meeting these filtering criteria were verified with IGV [89, 90]. *Sanger sequencing was only performed for protein truncating, splicing, and selected missense variants

Framework for the Identification of Potentially Pathogenic Variants. Integrated laboratory processing and bioinformatic analysis procedures for comprehensive complete gene variant determination and analysis. Intermediate datasets resulting from filtering are represented in yellow and final datasets in green. Non-bioinformatic steps, such as sample preparation are represented in blue and prediction programs in purple. Sequencing analysis yields base calls for all samples. CASAVA [85] and CRAC [86] were used to align these sequencing results to hg19. GATK [88] was used to call variants from this data against GRCh37 release of the reference human genome. Variants with a quality score < 50 and/or call confidence score < 30 were eliminated along with variants falling outside of our target regions. SNPnexus [112–114] was used to identify the genomic location of the variants. Nonsense and indels were noted and prediction tools were used to assess the potential pathogenicity of missense variants. The Shannon Pipeline [91] evaluated the effect of a variant on natural and cryptic SSs, as well as SRFBSs. ASSEDA [38] was used to predict the potential isoforms as a result of these variants. PWMs for 83 TFs were built using an information weight matrix generator based on Bipad [106]. Mutation Analyzer evaluated the effect of variants found 10 kb upstream up to the first intron on protein binding. Bit thresholds (R
i values) for filtering variants on software program outputs are indicated. Variants falling within the UTR sequences were assessed using SNPfold [20], and the most probable variants that alter mRNA structure (p < 0.1) were then processed using mFold to predict the effect on stability [83]. All UTR variants were scanned with a modified version of the Shannon Pipeline, which uses PWMs computed from nucleotide frequencies for 28 RBPs in RBPDB [109] and 76 RBPs in CISBP-RNA [110]. All variants meeting these filtering criteria were verified with IGV [89, 90]. *Sanger sequencing was only performed for protein truncating, splicing, and selected missense variants

All variants were analyzed using the Shannon Human Splicing Mutation Pipeline, a genome-scale variant analysis program that predicts the effects of variants on mRNA splicing [91, 92]. Variants were flagged based on criteria reported in Shirley et al. (2013): weakened natural site ≥ 1.0 bits, or strengthened cryptic site (within 300 nt of the nearest exon) where cryptic site strength is equivalent or greater than the nearest natural site of the same phase [91]. The effects of flagged variants were further analyzed in detail using the Automated Splice Site and Exon Definition Analysis (ASSEDA) server [38].

Exonic variants and those found within 500 nt of an exon were assessed for their effects, if any, on SRFBSs [38]. Sequence logos for splicing regulatory factors (SRFs) (SRSF1, SRSF2, SRSF5, SRSF6, hnRNPH, hnRNPA1, ELAVL1, TIA1, and PTB) and their Rsequence values (the mean information content [93]) are provided in Caminsky et al. (2015) [36]. Because these motifs occur frequently in unspliced transcripts, only variants with large information changes were flagged, notably those with a) ≥ 4.0 bit decrease, i.e. at least a 16-fold reduction in binding site affinity, with Ri,initial ≥ Rsequence for the particular factor analyzed, or b) ≥ 4.0 bit increase in a site where Ri,final ≥ 0 bits. ASSEDA was used to calculate Ri,total, with the option selected to include the given SRF in the calculation. Variants decreasing Ri,total by < 3.0 bits (i.e. 8-fold) were predicted to potentially have benign effects on expression, and were not considered further.

Activation of pseudoexons through creating/strengthening of an intronic cryptic SS was also assessed [94]. Changes in intronic cryptic sites, where ΔRi > 1 bit and Ri,final ≥ (Rsequence – 1 standard deviation [S.D.] of Rsequence), were identified. A pseudoexon was predicted if a pre-existing cryptic site of opposite polarity (with Ri > [Rsequence - 1 S.D.]) and in the proper orientation for formation of exons between 10–250 nt in length was present. In addition, the minimum intronic distance between the pseudoexon and either adjacent natural exon was 100 nt. The acceptor site of the pseudoexon was also required to have a strong hnRNPA1 site located within 10 nt (Ri ≥ Rsequence) [38] to ensure accurate proofreading of the exon [37].

Next, variants affecting the strength of SRFs were analyzed by a contextual exon definition analysis of ΔRi,total. The context refers to the documented splicing activity of an SRF. For example, TIA1 has been shown to be an intronic enhancer of exon definition, so only intronic sites were considered. Similarly, hnRNPA1 proofreads the 3’ SS (acceptor) and inhibits exon recognition elsewhere [95]. Variants that lead to redundant SRFBS changes (i.e. one site is abolished and another proximate site [≤2 nt] of equivalent strength is activated) were assumed to have a neutral effect on splicing. If the strength of a site bound by PTB (polypyrimidine tract binding protein) was affected, its impact on binding by other factors was analyzed, as PTB impedes binding of other factors with overlapping recognition sites, but does not directly enhance or inhibit splicing itself [96].

To determine effects of variants on transcription factor (TF) binding, we first established which TFs bound to the sequenced regions of the gene promoters (and first exons) in this study by using ChIP-seq data from 125 cell types (Additional file 1: Methods) [97]. We identified 141 TFs with evidence for binding to the promoters of the genes we sequenced, including c-Myc, C/EBPβ, and Sp1, shown to transcriptionally regulate BRCA1, TP53, and ATM, respectively [98–100]. Furthermore, polymorphisms in TCF7L2, known to bind enhancer regions of a wide variety of genes in a tissue-specific manner [101], have been shown to increase risk of sporadic [102] and hereditary breast [103], as well as other types of cancer [104, 105].

IT-based models of the 141 TFs of interest were derived by entropy minimization of the DNase accessible ChIP-seq subsets [106]. Details are provided in Lu R, Mucaki E, and Rogan PK (BioRxiv; http://dx.doi.org/10.1101/042853). While some data sets would only yield noise or co-factor motifs (i.e. co-factors that bind via tethering, or histone modifying proteins [107]), techniques such as motif masking and increasing the number of Monte Carlo cycles yielded models for 83 TFs resembling each factor’s published motif. Additional file 3: Table S1 contains the final list of TFs and the models we built (described below) [108].

These TFBS models (N = 83) were used to scan all variants called in the promoter regions (10 kb upstream of transcriptional start site to the end of IVS1) of HBOC genes for changes in Ri. Binding site changes that weaken interactions with the corresponding TF (to Ri ≤ Rsequence) are likely to affect regulation of the adjacent target gene. Stringent criteria were used to prioritize the most likely variants and thus only changes to strong TFBSs (Ri,initial ≥ Rsequence), where reduction in strength was significant (ΔRi ≥ 4.0 bits), were considered. Alternatively, novel or strengthened TFBSs were also considered sources of dysregulated transcription. These sites were defined as having Ri,final ≥ Rsequence and as being the strongest predicted site in the corresponding genomic interval (i.e. exceeding the Ri values of adjacent sites unaltered by the variant). Variants were not prioritized if the TF was known to a) enhance transcription and IT analysis predicted stronger binding, or b) repress transcription and IT analysis predicted weaker binding.

Two complementary strategies were used to assess the possible impact of variants within UTRs. First, SNPfold software was used to assess the effect of a variant on 2° structure of the UTR (Additional file 1: Methods) [20]. Variants flagged by SNPfold with the highest probability of altering stable 2° structures in mRNA (where p-value < 0.1) were prioritized. To evaluate these predictions, oligonucleotides containing complete wild-type and variant UTR sequences (Additional file 4: Table S2) were transcribed in vitro and followed by SHAPE analysis, a method that can confirm structural changes in mRNA [44].

Second, the effects of variants on the strength of RBBSs were predicted. Frequency-based, position weight matrices (PWMs) for 156 RNA-binding proteins (RBPs) were obtained from the RNA-Binding Protein DataBase (RBPDB) [109] and the Catalog of Inferred Sequence Binding Preferences of RNA binding proteins (CISBP-RNA) [110, 111]. These were used to compute information weight matrices (based on the method described by Schneider et al. 1984; N = 147) (see Additional file 1: Methods) [40]. All UTR variants were assessed using a modified version of the Shannon Pipeline [91] containing the RBPDB and CISBP-RNA models. Results were filtered to include a) variants with |ΔRi| ≥ 4.0 bits, b) variants creating or strengthening sites (Ri,final ≥ Rsequence and the Ri,initial < Rsequence), and c) RBBSs not overlapping or occurring within 10 nt of a stronger, pre-existing site of another RBP.

The predicted effects of all coding variants were assessed with SNPnexus [112–114], an annotation tool that can be applied to known and novel variants using up-to-date dbSNP and UCSC human genome annotations. Variants predicted to cause premature protein truncation were given higher priority than those resulting in missense (or synonymous) coding changes. Missense variants were first cross referenced with dbSNP 142 [115]. Population frequencies from the Exome Variant Server [116] and 1000Genomes [117] are also provided. The predicted effects on protein conservation and function of the remaining variants were evaluated by in silico tools: PolyPhen-2 [118], Mutation Assessor (release 2) [119, 120], and PROVEAN (v1.1.3) [121, 122]. Default settings were applied and in the case of PROVEAN, the “PROVEAN Human Genome Variants Tool” was used, which includes SIFT predictions as a part of its output. Variants predicted by all four programs to be benign were less likely to have a deleterious impact on protein activity; however this did not exclude them from mRNA splicing analysis (described above in IT-Based Variant Analysis). All rare and novel variants were cross-referenced with general mutation databases (ClinVar [123, 124], Human Gene Mutation Database [HGMD] [125, 126], Leiden Open Variant Database [LOVD] [127–134], Domain Mapping of Disease Mutations [DM2] [135], Expert Protein Analysis System [ExPASy] [136] and UniProt [137, 138]), and gene-specific databases (BRCA1/2: the Breast Cancer Information Core database [BIC] [139] and Evidence-based Network for the Interpretation of Germline Mutant Alleles [ENIGMA] [140]; TP53: International Agency for Research on Cancer [IARC] [141]), as well as published reports to prioritize them for further workup.

Flagged variants were prioritized if they were likely to encode a dysfunctional protein (indels, nonsense codon > 50 amino acids from the C-terminus, or abolition of a natural SS resulting in out-of-frame exon skipping) or if they exceeded established thresholds for fold changes in binding affinity based on IT (see Methods above). In several instances, our classification was superseded by previous functional or pedigree analyses (reported in published literature or databases) that categorized these variants as pathogenic or benign.

We identified the BRCA1 exon 17 nonsense variant c.5136G > A (chr17:41215907C > T; rs80357418; 2-5A) [142] in the sample that was provided as a positive control. This was the same mutation identified by the MGL as pathogenic for this patient. We also prioritized another variant in this patient (Table 1) [143].Table 1Prioritized variants in the positive controlGenemRNA ProteinrsID (dbSNP 142)CategoryConsequenceRef
BRCA1
c.5136G > Ars80357418Nonsense151 AA short[142]p.Trp1712Ter
BRCA2
c.3218A > Grs80358566SRFBSRepressor action of hnRNPA1 at this site abolished (5.2 to 0.4 bits). Blocking action of PTB removed as site is abolished (5.5 to -7.5 bits) and may uncover binding sites of other SRFs.p.Gln1073ArgMissenseListed in ClinVar as conflicting interpretations (likely benign, unknown) and in BIC as unknown clinical importance. 2 in silico programs called deleterious. The variant occurs between repeat motifs BRC1 and BRC2 of BRCA2, a region in which pathogenic missense mutations have not yet been identified.[143]SRFBSRepressor action of hnRNPA1 at this site abolished (5.2 to 0.4 bits). Blocking action of PTB removed as site is abolished (5.5 to -7.5 bits) and may uncover binding sites of other SRFs.

Prioritized variants in the positive control

Protein-truncating, prioritized splicing, and selected prioritized missense variants were verified by Sanger sequencing. Primers of PCR amplicons are indicated in Additional file 5: Table S3.

Potential large rearrangements were detected with BreakDancer software [144], which identifies novel genomic rearrangements based on the respective orientation and distance between ends of the same read (and exceeding the lengths of NGS library inserts). This approach can, in theory, approximately localize deletions, duplications, or other types of breakpoints within exons, introns, and regulatory regions (eg. promoters) that could affect gene expression and function. We required at least 4 reads per suspected rearrangement in a sample separated by >700 nt, with each end mapping to proximate genomic reference coordinates to infer a potential deletion. Synthetic and cost limitations in the maximum genomic real estate covered by the capture reagent led to a tradeoff between extending the span of captured genomic intervals and higher tiling densities over shorter sequences, ie. exons, to achieve the level of coverage to reliably detect deletions based on read counts alone.

Our complete gene enrichment strategy with independent capture of both genomic strands enabled and facilitated development of a new algorithm to identify potential hemizygous genomic intervals in these individuals. In each subject, we first searched for contiguous long stretches (usually > > 1 kb) of non-polymorphic segments with diminished repetitive element content (<10 %), which is consistent with the possibility of these regions harboring a deletion. Then, we determined the likelihood of homo- or hemizygosity by comparing the degree of heterozygosity of variants in each of these intervals in for an individual with all of the other individuals sequenced with this protocol in this population. Regions containing haplotype blocks in strong linkage disequilibrium (LD; from HapMap [145]) were then excluded as candidate deletion intervals. Some individuals without a deletion are expected to be non-polymorphic, because detection of heterozygosity depends on genomic length of the region, marker informativeness, and the level of LD for those markers. We required that > 80 % of the control individuals be heterozyogous for at least two well-distributed loci within these intervals. Highly informative SNPs with a random genomic distribution in the controls (and other public databases) and which were non-polymorphic in the individual with the suspected deletion were weighted more heavily in inferring potential hemizygosity. This analysis was implemented using a Perl script that identified the most likely intervals of hemizygosity, which were then crossreferenced with the corresponding genomic intervals in HapMap.

### Design of tiled capture array for HBOC gene panel

Nucleic acid hybridization capture reagents designed from genomic sequences generally avoid repetitive sequence content to avoid cross hybridization [77]. Complete gene sequences harbor numerous repetitive sequences, and an excess of denatured C0t-1 DNA is usually added to hybridization to prevent inclusion of these sequences [78]. RepeatMasker software completely masks all repetitive and low-complexity sequences [79]. We increased sequence coverage in complete genes with capture probes by enriching for both single-copy and divergent repeat (>30 % divergence) regions, such that, under the correct hybridization and wash conditions, all probes hybridize only to their correct genomic locations [77]. This step was incorporated into a modified version of Gnirke and colleagues’ (2009) in-solution hybridization enrichment protocol, in which the majority of library preparation, pull-down, and wash steps were automated using a BioMek® FXP Automation Workstation (Beckman Coulter, Mississauga, Canada) [80].

Genes ATM (RefSeq: NM_000051.3, NP_000042.3), BRCA1 (RefSeq: NM_007294.3, NP_009225.1), BRCA2 (RefSeq: NM_000059.3, NP_000050.2), CDH1 (RefSeq: NM_004360.3, NP_004351.1), CHEK2 (RefSeq: NM_145862.2, NP_665861.1), PALB2 (RefSeq: NM_024675.3, NP_078951.2), and TP53 (RefSeq: NM_000546.5, NP_000537.3) were selected for capture probe design by targeting single copy or highly divergent repeat regions (spanning 10 kb up- and downstream of each gene relative to the most upstream first exon and most downstream final exon in RefSeq) using an ab initio approach [77]. If a region was excluded by ab initio but lacked a conserved repeat element (i.e. divergence > 30 %) [79], the region was added back into the probe-design sequence file. Probe sequences were selected using PICKY 2.2 software [81]. These probes were used in solution hybridization to capture our target sequences, followed by NGS on an Illumina Genome Analyzer IIx (Additional file 1: Methods).

Genomic sequences from both strands were captured using overlapping oligonucleotide sequence designs covering 342,075 nt among the 7 genes (Fig. 1). In total, 11,841 oligonucleotides were synthesized from the transcribed strand consisting of the complete, single copy coding, and flanking regions of ATM (3513), BRCA1 (1587), BRCA2 (2386), CDH1 (1867), CHEK2 (889), PALB2 (811), and TP53 (788). Additionally, 11,828 antisense strand oligos were synthesized (3497 ATM, 1591 BRCA1, 2395 BRCA2, 1860 CDH1, 883 CHEK2, 826 PALB2, and 776 TP53). Any intronic or intergenic regions without probe coverage are most likely due to the presence of conserved repetitive elements or other paralogous sequences.Fig. 1Capture Probe Coverage over Sequenced Genes. The genomic structure of the 7 genes chosen are displayed with the UCSC Genome Browser. Top row for each gene is a custom track with the “dense” visualization modality selected with black regions indicating the intervals covered by the oligonucleotide capture reagent. Regions without probe coverage contain conserved repetitive sequences or correspond to paralogous sequences that are unsuitable for probe design

Capture Probe Coverage over Sequenced Genes. The genomic structure of the 7 genes chosen are displayed with the UCSC Genome Browser. Top row for each gene is a custom track with the “dense” visualization modality selected with black regions indicating the intervals covered by the oligonucleotide capture reagent. Regions without probe coverage contain conserved repetitive sequences or correspond to paralogous sequences that are unsuitable for probe design

For regions lacking probe coverage (of ≥ 10 nt, N = 141; 8 in ATM, 26 in BRCA1, 10 in BRCA2, 29 in CDH1, 36 in CHEK2, 15 in PALB2, and 17 in TP53), probes were selected based on predicted Tms similar to other probes, limited alignment to other sequences in the transcriptome (<10 times), and avoidance of stable, base-paired 2° structures (with unaFOLD) [82, 83]. The average coverage of these sequenced regions was 14.1–24.9 % lower than other probe sets, indicating that capture was less efficient, though still successful.

### HBOC samples for oligo capture and high-throughput sequencing

Genomic DNA from 102 patients previously tested for inherited breast/ovarian cancer without evidence of a predisposing genetic mutation, was obtained from the Molecular Genetics Laboratory (MGL) at the London Health Sciences Centre in London, Ontario, Canada. Patients qualified for genetic susceptibility testing as determined by the Ontario Ministry of Health and Long-Term Care BRCA1 and BRCA2 genetic testing criteria [84] (see Additional file 2). The University of Western Ontario research ethics board (REB) approved this anonymized study of these individuals to evaluate the analytical methods presented here. BRCA1 and BRCA2 were previously analyzed by Protein Truncation Test (PTT) and Multiplex Ligation-dependent Probe Amplification (MLPA). The exons of several patients (N = 14) had also been Sanger sequenced. No pathogenic sequence change was found in any of these individuals. In addition, one patient with a known pathogenic BRCA variant was re-sequenced by NGS as a positive control.

### Sequence alignment and variant calling

Variant analysis involved the steps of detection, filtering, IT-based and coding sequence analysis, and prioritization (Fig. 2). Sequencing data were demultiplexed and aligned to the specific chromosomes of our sequenced genes (hg19) using both CASAVA (Consensus Assessment of Sequencing and Variation; v1.8.2) [85] and CRAC (Complex Reads Analysis and Classification; v1.3.0) [86] software. Alignments were prepared for variant calling using Picard [87] and variant calling was performed on both versions of the aligned sequences using the UnifiedGenotyper tool in the Genome Analysis Toolkit (GATK) [88]. We used the recommended minimum phred base quality score of 30, and results were exported in variant call format (VCF; v4.1). A software program was developed to exclude variants called outside of targeted capture regions and those with quality scores < 50. Variants flagged by bioinformatic analysis (described below) were also assessed by manually inspecting the reads in the region using the Integrative Genomics Viewer (IGV; version 2.3) [89, 90] to note and eliminate obvious false positives (i.e. variant called due to polyhomonucleotide run dephasing, or PCR duplicates that were not eliminated by Picard). Finally, common variants (≥1 % allele frequency based on dbSNP 142 or > 5 individuals in our study cohort) were not prioritized.Fig. 2Framework for the Identification of Potentially Pathogenic Variants. Integrated laboratory processing and bioinformatic analysis procedures for comprehensive complete gene variant determination and analysis. Intermediate datasets resulting from filtering are represented in yellow and final datasets in green. Non-bioinformatic steps, such as sample preparation are represented in blue and prediction programs in purple. Sequencing analysis yields base calls for all samples. CASAVA [85] and CRAC [86] were used to align these sequencing results to hg19. GATK [88] was used to call variants from this data against GRCh37 release of the reference human genome. Variants with a quality score < 50 and/or call confidence score < 30 were eliminated along with variants falling outside of our target regions. SNPnexus [112–114] was used to identify the genomic location of the variants. Nonsense and indels were noted and prediction tools were used to assess the potential pathogenicity of missense variants. The Shannon Pipeline [91] evaluated the effect of a variant on natural and cryptic SSs, as well as SRFBSs. ASSEDA [38] was used to predict the potential isoforms as a result of these variants. PWMs for 83 TFs were built using an information weight matrix generator based on Bipad [106]. Mutation Analyzer evaluated the effect of variants found 10 kb upstream up to the first intron on protein binding. Bit thresholds (R
i values) for filtering variants on software program outputs are indicated. Variants falling within the UTR sequences were assessed using SNPfold [20], and the most probable variants that alter mRNA structure (p < 0.1) were then processed using mFold to predict the effect on stability [83]. All UTR variants were scanned with a modified version of the Shannon Pipeline, which uses PWMs computed from nucleotide frequencies for 28 RBPs in RBPDB [109] and 76 RBPs in CISBP-RNA [110]. All variants meeting these filtering criteria were verified with IGV [89, 90]. *Sanger sequencing was only performed for protein truncating, splicing, and selected missense variants

Framework for the Identification of Potentially Pathogenic Variants. Integrated laboratory processing and bioinformatic analysis procedures for comprehensive complete gene variant determination and analysis. Intermediate datasets resulting from filtering are represented in yellow and final datasets in green. Non-bioinformatic steps, such as sample preparation are represented in blue and prediction programs in purple. Sequencing analysis yields base calls for all samples. CASAVA [85] and CRAC [86] were used to align these sequencing results to hg19. GATK [88] was used to call variants from this data against GRCh37 release of the reference human genome. Variants with a quality score < 50 and/or call confidence score < 30 were eliminated along with variants falling outside of our target regions. SNPnexus [112–114] was used to identify the genomic location of the variants. Nonsense and indels were noted and prediction tools were used to assess the potential pathogenicity of missense variants. The Shannon Pipeline [91] evaluated the effect of a variant on natural and cryptic SSs, as well as SRFBSs. ASSEDA [38] was used to predict the potential isoforms as a result of these variants. PWMs for 83 TFs were built using an information weight matrix generator based on Bipad [106]. Mutation Analyzer evaluated the effect of variants found 10 kb upstream up to the first intron on protein binding. Bit thresholds (R
i values) for filtering variants on software program outputs are indicated. Variants falling within the UTR sequences were assessed using SNPfold [20], and the most probable variants that alter mRNA structure (p < 0.1) were then processed using mFold to predict the effect on stability [83]. All UTR variants were scanned with a modified version of the Shannon Pipeline, which uses PWMs computed from nucleotide frequencies for 28 RBPs in RBPDB [109] and 76 RBPs in CISBP-RNA [110]. All variants meeting these filtering criteria were verified with IGV [89, 90]. *Sanger sequencing was only performed for protein truncating, splicing, and selected missense variants

### IT-based variant analysis

All variants were analyzed using the Shannon Human Splicing Mutation Pipeline, a genome-scale variant analysis program that predicts the effects of variants on mRNA splicing [91, 92]. Variants were flagged based on criteria reported in Shirley et al. (2013): weakened natural site ≥ 1.0 bits, or strengthened cryptic site (within 300 nt of the nearest exon) where cryptic site strength is equivalent or greater than the nearest natural site of the same phase [91]. The effects of flagged variants were further analyzed in detail using the Automated Splice Site and Exon Definition Analysis (ASSEDA) server [38].

Exonic variants and those found within 500 nt of an exon were assessed for their effects, if any, on SRFBSs [38]. Sequence logos for splicing regulatory factors (SRFs) (SRSF1, SRSF2, SRSF5, SRSF6, hnRNPH, hnRNPA1, ELAVL1, TIA1, and PTB) and their Rsequence values (the mean information content [93]) are provided in Caminsky et al. (2015) [36]. Because these motifs occur frequently in unspliced transcripts, only variants with large information changes were flagged, notably those with a) ≥ 4.0 bit decrease, i.e. at least a 16-fold reduction in binding site affinity, with Ri,initial ≥ Rsequence for the particular factor analyzed, or b) ≥ 4.0 bit increase in a site where Ri,final ≥ 0 bits. ASSEDA was used to calculate Ri,total, with the option selected to include the given SRF in the calculation. Variants decreasing Ri,total by < 3.0 bits (i.e. 8-fold) were predicted to potentially have benign effects on expression, and were not considered further.

Activation of pseudoexons through creating/strengthening of an intronic cryptic SS was also assessed [94]. Changes in intronic cryptic sites, where ΔRi > 1 bit and Ri,final ≥ (Rsequence – 1 standard deviation [S.D.] of Rsequence), were identified. A pseudoexon was predicted if a pre-existing cryptic site of opposite polarity (with Ri > [Rsequence - 1 S.D.]) and in the proper orientation for formation of exons between 10–250 nt in length was present. In addition, the minimum intronic distance between the pseudoexon and either adjacent natural exon was 100 nt. The acceptor site of the pseudoexon was also required to have a strong hnRNPA1 site located within 10 nt (Ri ≥ Rsequence) [38] to ensure accurate proofreading of the exon [37].

Next, variants affecting the strength of SRFs were analyzed by a contextual exon definition analysis of ΔRi,total. The context refers to the documented splicing activity of an SRF. For example, TIA1 has been shown to be an intronic enhancer of exon definition, so only intronic sites were considered. Similarly, hnRNPA1 proofreads the 3’ SS (acceptor) and inhibits exon recognition elsewhere [95]. Variants that lead to redundant SRFBS changes (i.e. one site is abolished and another proximate site [≤2 nt] of equivalent strength is activated) were assumed to have a neutral effect on splicing. If the strength of a site bound by PTB (polypyrimidine tract binding protein) was affected, its impact on binding by other factors was analyzed, as PTB impedes binding of other factors with overlapping recognition sites, but does not directly enhance or inhibit splicing itself [96].

To determine effects of variants on transcription factor (TF) binding, we first established which TFs bound to the sequenced regions of the gene promoters (and first exons) in this study by using ChIP-seq data from 125 cell types (Additional file 1: Methods) [97]. We identified 141 TFs with evidence for binding to the promoters of the genes we sequenced, including c-Myc, C/EBPβ, and Sp1, shown to transcriptionally regulate BRCA1, TP53, and ATM, respectively [98–100]. Furthermore, polymorphisms in TCF7L2, known to bind enhancer regions of a wide variety of genes in a tissue-specific manner [101], have been shown to increase risk of sporadic [102] and hereditary breast [103], as well as other types of cancer [104, 105].

IT-based models of the 141 TFs of interest were derived by entropy minimization of the DNase accessible ChIP-seq subsets [106]. Details are provided in Lu R, Mucaki E, and Rogan PK (BioRxiv; http://dx.doi.org/10.1101/042853). While some data sets would only yield noise or co-factor motifs (i.e. co-factors that bind via tethering, or histone modifying proteins [107]), techniques such as motif masking and increasing the number of Monte Carlo cycles yielded models for 83 TFs resembling each factor’s published motif. Additional file 3: Table S1 contains the final list of TFs and the models we built (described below) [108].

These TFBS models (N = 83) were used to scan all variants called in the promoter regions (10 kb upstream of transcriptional start site to the end of IVS1) of HBOC genes for changes in Ri. Binding site changes that weaken interactions with the corresponding TF (to Ri ≤ Rsequence) are likely to affect regulation of the adjacent target gene. Stringent criteria were used to prioritize the most likely variants and thus only changes to strong TFBSs (Ri,initial ≥ Rsequence), where reduction in strength was significant (ΔRi ≥ 4.0 bits), were considered. Alternatively, novel or strengthened TFBSs were also considered sources of dysregulated transcription. These sites were defined as having Ri,final ≥ Rsequence and as being the strongest predicted site in the corresponding genomic interval (i.e. exceeding the Ri values of adjacent sites unaltered by the variant). Variants were not prioritized if the TF was known to a) enhance transcription and IT analysis predicted stronger binding, or b) repress transcription and IT analysis predicted weaker binding.

Two complementary strategies were used to assess the possible impact of variants within UTRs. First, SNPfold software was used to assess the effect of a variant on 2° structure of the UTR (Additional file 1: Methods) [20]. Variants flagged by SNPfold with the highest probability of altering stable 2° structures in mRNA (where p-value < 0.1) were prioritized. To evaluate these predictions, oligonucleotides containing complete wild-type and variant UTR sequences (Additional file 4: Table S2) were transcribed in vitro and followed by SHAPE analysis, a method that can confirm structural changes in mRNA [44].

Second, the effects of variants on the strength of RBBSs were predicted. Frequency-based, position weight matrices (PWMs) for 156 RNA-binding proteins (RBPs) were obtained from the RNA-Binding Protein DataBase (RBPDB) [109] and the Catalog of Inferred Sequence Binding Preferences of RNA binding proteins (CISBP-RNA) [110, 111]. These were used to compute information weight matrices (based on the method described by Schneider et al. 1984; N = 147) (see Additional file 1: Methods) [40]. All UTR variants were assessed using a modified version of the Shannon Pipeline [91] containing the RBPDB and CISBP-RNA models. Results were filtered to include a) variants with |ΔRi| ≥ 4.0 bits, b) variants creating or strengthening sites (Ri,final ≥ Rsequence and the Ri,initial < Rsequence), and c) RBBSs not overlapping or occurring within 10 nt of a stronger, pre-existing site of another RBP.

### Exonic protein-altering variant analysis

The predicted effects of all coding variants were assessed with SNPnexus [112–114], an annotation tool that can be applied to known and novel variants using up-to-date dbSNP and UCSC human genome annotations. Variants predicted to cause premature protein truncation were given higher priority than those resulting in missense (or synonymous) coding changes. Missense variants were first cross referenced with dbSNP 142 [115]. Population frequencies from the Exome Variant Server [116] and 1000Genomes [117] are also provided. The predicted effects on protein conservation and function of the remaining variants were evaluated by in silico tools: PolyPhen-2 [118], Mutation Assessor (release 2) [119, 120], and PROVEAN (v1.1.3) [121, 122]. Default settings were applied and in the case of PROVEAN, the “PROVEAN Human Genome Variants Tool” was used, which includes SIFT predictions as a part of its output. Variants predicted by all four programs to be benign were less likely to have a deleterious impact on protein activity; however this did not exclude them from mRNA splicing analysis (described above in IT-Based Variant Analysis). All rare and novel variants were cross-referenced with general mutation databases (ClinVar [123, 124], Human Gene Mutation Database [HGMD] [125, 126], Leiden Open Variant Database [LOVD] [127–134], Domain Mapping of Disease Mutations [DM2] [135], Expert Protein Analysis System [ExPASy] [136] and UniProt [137, 138]), and gene-specific databases (BRCA1/2: the Breast Cancer Information Core database [BIC] [139] and Evidence-based Network for the Interpretation of Germline Mutant Alleles [ENIGMA] [140]; TP53: International Agency for Research on Cancer [IARC] [141]), as well as published reports to prioritize them for further workup.

### Variant classification

Flagged variants were prioritized if they were likely to encode a dysfunctional protein (indels, nonsense codon > 50 amino acids from the C-terminus, or abolition of a natural SS resulting in out-of-frame exon skipping) or if they exceeded established thresholds for fold changes in binding affinity based on IT (see Methods above). In several instances, our classification was superseded by previous functional or pedigree analyses (reported in published literature or databases) that categorized these variants as pathogenic or benign.

### Positive control

We identified the BRCA1 exon 17 nonsense variant c.5136G > A (chr17:41215907C > T; rs80357418; 2-5A) [142] in the sample that was provided as a positive control. This was the same mutation identified by the MGL as pathogenic for this patient. We also prioritized another variant in this patient (Table 1) [143].Table 1Prioritized variants in the positive controlGenemRNA ProteinrsID (dbSNP 142)CategoryConsequenceRef
BRCA1
c.5136G > Ars80357418Nonsense151 AA short[142]p.Trp1712Ter
BRCA2
c.3218A > Grs80358566SRFBSRepressor action of hnRNPA1 at this site abolished (5.2 to 0.4 bits). Blocking action of PTB removed as site is abolished (5.5 to -7.5 bits) and may uncover binding sites of other SRFs.p.Gln1073ArgMissenseListed in ClinVar as conflicting interpretations (likely benign, unknown) and in BIC as unknown clinical importance. 2 in silico programs called deleterious. The variant occurs between repeat motifs BRC1 and BRC2 of BRCA2, a region in which pathogenic missense mutations have not yet been identified.[143]SRFBSRepressor action of hnRNPA1 at this site abolished (5.2 to 0.4 bits). Blocking action of PTB removed as site is abolished (5.5 to -7.5 bits) and may uncover binding sites of other SRFs.

Prioritized variants in the positive control

### Variant validation

Protein-truncating, prioritized splicing, and selected prioritized missense variants were verified by Sanger sequencing. Primers of PCR amplicons are indicated in Additional file 5: Table S3.

### Deletion analysis

Potential large rearrangements were detected with BreakDancer software [144], which identifies novel genomic rearrangements based on the respective orientation and distance between ends of the same read (and exceeding the lengths of NGS library inserts). This approach can, in theory, approximately localize deletions, duplications, or other types of breakpoints within exons, introns, and regulatory regions (eg. promoters) that could affect gene expression and function. We required at least 4 reads per suspected rearrangement in a sample separated by >700 nt, with each end mapping to proximate genomic reference coordinates to infer a potential deletion. Synthetic and cost limitations in the maximum genomic real estate covered by the capture reagent led to a tradeoff between extending the span of captured genomic intervals and higher tiling densities over shorter sequences, ie. exons, to achieve the level of coverage to reliably detect deletions based on read counts alone.

Our complete gene enrichment strategy with independent capture of both genomic strands enabled and facilitated development of a new algorithm to identify potential hemizygous genomic intervals in these individuals. In each subject, we first searched for contiguous long stretches (usually > > 1 kb) of non-polymorphic segments with diminished repetitive element content (<10 %), which is consistent with the possibility of these regions harboring a deletion. Then, we determined the likelihood of homo- or hemizygosity by comparing the degree of heterozygosity of variants in each of these intervals in for an individual with all of the other individuals sequenced with this protocol in this population. Regions containing haplotype blocks in strong linkage disequilibrium (LD; from HapMap [145]) were then excluded as candidate deletion intervals. Some individuals without a deletion are expected to be non-polymorphic, because detection of heterozygosity depends on genomic length of the region, marker informativeness, and the level of LD for those markers. We required that > 80 % of the control individuals be heterozyogous for at least two well-distributed loci within these intervals. Highly informative SNPs with a random genomic distribution in the controls (and other public databases) and which were non-polymorphic in the individual with the suspected deletion were weighted more heavily in inferring potential hemizygosity. This analysis was implemented using a Perl script that identified the most likely intervals of hemizygosity, which were then crossreferenced with the corresponding genomic intervals in HapMap.

### Junctional read detection

Potential large rearrangements were detected with BreakDancer software [144], which identifies novel genomic rearrangements based on the respective orientation and distance between ends of the same read (and exceeding the lengths of NGS library inserts). This approach can, in theory, approximately localize deletions, duplications, or other types of breakpoints within exons, introns, and regulatory regions (eg. promoters) that could affect gene expression and function. We required at least 4 reads per suspected rearrangement in a sample separated by >700 nt, with each end mapping to proximate genomic reference coordinates to infer a potential deletion. Synthetic and cost limitations in the maximum genomic real estate covered by the capture reagent led to a tradeoff between extending the span of captured genomic intervals and higher tiling densities over shorter sequences, ie. exons, to achieve the level of coverage to reliably detect deletions based on read counts alone.

### Prioritization based on potential hemizygosity

Our complete gene enrichment strategy with independent capture of both genomic strands enabled and facilitated development of a new algorithm to identify potential hemizygous genomic intervals in these individuals. In each subject, we first searched for contiguous long stretches (usually > > 1 kb) of non-polymorphic segments with diminished repetitive element content (<10 %), which is consistent with the possibility of these regions harboring a deletion. Then, we determined the likelihood of homo- or hemizygosity by comparing the degree of heterozygosity of variants in each of these intervals in for an individual with all of the other individuals sequenced with this protocol in this population. Regions containing haplotype blocks in strong linkage disequilibrium (LD; from HapMap [145]) were then excluded as candidate deletion intervals. Some individuals without a deletion are expected to be non-polymorphic, because detection of heterozygosity depends on genomic length of the region, marker informativeness, and the level of LD for those markers. We required that > 80 % of the control individuals be heterozyogous for at least two well-distributed loci within these intervals. Highly informative SNPs with a random genomic distribution in the controls (and other public databases) and which were non-polymorphic in the individual with the suspected deletion were weighted more heavily in inferring potential hemizygosity. This analysis was implemented using a Perl script that identified the most likely intervals of hemizygosity, which were then crossreferenced with the corresponding genomic intervals in HapMap.

### Results

The average coverage of capture region per individual was 90.8x (range of 53.8 to 118.2x between 32 samples) with 98.8 % of the probe-covered nucleotides having ≥ 10 reads. Samples with fewer than 10 reads per nucleotide were re-sequenced and the results of both runs were combined. The combined coverage of these samples was, on average, 48.2x (±36.2).

The consistency of both library preparation and capture protocols was improved from initial runs, which significantly impacted sequence coverage (Additional file 1: Methods). Of the 102 patients tested, 14 had been previously Sanger sequenced for BRCA1 and BRCA2 exons. Confirmation of previously discovered SNVs served to assess the methodological improvements introduced during NGS and ultimately, to increase confidence in variant calling. Initially, only 15 of 49 SNVs in 3 samples were detected. The detection rate of SNVs was improved to 100 % as the protocol progressed. All known SNVs (N = 157) were called in subsequent sequencing runs where purification steps were replaced with solid phase reversible immobilization beads and where RNA bait was transcribed the same day as capture. To minimize false positive variant calls, sequence read data were aligned with CASAVA and CRAC, variants were called for each alignment with GATK, and discrepancies were then resolved by manual review.

GATK called 14,164 unique SNVs and 1147 indels. Only 3777 (15.3 %) SNVs were present in both CASAVA and CRAC-alignments for at least one patient, and even fewer indel calls were concordant between both methods (N = 110; 6.2 %). For all other SNVs and indels, CASAVA called 6871 and 1566, respectively, whereas CRAC called 13,958 and 110, respectively. Some variants were counted more than once if they were called by different alignment programs in two or more patients. Intronic and intergenic variants proximate to low complexity sequences tend to generate false positive variants due to ambiguous alignment, a well known technical issue in short read sequence analysis [146, 147], contributing to this discrepancy. For example, CRAC correctly called a 19 nt deletion of BRCA1 (rs80359876; also confirmed by Sanger sequencing) but CASAVA flagged the deleted segment as a series of false-positives (Additional file 6: Figure S1). For these reasons, all variants were manually reviewed.

The Shannon Pipeline reported 99 unique variants in natural donor or acceptor SSs. After technical and frequency filtering criteria were applied, 12 variants remained (Additional file 7: Table S4). IT analysis allowed for the prioritization of 3 variants, summarized in Table 2.Table 2Variants prioritized by IT analysisPatient IDGenemRNArsID (dbSNP 142)Information ChangeConsequencef or Binding Factor Affected
R
i,initial

R
i,final
ΔR
i or R
i
e
Allele Frequency (%)d
(bits)(bits)(bits)Abolished Natural SS7-4 F
ATM
c.3747-1G > Aa
Novel11.00.1−10.9Exon skipping and use of alternative splice forms4-1 F
ATM
c.6347 + 1G > Tb
Novel10.4−8.3−18.6Exon skippingLeaky Natural SS4-2B
CHEK2
c.320-5 T > Aa
rs1219087006.84.1−2.7Leaky splicing with intron inclusion0.08Activated Cryptic SS7-3E
BRCA1
c.548-293G > Ars117281398−12.12.614.7Cryptic site not expected to be used. Total information for natural exon is stronger than cryptic exon.0.747-4A
BRCA2
c.7618-269_7618-260del10Novel3.99.45.5Cryptic site not expected to be used. Total information for natural exon is stronger than cryptic exon.Pseudoexon formation due to activated acceptor SS7-3 F
BRCA2
c.8332-805G > ANovel−9.35.45.6e
6065/211/592f
7-3D
CDH1
c.164-2023A > Grs184740925−6.64.36.5e
61,236/224/1798f
0.35-3H
CDH1
c.2296-174 T > Ars5654888667.38.55.0e
1175/50/124f
0.02Pseudoexon formation due to activated donor SS3-6A
BRCA1
c.212 + 253G > Ars1893521914.16.75.2e
186/63/1250f
0.085-2G
BRCA2
c.7007 + 2691G > Ars3678905774.77.27.7e
2589/103/5272f
0.02Affected TFBSs7-4B
BRCA1
c.-8895G > ANovel10.9−0.2−11.1GATA-3 (GATA3)5-3E
CDH1
c.-54G > Crs50308741.712.010.4E2F-4 (E2F4)7-4E0.165-2B
PALB2
c.-291C > Grs55282422712.1−1.3−13.4GABPα (GABPA)0.17-2 F
TP53
c.-28-3132 T > Crs17882863−6.310.917.2RUNX3 (RUNX3)0.34-1A
TP53
c.-28-1102 T > Crs1134516735.112.37.2E2F-4 (E2F4)0.48.012.94.8Sp1 (SP1)Affected RBBSs7-4G
ATM
c.-244 T > Ars5399482189.8−19.9−29.7RBFOXc.-744 T > A0.04c.-1929 T > Ac.-3515 T > A5-3C
CDH1
c.*424 T > ANovel−20.39.629.9SF3B48.21.8−6.4CELF47-2E
CHEK2
c.-588G > Ars14156834210.93.7−7.2BX511012.14-3C.5-4G
CHEK2
c.-345C > Tc
rs1378530073.311.48.2SF3B43-1A
TP53
c.-107 T > Crs11353009010.54.5−6.0ELAVL14-1Hc.-188 T > C0.724-2H
TP53
c.*1175A > Crs7837822210.74.1−6.6KHDRBS17-2 Fc.*1376A > C0.26c.*1464A > C
aConfirmed by Sanger sequencing
bAmbiguous Sanger sequencing results
cPrioritized under missense change and was therefore verified with Sanger sequencing. Variant was confirmed
dIf available
e
R
i of site of opposite polarity in the pseudoexon
fConsequences for pseudoexon formation describe how the intron is divided: “new intron A length/pseudoexon length/new exon B lengthNone of the variants have been previously reported by other groups with the exception of CHEK2 c.320-5T>A [148]

Variants prioritized by IT analysis

aConfirmed by Sanger sequencing

bAmbiguous Sanger sequencing results

cPrioritized under missense change and was therefore verified with Sanger sequencing. Variant was confirmed

dIf available

e
R
i of site of opposite polarity in the pseudoexon

fConsequences for pseudoexon formation describe how the intron is divided: “new intron A length/pseudoexon length/new exon B length

None of the variants have been previously reported by other groups with the exception of CHEK2 c.320-5T>A [148]

First, the novel ATM variant c.3747-1G > A (chr11:108,154,953G > A; sample number 7-4 F) abolishes the natural acceptor of exon 26 (11.0 to 0.1 bits). ASSEDA reports the presence of a 5.3 bit cryptic acceptor site 13 nt downstream of the natural site, but the effect of the variant on a pre-existing cryptic site is negligible (~0.1 bits). The cryptic exon would lead to exon deletion and frameshift (Fig. 3a). ASSEDA also predicts skipping of the 246 nt exon, as the Ri,final of the natural acceptor is now below Ri,minimum (1.6 bits), altering the reading frame. Second, the novel ATM c.6347 + 1G > T (chr11:108188249G > T; 4-1 F) abolishes the 10.4 bit natural donor site of exon 44 (ΔRi = -18.6 bits), and is predicted to cause exon skipping. Finally, the previously reported CHEK2 variant, c.320-5A > T (chr22:29,121,360 T > A; rs121908700; 4-2B) [148] weakens the natural acceptor of exon 3 (6.8 to 4.1 bits), and may activate a cryptic acceptor (7.4 bits) 92 nt upstream of the natural acceptor site which would shift the reading frame (Fig. 4). A constitutive, frameshifted alternative isoform of CHEK2 lacking exons 3 and 4 has been reported, but skipping of exon 3 alone is not normally observed.Fig. 3Predicted Isoforms and Relative Abundances as a Consequence of ATM splice variant c.3747-1G > A. Intronic ATM variant c.3747-1G > A abolishes (11.0 to 0.1 bits) the natural acceptor of exon 26 (total of 63 exons). a ASSEDA predicts skipping of the natural exon (R
i,total from 14.5 to 3.6 bits [a 1910 fold decrease in exon strength]; isoform 7) and/or activation of a pre-existing cryptic acceptor site 13 nt downstream (R
i,total for cryptic exon = 9.0 bits; isoform 1) of the natural site leading to exon truncation. The reading frame is altered in both mutant isoforms. The other isoforms use weak, alternate acceptor/donor sites leading to cryptic exons with much lower total information. b Before the mutation, isoform 7 is expected to be the most abundant splice form. c After the mutation, isoform 1 is predicted to become the most abundant splice form and the wild-type isoform is not expected to be expressedFig. 4Predicted Isoforms and Relative Abundances as a Consequence of CHEK2 splice variant c.320-5 T > A. Intronic CHEK2 variant c.320-5 T > A weakens (6.8 to 4.1 bits) the natural acceptor of exon 3 (total of 15 exons). a ASSEDA reports the weakening of the natural exon strength (R
i,total reduced from 13.2 to 10.5 bits), which would result in reduced splicing of the exon otherwise known as leaky splicing. A pre-existing cryptic acceptor exists 92 nt upstream of the natural site, leading to a cryptic exon with similar strength to the mutated exon (R
i,total = 10.0 bits). This cryptic exon would contain 92 nt of the intron. b Before the mutation, isoform 1 is expected to be the only isoform expressed. c After the mutation, isoform 1 (wild-type) is predicted to become relatively less abundant and isoform 2 is expected to be expressed, although less abundant in relation to isoform 1

Predicted Isoforms and Relative Abundances as a Consequence of ATM splice variant c.3747-1G > A. Intronic ATM variant c.3747-1G > A abolishes (11.0 to 0.1 bits) the natural acceptor of exon 26 (total of 63 exons). a ASSEDA predicts skipping of the natural exon (R
i,total from 14.5 to 3.6 bits [a 1910 fold decrease in exon strength]; isoform 7) and/or activation of a pre-existing cryptic acceptor site 13 nt downstream (R
i,total for cryptic exon = 9.0 bits; isoform 1) of the natural site leading to exon truncation. The reading frame is altered in both mutant isoforms. The other isoforms use weak, alternate acceptor/donor sites leading to cryptic exons with much lower total information. b Before the mutation, isoform 7 is expected to be the most abundant splice form. c After the mutation, isoform 1 is predicted to become the most abundant splice form and the wild-type isoform is not expected to be expressed

Predicted Isoforms and Relative Abundances as a Consequence of CHEK2 splice variant c.320-5 T > A. Intronic CHEK2 variant c.320-5 T > A weakens (6.8 to 4.1 bits) the natural acceptor of exon 3 (total of 15 exons). a ASSEDA reports the weakening of the natural exon strength (R
i,total reduced from 13.2 to 10.5 bits), which would result in reduced splicing of the exon otherwise known as leaky splicing. A pre-existing cryptic acceptor exists 92 nt upstream of the natural site, leading to a cryptic exon with similar strength to the mutated exon (R
i,total = 10.0 bits). This cryptic exon would contain 92 nt of the intron. b Before the mutation, isoform 1 is expected to be the only isoform expressed. c After the mutation, isoform 1 (wild-type) is predicted to become relatively less abundant and isoform 2 is expected to be expressed, although less abundant in relation to isoform 1

Variants either strengthening (N = 4) or slightly weakening (ΔRi < 1.0 bits; N = 4) a natural site were not prioritized. In addition, we rejected the ATM variant (c.1066-6 T > G; chr11:108,119,654 T > G; 4-1E and 7-2B), which slightly weakens the natural acceptor of exon 9 (11.0 to 8.1 bits). Although other studies have shown leaky expression as a result of this variant [149], a more recent meta-analysis concluded that this variant is not associated with increased breast cancer risk [150].

Two variants produced information changes that could potentially impact cryptic splicing, but were not prioritized for the following reasons (Table 2). The first variant, novel BRCA2 deletion c.7618-269_7618-260del10 (chr13:32931610_32931619del10; 7-4A) strengthens a cryptic acceptor site 245 nt upstream from the natural acceptor of exon 16 (Ri,final = 9.4 bits, ΔRi = 5.5 bits). Being 5.7-fold stronger than the natural site (6.9 bits), two potential cryptic isoforms were predicted, however the exon strengths of both are weaker than the unaffected natural exon (Ri,total = 6.6 bits) and thus neither were prioritized. The larger gap surprisal penalties explain the differences in exon strength. The natural donor SS may still be used in conjunction with the abovementioned cryptic SS, resulting in an exon with Ri,total = 3.5 bits. Alternatively, the cryptic site and a weak donor site 180 nt upstream of the natural donor (Ri = 0.7 vs 1.4, cryptic and natural donors, respectively) result in an exon with Ri,total = 6.5 bits. The second variant, BRCA1 c.548-293G > A (chr17:41249599C > T; 7-3E), creates a weak cryptic acceptor (Ri,final = 2.6 bits, ΔRi = 6.2 bits) 291 nt upstream of the natural acceptor for exon 8 (Ri = 0.5). Although the cryptic exon is strengthened (final Ri,total = 6.9 bits, ΔRi = 14.7 bits), ASSEDA predicts the level of expression of this exon to be negligible, as it is weaker than the natural exon (Ri,total = 8.4 bits) due to the increased length of the predicted exon (+291 nt) [38].

The Shannon Pipeline initially reported 1583 unique variants creating or strengthening intronic cryptic sites. We prioritized 5 variants, 1 of which is novel (BRCA2 c.8332-805G > A; 7-3 F), that were within 250 nt of a pre-existing complementary cryptic site and have an hnRNPA1 site within 5 nt of the acceptor (Table 2). If used, 3 of these pseudoexons would lead to a frameshifted transcript.

Variants within 500 nt of an exon junction and all exonic variants (N = 4015) were investigated for their potential effects on affinity of sites to corresponding SRFs [38]. IT analysis flagged 54 variants significantly altering the strength of at least one binding site (Additional file 8: Table S5). A careful review of the variants, the factor affected, and the position of the binding site relative to the natural SS, prioritized 36 variants (21 novel), of which 4 are in exons and 32 are in introns. As an example, a novel CHEK2 exon 2 variant c.69C > A (p.Gly23=) is predicted to increase the strength of an hnRNP A1 site (0.7 to 5.3 bits) and decrease total exon strength (ΔRi,total = -5.7 bits). A similar type of exonic variant in FANCM, which was predicted to create an exonic hnRNP A1 site by IT, has been shown to bind this exonic repressor and induce exon skipping [37].

We assessed SNVs with models of 83 TFs experimentally shown to bind (Additional file 3: Table S1) upstream or within the first exon and intron of our sequenced genes (N = 2177). Thirteen variants expected to significantly affect TF binding were flagged (Additional file 9: Table S6). The final filtering step considered the known function of the TF in transcription, resulting in 5 variants (Table 2) in 6 patients (one variant was identified in two patients). Four of these variants have been previously reported (rs5030874, rs552824227, rs17882863, rs113451673) and one is novel (c.-8895G > A; 7-4B).

There were 364 unique UTR variants found by sequencing. These variants were evaluated for their effects on mRNA 2° structure (including that of splice forms with alternate UTRs in the cases of BRCA1 and TP53) through SNPfold, resulting in 5 flagged variants (Table 3), all of which have been previously reported.Table 3Variants predicted by SNPfold to affect UTR structureClassa
Patient IDGenemRNAUTR positionrsID (dbSNP 142)Ranke

p-valueAllele Frequency (%)d
FIn 26 patients
BRCA2
b
c.-52A > G5’ UTRrs2061182/9000.00214.86FIn 40 patients
BRCA2
b
c.*532A > G3’ UTRrs11571836239/27000.08919.75P7-4C
CDH1
c
c.-71C > G5’ UTRrs3403377169/6000.1150.56F4-2E
TP53
b
c.*485G > A3’ UTRrs4968187169/45000.0385-4A5.11F2-1A, 7-1B, 5-2A.7-1D, 7-2B, 7-2F
TP53
b
c.*826G > A3’ UTRrs17884306371/45000.0827-4C5.71
aF:Flagged; P:Prioritized
bLong Range UTR SNPfold Analysis
cLocal Range SNPfold Analysis
dIf available
eRank of the SNP, in terms of how much it changes the mRNA structure compared to all other possible mutations

Variants predicted by SNPfold to affect UTR structure

aF:Flagged; P:Prioritized

bLong Range UTR SNPfold Analysis

cLocal Range SNPfold Analysis

dIf available

eRank of the SNP, in terms of how much it changes the mRNA structure compared to all other possible mutations

Analysis of three variants using mFOLD [83] revealed likely changes to the UTR structure (Fig. 5). Two variants with possible 2° structure effects were common (BRCA2 c.-52A > G [N = 26 samples] and c.*532A > G [N = 40]) and not prioritized. The 5’ UTR CDH1 variant c.-71C > G (chr16:68771248C > G; rs34033771; 7-4C) disrupts a double-stranded hairpin region to create a larger loop structure, thus increasing binding accessibility (Fig. 5a and b). Analysis using RBPDB and CISBP-RNA-derived IT models suggests this variant affects binding by NCL (Nucleolin, a transcription coactivator) by decreasing binding affinity 14-fold (Ri,initial = 6.6 bits, ΔRi = -3.8 bits) (Additional file 10: Table S7). This RBP has been shown to bind to the 5’ and 3’ UTR of p53 mRNA and plays a role in repressing its translation [151].Fig. 5Predicted Alteration in UTR Structure Using mFOLD for Variants Flagged by SNPfold. Wild-type and variant structures are displayed, with the variant indicated by a red arrow. a Predicted wild-type structure of CDH1 5’UTR surrounding c.-71. b Predicted CDH1 5’UTR structure due to c.-71C > G variant. c Predicted wild-type TP53 3’UTR structure surrounding c.*485. d Predicted TP53 5’UTR structure due to c.*485G > A variant. e Predicted wild-type TP53 3’UTR structure surrounding c.*826. f Predicted TP53 5’UTR structure due to c.*826G > A variant. §SHAPE analysis revealed differences in reactivity between mutant and variant mRNAs, confirming alterations to 2° structure

Predicted Alteration in UTR Structure Using mFOLD for Variants Flagged by SNPfold. Wild-type and variant structures are displayed, with the variant indicated by a red arrow. a Predicted wild-type structure of CDH1 5’UTR surrounding c.-71. b Predicted CDH1 5’UTR structure due to c.-71C > G variant. c Predicted wild-type TP53 3’UTR structure surrounding c.*485. d Predicted TP53 5’UTR structure due to c.*485G > A variant. e Predicted wild-type TP53 3’UTR structure surrounding c.*826. f Predicted TP53 5’UTR structure due to c.*826G > A variant. §SHAPE analysis revealed differences in reactivity between mutant and variant mRNAs, confirming alterations to 2° structure

In addition, the TP53 variant c.*485G > A (NM_000546.5: chr17:7572442C > T; rs4968187) is found at the 3’ UTR and was identified in two patients (4-2E and 5-4A). In silico mRNA folding analysis demonstrated this variant disrupts a G/C bond of a loop in the highest ranked potential mRNA structure (Fig. 5c and d). Also, SHAPE analysis showed a difference in 2° structure between the wild-type and mutant (data not shown). IT analysis with RBBS models indicated that this variant significantly increases the binding affinity of SF3B4 by > 48-fold (Ri,final = 11.0 bits, ΔRi = 5.6 bits) (Additional file 10: Table S7). This RBP is one of four subunits comprising the splice factor 3B, which binds upstream of the branch-point sequence in pre-mRNA [152].

The third flagged variant also occurs in the 3’ UTR of TP53 (c.*826G > A; chr17:7572,101C > T; rs17884306), and was identified in 6 patients (2-1A, 7-1B, 5-2A.7-1D, 7-2B, 7-2 F, and 7-4C). It disrupts a potential loop structure, stabilizing a double-stranded hairpin, and possibly making it less accessible (Fig. 5e and f). Analysis using RBPDB-derived models suggests this variant could affect the binding of both RBFOX2 and SF3B4 (Additional file 10: Table S7). A binding site for RBFOX2, which acts as a promoter of alternative splicing by favoring the inclusion of alternative exons [153], is created (Ri,final = 9.8 bits; ΔRi = -6.5 bits). This variant is also expected to simultaneously abolish a SF3B4 binding site (Ri,final = -20.3 bits; ΔRi = -29.9 bits).

RBPDB- and CISBP-RNA-derived information model analysis of all UTR variants resulted in the prioritization of 1 novel, and 5 previously-reported variants (Table 2). No patient within the cohort exhibited more than one prioritized RBBS variant.

To evaluate the background rate of prioritizing variants flagged by this method, all 5’ and 3’ UTR SNVs in dbSNP144 for the 7 genes sequenced (excluding those already flagged in Table 3) were evaluated by SNPfold and our RBP information models. Of 1207 SNVs, only 10 were prioritized with both methods, which results in a background rate of 0.83 %.

Exonic variants called by GATK (N = 245) included insertions, deletions, nonsense, missense, and synonymous changes.

We identified 3 patients with different indels (Table 4). One was a PALB2 insertion c.1617_1618insTT (chr16:23646249_23646250insAA; 5-3A) in exon 4, previously reported in ClinVar as pathogenic. This mutation results in a frameshift and premature translation termination by 626 residues, abolishing domain interactions with RAD51, BRCA2, and POLH [137]. We also identified two known frameshift mutations in BRCA1: c.4964_4982del19 in exon 15 (chr17:41222949_41222967del19; rs80359876; 5-1B) and c.5266_5267insC in exon 19 (chr17:41209079_41209080insG; rs397507247; 5-3C) [148, 154]. Both are indicated as pathogenic and common in the BIC Database due to the loss of one or both C-terminal BRCT repeat domains [137]. Truncation of these domains produces instability and impairs nuclear transcript localization [155], and this bipartite domain is responsible for binding phosphoproteins that are phosphorylated in response to DNA damage [156, 157].Table 4Variants resulting in premature protein truncationPatient IDGeneExonmRNA ProteinrsID (dbSNP 142)ClinVard,e,f
DetailsRefAllele Frequency (%)c
Insertions/Deletions5-1B
BRCA1
15 of 23c.4964_4982del19a
rs803598766d; Pathogenic/likely pathogenice; Familial breast and breast-ovarian cancer, Hereditary cancer-predisposing syndromef.STOP at p.1670-p.Ser1655Tyrfs193 AA short5-3C
BRCA1
19 of 23c.5266_5267insCa
rs39750724713d; Pathogenic, risk factore; Familial breast, breast-ovarian, and pancreatic cancer, Hereditary cancer-predisposing syndromef.STOP at p.1788[148, 154]p.Gln1756Profs75 AA short5-3A
PALB2
4 of 13c.1617_1618insTTa
-1d; Pathogenice; Hereditary cancer-predisposing syndromef.STOP at p.561-p.Asn540Leufs626 AA shortStop Codons7-1G
BRCA2
15 of 27c.7558C > Tb
rs803589815d; Pathogenice; Familial breast, and breast-ovarian cancer, Hereditary cancer-predisposing syndromef.899 AA short[158]p.Arg2520Ter4-4A
BRCA2
25 of 27c.9294C > Ga
rs803592003d; Pathogenice; Familial breast and breast-ovarian cancerf.321 AA short[159]p.Tyr3098Ter7-3A
PALB2
4 of 13c.1240C > Ta
rs1801771003d; Pathogenice; Familial breast cancer, Hereditary cancer-predisposing syndromef.773 AA short[58]p.Arg414Ter4-4D
PALB2
4 of 13c.1042C > Ta
Novel-839 AA short-p.Gln348Ter
aConfirmed by Sanger sequencing
bNot confirmed by Sanger sequencing
cIf available
dNumber of submissions
eClinical significance
fCondition(s)

Variants resulting in premature protein truncation

aConfirmed by Sanger sequencing

bNot confirmed by Sanger sequencing

cIf available

dNumber of submissions

eClinical significance

fCondition(s)

We also identified 4 nonsense mutations, one of which was novel in exon 4 of PALB2 (c.1042C > T; chr16:23646825G > A; 4-4D). Another in PALB2 has been previously reported (c.1240C > T; chr16:23646627G > A; rs180177100; 7-3A) [58]. As a consequence, functional domains of PALB2 that interact with BRCA1, RAD51, BRCA2, and POLH are lost [137]. Two known nonsense mutations were found in BRCA2, c.7558C > T in exon 15 [158] and c.9294C > G in exon 25 [159]. The first (chr13:32930687C > T; rs80358981; 7-1G) causes the loss of the BRCA2 region that binds FANCD2, responsible for loading BRCA2 onto damaged chromatin [160]. The second (chr13:32968863C > G, rs80359200; 4-4A) does not occur within a known functional domain, however the transcript is likely to be degraded by nonsense mediated decay [161].

GATK called 61 missense variants, of which 18 were identified in 6 patients or more and 19 had allele frequencies > 1.0 % (Additional file 11: Table S8). The 40 remaining variants (15 ATM, 8 BRCA1, 9 BRCA2, 2 CDH1, 2 CHEK2, 3 PALB2, and 1 TP53) were assessed using a combination of gene specific databases, published classifications, and 4 in silico tools (Additional file 12: Table S9). We prioritized 27 variants, 2 of which were novel. None of the non-prioritized variants were predicted to be damaging by more than 2 of 4 conservation-based software programs.

Initially, 15,311 unique variants were identified by complete gene sequencing of 7 HBOC genes. Of these, 132 were flagged after filtering, and further reduced by IT-based variant analysis and consultation of the published literature to 87 prioritized variants. Figure 6 illustrates the decrease in the number of unique variants per patient at each step of our identification and prioritization process. The distribution of prioritized variants by gene is 34 in ATM, 13 in BRCA1, 11 in BRCA2, 8 in CDH1, 6 in CHEK2, 10 in PALB2, and 5 in TP53 (Additional file 13: Table S10), which are categorized by type in Table 5.Fig. 6Ladder Plot Representing Variant Identification and Prioritization. Each line is representative of a different sample in each sequencing run (a-e), illustrating the number of unique variants at important steps throughout the variant prioritization process. The left-most point indicates the total number of unique variants. The second point represents the number of unique variants remaining after common (>5 patients within cohort and/or ≥ 1.0 % allele frequency) and false-positive variants were removed. The right-most point represents the final number of unique. No variants were prioritized in the following patients: 2-1A, 2-5A, 2-6A, 3-2A, 3-3A, 3-4A, 3-5A, 3-8A, 4-1B, 4-2C, 4-2 F, 4-3B, 4-3D, 4-4B, 4-4E, 5-1G, 5-1H, 5-3D, 5-4C, 5-4D, 5-4 F, 5-4G, 5-4H, 7-1B, 7-1C, 7-1D, 7-1H, 7-2B, 7-2C, 7-2H, 7-3H, 7-4A, 7-4D, 7-4H. The average number of variants per patient at each step is indicated in a table below each plot, along with the percent reduction in variants from one step to anotherTable 5Summary of prioritized variants by geneIndelNonsenseMissenseNatural SplicingCryptic SplicingPseudoexonSR FactorTFUTR StructureUTR BindingTotal
ATM
00142001800134 a

BRCA1
202001710013
BRCA2
023002400011
CDH1
00200211118
CHEK2
00210030026 a

PALB2
123000310010
TP53
00100002025Total342730536516Three variants were prioritized under multiple categories: ATM chr11:108121730A > G (missense and SRFBS), CHEK2 chr22:29121242G > A (missense, UTR binding), and CHEK2 chr22:29130520C > T (missense, UTR binding)
a Counts represent the number of unique variants identified (i.e. a variant is not counted twice if it appeared in multiple individuals)

Ladder Plot Representing Variant Identification and Prioritization. Each line is representative of a different sample in each sequencing run (a-e), illustrating the number of unique variants at important steps throughout the variant prioritization process. The left-most point indicates the total number of unique variants. The second point represents the number of unique variants remaining after common (>5 patients within cohort and/or ≥ 1.0 % allele frequency) and false-positive variants were removed. The right-most point represents the final number of unique. No variants were prioritized in the following patients: 2-1A, 2-5A, 2-6A, 3-2A, 3-3A, 3-4A, 3-5A, 3-8A, 4-1B, 4-2C, 4-2 F, 4-3B, 4-3D, 4-4B, 4-4E, 5-1G, 5-1H, 5-3D, 5-4C, 5-4D, 5-4 F, 5-4G, 5-4H, 7-1B, 7-1C, 7-1D, 7-1H, 7-2B, 7-2C, 7-2H, 7-3H, 7-4A, 7-4D, 7-4H. The average number of variants per patient at each step is indicated in a table below each plot, along with the percent reduction in variants from one step to another

Summary of prioritized variants by gene

Three variants were prioritized under multiple categories: ATM chr11:108121730A > G (missense and SRFBS), CHEK2 chr22:29121242G > A (missense, UTR binding), and CHEK2 chr22:29130520C > T (missense, UTR binding)

a Counts represent the number of unique variants identified (i.e. a variant is not counted twice if it appeared in multiple individuals)

Three prioritized variants have multiple predicted roles: ATM c.1538A > G in missense and SRFBS, CHEK2 c.190G > A in missense and UTR binding, and CHEK2 c.433C > T in missense and UTR binding. Of the 102 patients that were sequenced, 72 (70.6 %) exhibited at least one prioritized variant, and some patients harbored more than one prioritized variant (N = 33; 32 %). Additional file 14: Table S11 presents a summary of all flagged and prioritized variants for patients with at least one prioritized variant.

Using BreakDancer, none of the individuals analyzed exhibited large rearrangements that met the level of stringency required, but a small intragenic rearrangement in BRCA1 was identified and confirmed by Sanger sequencing. Attempts to detect deletions with BreakDancer only flagged single, non-contiguous paired-end reads, rather than a series of reads clustered within the same region within the same individual, which would be necessary to indicate the presence of a true deletion or structural rearrangement.

After prioritizing individuals for potential hemizygosity in the sequenced regions, potential deletions were detected in BRCA2 and CDH1. Patient UWO5-4D exhibited a non-polymorphic 32.1 kb interval in BRCA2, spanning introns 1 to 13, that was absent from all of the other individuals (chr13:32890227-32922331). Haploview (hapmap.org) showed very low levels of LD in this region. The potential deletion may extend further downstream, however the presence of a haploblock, covering the entire sequenced interval beyond exon 11, with significant LD precludes delineation of the telomeric breakpoint. We also flagged a non-polymorphic 2.6 kb interval near the 3’ end of CDH1 in 6 individuals (UWO3-5, UWO4-2C, UWO4-4E, UWO4-4 F, UWO4-2G, UWO5-2H). This is a low LD region spanning chr16:68861286-68863887 that includes exons 14 and 15, and is polymorphic in all of the other individuals sequenced. CDH1 mutations are characteristically present in families with predisposition to gastric cancer, however breast cancer frequently co-occurs [69]. A study of CDH1 deletions in inherited gastric cancer identified two families with deletions that overlap the intervals prioritized in the present study [162].

The analysis and prioritization of non-coding variants can also be accomplished using Combined Annotation Dependent Depletion (CADD; [163]), which uses known and simulated variants to compute a C-score, an ad hoc measure of how deleterious is likely to be. The suggested C-score cutoff is between 10 and 20, though it is stated that any selected cutoff value would be arbitrary (http://cadd.gs.washington.edu/info). This contrasts with information-based methods, which are based on thermodynamically-defined thresholds. To directly compare methods, CADD scores were obtained for all prioritized or flagged SNVs. Half of prioritized variants met this cutoff (C > 10), while only 28.6 % of flagged variants did the same. All prioritized nonsense variants (4/4) and 26/27 missense variants had strong C-scores. Prioritized non-coding variant categories that correlated well with CADD include natural splicing variants (4/4), UTR structure variants (1/1), and RBPs (4/6). Weakly correlated variants included those affecting SRFBPs (5/36), TFBS (2/5), and pseudoexon activating variants (0/5). Missense mutations comprised 75 % of the flagged variants with C > 10. The aforementioned flagged splicing variant ATM c.1066-6 T > G also exceeded the threshold C value (C = 11.9). Meanwhile, the flagged TP53 variant, shown by SHAPE analysis to alter UTR structure, did not (C = 5.3). Despite consistency between some variant categories, the underlying assumptions of each approach probably explain why these results differ for non-coding variants. The limited numbers validated, deleterious non-coding variants also contributes to the accuracy of these predictions [163].

We verified prioritized protein-truncating (N = 7) and splicing (N = 4) variants by Sanger sequencing (Tables 2 and 4, respectively). In addition, two missense variants (BRCA2 c.7958 T > C and CHEK2 c.433C > T) were re-sequenced, since they are indicated as likely pathogenic/pathogenic in ClinVar (Additional file 12: Table S9). All protein-truncating variants were confirmed, with one exception (BRCA2 c.7558C > T, no evidence for the variant was present for either strand). Two of the mRNA splicing mutations were confirmed on both strands, while the other two were confirmed on a single strand (ATM c.6347 + 1G > T and ATM c.1066-6 T > G). Both documented pathogenic missense variants were also confirmed.

### Capture, sequencing, and alignment

The average coverage of capture region per individual was 90.8x (range of 53.8 to 118.2x between 32 samples) with 98.8 % of the probe-covered nucleotides having ≥ 10 reads. Samples with fewer than 10 reads per nucleotide were re-sequenced and the results of both runs were combined. The combined coverage of these samples was, on average, 48.2x (±36.2).

The consistency of both library preparation and capture protocols was improved from initial runs, which significantly impacted sequence coverage (Additional file 1: Methods). Of the 102 patients tested, 14 had been previously Sanger sequenced for BRCA1 and BRCA2 exons. Confirmation of previously discovered SNVs served to assess the methodological improvements introduced during NGS and ultimately, to increase confidence in variant calling. Initially, only 15 of 49 SNVs in 3 samples were detected. The detection rate of SNVs was improved to 100 % as the protocol progressed. All known SNVs (N = 157) were called in subsequent sequencing runs where purification steps were replaced with solid phase reversible immobilization beads and where RNA bait was transcribed the same day as capture. To minimize false positive variant calls, sequence read data were aligned with CASAVA and CRAC, variants were called for each alignment with GATK, and discrepancies were then resolved by manual review.

GATK called 14,164 unique SNVs and 1147 indels. Only 3777 (15.3 %) SNVs were present in both CASAVA and CRAC-alignments for at least one patient, and even fewer indel calls were concordant between both methods (N = 110; 6.2 %). For all other SNVs and indels, CASAVA called 6871 and 1566, respectively, whereas CRAC called 13,958 and 110, respectively. Some variants were counted more than once if they were called by different alignment programs in two or more patients. Intronic and intergenic variants proximate to low complexity sequences tend to generate false positive variants due to ambiguous alignment, a well known technical issue in short read sequence analysis [146, 147], contributing to this discrepancy. For example, CRAC correctly called a 19 nt deletion of BRCA1 (rs80359876; also confirmed by Sanger sequencing) but CASAVA flagged the deleted segment as a series of false-positives (Additional file 6: Figure S1). For these reasons, all variants were manually reviewed.

### IT-based variant identification and prioritization

The Shannon Pipeline reported 99 unique variants in natural donor or acceptor SSs. After technical and frequency filtering criteria were applied, 12 variants remained (Additional file 7: Table S4). IT analysis allowed for the prioritization of 3 variants, summarized in Table 2.Table 2Variants prioritized by IT analysisPatient IDGenemRNArsID (dbSNP 142)Information ChangeConsequencef or Binding Factor Affected
R
i,initial

R
i,final
ΔR
i or R
i
e
Allele Frequency (%)d
(bits)(bits)(bits)Abolished Natural SS7-4 F
ATM
c.3747-1G > Aa
Novel11.00.1−10.9Exon skipping and use of alternative splice forms4-1 F
ATM
c.6347 + 1G > Tb
Novel10.4−8.3−18.6Exon skippingLeaky Natural SS4-2B
CHEK2
c.320-5 T > Aa
rs1219087006.84.1−2.7Leaky splicing with intron inclusion0.08Activated Cryptic SS7-3E
BRCA1
c.548-293G > Ars117281398−12.12.614.7Cryptic site not expected to be used. Total information for natural exon is stronger than cryptic exon.0.747-4A
BRCA2
c.7618-269_7618-260del10Novel3.99.45.5Cryptic site not expected to be used. Total information for natural exon is stronger than cryptic exon.Pseudoexon formation due to activated acceptor SS7-3 F
BRCA2
c.8332-805G > ANovel−9.35.45.6e
6065/211/592f
7-3D
CDH1
c.164-2023A > Grs184740925−6.64.36.5e
61,236/224/1798f
0.35-3H
CDH1
c.2296-174 T > Ars5654888667.38.55.0e
1175/50/124f
0.02Pseudoexon formation due to activated donor SS3-6A
BRCA1
c.212 + 253G > Ars1893521914.16.75.2e
186/63/1250f
0.085-2G
BRCA2
c.7007 + 2691G > Ars3678905774.77.27.7e
2589/103/5272f
0.02Affected TFBSs7-4B
BRCA1
c.-8895G > ANovel10.9−0.2−11.1GATA-3 (GATA3)5-3E
CDH1
c.-54G > Crs50308741.712.010.4E2F-4 (E2F4)7-4E0.165-2B
PALB2
c.-291C > Grs55282422712.1−1.3−13.4GABPα (GABPA)0.17-2 F
TP53
c.-28-3132 T > Crs17882863−6.310.917.2RUNX3 (RUNX3)0.34-1A
TP53
c.-28-1102 T > Crs1134516735.112.37.2E2F-4 (E2F4)0.48.012.94.8Sp1 (SP1)Affected RBBSs7-4G
ATM
c.-244 T > Ars5399482189.8−19.9−29.7RBFOXc.-744 T > A0.04c.-1929 T > Ac.-3515 T > A5-3C
CDH1
c.*424 T > ANovel−20.39.629.9SF3B48.21.8−6.4CELF47-2E
CHEK2
c.-588G > Ars14156834210.93.7−7.2BX511012.14-3C.5-4G
CHEK2
c.-345C > Tc
rs1378530073.311.48.2SF3B43-1A
TP53
c.-107 T > Crs11353009010.54.5−6.0ELAVL14-1Hc.-188 T > C0.724-2H
TP53
c.*1175A > Crs7837822210.74.1−6.6KHDRBS17-2 Fc.*1376A > C0.26c.*1464A > C
aConfirmed by Sanger sequencing
bAmbiguous Sanger sequencing results
cPrioritized under missense change and was therefore verified with Sanger sequencing. Variant was confirmed
dIf available
e
R
i of site of opposite polarity in the pseudoexon
fConsequences for pseudoexon formation describe how the intron is divided: “new intron A length/pseudoexon length/new exon B lengthNone of the variants have been previously reported by other groups with the exception of CHEK2 c.320-5T>A [148]

Variants prioritized by IT analysis

aConfirmed by Sanger sequencing

bAmbiguous Sanger sequencing results

cPrioritized under missense change and was therefore verified with Sanger sequencing. Variant was confirmed

dIf available

e
R
i of site of opposite polarity in the pseudoexon

fConsequences for pseudoexon formation describe how the intron is divided: “new intron A length/pseudoexon length/new exon B length

None of the variants have been previously reported by other groups with the exception of CHEK2 c.320-5T>A [148]

First, the novel ATM variant c.3747-1G > A (chr11:108,154,953G > A; sample number 7-4 F) abolishes the natural acceptor of exon 26 (11.0 to 0.1 bits). ASSEDA reports the presence of a 5.3 bit cryptic acceptor site 13 nt downstream of the natural site, but the effect of the variant on a pre-existing cryptic site is negligible (~0.1 bits). The cryptic exon would lead to exon deletion and frameshift (Fig. 3a). ASSEDA also predicts skipping of the 246 nt exon, as the Ri,final of the natural acceptor is now below Ri,minimum (1.6 bits), altering the reading frame. Second, the novel ATM c.6347 + 1G > T (chr11:108188249G > T; 4-1 F) abolishes the 10.4 bit natural donor site of exon 44 (ΔRi = -18.6 bits), and is predicted to cause exon skipping. Finally, the previously reported CHEK2 variant, c.320-5A > T (chr22:29,121,360 T > A; rs121908700; 4-2B) [148] weakens the natural acceptor of exon 3 (6.8 to 4.1 bits), and may activate a cryptic acceptor (7.4 bits) 92 nt upstream of the natural acceptor site which would shift the reading frame (Fig. 4). A constitutive, frameshifted alternative isoform of CHEK2 lacking exons 3 and 4 has been reported, but skipping of exon 3 alone is not normally observed.Fig. 3Predicted Isoforms and Relative Abundances as a Consequence of ATM splice variant c.3747-1G > A. Intronic ATM variant c.3747-1G > A abolishes (11.0 to 0.1 bits) the natural acceptor of exon 26 (total of 63 exons). a ASSEDA predicts skipping of the natural exon (R
i,total from 14.5 to 3.6 bits [a 1910 fold decrease in exon strength]; isoform 7) and/or activation of a pre-existing cryptic acceptor site 13 nt downstream (R
i,total for cryptic exon = 9.0 bits; isoform 1) of the natural site leading to exon truncation. The reading frame is altered in both mutant isoforms. The other isoforms use weak, alternate acceptor/donor sites leading to cryptic exons with much lower total information. b Before the mutation, isoform 7 is expected to be the most abundant splice form. c After the mutation, isoform 1 is predicted to become the most abundant splice form and the wild-type isoform is not expected to be expressedFig. 4Predicted Isoforms and Relative Abundances as a Consequence of CHEK2 splice variant c.320-5 T > A. Intronic CHEK2 variant c.320-5 T > A weakens (6.8 to 4.1 bits) the natural acceptor of exon 3 (total of 15 exons). a ASSEDA reports the weakening of the natural exon strength (R
i,total reduced from 13.2 to 10.5 bits), which would result in reduced splicing of the exon otherwise known as leaky splicing. A pre-existing cryptic acceptor exists 92 nt upstream of the natural site, leading to a cryptic exon with similar strength to the mutated exon (R
i,total = 10.0 bits). This cryptic exon would contain 92 nt of the intron. b Before the mutation, isoform 1 is expected to be the only isoform expressed. c After the mutation, isoform 1 (wild-type) is predicted to become relatively less abundant and isoform 2 is expected to be expressed, although less abundant in relation to isoform 1

Predicted Isoforms and Relative Abundances as a Consequence of ATM splice variant c.3747-1G > A. Intronic ATM variant c.3747-1G > A abolishes (11.0 to 0.1 bits) the natural acceptor of exon 26 (total of 63 exons). a ASSEDA predicts skipping of the natural exon (R
i,total from 14.5 to 3.6 bits [a 1910 fold decrease in exon strength]; isoform 7) and/or activation of a pre-existing cryptic acceptor site 13 nt downstream (R
i,total for cryptic exon = 9.0 bits; isoform 1) of the natural site leading to exon truncation. The reading frame is altered in both mutant isoforms. The other isoforms use weak, alternate acceptor/donor sites leading to cryptic exons with much lower total information. b Before the mutation, isoform 7 is expected to be the most abundant splice form. c After the mutation, isoform 1 is predicted to become the most abundant splice form and the wild-type isoform is not expected to be expressed

Predicted Isoforms and Relative Abundances as a Consequence of CHEK2 splice variant c.320-5 T > A. Intronic CHEK2 variant c.320-5 T > A weakens (6.8 to 4.1 bits) the natural acceptor of exon 3 (total of 15 exons). a ASSEDA reports the weakening of the natural exon strength (R
i,total reduced from 13.2 to 10.5 bits), which would result in reduced splicing of the exon otherwise known as leaky splicing. A pre-existing cryptic acceptor exists 92 nt upstream of the natural site, leading to a cryptic exon with similar strength to the mutated exon (R
i,total = 10.0 bits). This cryptic exon would contain 92 nt of the intron. b Before the mutation, isoform 1 is expected to be the only isoform expressed. c After the mutation, isoform 1 (wild-type) is predicted to become relatively less abundant and isoform 2 is expected to be expressed, although less abundant in relation to isoform 1

Variants either strengthening (N = 4) or slightly weakening (ΔRi < 1.0 bits; N = 4) a natural site were not prioritized. In addition, we rejected the ATM variant (c.1066-6 T > G; chr11:108,119,654 T > G; 4-1E and 7-2B), which slightly weakens the natural acceptor of exon 9 (11.0 to 8.1 bits). Although other studies have shown leaky expression as a result of this variant [149], a more recent meta-analysis concluded that this variant is not associated with increased breast cancer risk [150].

Two variants produced information changes that could potentially impact cryptic splicing, but were not prioritized for the following reasons (Table 2). The first variant, novel BRCA2 deletion c.7618-269_7618-260del10 (chr13:32931610_32931619del10; 7-4A) strengthens a cryptic acceptor site 245 nt upstream from the natural acceptor of exon 16 (Ri,final = 9.4 bits, ΔRi = 5.5 bits). Being 5.7-fold stronger than the natural site (6.9 bits), two potential cryptic isoforms were predicted, however the exon strengths of both are weaker than the unaffected natural exon (Ri,total = 6.6 bits) and thus neither were prioritized. The larger gap surprisal penalties explain the differences in exon strength. The natural donor SS may still be used in conjunction with the abovementioned cryptic SS, resulting in an exon with Ri,total = 3.5 bits. Alternatively, the cryptic site and a weak donor site 180 nt upstream of the natural donor (Ri = 0.7 vs 1.4, cryptic and natural donors, respectively) result in an exon with Ri,total = 6.5 bits. The second variant, BRCA1 c.548-293G > A (chr17:41249599C > T; 7-3E), creates a weak cryptic acceptor (Ri,final = 2.6 bits, ΔRi = 6.2 bits) 291 nt upstream of the natural acceptor for exon 8 (Ri = 0.5). Although the cryptic exon is strengthened (final Ri,total = 6.9 bits, ΔRi = 14.7 bits), ASSEDA predicts the level of expression of this exon to be negligible, as it is weaker than the natural exon (Ri,total = 8.4 bits) due to the increased length of the predicted exon (+291 nt) [38].

The Shannon Pipeline initially reported 1583 unique variants creating or strengthening intronic cryptic sites. We prioritized 5 variants, 1 of which is novel (BRCA2 c.8332-805G > A; 7-3 F), that were within 250 nt of a pre-existing complementary cryptic site and have an hnRNPA1 site within 5 nt of the acceptor (Table 2). If used, 3 of these pseudoexons would lead to a frameshifted transcript.

Variants within 500 nt of an exon junction and all exonic variants (N = 4015) were investigated for their potential effects on affinity of sites to corresponding SRFs [38]. IT analysis flagged 54 variants significantly altering the strength of at least one binding site (Additional file 8: Table S5). A careful review of the variants, the factor affected, and the position of the binding site relative to the natural SS, prioritized 36 variants (21 novel), of which 4 are in exons and 32 are in introns. As an example, a novel CHEK2 exon 2 variant c.69C > A (p.Gly23=) is predicted to increase the strength of an hnRNP A1 site (0.7 to 5.3 bits) and decrease total exon strength (ΔRi,total = -5.7 bits). A similar type of exonic variant in FANCM, which was predicted to create an exonic hnRNP A1 site by IT, has been shown to bind this exonic repressor and induce exon skipping [37].

We assessed SNVs with models of 83 TFs experimentally shown to bind (Additional file 3: Table S1) upstream or within the first exon and intron of our sequenced genes (N = 2177). Thirteen variants expected to significantly affect TF binding were flagged (Additional file 9: Table S6). The final filtering step considered the known function of the TF in transcription, resulting in 5 variants (Table 2) in 6 patients (one variant was identified in two patients). Four of these variants have been previously reported (rs5030874, rs552824227, rs17882863, rs113451673) and one is novel (c.-8895G > A; 7-4B).

There were 364 unique UTR variants found by sequencing. These variants were evaluated for their effects on mRNA 2° structure (including that of splice forms with alternate UTRs in the cases of BRCA1 and TP53) through SNPfold, resulting in 5 flagged variants (Table 3), all of which have been previously reported.Table 3Variants predicted by SNPfold to affect UTR structureClassa
Patient IDGenemRNAUTR positionrsID (dbSNP 142)Ranke

p-valueAllele Frequency (%)d
FIn 26 patients
BRCA2
b
c.-52A > G5’ UTRrs2061182/9000.00214.86FIn 40 patients
BRCA2
b
c.*532A > G3’ UTRrs11571836239/27000.08919.75P7-4C
CDH1
c
c.-71C > G5’ UTRrs3403377169/6000.1150.56F4-2E
TP53
b
c.*485G > A3’ UTRrs4968187169/45000.0385-4A5.11F2-1A, 7-1B, 5-2A.7-1D, 7-2B, 7-2F
TP53
b
c.*826G > A3’ UTRrs17884306371/45000.0827-4C5.71
aF:Flagged; P:Prioritized
bLong Range UTR SNPfold Analysis
cLocal Range SNPfold Analysis
dIf available
eRank of the SNP, in terms of how much it changes the mRNA structure compared to all other possible mutations

Variants predicted by SNPfold to affect UTR structure

aF:Flagged; P:Prioritized

bLong Range UTR SNPfold Analysis

cLocal Range SNPfold Analysis

dIf available

eRank of the SNP, in terms of how much it changes the mRNA structure compared to all other possible mutations

Analysis of three variants using mFOLD [83] revealed likely changes to the UTR structure (Fig. 5). Two variants with possible 2° structure effects were common (BRCA2 c.-52A > G [N = 26 samples] and c.*532A > G [N = 40]) and not prioritized. The 5’ UTR CDH1 variant c.-71C > G (chr16:68771248C > G; rs34033771; 7-4C) disrupts a double-stranded hairpin region to create a larger loop structure, thus increasing binding accessibility (Fig. 5a and b). Analysis using RBPDB and CISBP-RNA-derived IT models suggests this variant affects binding by NCL (Nucleolin, a transcription coactivator) by decreasing binding affinity 14-fold (Ri,initial = 6.6 bits, ΔRi = -3.8 bits) (Additional file 10: Table S7). This RBP has been shown to bind to the 5’ and 3’ UTR of p53 mRNA and plays a role in repressing its translation [151].Fig. 5Predicted Alteration in UTR Structure Using mFOLD for Variants Flagged by SNPfold. Wild-type and variant structures are displayed, with the variant indicated by a red arrow. a Predicted wild-type structure of CDH1 5’UTR surrounding c.-71. b Predicted CDH1 5’UTR structure due to c.-71C > G variant. c Predicted wild-type TP53 3’UTR structure surrounding c.*485. d Predicted TP53 5’UTR structure due to c.*485G > A variant. e Predicted wild-type TP53 3’UTR structure surrounding c.*826. f Predicted TP53 5’UTR structure due to c.*826G > A variant. §SHAPE analysis revealed differences in reactivity between mutant and variant mRNAs, confirming alterations to 2° structure

Predicted Alteration in UTR Structure Using mFOLD for Variants Flagged by SNPfold. Wild-type and variant structures are displayed, with the variant indicated by a red arrow. a Predicted wild-type structure of CDH1 5’UTR surrounding c.-71. b Predicted CDH1 5’UTR structure due to c.-71C > G variant. c Predicted wild-type TP53 3’UTR structure surrounding c.*485. d Predicted TP53 5’UTR structure due to c.*485G > A variant. e Predicted wild-type TP53 3’UTR structure surrounding c.*826. f Predicted TP53 5’UTR structure due to c.*826G > A variant. §SHAPE analysis revealed differences in reactivity between mutant and variant mRNAs, confirming alterations to 2° structure

In addition, the TP53 variant c.*485G > A (NM_000546.5: chr17:7572442C > T; rs4968187) is found at the 3’ UTR and was identified in two patients (4-2E and 5-4A). In silico mRNA folding analysis demonstrated this variant disrupts a G/C bond of a loop in the highest ranked potential mRNA structure (Fig. 5c and d). Also, SHAPE analysis showed a difference in 2° structure between the wild-type and mutant (data not shown). IT analysis with RBBS models indicated that this variant significantly increases the binding affinity of SF3B4 by > 48-fold (Ri,final = 11.0 bits, ΔRi = 5.6 bits) (Additional file 10: Table S7). This RBP is one of four subunits comprising the splice factor 3B, which binds upstream of the branch-point sequence in pre-mRNA [152].

The third flagged variant also occurs in the 3’ UTR of TP53 (c.*826G > A; chr17:7572,101C > T; rs17884306), and was identified in 6 patients (2-1A, 7-1B, 5-2A.7-1D, 7-2B, 7-2 F, and 7-4C). It disrupts a potential loop structure, stabilizing a double-stranded hairpin, and possibly making it less accessible (Fig. 5e and f). Analysis using RBPDB-derived models suggests this variant could affect the binding of both RBFOX2 and SF3B4 (Additional file 10: Table S7). A binding site for RBFOX2, which acts as a promoter of alternative splicing by favoring the inclusion of alternative exons [153], is created (Ri,final = 9.8 bits; ΔRi = -6.5 bits). This variant is also expected to simultaneously abolish a SF3B4 binding site (Ri,final = -20.3 bits; ΔRi = -29.9 bits).

RBPDB- and CISBP-RNA-derived information model analysis of all UTR variants resulted in the prioritization of 1 novel, and 5 previously-reported variants (Table 2). No patient within the cohort exhibited more than one prioritized RBBS variant.

To evaluate the background rate of prioritizing variants flagged by this method, all 5’ and 3’ UTR SNVs in dbSNP144 for the 7 genes sequenced (excluding those already flagged in Table 3) were evaluated by SNPfold and our RBP information models. Of 1207 SNVs, only 10 were prioritized with both methods, which results in a background rate of 0.83 %.

### Natural SS variants

The Shannon Pipeline reported 99 unique variants in natural donor or acceptor SSs. After technical and frequency filtering criteria were applied, 12 variants remained (Additional file 7: Table S4). IT analysis allowed for the prioritization of 3 variants, summarized in Table 2.Table 2Variants prioritized by IT analysisPatient IDGenemRNArsID (dbSNP 142)Information ChangeConsequencef or Binding Factor Affected
R
i,initial

R
i,final
ΔR
i or R
i
e
Allele Frequency (%)d
(bits)(bits)(bits)Abolished Natural SS7-4 F
ATM
c.3747-1G > Aa
Novel11.00.1−10.9Exon skipping and use of alternative splice forms4-1 F
ATM
c.6347 + 1G > Tb
Novel10.4−8.3−18.6Exon skippingLeaky Natural SS4-2B
CHEK2
c.320-5 T > Aa
rs1219087006.84.1−2.7Leaky splicing with intron inclusion0.08Activated Cryptic SS7-3E
BRCA1
c.548-293G > Ars117281398−12.12.614.7Cryptic site not expected to be used. Total information for natural exon is stronger than cryptic exon.0.747-4A
BRCA2
c.7618-269_7618-260del10Novel3.99.45.5Cryptic site not expected to be used. Total information for natural exon is stronger than cryptic exon.Pseudoexon formation due to activated acceptor SS7-3 F
BRCA2
c.8332-805G > ANovel−9.35.45.6e
6065/211/592f
7-3D
CDH1
c.164-2023A > Grs184740925−6.64.36.5e
61,236/224/1798f
0.35-3H
CDH1
c.2296-174 T > Ars5654888667.38.55.0e
1175/50/124f
0.02Pseudoexon formation due to activated donor SS3-6A
BRCA1
c.212 + 253G > Ars1893521914.16.75.2e
186/63/1250f
0.085-2G
BRCA2
c.7007 + 2691G > Ars3678905774.77.27.7e
2589/103/5272f
0.02Affected TFBSs7-4B
BRCA1
c.-8895G > ANovel10.9−0.2−11.1GATA-3 (GATA3)5-3E
CDH1
c.-54G > Crs50308741.712.010.4E2F-4 (E2F4)7-4E0.165-2B
PALB2
c.-291C > Grs55282422712.1−1.3−13.4GABPα (GABPA)0.17-2 F
TP53
c.-28-3132 T > Crs17882863−6.310.917.2RUNX3 (RUNX3)0.34-1A
TP53
c.-28-1102 T > Crs1134516735.112.37.2E2F-4 (E2F4)0.48.012.94.8Sp1 (SP1)Affected RBBSs7-4G
ATM
c.-244 T > Ars5399482189.8−19.9−29.7RBFOXc.-744 T > A0.04c.-1929 T > Ac.-3515 T > A5-3C
CDH1
c.*424 T > ANovel−20.39.629.9SF3B48.21.8−6.4CELF47-2E
CHEK2
c.-588G > Ars14156834210.93.7−7.2BX511012.14-3C.5-4G
CHEK2
c.-345C > Tc
rs1378530073.311.48.2SF3B43-1A
TP53
c.-107 T > Crs11353009010.54.5−6.0ELAVL14-1Hc.-188 T > C0.724-2H
TP53
c.*1175A > Crs7837822210.74.1−6.6KHDRBS17-2 Fc.*1376A > C0.26c.*1464A > C
aConfirmed by Sanger sequencing
bAmbiguous Sanger sequencing results
cPrioritized under missense change and was therefore verified with Sanger sequencing. Variant was confirmed
dIf available
e
R
i of site of opposite polarity in the pseudoexon
fConsequences for pseudoexon formation describe how the intron is divided: “new intron A length/pseudoexon length/new exon B lengthNone of the variants have been previously reported by other groups with the exception of CHEK2 c.320-5T>A [148]

Variants prioritized by IT analysis

aConfirmed by Sanger sequencing

bAmbiguous Sanger sequencing results

cPrioritized under missense change and was therefore verified with Sanger sequencing. Variant was confirmed

dIf available

e
R
i of site of opposite polarity in the pseudoexon

fConsequences for pseudoexon formation describe how the intron is divided: “new intron A length/pseudoexon length/new exon B length

None of the variants have been previously reported by other groups with the exception of CHEK2 c.320-5T>A [148]

First, the novel ATM variant c.3747-1G > A (chr11:108,154,953G > A; sample number 7-4 F) abolishes the natural acceptor of exon 26 (11.0 to 0.1 bits). ASSEDA reports the presence of a 5.3 bit cryptic acceptor site 13 nt downstream of the natural site, but the effect of the variant on a pre-existing cryptic site is negligible (~0.1 bits). The cryptic exon would lead to exon deletion and frameshift (Fig. 3a). ASSEDA also predicts skipping of the 246 nt exon, as the Ri,final of the natural acceptor is now below Ri,minimum (1.6 bits), altering the reading frame. Second, the novel ATM c.6347 + 1G > T (chr11:108188249G > T; 4-1 F) abolishes the 10.4 bit natural donor site of exon 44 (ΔRi = -18.6 bits), and is predicted to cause exon skipping. Finally, the previously reported CHEK2 variant, c.320-5A > T (chr22:29,121,360 T > A; rs121908700; 4-2B) [148] weakens the natural acceptor of exon 3 (6.8 to 4.1 bits), and may activate a cryptic acceptor (7.4 bits) 92 nt upstream of the natural acceptor site which would shift the reading frame (Fig. 4). A constitutive, frameshifted alternative isoform of CHEK2 lacking exons 3 and 4 has been reported, but skipping of exon 3 alone is not normally observed.Fig. 3Predicted Isoforms and Relative Abundances as a Consequence of ATM splice variant c.3747-1G > A. Intronic ATM variant c.3747-1G > A abolishes (11.0 to 0.1 bits) the natural acceptor of exon 26 (total of 63 exons). a ASSEDA predicts skipping of the natural exon (R
i,total from 14.5 to 3.6 bits [a 1910 fold decrease in exon strength]; isoform 7) and/or activation of a pre-existing cryptic acceptor site 13 nt downstream (R
i,total for cryptic exon = 9.0 bits; isoform 1) of the natural site leading to exon truncation. The reading frame is altered in both mutant isoforms. The other isoforms use weak, alternate acceptor/donor sites leading to cryptic exons with much lower total information. b Before the mutation, isoform 7 is expected to be the most abundant splice form. c After the mutation, isoform 1 is predicted to become the most abundant splice form and the wild-type isoform is not expected to be expressedFig. 4Predicted Isoforms and Relative Abundances as a Consequence of CHEK2 splice variant c.320-5 T > A. Intronic CHEK2 variant c.320-5 T > A weakens (6.8 to 4.1 bits) the natural acceptor of exon 3 (total of 15 exons). a ASSEDA reports the weakening of the natural exon strength (R
i,total reduced from 13.2 to 10.5 bits), which would result in reduced splicing of the exon otherwise known as leaky splicing. A pre-existing cryptic acceptor exists 92 nt upstream of the natural site, leading to a cryptic exon with similar strength to the mutated exon (R
i,total = 10.0 bits). This cryptic exon would contain 92 nt of the intron. b Before the mutation, isoform 1 is expected to be the only isoform expressed. c After the mutation, isoform 1 (wild-type) is predicted to become relatively less abundant and isoform 2 is expected to be expressed, although less abundant in relation to isoform 1

Predicted Isoforms and Relative Abundances as a Consequence of ATM splice variant c.3747-1G > A. Intronic ATM variant c.3747-1G > A abolishes (11.0 to 0.1 bits) the natural acceptor of exon 26 (total of 63 exons). a ASSEDA predicts skipping of the natural exon (R
i,total from 14.5 to 3.6 bits [a 1910 fold decrease in exon strength]; isoform 7) and/or activation of a pre-existing cryptic acceptor site 13 nt downstream (R
i,total for cryptic exon = 9.0 bits; isoform 1) of the natural site leading to exon truncation. The reading frame is altered in both mutant isoforms. The other isoforms use weak, alternate acceptor/donor sites leading to cryptic exons with much lower total information. b Before the mutation, isoform 7 is expected to be the most abundant splice form. c After the mutation, isoform 1 is predicted to become the most abundant splice form and the wild-type isoform is not expected to be expressed

Predicted Isoforms and Relative Abundances as a Consequence of CHEK2 splice variant c.320-5 T > A. Intronic CHEK2 variant c.320-5 T > A weakens (6.8 to 4.1 bits) the natural acceptor of exon 3 (total of 15 exons). a ASSEDA reports the weakening of the natural exon strength (R
i,total reduced from 13.2 to 10.5 bits), which would result in reduced splicing of the exon otherwise known as leaky splicing. A pre-existing cryptic acceptor exists 92 nt upstream of the natural site, leading to a cryptic exon with similar strength to the mutated exon (R
i,total = 10.0 bits). This cryptic exon would contain 92 nt of the intron. b Before the mutation, isoform 1 is expected to be the only isoform expressed. c After the mutation, isoform 1 (wild-type) is predicted to become relatively less abundant and isoform 2 is expected to be expressed, although less abundant in relation to isoform 1

Variants either strengthening (N = 4) or slightly weakening (ΔRi < 1.0 bits; N = 4) a natural site were not prioritized. In addition, we rejected the ATM variant (c.1066-6 T > G; chr11:108,119,654 T > G; 4-1E and 7-2B), which slightly weakens the natural acceptor of exon 9 (11.0 to 8.1 bits). Although other studies have shown leaky expression as a result of this variant [149], a more recent meta-analysis concluded that this variant is not associated with increased breast cancer risk [150].

### Cryptic SS activation

Two variants produced information changes that could potentially impact cryptic splicing, but were not prioritized for the following reasons (Table 2). The first variant, novel BRCA2 deletion c.7618-269_7618-260del10 (chr13:32931610_32931619del10; 7-4A) strengthens a cryptic acceptor site 245 nt upstream from the natural acceptor of exon 16 (Ri,final = 9.4 bits, ΔRi = 5.5 bits). Being 5.7-fold stronger than the natural site (6.9 bits), two potential cryptic isoforms were predicted, however the exon strengths of both are weaker than the unaffected natural exon (Ri,total = 6.6 bits) and thus neither were prioritized. The larger gap surprisal penalties explain the differences in exon strength. The natural donor SS may still be used in conjunction with the abovementioned cryptic SS, resulting in an exon with Ri,total = 3.5 bits. Alternatively, the cryptic site and a weak donor site 180 nt upstream of the natural donor (Ri = 0.7 vs 1.4, cryptic and natural donors, respectively) result in an exon with Ri,total = 6.5 bits. The second variant, BRCA1 c.548-293G > A (chr17:41249599C > T; 7-3E), creates a weak cryptic acceptor (Ri,final = 2.6 bits, ΔRi = 6.2 bits) 291 nt upstream of the natural acceptor for exon 8 (Ri = 0.5). Although the cryptic exon is strengthened (final Ri,total = 6.9 bits, ΔRi = 14.7 bits), ASSEDA predicts the level of expression of this exon to be negligible, as it is weaker than the natural exon (Ri,total = 8.4 bits) due to the increased length of the predicted exon (+291 nt) [38].

### Pseudoexon formation

The Shannon Pipeline initially reported 1583 unique variants creating or strengthening intronic cryptic sites. We prioritized 5 variants, 1 of which is novel (BRCA2 c.8332-805G > A; 7-3 F), that were within 250 nt of a pre-existing complementary cryptic site and have an hnRNPA1 site within 5 nt of the acceptor (Table 2). If used, 3 of these pseudoexons would lead to a frameshifted transcript.

### SRF binding

Variants within 500 nt of an exon junction and all exonic variants (N = 4015) were investigated for their potential effects on affinity of sites to corresponding SRFs [38]. IT analysis flagged 54 variants significantly altering the strength of at least one binding site (Additional file 8: Table S5). A careful review of the variants, the factor affected, and the position of the binding site relative to the natural SS, prioritized 36 variants (21 novel), of which 4 are in exons and 32 are in introns. As an example, a novel CHEK2 exon 2 variant c.69C > A (p.Gly23=) is predicted to increase the strength of an hnRNP A1 site (0.7 to 5.3 bits) and decrease total exon strength (ΔRi,total = -5.7 bits). A similar type of exonic variant in FANCM, which was predicted to create an exonic hnRNP A1 site by IT, has been shown to bind this exonic repressor and induce exon skipping [37].

### TF binding

We assessed SNVs with models of 83 TFs experimentally shown to bind (Additional file 3: Table S1) upstream or within the first exon and intron of our sequenced genes (N = 2177). Thirteen variants expected to significantly affect TF binding were flagged (Additional file 9: Table S6). The final filtering step considered the known function of the TF in transcription, resulting in 5 variants (Table 2) in 6 patients (one variant was identified in two patients). Four of these variants have been previously reported (rs5030874, rs552824227, rs17882863, rs113451673) and one is novel (c.-8895G > A; 7-4B).

### UTR structure and protein binding

There were 364 unique UTR variants found by sequencing. These variants were evaluated for their effects on mRNA 2° structure (including that of splice forms with alternate UTRs in the cases of BRCA1 and TP53) through SNPfold, resulting in 5 flagged variants (Table 3), all of which have been previously reported.Table 3Variants predicted by SNPfold to affect UTR structureClassa
Patient IDGenemRNAUTR positionrsID (dbSNP 142)Ranke

p-valueAllele Frequency (%)d
FIn 26 patients
BRCA2
b
c.-52A > G5’ UTRrs2061182/9000.00214.86FIn 40 patients
BRCA2
b
c.*532A > G3’ UTRrs11571836239/27000.08919.75P7-4C
CDH1
c
c.-71C > G5’ UTRrs3403377169/6000.1150.56F4-2E
TP53
b
c.*485G > A3’ UTRrs4968187169/45000.0385-4A5.11F2-1A, 7-1B, 5-2A.7-1D, 7-2B, 7-2F
TP53
b
c.*826G > A3’ UTRrs17884306371/45000.0827-4C5.71
aF:Flagged; P:Prioritized
bLong Range UTR SNPfold Analysis
cLocal Range SNPfold Analysis
dIf available
eRank of the SNP, in terms of how much it changes the mRNA structure compared to all other possible mutations

Variants predicted by SNPfold to affect UTR structure

aF:Flagged; P:Prioritized

bLong Range UTR SNPfold Analysis

cLocal Range SNPfold Analysis

dIf available

eRank of the SNP, in terms of how much it changes the mRNA structure compared to all other possible mutations

Analysis of three variants using mFOLD [83] revealed likely changes to the UTR structure (Fig. 5). Two variants with possible 2° structure effects were common (BRCA2 c.-52A > G [N = 26 samples] and c.*532A > G [N = 40]) and not prioritized. The 5’ UTR CDH1 variant c.-71C > G (chr16:68771248C > G; rs34033771; 7-4C) disrupts a double-stranded hairpin region to create a larger loop structure, thus increasing binding accessibility (Fig. 5a and b). Analysis using RBPDB and CISBP-RNA-derived IT models suggests this variant affects binding by NCL (Nucleolin, a transcription coactivator) by decreasing binding affinity 14-fold (Ri,initial = 6.6 bits, ΔRi = -3.8 bits) (Additional file 10: Table S7). This RBP has been shown to bind to the 5’ and 3’ UTR of p53 mRNA and plays a role in repressing its translation [151].Fig. 5Predicted Alteration in UTR Structure Using mFOLD for Variants Flagged by SNPfold. Wild-type and variant structures are displayed, with the variant indicated by a red arrow. a Predicted wild-type structure of CDH1 5’UTR surrounding c.-71. b Predicted CDH1 5’UTR structure due to c.-71C > G variant. c Predicted wild-type TP53 3’UTR structure surrounding c.*485. d Predicted TP53 5’UTR structure due to c.*485G > A variant. e Predicted wild-type TP53 3’UTR structure surrounding c.*826. f Predicted TP53 5’UTR structure due to c.*826G > A variant. §SHAPE analysis revealed differences in reactivity between mutant and variant mRNAs, confirming alterations to 2° structure

Predicted Alteration in UTR Structure Using mFOLD for Variants Flagged by SNPfold. Wild-type and variant structures are displayed, with the variant indicated by a red arrow. a Predicted wild-type structure of CDH1 5’UTR surrounding c.-71. b Predicted CDH1 5’UTR structure due to c.-71C > G variant. c Predicted wild-type TP53 3’UTR structure surrounding c.*485. d Predicted TP53 5’UTR structure due to c.*485G > A variant. e Predicted wild-type TP53 3’UTR structure surrounding c.*826. f Predicted TP53 5’UTR structure due to c.*826G > A variant. §SHAPE analysis revealed differences in reactivity between mutant and variant mRNAs, confirming alterations to 2° structure

In addition, the TP53 variant c.*485G > A (NM_000546.5: chr17:7572442C > T; rs4968187) is found at the 3’ UTR and was identified in two patients (4-2E and 5-4A). In silico mRNA folding analysis demonstrated this variant disrupts a G/C bond of a loop in the highest ranked potential mRNA structure (Fig. 5c and d). Also, SHAPE analysis showed a difference in 2° structure between the wild-type and mutant (data not shown). IT analysis with RBBS models indicated that this variant significantly increases the binding affinity of SF3B4 by > 48-fold (Ri,final = 11.0 bits, ΔRi = 5.6 bits) (Additional file 10: Table S7). This RBP is one of four subunits comprising the splice factor 3B, which binds upstream of the branch-point sequence in pre-mRNA [152].

The third flagged variant also occurs in the 3’ UTR of TP53 (c.*826G > A; chr17:7572,101C > T; rs17884306), and was identified in 6 patients (2-1A, 7-1B, 5-2A.7-1D, 7-2B, 7-2 F, and 7-4C). It disrupts a potential loop structure, stabilizing a double-stranded hairpin, and possibly making it less accessible (Fig. 5e and f). Analysis using RBPDB-derived models suggests this variant could affect the binding of both RBFOX2 and SF3B4 (Additional file 10: Table S7). A binding site for RBFOX2, which acts as a promoter of alternative splicing by favoring the inclusion of alternative exons [153], is created (Ri,final = 9.8 bits; ΔRi = -6.5 bits). This variant is also expected to simultaneously abolish a SF3B4 binding site (Ri,final = -20.3 bits; ΔRi = -29.9 bits).

RBPDB- and CISBP-RNA-derived information model analysis of all UTR variants resulted in the prioritization of 1 novel, and 5 previously-reported variants (Table 2). No patient within the cohort exhibited more than one prioritized RBBS variant.

To evaluate the background rate of prioritizing variants flagged by this method, all 5’ and 3’ UTR SNVs in dbSNP144 for the 7 genes sequenced (excluding those already flagged in Table 3) were evaluated by SNPfold and our RBP information models. Of 1207 SNVs, only 10 were prioritized with both methods, which results in a background rate of 0.83 %.

### Exonic variants altering protein sequence

Exonic variants called by GATK (N = 245) included insertions, deletions, nonsense, missense, and synonymous changes.

We identified 3 patients with different indels (Table 4). One was a PALB2 insertion c.1617_1618insTT (chr16:23646249_23646250insAA; 5-3A) in exon 4, previously reported in ClinVar as pathogenic. This mutation results in a frameshift and premature translation termination by 626 residues, abolishing domain interactions with RAD51, BRCA2, and POLH [137]. We also identified two known frameshift mutations in BRCA1: c.4964_4982del19 in exon 15 (chr17:41222949_41222967del19; rs80359876; 5-1B) and c.5266_5267insC in exon 19 (chr17:41209079_41209080insG; rs397507247; 5-3C) [148, 154]. Both are indicated as pathogenic and common in the BIC Database due to the loss of one or both C-terminal BRCT repeat domains [137]. Truncation of these domains produces instability and impairs nuclear transcript localization [155], and this bipartite domain is responsible for binding phosphoproteins that are phosphorylated in response to DNA damage [156, 157].Table 4Variants resulting in premature protein truncationPatient IDGeneExonmRNA ProteinrsID (dbSNP 142)ClinVard,e,f
DetailsRefAllele Frequency (%)c
Insertions/Deletions5-1B
BRCA1
15 of 23c.4964_4982del19a
rs803598766d; Pathogenic/likely pathogenice; Familial breast and breast-ovarian cancer, Hereditary cancer-predisposing syndromef.STOP at p.1670-p.Ser1655Tyrfs193 AA short5-3C
BRCA1
19 of 23c.5266_5267insCa
rs39750724713d; Pathogenic, risk factore; Familial breast, breast-ovarian, and pancreatic cancer, Hereditary cancer-predisposing syndromef.STOP at p.1788[148, 154]p.Gln1756Profs75 AA short5-3A
PALB2
4 of 13c.1617_1618insTTa
-1d; Pathogenice; Hereditary cancer-predisposing syndromef.STOP at p.561-p.Asn540Leufs626 AA shortStop Codons7-1G
BRCA2
15 of 27c.7558C > Tb
rs803589815d; Pathogenice; Familial breast, and breast-ovarian cancer, Hereditary cancer-predisposing syndromef.899 AA short[158]p.Arg2520Ter4-4A
BRCA2
25 of 27c.9294C > Ga
rs803592003d; Pathogenice; Familial breast and breast-ovarian cancerf.321 AA short[159]p.Tyr3098Ter7-3A
PALB2
4 of 13c.1240C > Ta
rs1801771003d; Pathogenice; Familial breast cancer, Hereditary cancer-predisposing syndromef.773 AA short[58]p.Arg414Ter4-4D
PALB2
4 of 13c.1042C > Ta
Novel-839 AA short-p.Gln348Ter
aConfirmed by Sanger sequencing
bNot confirmed by Sanger sequencing
cIf available
dNumber of submissions
eClinical significance
fCondition(s)

Variants resulting in premature protein truncation

aConfirmed by Sanger sequencing

bNot confirmed by Sanger sequencing

cIf available

dNumber of submissions

eClinical significance

fCondition(s)

We also identified 4 nonsense mutations, one of which was novel in exon 4 of PALB2 (c.1042C > T; chr16:23646825G > A; 4-4D). Another in PALB2 has been previously reported (c.1240C > T; chr16:23646627G > A; rs180177100; 7-3A) [58]. As a consequence, functional domains of PALB2 that interact with BRCA1, RAD51, BRCA2, and POLH are lost [137]. Two known nonsense mutations were found in BRCA2, c.7558C > T in exon 15 [158] and c.9294C > G in exon 25 [159]. The first (chr13:32930687C > T; rs80358981; 7-1G) causes the loss of the BRCA2 region that binds FANCD2, responsible for loading BRCA2 onto damaged chromatin [160]. The second (chr13:32968863C > G, rs80359200; 4-4A) does not occur within a known functional domain, however the transcript is likely to be degraded by nonsense mediated decay [161].

GATK called 61 missense variants, of which 18 were identified in 6 patients or more and 19 had allele frequencies > 1.0 % (Additional file 11: Table S8). The 40 remaining variants (15 ATM, 8 BRCA1, 9 BRCA2, 2 CDH1, 2 CHEK2, 3 PALB2, and 1 TP53) were assessed using a combination of gene specific databases, published classifications, and 4 in silico tools (Additional file 12: Table S9). We prioritized 27 variants, 2 of which were novel. None of the non-prioritized variants were predicted to be damaging by more than 2 of 4 conservation-based software programs.

### Protein-truncating variants

We identified 3 patients with different indels (Table 4). One was a PALB2 insertion c.1617_1618insTT (chr16:23646249_23646250insAA; 5-3A) in exon 4, previously reported in ClinVar as pathogenic. This mutation results in a frameshift and premature translation termination by 626 residues, abolishing domain interactions with RAD51, BRCA2, and POLH [137]. We also identified two known frameshift mutations in BRCA1: c.4964_4982del19 in exon 15 (chr17:41222949_41222967del19; rs80359876; 5-1B) and c.5266_5267insC in exon 19 (chr17:41209079_41209080insG; rs397507247; 5-3C) [148, 154]. Both are indicated as pathogenic and common in the BIC Database due to the loss of one or both C-terminal BRCT repeat domains [137]. Truncation of these domains produces instability and impairs nuclear transcript localization [155], and this bipartite domain is responsible for binding phosphoproteins that are phosphorylated in response to DNA damage [156, 157].Table 4Variants resulting in premature protein truncationPatient IDGeneExonmRNA ProteinrsID (dbSNP 142)ClinVard,e,f
DetailsRefAllele Frequency (%)c
Insertions/Deletions5-1B
BRCA1
15 of 23c.4964_4982del19a
rs803598766d; Pathogenic/likely pathogenice; Familial breast and breast-ovarian cancer, Hereditary cancer-predisposing syndromef.STOP at p.1670-p.Ser1655Tyrfs193 AA short5-3C
BRCA1
19 of 23c.5266_5267insCa
rs39750724713d; Pathogenic, risk factore; Familial breast, breast-ovarian, and pancreatic cancer, Hereditary cancer-predisposing syndromef.STOP at p.1788[148, 154]p.Gln1756Profs75 AA short5-3A
PALB2
4 of 13c.1617_1618insTTa
-1d; Pathogenice; Hereditary cancer-predisposing syndromef.STOP at p.561-p.Asn540Leufs626 AA shortStop Codons7-1G
BRCA2
15 of 27c.7558C > Tb
rs803589815d; Pathogenice; Familial breast, and breast-ovarian cancer, Hereditary cancer-predisposing syndromef.899 AA short[158]p.Arg2520Ter4-4A
BRCA2
25 of 27c.9294C > Ga
rs803592003d; Pathogenice; Familial breast and breast-ovarian cancerf.321 AA short[159]p.Tyr3098Ter7-3A
PALB2
4 of 13c.1240C > Ta
rs1801771003d; Pathogenice; Familial breast cancer, Hereditary cancer-predisposing syndromef.773 AA short[58]p.Arg414Ter4-4D
PALB2
4 of 13c.1042C > Ta
Novel-839 AA short-p.Gln348Ter
aConfirmed by Sanger sequencing
bNot confirmed by Sanger sequencing
cIf available
dNumber of submissions
eClinical significance
fCondition(s)

Variants resulting in premature protein truncation

aConfirmed by Sanger sequencing

bNot confirmed by Sanger sequencing

cIf available

dNumber of submissions

eClinical significance

fCondition(s)

We also identified 4 nonsense mutations, one of which was novel in exon 4 of PALB2 (c.1042C > T; chr16:23646825G > A; 4-4D). Another in PALB2 has been previously reported (c.1240C > T; chr16:23646627G > A; rs180177100; 7-3A) [58]. As a consequence, functional domains of PALB2 that interact with BRCA1, RAD51, BRCA2, and POLH are lost [137]. Two known nonsense mutations were found in BRCA2, c.7558C > T in exon 15 [158] and c.9294C > G in exon 25 [159]. The first (chr13:32930687C > T; rs80358981; 7-1G) causes the loss of the BRCA2 region that binds FANCD2, responsible for loading BRCA2 onto damaged chromatin [160]. The second (chr13:32968863C > G, rs80359200; 4-4A) does not occur within a known functional domain, however the transcript is likely to be degraded by nonsense mediated decay [161].

### Missense

GATK called 61 missense variants, of which 18 were identified in 6 patients or more and 19 had allele frequencies > 1.0 % (Additional file 11: Table S8). The 40 remaining variants (15 ATM, 8 BRCA1, 9 BRCA2, 2 CDH1, 2 CHEK2, 3 PALB2, and 1 TP53) were assessed using a combination of gene specific databases, published classifications, and 4 in silico tools (Additional file 12: Table S9). We prioritized 27 variants, 2 of which were novel. None of the non-prioritized variants were predicted to be damaging by more than 2 of 4 conservation-based software programs.

### Variant classification

Initially, 15,311 unique variants were identified by complete gene sequencing of 7 HBOC genes. Of these, 132 were flagged after filtering, and further reduced by IT-based variant analysis and consultation of the published literature to 87 prioritized variants. Figure 6 illustrates the decrease in the number of unique variants per patient at each step of our identification and prioritization process. The distribution of prioritized variants by gene is 34 in ATM, 13 in BRCA1, 11 in BRCA2, 8 in CDH1, 6 in CHEK2, 10 in PALB2, and 5 in TP53 (Additional file 13: Table S10), which are categorized by type in Table 5.Fig. 6Ladder Plot Representing Variant Identification and Prioritization. Each line is representative of a different sample in each sequencing run (a-e), illustrating the number of unique variants at important steps throughout the variant prioritization process. The left-most point indicates the total number of unique variants. The second point represents the number of unique variants remaining after common (>5 patients within cohort and/or ≥ 1.0 % allele frequency) and false-positive variants were removed. The right-most point represents the final number of unique. No variants were prioritized in the following patients: 2-1A, 2-5A, 2-6A, 3-2A, 3-3A, 3-4A, 3-5A, 3-8A, 4-1B, 4-2C, 4-2 F, 4-3B, 4-3D, 4-4B, 4-4E, 5-1G, 5-1H, 5-3D, 5-4C, 5-4D, 5-4 F, 5-4G, 5-4H, 7-1B, 7-1C, 7-1D, 7-1H, 7-2B, 7-2C, 7-2H, 7-3H, 7-4A, 7-4D, 7-4H. The average number of variants per patient at each step is indicated in a table below each plot, along with the percent reduction in variants from one step to anotherTable 5Summary of prioritized variants by geneIndelNonsenseMissenseNatural SplicingCryptic SplicingPseudoexonSR FactorTFUTR StructureUTR BindingTotal
ATM
00142001800134 a

BRCA1
202001710013
BRCA2
023002400011
CDH1
00200211118
CHEK2
00210030026 a

PALB2
123000310010
TP53
00100002025Total342730536516Three variants were prioritized under multiple categories: ATM chr11:108121730A > G (missense and SRFBS), CHEK2 chr22:29121242G > A (missense, UTR binding), and CHEK2 chr22:29130520C > T (missense, UTR binding)
a Counts represent the number of unique variants identified (i.e. a variant is not counted twice if it appeared in multiple individuals)

Ladder Plot Representing Variant Identification and Prioritization. Each line is representative of a different sample in each sequencing run (a-e), illustrating the number of unique variants at important steps throughout the variant prioritization process. The left-most point indicates the total number of unique variants. The second point represents the number of unique variants remaining after common (>5 patients within cohort and/or ≥ 1.0 % allele frequency) and false-positive variants were removed. The right-most point represents the final number of unique. No variants were prioritized in the following patients: 2-1A, 2-5A, 2-6A, 3-2A, 3-3A, 3-4A, 3-5A, 3-8A, 4-1B, 4-2C, 4-2 F, 4-3B, 4-3D, 4-4B, 4-4E, 5-1G, 5-1H, 5-3D, 5-4C, 5-4D, 5-4 F, 5-4G, 5-4H, 7-1B, 7-1C, 7-1D, 7-1H, 7-2B, 7-2C, 7-2H, 7-3H, 7-4A, 7-4D, 7-4H. The average number of variants per patient at each step is indicated in a table below each plot, along with the percent reduction in variants from one step to another

Summary of prioritized variants by gene

Three variants were prioritized under multiple categories: ATM chr11:108121730A > G (missense and SRFBS), CHEK2 chr22:29121242G > A (missense, UTR binding), and CHEK2 chr22:29130520C > T (missense, UTR binding)

a Counts represent the number of unique variants identified (i.e. a variant is not counted twice if it appeared in multiple individuals)

Three prioritized variants have multiple predicted roles: ATM c.1538A > G in missense and SRFBS, CHEK2 c.190G > A in missense and UTR binding, and CHEK2 c.433C > T in missense and UTR binding. Of the 102 patients that were sequenced, 72 (70.6 %) exhibited at least one prioritized variant, and some patients harbored more than one prioritized variant (N = 33; 32 %). Additional file 14: Table S11 presents a summary of all flagged and prioritized variants for patients with at least one prioritized variant.

### Prioritization of potential deletions

Using BreakDancer, none of the individuals analyzed exhibited large rearrangements that met the level of stringency required, but a small intragenic rearrangement in BRCA1 was identified and confirmed by Sanger sequencing. Attempts to detect deletions with BreakDancer only flagged single, non-contiguous paired-end reads, rather than a series of reads clustered within the same region within the same individual, which would be necessary to indicate the presence of a true deletion or structural rearrangement.

After prioritizing individuals for potential hemizygosity in the sequenced regions, potential deletions were detected in BRCA2 and CDH1. Patient UWO5-4D exhibited a non-polymorphic 32.1 kb interval in BRCA2, spanning introns 1 to 13, that was absent from all of the other individuals (chr13:32890227-32922331). Haploview (hapmap.org) showed very low levels of LD in this region. The potential deletion may extend further downstream, however the presence of a haploblock, covering the entire sequenced interval beyond exon 11, with significant LD precludes delineation of the telomeric breakpoint. We also flagged a non-polymorphic 2.6 kb interval near the 3’ end of CDH1 in 6 individuals (UWO3-5, UWO4-2C, UWO4-4E, UWO4-4 F, UWO4-2G, UWO5-2H). This is a low LD region spanning chr16:68861286-68863887 that includes exons 14 and 15, and is polymorphic in all of the other individuals sequenced. CDH1 mutations are characteristically present in families with predisposition to gastric cancer, however breast cancer frequently co-occurs [69]. A study of CDH1 deletions in inherited gastric cancer identified two families with deletions that overlap the intervals prioritized in the present study [162].

### Comparison to combined annotation dependent depletion

The analysis and prioritization of non-coding variants can also be accomplished using Combined Annotation Dependent Depletion (CADD; [163]), which uses known and simulated variants to compute a C-score, an ad hoc measure of how deleterious is likely to be. The suggested C-score cutoff is between 10 and 20, though it is stated that any selected cutoff value would be arbitrary (http://cadd.gs.washington.edu/info). This contrasts with information-based methods, which are based on thermodynamically-defined thresholds. To directly compare methods, CADD scores were obtained for all prioritized or flagged SNVs. Half of prioritized variants met this cutoff (C > 10), while only 28.6 % of flagged variants did the same. All prioritized nonsense variants (4/4) and 26/27 missense variants had strong C-scores. Prioritized non-coding variant categories that correlated well with CADD include natural splicing variants (4/4), UTR structure variants (1/1), and RBPs (4/6). Weakly correlated variants included those affecting SRFBPs (5/36), TFBS (2/5), and pseudoexon activating variants (0/5). Missense mutations comprised 75 % of the flagged variants with C > 10. The aforementioned flagged splicing variant ATM c.1066-6 T > G also exceeded the threshold C value (C = 11.9). Meanwhile, the flagged TP53 variant, shown by SHAPE analysis to alter UTR structure, did not (C = 5.3). Despite consistency between some variant categories, the underlying assumptions of each approach probably explain why these results differ for non-coding variants. The limited numbers validated, deleterious non-coding variants also contributes to the accuracy of these predictions [163].

### Variant verification

We verified prioritized protein-truncating (N = 7) and splicing (N = 4) variants by Sanger sequencing (Tables 2 and 4, respectively). In addition, two missense variants (BRCA2 c.7958 T > C and CHEK2 c.433C > T) were re-sequenced, since they are indicated as likely pathogenic/pathogenic in ClinVar (Additional file 12: Table S9). All protein-truncating variants were confirmed, with one exception (BRCA2 c.7558C > T, no evidence for the variant was present for either strand). Two of the mRNA splicing mutations were confirmed on both strands, while the other two were confirmed on a single strand (ATM c.6347 + 1G > T and ATM c.1066-6 T > G). Both documented pathogenic missense variants were also confirmed.

### Discussion

NGS technology offers advantages in throughput and variant detection [126], but the task of interpreting the sheer volume of variants in complete gene or genome data can be daunting. The whole genome of a Yoruban male contained approximately 4.2 million SNVs and 0.4 million structural variants [164]. The variant density in the present study (average 948 variants per patient) was 5.3-fold lower than the same regions in HapMap sample NA12878 in Illumina Platinum Genomes Project (5029 variants) [165]. The difference can be attributed primarily to the exclusion of polymorphisms in highly repetitive regions in our study.

Conventional coding sequence analysis, combined with an IT-based approach for regulatory and splicing-related variants, reduced the set to a manageable number of prioritized variants. Unification of non-coding analysis of diverse protein-nucleic acid interactions using the IT framework accomplishes this by applying thermodynamic-based thresholds to binding affinity changes and by selecting the most significant binding site information changes, regardless of whether the motifs of different factors overlap.

Previously, rule-based systems have been proposed for variant severity classification [166, 167]. Functional validation and risk analyses of these variants are a prerequisite for classification, but this would not be practical to accomplish without first limiting the subset of variants analyzed. With the exception of some (but not all [37]) protein truncating variants, classification is generally not achievable by sequence analysis alone. Only a minority of variants with extreme likelihoods of pathogenic or benign phenotypes are clearly delineated because only these types of variants are considered actionable [166, 167]. The proposed classification systems preferably require functional, co-segregation, and risk analyses to stratify patients. Nevertheless, the majority of variants are VUS, especially in the case of variants occurring beyond exon boundaries. Of the 5713 variants in the BIC database, the clinical significance of 4102 BRCA1 and BRCA2 variants are either unknown (1904) or pending (2198), and only 1535 have been classified as pathogenic (Class 5) [168]. Our results cannot be considered equivalent to validation, which usually include expression assays [36] or the use of RNA-seq data [169] (splicing), qRT-PCR [170] (transcription), SHAPE analysis (mRNA 2° structure) [44], or binding assays to determine functional effects of variants. Classification of VUS in BRCA1 and BRCA2 by the ENIGMA Consortium addresses mRNA splicing and missense variants. Criteria define risk based on whether the variant occurs within a protein structural domain, the impact on protein function, and the segregation pattern of variant with disease in pedigrees [171]. These guidelines cannot be fully implemented here for several reasons: a) patients were anonymized in this study, precluding segregation analysis, b) the splicing mutation guideline does not take into account predicted leaky or cryptic SS mutations, nor other non-canonical changes that have been demonstrated to alter the expression of these and numerous other genes, c) conserved domains have not been identified in regions of the proteins encoded by these genes, especially BRCA2, where many missense mutations reside, and d) the guidelines are currently silent as to the potential impact of regulatory variants affecting splicing, RNA stability, and transcriptional regulation.

While the miRNA variant prediction program mrSNP [172] was used to evaluate all of the 3’ UTR variants, 41.4 % of the variants were predicted to alter the stability of the miRNA-target mRNA duplex for at least one miRNA expressed in breast tissue. However, only 2 of these interactions could be confirmed using TarBase [173], and these variants could not be prioritized for disruption of miRNA regulation. Other post-transcriptional processes, including miRNA regulation, that were not addressed in this study, may also be amenable to such IT-based modeling. With the proposed approach, functional prediction of such variants could precede or at least inform the classification of VUS.

It is unrealistic to expect all variants to be functionally analyzed, just as it may not be feasible to assess family members for a suspected pathogenic variant detected in a proband. The prioritization procedure reduces the chance that significant variants have been overlooked. Capturing coding and non-coding regions of HBOC-related genes, combined with the framework for assessing variants, balances the need to comprehensively detect all variation in a gene panel with the goal of identifying variants likely to be phenotypically relevant.

The location of variants in relation to known protein domains was documented in this study, but was not directly incorporated into our prioritization method. The locations and impact of splicing mutations in BRCA1 and BRCA2 were mapped to the known functional domains of the encoded proteins [174]. A high concentration of variants predicted to result in splicing changes occurred in the BRCT, RING finger, and NLS domains of BRCA1. However, BRCA2 variants were generally concentrated outside of known functional domains (aside from the C-terminal domain). Because of these inconsistencies, domain-mapping was not integrated with IT based prioritization. However, where adequate information on structure-function relationships is available (eg. TP53), we suggest that such analysis be carried out subsequent to IT-based variant prioritization.

Although coding variants are typically the sole focus of a molecular diagnostic laboratory (with the exception of the canonical dinucleotide positions within SS), non-coding mutations have long been known to be disease causing [19, 36, 175–183]. In this study, variant density in non-coding regions significantly exceeded exonic variants by > 60-fold, which, in absolute terms, constituted 1.6 % of the 15,311 variants. This is comparable to whole genome sequencing studies, which typically result in 3-4 million variants per individual, with < 2 % occurring in protein coding regions [184]. IT analysis prioritized 3 natural SS, 36 SRFBS, 5 TFBS, and 6 RBBS variants and 5 predicted to create pseudoexons. Two SS variants in ATM (c.3747-1G > A and c.6347 + 1G > T) were predicted to completely abolish the natural site and cause exon skipping. A CHEK2 variant (c.320-5A > T) was predicted to result in leaky splicing.

The IT-based framework evaluates all variants on a common scale, based on bit values, the universal unit that predicts changes in binding affinity [185]. A variant can alter the strength of one or a “set” of binding sites; the magnitude and direction of these changes is used to rank their significance. The models used to derive information weight matrices take into account the frequency of all observed bases at a given position of a binding motif, making them more accurate than consensus sequence and conservation-based approaches [36].

IT has been widely used to analyze natural and cryptic SSs [36], but its use in SRFBS analysis was only introduced recently [38]. For this reason, we assigned conservative, minimum thresholds for reporting information changes. Although there are examples of disease-causing variants resulting in small changes in Ri [174, 186–192], the majority of deleterious splicing mutations that have been verified functionally, produce large information changes. Among 698 experimentally verified deleterious variants in 117 studies, only 1.96 % resulted in < 1.0 bit change [36]. For SRFBS variants, the absolute information changes for deleterious variants ranged from 0.2 to 17.1 bits (mean 4.7 ± 3.8). This first application of IT in TFBS and RBBS analysis, however, lacks a large reference set of validated mutations for the distribution of information changes associated with deleterious variants. The release of new ChIP-seq datasets will enable IT models to be derived for TFs currently unmodeled and will improve existing models [193].

Pseudoexon activation results in disease-causing mutations [194], however such consequences are not customarily screened for in mRNA splicing analysis. IT analysis was used to detect variants that predict pseudoexon formation and 5 variants were prioritized. Previously, we have predicted experimentally proven pseudoexons with IT (Ref 42: Table 2, No #2; and Ref 195: Table 2, No #7) [42, 195]. Although it was not possible to confirm prioritized variants in the current study predicted to activate pseudoexons because of their low allele frequencies, common intronic variants that were predicted to form pseudoexons were analyzed. We then searched for evidence of pseudoexon activation in mapped human EST and mRNA tracks [196] and RNA-seq data of breast normal and tumour tissue from the Cancer Genome Atlas project [15]. One of these variants (rs6005843) appeared to splice the human EST HY160109 [197] at the predicted cryptic SS and is expressed within the pseudoexon boundaries.

Variants that were common within our population sample (i.e. occurring in > 5 individuals) and/or common in the general population (>1.0 % allele frequency) reduced the list of flagged variants substantially. This is now a commonly accepted approach for reducing candidate disease variants [166], based on the principle that the disease-causing variants occur at lower population frequencies. Variants occurring in > 5 patients all either had allele frequencies above 1.0 % or, as shown previously, resulted in very small ΔRi values [198].

The genomic context of sequence changes can influence the interpretation of a particular variant [36]. For example, variants causing significant information changes may be interpreted as inconsequential if they are functionally redundant or enhancing existing binding site function (see IT-Based Variant Analysis for details). Our understanding of the roles and context of these cognate protein factors is incomplete, which affects confidence in interpretation of variants that alter binding. Also, certain factors with important roles in the regulation of these genes, but that do not bind DNA directly or in a sequence-specific manner (eg. CtBP2 [199]), could not be included. Therefore, some variants may have been incorrectly excluded.

Although individuals can be prioritized based on potential hemizygosity, this does not definitively identify deletions. Nevertheless, it should be possible to prioritize those individuals worthy of further detailed diagnostic workup. It has not escaped our attention that the weighted probabilities obtained from this analysis could be represented and formalized using the same units of Shannon information (in bits) as the other sequence changes we have described, analogous to single or multinucleotide gene variants predicted to affect nucleic acid binding sites. Full development and validation of this method is in progress.

We also identified 4 nonsense and 3 indels in this cohort. In one individual, a 19 nt BRCA1 deletion in exon 15 causes a frameshift leading to a stop codon within 14 codons downstream. This variant, rs80359876, is considered clinically relevant. Interestingly, this deletion overlaps two other published deletions in this exon (rs397509209 and rs80359884). This raises the question as to whether this region of the BRCA1 gene is a hotspot for replication errors. DNA folding analysis indicates a possible 15 nt long stem-loop spanning this interval as the most stable predicted structure (data not shown). This 15 nt structure occurs entirely within the rs80359876 and rs397509209 deletions and partially overlaps rs80359884 (13 of 15 nt of the stem loop). It is plausible that the 2° structure of this sequence predisposes to a replication error that leads to the observed deletion.

Missense coding variants were also assessed using multiple in silico tools and evaluated based on allele frequency, literature references, and gene-specific databases. Of the 27 prioritized missense variants, the previously reported CHEK2 variant c.433G > A (chr22:29121242G > A; rs137853007) stood out, as it was identified in one patient (4-3C.5-4G) and is predicted by all 4 in silico tools to have a damaging effect on protein function. Accordingly, Wu et al. (2001) demonstrated reduced in vitro kinase activity and phosphorylation by ATM kinase compared to the wild-type CHEK2 protein [200], presumably due to the variant’s occurrence within the forkhead homology-associated domain, involved in protein-phosphoprotein interactions [201]. Implicated in Li-Fraumeni syndrome, known to increase the risk of developing several types of cancer including breast [202, 203], the CHEK2: c.433G > A variant is expected to result in a misfolded protein that would be targeted for degradation via the ubiquitin-proteosome pathway [204]. Another important missense variant is c.7958 T > C (chr13:32,936,812 T > C; rs80359022; 4-4C) in exon 17 of BRCA2. Although classified as being of unknown clinical importance in both BIC and ClinVar, it has been classified as pathogenic based on posterior probability calculations [205].

It is unlikely that all prioritized variants are pathogenic in patients carrying more than one prioritized variant. Nevertheless, a polygenic model for breast cancer susceptibility, whereby multiple moderate and low-risk alleles contribute to increased risk of HBOC may also account for multiple prioritized variants [206, 207]. There was a significant fraction of patients (29.4 %) in whom no variants were prioritized. This could be due to a) the inability of the analysis to predict a variant affecting the binding sites analyzed, b) a pathogenic variant affecting a function that was not analyzed or in a gene that was not sequenced, c) a large rearrangement/deletion where both breakpoints occur beyond the captured genomic intervals (which is unlikely, as this would have been observed as an extended non-polymorphic sequence), or d) the significant family history was not due to heritable, but instead to shared environmental influences.

BRCA coding variants were found in individuals who were previously screened for lesions in these genes, suggesting this NGS protocol is a more sensitive approach for detecting coding changes. However, previous testing of a number of these patients had been predominantly based on PTT and MLPA, which have lower sensitivity for detecting mutations than sequence analysis. Nevertheless, we identified 2 BRCA1 and 2 BRCA2 variants predicted to encode prematurely truncated proteins. Fewer non-coding BRCA variants were prioritized (15.7 %) than expected by linkage analysis [49], however this presumes at least 4 affected breast cancer diagnoses per pedigree, and, in the present study, the number of affected individuals per family was not known.

Prioritization of a variant does not equate with pathogenicity. Some prioritized variants may not increase risk, but may simply modify a primary unrecognized pathogenic mutation. A patient with a known BRCA1 nonsense variant, used as a positive control, was also found to possess an additional prioritized variant in BRCA2 (missense variant chr13:32911710A > G), which was flagged by PROVEAN and SIFT as damaging, as well as flagged for changing an SRFBS for abolishing a PTB site (while simultaneously abolished an exonic hnRNPA1 site). This variant has been identified in cases of early onset prostate cancer and is considered a VUS in ClinVar [143]. Similarly, variants prioritized in multiple patients may act as risk modifiers rather than pathogenic mutations. A larger cohort of patients with known pathogenic mutations would be necessary to calculate a background/basal rate of falsely flagged variants.

Other groups have attempted to develop comprehensive approaches for variant analysis, analogous to the one proposed here [208–210]. While most employ high-throughput sequencing and classify variants, either the sequences analyzed or the types of variants assessed tend to be limited. In particular, non-coding sequences have not been sequenced or studied to the same extent, and none of these analytical approaches have adopted a common framework for mutation analysis.

Our published oligonucleotide design method [77] produced an average sequence coverage of 98.8 %. The capture reagent did not overlap conserved highly repetitive regions, but included divergent repetitive sequences. Nevertheless, neighboring probes generated reads with partial overlap of repetitive intervals. As previously reported [147], we noted that false positive variant calls within intronic and intergenic regions were the most common consequence of dephasing in low complexity, pyrimidine-enriched intervals. This was not alleviated by processing data with software programs based on different alignment or calling algorithms. Manual review of all intronic or intergenic variants became imperative. As these sequences can still affect functional binding elements detectable by IT analysis (i.e. 3’ SSs and SRFBSs), it may prove essential to adopt or develop alignment software that explicitly and correctly identifies variants in these regions [147]. Most variants were confirmed with Sanger sequencing (10/13), and those that could not be confirmed are not necessarily false positives. A recent study demonstrated that NGS can identify variants that Sanger sequencing cannot, and reproducing sequencing results by NGS may be worthwhile before eliminating such variants [211].

### Non-coding variants

Although coding variants are typically the sole focus of a molecular diagnostic laboratory (with the exception of the canonical dinucleotide positions within SS), non-coding mutations have long been known to be disease causing [19, 36, 175–183]. In this study, variant density in non-coding regions significantly exceeded exonic variants by > 60-fold, which, in absolute terms, constituted 1.6 % of the 15,311 variants. This is comparable to whole genome sequencing studies, which typically result in 3-4 million variants per individual, with < 2 % occurring in protein coding regions [184]. IT analysis prioritized 3 natural SS, 36 SRFBS, 5 TFBS, and 6 RBBS variants and 5 predicted to create pseudoexons. Two SS variants in ATM (c.3747-1G > A and c.6347 + 1G > T) were predicted to completely abolish the natural site and cause exon skipping. A CHEK2 variant (c.320-5A > T) was predicted to result in leaky splicing.

The IT-based framework evaluates all variants on a common scale, based on bit values, the universal unit that predicts changes in binding affinity [185]. A variant can alter the strength of one or a “set” of binding sites; the magnitude and direction of these changes is used to rank their significance. The models used to derive information weight matrices take into account the frequency of all observed bases at a given position of a binding motif, making them more accurate than consensus sequence and conservation-based approaches [36].

IT has been widely used to analyze natural and cryptic SSs [36], but its use in SRFBS analysis was only introduced recently [38]. For this reason, we assigned conservative, minimum thresholds for reporting information changes. Although there are examples of disease-causing variants resulting in small changes in Ri [174, 186–192], the majority of deleterious splicing mutations that have been verified functionally, produce large information changes. Among 698 experimentally verified deleterious variants in 117 studies, only 1.96 % resulted in < 1.0 bit change [36]. For SRFBS variants, the absolute information changes for deleterious variants ranged from 0.2 to 17.1 bits (mean 4.7 ± 3.8). This first application of IT in TFBS and RBBS analysis, however, lacks a large reference set of validated mutations for the distribution of information changes associated with deleterious variants. The release of new ChIP-seq datasets will enable IT models to be derived for TFs currently unmodeled and will improve existing models [193].

Pseudoexon activation results in disease-causing mutations [194], however such consequences are not customarily screened for in mRNA splicing analysis. IT analysis was used to detect variants that predict pseudoexon formation and 5 variants were prioritized. Previously, we have predicted experimentally proven pseudoexons with IT (Ref 42: Table 2, No #2; and Ref 195: Table 2, No #7) [42, 195]. Although it was not possible to confirm prioritized variants in the current study predicted to activate pseudoexons because of their low allele frequencies, common intronic variants that were predicted to form pseudoexons were analyzed. We then searched for evidence of pseudoexon activation in mapped human EST and mRNA tracks [196] and RNA-seq data of breast normal and tumour tissue from the Cancer Genome Atlas project [15]. One of these variants (rs6005843) appeared to splice the human EST HY160109 [197] at the predicted cryptic SS and is expressed within the pseudoexon boundaries.

Variants that were common within our population sample (i.e. occurring in > 5 individuals) and/or common in the general population (>1.0 % allele frequency) reduced the list of flagged variants substantially. This is now a commonly accepted approach for reducing candidate disease variants [166], based on the principle that the disease-causing variants occur at lower population frequencies. Variants occurring in > 5 patients all either had allele frequencies above 1.0 % or, as shown previously, resulted in very small ΔRi values [198].

The genomic context of sequence changes can influence the interpretation of a particular variant [36]. For example, variants causing significant information changes may be interpreted as inconsequential if they are functionally redundant or enhancing existing binding site function (see IT-Based Variant Analysis for details). Our understanding of the roles and context of these cognate protein factors is incomplete, which affects confidence in interpretation of variants that alter binding. Also, certain factors with important roles in the regulation of these genes, but that do not bind DNA directly or in a sequence-specific manner (eg. CtBP2 [199]), could not be included. Therefore, some variants may have been incorrectly excluded.

### Prioritization of potential deletions

Although individuals can be prioritized based on potential hemizygosity, this does not definitively identify deletions. Nevertheless, it should be possible to prioritize those individuals worthy of further detailed diagnostic workup. It has not escaped our attention that the weighted probabilities obtained from this analysis could be represented and formalized using the same units of Shannon information (in bits) as the other sequence changes we have described, analogous to single or multinucleotide gene variants predicted to affect nucleic acid binding sites. Full development and validation of this method is in progress.

### Coding sequence changes

We also identified 4 nonsense and 3 indels in this cohort. In one individual, a 19 nt BRCA1 deletion in exon 15 causes a frameshift leading to a stop codon within 14 codons downstream. This variant, rs80359876, is considered clinically relevant. Interestingly, this deletion overlaps two other published deletions in this exon (rs397509209 and rs80359884). This raises the question as to whether this region of the BRCA1 gene is a hotspot for replication errors. DNA folding analysis indicates a possible 15 nt long stem-loop spanning this interval as the most stable predicted structure (data not shown). This 15 nt structure occurs entirely within the rs80359876 and rs397509209 deletions and partially overlaps rs80359884 (13 of 15 nt of the stem loop). It is plausible that the 2° structure of this sequence predisposes to a replication error that leads to the observed deletion.

Missense coding variants were also assessed using multiple in silico tools and evaluated based on allele frequency, literature references, and gene-specific databases. Of the 27 prioritized missense variants, the previously reported CHEK2 variant c.433G > A (chr22:29121242G > A; rs137853007) stood out, as it was identified in one patient (4-3C.5-4G) and is predicted by all 4 in silico tools to have a damaging effect on protein function. Accordingly, Wu et al. (2001) demonstrated reduced in vitro kinase activity and phosphorylation by ATM kinase compared to the wild-type CHEK2 protein [200], presumably due to the variant’s occurrence within the forkhead homology-associated domain, involved in protein-phosphoprotein interactions [201]. Implicated in Li-Fraumeni syndrome, known to increase the risk of developing several types of cancer including breast [202, 203], the CHEK2: c.433G > A variant is expected to result in a misfolded protein that would be targeted for degradation via the ubiquitin-proteosome pathway [204]. Another important missense variant is c.7958 T > C (chr13:32,936,812 T > C; rs80359022; 4-4C) in exon 17 of BRCA2. Although classified as being of unknown clinical importance in both BIC and ClinVar, it has been classified as pathogenic based on posterior probability calculations [205].

It is unlikely that all prioritized variants are pathogenic in patients carrying more than one prioritized variant. Nevertheless, a polygenic model for breast cancer susceptibility, whereby multiple moderate and low-risk alleles contribute to increased risk of HBOC may also account for multiple prioritized variants [206, 207]. There was a significant fraction of patients (29.4 %) in whom no variants were prioritized. This could be due to a) the inability of the analysis to predict a variant affecting the binding sites analyzed, b) a pathogenic variant affecting a function that was not analyzed or in a gene that was not sequenced, c) a large rearrangement/deletion where both breakpoints occur beyond the captured genomic intervals (which is unlikely, as this would have been observed as an extended non-polymorphic sequence), or d) the significant family history was not due to heritable, but instead to shared environmental influences.

BRCA coding variants were found in individuals who were previously screened for lesions in these genes, suggesting this NGS protocol is a more sensitive approach for detecting coding changes. However, previous testing of a number of these patients had been predominantly based on PTT and MLPA, which have lower sensitivity for detecting mutations than sequence analysis. Nevertheless, we identified 2 BRCA1 and 2 BRCA2 variants predicted to encode prematurely truncated proteins. Fewer non-coding BRCA variants were prioritized (15.7 %) than expected by linkage analysis [49], however this presumes at least 4 affected breast cancer diagnoses per pedigree, and, in the present study, the number of affected individuals per family was not known.

Prioritization of a variant does not equate with pathogenicity. Some prioritized variants may not increase risk, but may simply modify a primary unrecognized pathogenic mutation. A patient with a known BRCA1 nonsense variant, used as a positive control, was also found to possess an additional prioritized variant in BRCA2 (missense variant chr13:32911710A > G), which was flagged by PROVEAN and SIFT as damaging, as well as flagged for changing an SRFBS for abolishing a PTB site (while simultaneously abolished an exonic hnRNPA1 site). This variant has been identified in cases of early onset prostate cancer and is considered a VUS in ClinVar [143]. Similarly, variants prioritized in multiple patients may act as risk modifiers rather than pathogenic mutations. A larger cohort of patients with known pathogenic mutations would be necessary to calculate a background/basal rate of falsely flagged variants.

Other groups have attempted to develop comprehensive approaches for variant analysis, analogous to the one proposed here [208–210]. While most employ high-throughput sequencing and classify variants, either the sequences analyzed or the types of variants assessed tend to be limited. In particular, non-coding sequences have not been sequenced or studied to the same extent, and none of these analytical approaches have adopted a common framework for mutation analysis.

Our published oligonucleotide design method [77] produced an average sequence coverage of 98.8 %. The capture reagent did not overlap conserved highly repetitive regions, but included divergent repetitive sequences. Nevertheless, neighboring probes generated reads with partial overlap of repetitive intervals. As previously reported [147], we noted that false positive variant calls within intronic and intergenic regions were the most common consequence of dephasing in low complexity, pyrimidine-enriched intervals. This was not alleviated by processing data with software programs based on different alignment or calling algorithms. Manual review of all intronic or intergenic variants became imperative. As these sequences can still affect functional binding elements detectable by IT analysis (i.e. 3’ SSs and SRFBSs), it may prove essential to adopt or develop alignment software that explicitly and correctly identifies variants in these regions [147]. Most variants were confirmed with Sanger sequencing (10/13), and those that could not be confirmed are not necessarily false positives. A recent study demonstrated that NGS can identify variants that Sanger sequencing cannot, and reproducing sequencing results by NGS may be worthwhile before eliminating such variants [211].

### Conclusions

Through a comprehensive protocol based on high-throughput, IT-based and complementary coding sequence analyses, the numbers of VUS can be reduced to a manageable quantity of variants, prioritized by predicted function. While exonic variants corresponded to a small fraction of prioritized variants, there is considerably more evidence for their pathogenicity because clinical sequencing has concentrated in these regions. Our sequencing approach illustrates the importance of sequencing non-coding regions of genes to establish pathogenic mutations not already evident from changes in the amino acid based genetic code. We suggest our approach for variant flagging and prioritization bridges the phase between high-throughput sequencing, variant detection with the time-consuming process of variant classification, including pedigree analysis and functional validation. Subsequent to completion of the present study, ethics approval was obtained for a similar analysis of consented patients with clinical information. This work has since been described elsewhere [212].

Variants will be deposited with the ENIGMA Consortium (www.enigmaconsortium.org), which is a designated organization for curation of HBOC mutations and which is charged with protection of genetic privacy of participants. Additionally, all likely pathogenic variants were submitted to ClinVar (submission ID: SUB1332053) while other novel variants were submitted to dbSNP (ss1966658584-1966658622).

### Availability of supporting data

Variants will be deposited with the ENIGMA Consortium (www.enigmaconsortium.org), which is a designated organization for curation of HBOC mutations and which is charged with protection of genetic privacy of participants. Additionally, all likely pathogenic variants were submitted to ClinVar (submission ID: SUB1332053) while other novel variants were submitted to dbSNP (ss1966658584-1966658622).



# SUPPLEMENTAL FILE 1: 12920_2016_Article_178.pdf

# Preparing to download ...

[HHS Vulnerability Disclosure](https://www.hhs.gov/vulnerability-disclosure-policy/index.html)