# MAIN TEXT

## Exon splicing analysis of intronic variants in multigene cancer panel testing for hereditary breast/ovarian cancer

### Abstract

AbstractThe use of multigene panel testing for patients with a predisposition to breast/ovarian cancer is increasing as the identification of variants is useful for diagnosis and disease management. We identified pathogenic and likely pathogenic (P/LP) variants of high‐and moderate‐risk genes using a 23‐gene germline cancer panel in 518 patients with hereditary breast and ovarian cancers (HBOC). The frequency of P/LP variants was 12.4% (64/518) for high‐ and moderate‐penetrant genes, namely, BRCA2 (5.6%), BRCA1 (3.3%), CHEK2 (1.2%), MUTYH (0.8%), PALB2 (0.8%), MLH1 (0.4%), ATM (0.4%), BRIP1 (0.4%), TP53 (0.2%), and PMS2 (0.2%). Five patients possessed two P/LP variants in BRCA1/2 and other genes. We also compared the results from in silico splicing predictive tools and exon splicing patterns from patient samples by analyzing RT‐PCR product sequences in six P/LP intronic variants and two intronic variants of unknown significance (VUS). Altered transcriptional fragments were detected for P/LP intronic variants in BRCA1, BRIP1, CHEK2, PARB2, and PMS2. Notably, we identified an in‐frame deletion of the BRCA1 C‐terminal (BRCT) domain by exon skipping in BRCA1 c.5152+6T>C—as known VUS—indicating a risk for HBOC. Thus, exon splicing analysis can improve the identification of veiled intronic variants that would aid decision making and determination of hereditary cancer risk.

### INTRODUCTION

Breast cancer is a multifactorial disease caused by a combination of environmental and genetic factors.
1
, 
2
 Hereditary breast and ovarian cancers (HBOCs) account for approximately 5‐10% of breast cancer and 10%‐15% of ovarian cancer and primarily involve BRCA1 and BRCA2 variants.
2
, 
3
, 
4
 The identification of BRCA1 and BRCA2 germline variants significantly improve HBOC diagnosis, as they are predictors of cancer susceptibility for patients as well as their families.
5
, 
6

Recent advances in genetic sequencing technology have led to the discovery of novel genes that increase the risk of cancer in patients with familial predisposition.
5
, 
6
, 
7
, 
8
 However, the rapid introduction of multigene panel testing has raised several issues to be addressed for implementation in clinical settings.
9
 First, many of the tested genes are low‐ to moderate‐risk genes for which consensus management guidelines have not been established.
9
, 
10
 In the absence of identified variants, recommendations for cancer‐specific screening and prevention approaches for patients and family members are typically based on personal and/or family cancer history.
11
 Second, it is uncertain whether identifying such low‐ to moderate‐risk gene variants would influence the individual clinical management of patients referred for genetic testing.
11
 Although several studies have identified variants in moderate‐risk genes, such as ATM, BRIP1, CHEK2, BARD1, MRE11A, NBN, RAD50, RAD51, and XRCC2, as well as in high‐penetrant genes, including BRCA1/2, TP53, PTEN, STK11, CDH1, and PALB2,
12
, 
13
, 
14
 establishing clinical relevance and analyzing these variants across diverse ethnic populations is warranted.

Correct exon splicing is important for appropriate protein translation as alterations in this process can lead to aberrant cellular metabolism or functions. Abnormal splicing caused by mutation events may alter consensus splicing regulator sequences, leading to hereditary disorders.
15
 Although in silico bioinformatics algorithms were developed for evaluating the possible exon splicing effects of identified variants, the exact effects of variants should be demonstrated in functional assays.

In this study, we employed a comprehensive multigene panel that included 23 known or suspected cancer susceptibility genes to test Korean patients suspected of HBOC. We aimed to identify possible pathogenic or likely pathogenic (P/LP) variants as well as variants of unknown significance (VUS) for various genes including BRCA1. We also analyzed exon splicing patterns in intronic variants to evaluate their deleterious effects.

### MATERIAL AND METHODS

A total of 700 patients who were suspected of a familial predisposition to cancer were referred to a genetic counseling clinic in the Korea National Cancer Center and underwent BRCA1/2 testing between January 1, 2017 and December 31, 2018. Suspected clinical characteristics of HBOC were defined as follows: (a) at least one case of breast or ovarian cancer with a family history of breast and ovarian cancer; (b) first diagnosis of breast cancer onset ≤ 40 years old; (c) bilateral breast cancer or other primary cancer with other primary malignancy; and (d) simultaneous diagnosis of breast and ovarian cancers. All were in accordance with the criteria of HBOC testing according to the NCCN guidelines on genetic/familial high‐risk assessment: breast and ovarian (version 2, 2017).
13
 Of these, 518 patients who had agreed to the study underwent further evaluation with a customized 23‐gene hereditary cancer panel (Figure 1). Data on demographics, personal and familial history of cancer, and panel testing results were retrospectively collected for patients harboring P/LP variants. Clinicopathological characteristics of cancer, such as the stage and presence of the hormone receptor (HR) and human epidermal growth factor receptor 2 (HER2) states, were assessed by reviewing medical records. As a control group, 393 healthy female controls were recruited among individuals who had visited the National Cancer Center as part of a cancer‐screening program.

Germline variants in patients with hereditary breast and ovarian cancer detected using a high‐ and moderate‐penetrance hereditary cancer gene panel of 23 genes. Schematic representation of the patients and study workflow. A total of 700 breast/ovarian cancer patients visited the genetic counseling clinic between January 2017 and December 2018 at the National Cancer Center (Republic of Korea) and underwent BRCA1/2 testing. Of these, 518 patients were enrolled in the study and tested using a customized 23‐gene hereditary cancer panel. The frequency of pathogenic (P) or likely pathogenic (LP) variants was 12.4% (64/518)

BRCA1/2 genetic testing was performed by the Green Cross Company (Yongin, Republic of Korea) via direct sequencing. Briefly, genomic DNA was extracted from peripheral blood samples with a QIAamp DNA Blood Mini Kit (Qiagen) or Chemagic DNA Blood 200 Kit (Chemagen). Amplified products were sequenced on an ABI 3500xl Analyzer (Applied Biosystems) using the Bigdye Terminator v3.1 Cycle Sequencing Kit. Sequences were analyzed using Sequencher v5.0 software (Gene Codes). All variants are described according to HUGO‐approved systematic nomenclature (http://www.hgvs.org/mutnomen/).

Genomic DNA was extracted from peripheral blood of each patient. We employed a customized hereditary cancer panel (Celemics) that included all coding sequences and intron‐exon boundaries of the coding exon from 23 cancer predisposition genes (APC, ATM, BARD1, BRCA1, BRCA2, BRIP1, CDH1, CHEK2, EPCAM, MLH1, MSH2, MSH6, MUTYH, MEN1, NBN, PALB2, PMS2, PTEN, RAD50, RAD51C, RET, STK11, and TP53; Table 2). Products with each capture reaction were sequenced on an Illumina MiSeqDX (Illumina Inc) generating 2 × 150 bp paired‐end reads. Alignment of sequence reads, indexing of the reference genome (hg19), and variant calling were performed with a pipeline based on Genome Analysis Tool Kit (GATK) Best Practices.
16
 Alignment was performed with BWA‐mem (version 0.7.10),
17
 duplicated reads were marked with Picard (version 1.138; http://picard.sourceforge.net), and local alignment, base quality recalibration, and variant calling were performed using the GATK (version 3.5),
18
 samtools (version 0.1.19),
19
 FreeBayes (v0.9.21‐26‐gbfd9832), and Scalpel (version 0.5.3).
20

Variant annotation was performed with VEP (Ensembl Variant Effect Predictor)
21
 and dbNSFP v 3.0.
22
 We obtained all single‐base pair substitutions, insertions and/or deletions for each gene. Genetic variants were classified using a five‐tier system following the American College of Medical Genetics and Genomics (ACMG) guidelines as follows
23
: P, LP, VUS, likely benign, or benign.

We used the following five splice site prediction programs to predict the effect of intronic variants on the efficiency of splicing: Splice Site Finder (http://www.interactive‐biosoftware.com), GeneSplicer (http://www.cbcb.umd.edu/software/GeneSplicer), Splice Site Prediction by Neural Network (http://www.fruitfly.org/seq_tools/splice.html), MaxEntScan (http://genes.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html), and Human Splicing Finder (http://www.umd.be/HSF/). Analysis was conducted by the integrated software Alamut Visual (version 2.12; http://www.interactive‐biosoftware.com) using default settings for all predictions. A variation of more than 10% in at least two algorithms was considered to have an effect on splicing. All intronic variants classified as P/LP were analyzed and the two intronic VUS of BRCA1/2 were evaluated using the available samples.

We analyzed RNA transcripts from patient samples to compare outcomes via in silico splicing analysis. Total RNA was extracted from peripheral blood lymphocytes using the NucleoSpin RNA Blood Kit (Macherey‐Nagel) or from normal tissues using the AllPrep DNA/RNA Mini Kit (Qiagen) according to manufacturer's instructions. Using total RNA (1 µg) as template, cDNA was reverse‐transcribed using the Transcriptor First Strand cDNA Synthesis Kit (Roche Life Science) or ReverTra Ace qPCR RT Master Mix (Toyobo), followed by amplification via reverse transcription PCR (RT‐PCR) as previously reported.
24
, 
25
, 
26
, 
27
 Transcriptional products for intronic variants in BRCA1, BRCA2, BRIP1, CHEK2, PALB2, and PMS2 were obtained and validated by Sanger sequencing. Primer sequences are listed in Table S1.

Intronic variants were further genotyped in healthy female controls. Variants were identified by TaqMan SNP Genotyping Assays (Applied Biosystems) using the QuantStudio 7 Flex Real‐Time PCR System (Applied Biosystems).

### Study population

A total of 700 patients who were suspected of a familial predisposition to cancer were referred to a genetic counseling clinic in the Korea National Cancer Center and underwent BRCA1/2 testing between January 1, 2017 and December 31, 2018. Suspected clinical characteristics of HBOC were defined as follows: (a) at least one case of breast or ovarian cancer with a family history of breast and ovarian cancer; (b) first diagnosis of breast cancer onset ≤ 40 years old; (c) bilateral breast cancer or other primary cancer with other primary malignancy; and (d) simultaneous diagnosis of breast and ovarian cancers. All were in accordance with the criteria of HBOC testing according to the NCCN guidelines on genetic/familial high‐risk assessment: breast and ovarian (version 2, 2017).
13
 Of these, 518 patients who had agreed to the study underwent further evaluation with a customized 23‐gene hereditary cancer panel (Figure 1). Data on demographics, personal and familial history of cancer, and panel testing results were retrospectively collected for patients harboring P/LP variants. Clinicopathological characteristics of cancer, such as the stage and presence of the hormone receptor (HR) and human epidermal growth factor receptor 2 (HER2) states, were assessed by reviewing medical records. As a control group, 393 healthy female controls were recruited among individuals who had visited the National Cancer Center as part of a cancer‐screening program.

Germline variants in patients with hereditary breast and ovarian cancer detected using a high‐ and moderate‐penetrance hereditary cancer gene panel of 23 genes. Schematic representation of the patients and study workflow. A total of 700 breast/ovarian cancer patients visited the genetic counseling clinic between January 2017 and December 2018 at the National Cancer Center (Republic of Korea) and underwent BRCA1/2 testing. Of these, 518 patients were enrolled in the study and tested using a customized 23‐gene hereditary cancer panel. The frequency of pathogenic (P) or likely pathogenic (LP) variants was 12.4% (64/518)

### BRCA1 and BRCA2 direct sequencing

BRCA1/2 genetic testing was performed by the Green Cross Company (Yongin, Republic of Korea) via direct sequencing. Briefly, genomic DNA was extracted from peripheral blood samples with a QIAamp DNA Blood Mini Kit (Qiagen) or Chemagic DNA Blood 200 Kit (Chemagen). Amplified products were sequenced on an ABI 3500xl Analyzer (Applied Biosystems) using the Bigdye Terminator v3.1 Cycle Sequencing Kit. Sequences were analyzed using Sequencher v5.0 software (Gene Codes). All variants are described according to HUGO‐approved systematic nomenclature (http://www.hgvs.org/mutnomen/).

### Panel‐based sequencing assay

Genomic DNA was extracted from peripheral blood of each patient. We employed a customized hereditary cancer panel (Celemics) that included all coding sequences and intron‐exon boundaries of the coding exon from 23 cancer predisposition genes (APC, ATM, BARD1, BRCA1, BRCA2, BRIP1, CDH1, CHEK2, EPCAM, MLH1, MSH2, MSH6, MUTYH, MEN1, NBN, PALB2, PMS2, PTEN, RAD50, RAD51C, RET, STK11, and TP53; Table 2). Products with each capture reaction were sequenced on an Illumina MiSeqDX (Illumina Inc) generating 2 × 150 bp paired‐end reads. Alignment of sequence reads, indexing of the reference genome (hg19), and variant calling were performed with a pipeline based on Genome Analysis Tool Kit (GATK) Best Practices.
16
 Alignment was performed with BWA‐mem (version 0.7.10),
17
 duplicated reads were marked with Picard (version 1.138; http://picard.sourceforge.net), and local alignment, base quality recalibration, and variant calling were performed using the GATK (version 3.5),
18
 samtools (version 0.1.19),
19
 FreeBayes (v0.9.21‐26‐gbfd9832), and Scalpel (version 0.5.3).
20

### Variant classification

Variant annotation was performed with VEP (Ensembl Variant Effect Predictor)
21
 and dbNSFP v 3.0.
22
 We obtained all single‐base pair substitutions, insertions and/or deletions for each gene. Genetic variants were classified using a five‐tier system following the American College of Medical Genetics and Genomics (ACMG) guidelines as follows
23
: P, LP, VUS, likely benign, or benign.

### In silico exon splicing analysis of intronic variants

We used the following five splice site prediction programs to predict the effect of intronic variants on the efficiency of splicing: Splice Site Finder (http://www.interactive‐biosoftware.com), GeneSplicer (http://www.cbcb.umd.edu/software/GeneSplicer), Splice Site Prediction by Neural Network (http://www.fruitfly.org/seq_tools/splice.html), MaxEntScan (http://genes.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html), and Human Splicing Finder (http://www.umd.be/HSF/). Analysis was conducted by the integrated software Alamut Visual (version 2.12; http://www.interactive‐biosoftware.com) using default settings for all predictions. A variation of more than 10% in at least two algorithms was considered to have an effect on splicing. All intronic variants classified as P/LP were analyzed and the two intronic VUS of BRCA1/2 were evaluated using the available samples.

### Functional analysis and genotyping of splice acceptor variants

We analyzed RNA transcripts from patient samples to compare outcomes via in silico splicing analysis. Total RNA was extracted from peripheral blood lymphocytes using the NucleoSpin RNA Blood Kit (Macherey‐Nagel) or from normal tissues using the AllPrep DNA/RNA Mini Kit (Qiagen) according to manufacturer's instructions. Using total RNA (1 µg) as template, cDNA was reverse‐transcribed using the Transcriptor First Strand cDNA Synthesis Kit (Roche Life Science) or ReverTra Ace qPCR RT Master Mix (Toyobo), followed by amplification via reverse transcription PCR (RT‐PCR) as previously reported.
24
, 
25
, 
26
, 
27
 Transcriptional products for intronic variants in BRCA1, BRCA2, BRIP1, CHEK2, PALB2, and PMS2 were obtained and validated by Sanger sequencing. Primer sequences are listed in Table S1.

Intronic variants were further genotyped in healthy female controls. Variants were identified by TaqMan SNP Genotyping Assays (Applied Biosystems) using the QuantStudio 7 Flex Real‐Time PCR System (Applied Biosystems).

### RESULTS

Clinical characteristics of patients with breast/ovarian cancer subjected to multigene panel testing are listed in Table 1. A total of 507 patients (507/518; 97.9%) were diagnosed with breast cancer, whereas eight (8/518; 1.5%) were diagnosed with ovarian cancer; three patients were diagnosed with both cancers. Among patients with breast cancer, 35.8% were diagnosed before 40 years of age. Histologically, invasive ductal carcinoma (IDC) was predominant, and 51.8% of patients were categorized as stage I. Initial diagnosis of 54.5% of patients with ovarian cancer revealed that they were in their 50s; notably, these patients had a familial predisposition to cancer.

Demographic characteristics of patients with breast/ovarian cancer

Abbreviations: HER2, human epidermal growth factor receptor 2; HR, hormone receptors (estrogen receptor, progesterone receptor); TNBC, triple‐negative breast cancer.

P/LP variants were detected in 12.4% (64/518) of the 518 patients with breast and ovarian cancer tested with the multigene panel (Figure 1). As shown in Table 2, we divided the genes included in the panel into three groups based on hereditary penetrance. Of the P/LP variants, 10.0% were in group A, 1.5% in group B, and 1.7% in group C. In our cohort, group A variants were detected in BRCA2 (5.6%), BRCA1 (3.3%), PALB2 (0.8%), TP53 (0.2%), and PTEN (0.2%) genes. Group B variants were found in CHEK2 (1.2%) and ATM (0.4%), whereas group C variants were in MUTYH (0.8%), MLH1 (0.4%), BRIP1 (0.4%), and PMS2 (0.2%) (Table 2).

Pathogenic or likely pathogenic variants distributed according to hereditary penetrance

Of the 19 patients who harbored cancer‐predisposing germline variants except for BRCA1/2, 15 patients had a family history of cancer (Table 3). These patients were mainly diagnosed with IDC tumor type, and some also had bilateral breast cancer or multiple cancers (Table 3). The median age of diagnosis in these patients was 48.3 ± 11.8 years (ranging from 33 to 72 years; Table 2). Five patients possessed two P/LP variants in BRCA1/2 and other genes in the panel (Table S2). Variants in genes other than BRCA1/2 were as follows: one BRIP1 c.1794+1 G>A, two CHEK2 c.1555 C>T (p.Arg159*), and one MUTYH c.544 C>T (p.Arg182Cys); one patient had two pathogenic variants in the BRCA2. These patients were diagnosed with IDC tumor type HR+ and HER2‐ and had at least one clinical feature of suspected HBOC, ie, disease onset at less than 40 years of age, bilateral cancer, or family history of cancer.

Characteristics of BRCA 1/2‐variant‐negative patients with pathogenic/likely pathogenic variants in other cancer‐associated genes

Abbreviations: ACMG, American College of Medical Genetics and Genomics; DCIS, ductal carcinoma in situ; FDR, first‐degree relative; HER2, human epidermal growth factor receptor 2; HR, hormone receptors (estrogen receptor, progesterone receptor); IDC, invasive ductal carcinoma; N/A, not available; SDR, second‐degree relative; TNBC, triple‐negative breast cancer.

Reported by Nakagomi H et al
26

Via multigene panel analysis, we identified six P/LP intronic variants: BRCA1 c.302‐2A>C, BRCA1 c.5277+1G>A, CHEK2 c.846+1G>T, PALB2 c.2834+2T>C, BRIP1 c.1794+1 G>A, and PMS2 c.164‐1G>A (Table 4). To predict the exon splicing patterns of these intronic variants, we employed five in silico splice site prediction programs. The natural splicing sites of the variants were predicted to be affected (Table 4). To evaluate the predicted exon splicing effects by in silico programs, we analyzed mRNA transcripts from patient samples via RT‐PCR (Table 4 and Figure S1‐S6).

In silico exon splicing analysis and RT‐PCR results of intronic variants detected in patients with hereditary breast/ovarian cancer

Splice Site Finder

(0‐100)

Max ent scan

(0‐16)

NN

SPICE (0‐1)

Gene splicer

(0‐15)

Classification

(ACMG guidelines)

db

SNP

147

gnom

AD

_exome

_ALL

KRGDB_

1100

Pathogenic

(PVS1, PM2, PP5)

Abbreviations: A, acceptor site; ACMG, American College of Medical Genetics and Genomics; D, donor site; ExAC, Exome Aggregation Consortium; gnomAD, The Genome Aggregation Database; MUT, mutation; NE, splice site not evaluated by the algorithm; WT, wild type; 1000G, 1000 Genomes Project; ㅡ, splice site not detected.

CHEK2 c.846+1G>T was detected in patient PT13, who was diagnosed with HR‐, HER2/neu+ bilateral breast cancer (IDC); the proband's sister died from breast cancer at age 51 after being diagnosed at 47, but other family members were cancer free (Table 3 and Figure S1). CHEK2 c.846+1G>T was predicted to affect the donor splice site (Table 4). We revealed aberrant mRNA transcripts collected from the patient's lymphocytes that corresponded to the skipping of exon 7. This variant was predicted to cause an in‐frame deletion (266‐284 amino acids) in the kinase domain of CHEK2 protein (Figure S1).

PALB2 c.2834+2T>C was detected in patient PT39, who was diagnosed with HR+, HER2/neu‐ IDC at age 35; the proband had a strong family history of reproductive cancer with one paternal aunt being diagnosed with breast cancer at 45 years of age, and another diagnosed with lung cancer at 60 years of age (Table 3 and Figure S2). The proband's paternal grandfather was diagnosed with lung cancer in his late 70s, while the rest of the known family was cancer free. Via RNA analysis, PALB2 c.2834+2T>C was revealed to be the combined product of exon 7 and 9 by exon 8 skipping in PALB2 (Figure S2). This variant was predicted to truncate the PALB2 protein. Generation of truncated protein products by alternative splicing may be associated with an increased risk of breast cancer.
28

PT1 (a 50‐year‐old woman) was referred for genetic counseling after a diagnosis of bilateral IDC and was found to be HR+ and HER2/neu‐. Her family history was only notable with regards to her uncle, who was diagnosed with gastric cancer. Subsequent multigene panel testing revealed the presence of pathogenic BRCA1 c.923_924del, and BRIP1 c.1794+1G>A variants (Table S2 and Figure S3). BRIP1 c.1794+1G>A was predicted to affect the donor splice site, as evidenced by in silico analysis (Table 4). Functional analysis of the intronic variant identified a deletion within exon 12 in the BRIP1 gene, which caused abnormal transcriptional production (Figure S3). The proband's cancer‐related family history and our analysis suggest that this variant confers an increased risk of breast cancer.

PMS2 c.164‐1G>C was detected in patient PT43, who was diagnosed with breast and gastric cancers before age 50. Her father was diagnosed with colon cancer at 73 years old, while her uncle died of gastric cancer at age 60. The rest of the family was healthy (Table 3 and Figure S4). In silico analysis predicted that this variant may be problematic at the acceptor splice site or is activated at the cryptic site (Table 4). Splicing functional assays revealed that this variant induces aberrant splicing via partial exon 3 deletion (8 bp) and leads to the subsequent production of a truncated protein, as evidenced by RT‐PCR (Figure S4).

We identified 14 BRCA1 and 24 BRCA2 P/LP variants in our cohort (Table S3). Among these variants, we analyzed exon splicing in two BRCA1 intronic variants classified as P/LP BRCA1 c.302‐2A>C and c.5277+1G>A, which were predicted to affect the donor or acceptor splice sites (Table 4). RT‐PCR identified partial exon 7 deletion (10 bp), and this BRCA1 c.302‐2A>C variant was predicted to truncate BRCA1 protein (Table 4 and Figure S5). The BRCA1 c.302‐2A>C variant was detected in patient PT49, who was diagnosed with bilateral breast cancer. The patient's father died from lung cancer at 61 years old, and her sister was also diagnosed with bilateral breast cancer (Figure S5). BRCA1 c.5277+1G>A was identified in patient PT50, who was diagnosed with ovarian cancer at the age of 44; her family is cancer free. Two abnormal mRNA transcripts were detected in the patient's lymphocytes with BRCA1 c.5277+1G>A that were identified to have an 87‐bp insertion of intron 20 and exon 20 skipping, as seen in a previous study using mini‐gene splicing assays or blood samples (Figure S6).
26

We also examined the effect of one BRCA1 and one BRCA2 VUS on splicing (Table 4, Figure 2, and Figure S7). BRCA1 c.5152+6T>C and BRCA2 c.317‐10A>G are classified with uncertain significance (Table 4). BRCA1 c.5152+6T>C was expected to affect exon splicing, whereas no splicing alterations were predicted for BRCA2 c.317‐10A>G (Table 4). BRCA1 c.5152+6T>C mRNA transcripts were abnormal compared with their corresponding wild‐type transcripts; this variant produced a combined exon 17 and 19 by exon 18 skipping in BRCA1 and was predicted to be an in‐frame deletion of the BRCA1 C‐terminal (BRCT) domain of BRCA1 (Figure 2). However, mRNA transcript indicated that BRCA2 c.317‐10A>G had no effect on exon splicing (Figure S7). Although in silico predictions and RNA analysis revealed the pathogenicity of VUS variants, additional analysis is required to classify them as pathogenic.

Exon splicing analysis of the BRCA1 c.5152+6T>C variant of patient PT51. A, Schematic view of variant c.5152+6T>C localization in the BRCA1 gene. PCR primer alignment is indicated with the red and blue bars. Sequencing analysis for genomic DNA is presented below. B, RT‐PCR of lymphocyte‐derived RNA. Predicted scheme of mRNA transcript in control or patient samples (upper right panel). Agarose gel (2%) electrophoresis; lane 1: control sample; lane 2: patient sample. Two PCR products were detected in the patient sample (upper middle panel). Chromatogram sequences of the control and abnormal transcripts. Vertical line in the chromatogram indicates the exonic junction in transcripts. Exon 18 (78 bp) skipping between exon 17 and exon 19 was identified (upper left panel). Functional domains of BRCA1 and sequence alignment of the BRCA1 abnormal transcript (lower panel). Amino acid sequences of the splice variant (c.5152+6T>C) were aligned using a reference sequence (NP_009225.1) via NCBI BLAST (https://blast.ncbi.nlm.nih.gov/Blast.cgi). BRCA1 c.5152+6T>C was identified to encode a BRCA1 protein with an in‐frame deletion (26 amino acids) in the BRCA1 C‐terminal (BRCT) domain; this may affect the function of the BRCA1 BRCT domain. The red line indicates the location of the in‐frame deletion residues. C, Pedigree of patient PT51

We conducted genotyping of 393 healthy female Korean controls to further define P/LP or VUS intronic variants by comparing their prevalence. The eight intronic variants, analyzed in exon splicing assays, were not detected in healthy controls (Table 4). These results suggest that these intronic variants may affect the susceptibility to inherit breast and ovarian cancer.

### Baseline characteristics of the study population

Clinical characteristics of patients with breast/ovarian cancer subjected to multigene panel testing are listed in Table 1. A total of 507 patients (507/518; 97.9%) were diagnosed with breast cancer, whereas eight (8/518; 1.5%) were diagnosed with ovarian cancer; three patients were diagnosed with both cancers. Among patients with breast cancer, 35.8% were diagnosed before 40 years of age. Histologically, invasive ductal carcinoma (IDC) was predominant, and 51.8% of patients were categorized as stage I. Initial diagnosis of 54.5% of patients with ovarian cancer revealed that they were in their 50s; notably, these patients had a familial predisposition to cancer.

Demographic characteristics of patients with breast/ovarian cancer

Abbreviations: HER2, human epidermal growth factor receptor 2; HR, hormone receptors (estrogen receptor, progesterone receptor); TNBC, triple‐negative breast cancer.

### Cancer‐predisposing germline variants including the BRCA1/2

P/LP variants were detected in 12.4% (64/518) of the 518 patients with breast and ovarian cancer tested with the multigene panel (Figure 1). As shown in Table 2, we divided the genes included in the panel into three groups based on hereditary penetrance. Of the P/LP variants, 10.0% were in group A, 1.5% in group B, and 1.7% in group C. In our cohort, group A variants were detected in BRCA2 (5.6%), BRCA1 (3.3%), PALB2 (0.8%), TP53 (0.2%), and PTEN (0.2%) genes. Group B variants were found in CHEK2 (1.2%) and ATM (0.4%), whereas group C variants were in MUTYH (0.8%), MLH1 (0.4%), BRIP1 (0.4%), and PMS2 (0.2%) (Table 2).

Pathogenic or likely pathogenic variants distributed according to hereditary penetrance

Of the 19 patients who harbored cancer‐predisposing germline variants except for BRCA1/2, 15 patients had a family history of cancer (Table 3). These patients were mainly diagnosed with IDC tumor type, and some also had bilateral breast cancer or multiple cancers (Table 3). The median age of diagnosis in these patients was 48.3 ± 11.8 years (ranging from 33 to 72 years; Table 2). Five patients possessed two P/LP variants in BRCA1/2 and other genes in the panel (Table S2). Variants in genes other than BRCA1/2 were as follows: one BRIP1 c.1794+1 G>A, two CHEK2 c.1555 C>T (p.Arg159*), and one MUTYH c.544 C>T (p.Arg182Cys); one patient had two pathogenic variants in the BRCA2. These patients were diagnosed with IDC tumor type HR+ and HER2‐ and had at least one clinical feature of suspected HBOC, ie, disease onset at less than 40 years of age, bilateral cancer, or family history of cancer.

Characteristics of BRCA 1/2‐variant‐negative patients with pathogenic/likely pathogenic variants in other cancer‐associated genes

Abbreviations: ACMG, American College of Medical Genetics and Genomics; DCIS, ductal carcinoma in situ; FDR, first‐degree relative; HER2, human epidermal growth factor receptor 2; HR, hormone receptors (estrogen receptor, progesterone receptor); IDC, invasive ductal carcinoma; N/A, not available; SDR, second‐degree relative; TNBC, triple‐negative breast cancer.

Reported by Nakagomi H et al
26

### mRNA transcript and pedigree analysis of patients with intronic P/LP variants

Via multigene panel analysis, we identified six P/LP intronic variants: BRCA1 c.302‐2A>C, BRCA1 c.5277+1G>A, CHEK2 c.846+1G>T, PALB2 c.2834+2T>C, BRIP1 c.1794+1 G>A, and PMS2 c.164‐1G>A (Table 4). To predict the exon splicing patterns of these intronic variants, we employed five in silico splice site prediction programs. The natural splicing sites of the variants were predicted to be affected (Table 4). To evaluate the predicted exon splicing effects by in silico programs, we analyzed mRNA transcripts from patient samples via RT‐PCR (Table 4 and Figure S1‐S6).

In silico exon splicing analysis and RT‐PCR results of intronic variants detected in patients with hereditary breast/ovarian cancer

Splice Site Finder

(0‐100)

Max ent scan

(0‐16)

NN

SPICE (0‐1)

Gene splicer

(0‐15)

Classification

(ACMG guidelines)

db

SNP

147

gnom

AD

_exome

_ALL

KRGDB_

1100

Pathogenic

(PVS1, PM2, PP5)

Abbreviations: A, acceptor site; ACMG, American College of Medical Genetics and Genomics; D, donor site; ExAC, Exome Aggregation Consortium; gnomAD, The Genome Aggregation Database; MUT, mutation; NE, splice site not evaluated by the algorithm; WT, wild type; 1000G, 1000 Genomes Project; ㅡ, splice site not detected.

CHEK2 c.846+1G>T was detected in patient PT13, who was diagnosed with HR‐, HER2/neu+ bilateral breast cancer (IDC); the proband's sister died from breast cancer at age 51 after being diagnosed at 47, but other family members were cancer free (Table 3 and Figure S1). CHEK2 c.846+1G>T was predicted to affect the donor splice site (Table 4). We revealed aberrant mRNA transcripts collected from the patient's lymphocytes that corresponded to the skipping of exon 7. This variant was predicted to cause an in‐frame deletion (266‐284 amino acids) in the kinase domain of CHEK2 protein (Figure S1).

PALB2 c.2834+2T>C was detected in patient PT39, who was diagnosed with HR+, HER2/neu‐ IDC at age 35; the proband had a strong family history of reproductive cancer with one paternal aunt being diagnosed with breast cancer at 45 years of age, and another diagnosed with lung cancer at 60 years of age (Table 3 and Figure S2). The proband's paternal grandfather was diagnosed with lung cancer in his late 70s, while the rest of the known family was cancer free. Via RNA analysis, PALB2 c.2834+2T>C was revealed to be the combined product of exon 7 and 9 by exon 8 skipping in PALB2 (Figure S2). This variant was predicted to truncate the PALB2 protein. Generation of truncated protein products by alternative splicing may be associated with an increased risk of breast cancer.
28

PT1 (a 50‐year‐old woman) was referred for genetic counseling after a diagnosis of bilateral IDC and was found to be HR+ and HER2/neu‐. Her family history was only notable with regards to her uncle, who was diagnosed with gastric cancer. Subsequent multigene panel testing revealed the presence of pathogenic BRCA1 c.923_924del, and BRIP1 c.1794+1G>A variants (Table S2 and Figure S3). BRIP1 c.1794+1G>A was predicted to affect the donor splice site, as evidenced by in silico analysis (Table 4). Functional analysis of the intronic variant identified a deletion within exon 12 in the BRIP1 gene, which caused abnormal transcriptional production (Figure S3). The proband's cancer‐related family history and our analysis suggest that this variant confers an increased risk of breast cancer.

PMS2 c.164‐1G>C was detected in patient PT43, who was diagnosed with breast and gastric cancers before age 50. Her father was diagnosed with colon cancer at 73 years old, while her uncle died of gastric cancer at age 60. The rest of the family was healthy (Table 3 and Figure S4). In silico analysis predicted that this variant may be problematic at the acceptor splice site or is activated at the cryptic site (Table 4). Splicing functional assays revealed that this variant induces aberrant splicing via partial exon 3 deletion (8 bp) and leads to the subsequent production of a truncated protein, as evidenced by RT‐PCR (Figure S4).

### Exon splicing analysis for intronic VUS in BRCA1/2

We identified 14 BRCA1 and 24 BRCA2 P/LP variants in our cohort (Table S3). Among these variants, we analyzed exon splicing in two BRCA1 intronic variants classified as P/LP BRCA1 c.302‐2A>C and c.5277+1G>A, which were predicted to affect the donor or acceptor splice sites (Table 4). RT‐PCR identified partial exon 7 deletion (10 bp), and this BRCA1 c.302‐2A>C variant was predicted to truncate BRCA1 protein (Table 4 and Figure S5). The BRCA1 c.302‐2A>C variant was detected in patient PT49, who was diagnosed with bilateral breast cancer. The patient's father died from lung cancer at 61 years old, and her sister was also diagnosed with bilateral breast cancer (Figure S5). BRCA1 c.5277+1G>A was identified in patient PT50, who was diagnosed with ovarian cancer at the age of 44; her family is cancer free. Two abnormal mRNA transcripts were detected in the patient's lymphocytes with BRCA1 c.5277+1G>A that were identified to have an 87‐bp insertion of intron 20 and exon 20 skipping, as seen in a previous study using mini‐gene splicing assays or blood samples (Figure S6).
26

We also examined the effect of one BRCA1 and one BRCA2 VUS on splicing (Table 4, Figure 2, and Figure S7). BRCA1 c.5152+6T>C and BRCA2 c.317‐10A>G are classified with uncertain significance (Table 4). BRCA1 c.5152+6T>C was expected to affect exon splicing, whereas no splicing alterations were predicted for BRCA2 c.317‐10A>G (Table 4). BRCA1 c.5152+6T>C mRNA transcripts were abnormal compared with their corresponding wild‐type transcripts; this variant produced a combined exon 17 and 19 by exon 18 skipping in BRCA1 and was predicted to be an in‐frame deletion of the BRCA1 C‐terminal (BRCT) domain of BRCA1 (Figure 2). However, mRNA transcript indicated that BRCA2 c.317‐10A>G had no effect on exon splicing (Figure S7). Although in silico predictions and RNA analysis revealed the pathogenicity of VUS variants, additional analysis is required to classify them as pathogenic.

Exon splicing analysis of the BRCA1 c.5152+6T>C variant of patient PT51. A, Schematic view of variant c.5152+6T>C localization in the BRCA1 gene. PCR primer alignment is indicated with the red and blue bars. Sequencing analysis for genomic DNA is presented below. B, RT‐PCR of lymphocyte‐derived RNA. Predicted scheme of mRNA transcript in control or patient samples (upper right panel). Agarose gel (2%) electrophoresis; lane 1: control sample; lane 2: patient sample. Two PCR products were detected in the patient sample (upper middle panel). Chromatogram sequences of the control and abnormal transcripts. Vertical line in the chromatogram indicates the exonic junction in transcripts. Exon 18 (78 bp) skipping between exon 17 and exon 19 was identified (upper left panel). Functional domains of BRCA1 and sequence alignment of the BRCA1 abnormal transcript (lower panel). Amino acid sequences of the splice variant (c.5152+6T>C) were aligned using a reference sequence (NP_009225.1) via NCBI BLAST (https://blast.ncbi.nlm.nih.gov/Blast.cgi). BRCA1 c.5152+6T>C was identified to encode a BRCA1 protein with an in‐frame deletion (26 amino acids) in the BRCA1 C‐terminal (BRCT) domain; this may affect the function of the BRCA1 BRCT domain. The red line indicates the location of the in‐frame deletion residues. C, Pedigree of patient PT51

### Genotyping for intronic P/LP and VUS variants

We conducted genotyping of 393 healthy female Korean controls to further define P/LP or VUS intronic variants by comparing their prevalence. The eight intronic variants, analyzed in exon splicing assays, were not detected in healthy controls (Table 4). These results suggest that these intronic variants may affect the susceptibility to inherit breast and ovarian cancer.

### DISCUSSION

Targeted multigene panel analysis can provide detailed genetic information for the identification or management of patients with hereditary cancer.
29
, 
30
 Previous studies showed that expanded panel testing improves the identification of hereditary cancer risk for patients and their family members, as cancer susceptibility genes were identified in 1.9%‐8.1% of patients with BRCA1/2 variant‐negative breast/ovarian cancer (Table S4).
31
, 
32
, 
33
, 
34
, 
35
 By testing other genes besides BRCA1/2, we identified a frequency of 4.4% P/LP variants. These variants were identified in 5 of 19 patients (26.3%) with early‐onset breast cancer (<40 years old at onset). All patients included in the study met the criteria for HBOC genetic testing according to the NCCN 2017 guidelines;
13
 however, 31.6% (6/19) also had a family history of cancers other than HBOC (Table 3). This indicates that a multigene panel study is more effective than a stepwise single‐gene approach for HBOC genetic assessment, as is advised by the NCCN guidelines.
13

Nevertheless, multigene panel testing in clinical settings represents a considerable challenge as these panels include moderate or less well‐defined genes as well as high‐penetrant genes.
36
, 
37
, 
38
 Lack of clear management guidelines for variants in genes with undefined cancer risks or P/LP variants in genes can be problematic. A variant cannot be classified as a positive pathogenic result without an experimental study. Another concern is that the risk of overestimating the clinical interpretation of VUS results in low‐ to moderate‐risk genes. In the present study, we identified that 49.0% of patients had VUS within 23 genes including BRCA1/2 (data not shown). As the number of genes tested and the frequency of multigene panel testing continue to increase, the rate of VUS detection would also increase.
30
, 
33
, 
39

In this study, we observed MUTYH heterozygote c.857G>A (p.Gly286Glu) (three cases) and c.544C>T (p.Arg182Cys) variants (one case) with P/LP findings based on ClinVar data. It is known that biallelic (homozygous or compound heterozygous) MUTYH variants are related to MUTYH‐associated polyposis syndrome, which results in colorectal polyps and colorectal cancer; however, their association with malignancies other than colon cancer is less robust.
40
, 
41
 Previous studies have reported an increased risk of breast cancer, without statistical evidence in monoallelic MUTYH variants.
42
, 
43
, 
44
 In a study that enrolled Sephardi Jews of North African descent, homozygote or heterozygote carriers of p.Gly396Asp in MYUTH were found to be significantly increased in breast cancer patients (6.7%) compared with controls (3.7%) (OR, 1.39; 95% CI, 0.26‐7.53).
44
 Although a higher frequency of monoallelic MUTYH variants in families with both breast and colorectal cancer compared with those in the general population is increasingly being reported,
40
, 
42
, 
43
, 
45
 more evidence regarding the association between MUTYH variants and other cancers should be elucidated.

We performed mRNA transcript analysis of eight intronic variants in BRCA1, BRCA2, BRIP1, CHEK2, PALB2, and PMS2, which were classified as P/LP or VUS. The P/LP variants showed abnormal transcriptional fragments. CHEK2 (cell cycle checkpoint kinase 2) is a well‐established moderate‐penetrance breast cancer gene, but it lacks treatment and follow‐up guidelines.
1
, 
46

PALB2 (partner and localizer of BRCA2) serves a crucial role in the localization and stabilization of BRCA2 in nuclear chromatin, which is essential for BRCA2 to function in the homologous recombination‐mediated repair of double‐strand DNA breaks (DSBs)
47
, 
48
; PALB2 variants have been reported to be associated with pancreatic cancer development.
49

BRIP1 (BRCA1‐interacting protein C‐terminal helicase 1) encodes proteins that interact with BRCA1 during the repair of DSBs, and pathogenic variants of this gene have been investigated.
13
, 
50
 Germline pathogenic variants of PMS2 (PMS1 homolog 2) are implicated in Lynch syndrome and are associated with a significantly increased risk of breast cancer.
51
 Previous studies have attempted to identify some genes associated with DNA repair, such as ATM and CHEK2, which have also been added to breast–cancer‐specific gene panels.
13
, 
14
, 
52
, 
53
 However, there is controversy over whether these rare variants are clinically associated with a risk of breast cancer
36
, 
37
, 
54
; nonetheless, evidence regarding breast cancer incidence is limited.

In this study, the BRIP1 c.1794+1G>A (PT1) carrier was found to possess a c.923_924del variant in the BRCA1 gene (Table 4, Table S2, and Figure S3). Recently, BRIP1 c.1794+1G>A was registered as likely pathogenic in the ClinVar database, but its effects have not been reported in the literature. In our study, exon splicing analysis indicated that BRIP1 c.1794+1G>A results in exon 12 deletion, leading to a frameshift mutation that creates a premature stop codon in the BRIP1 protein (stop codon gained at 557 a.a.; reference sequence NP_114432.2; Table S2 and Figure S3). BRCA1 c.923_924del (p.Ser308Lysfs*11) is another frameshift mutation resulting from a deletion (Table S2). BRIP1 is a DNA helicase which interacts with the C‐terminal BRCT domain (1646‐1736, 1760‐1855 a.a.) of BRCA1 through its C‐terminal‐BRCA1‐binding domain (888‐1063 a.a.) and functions in BRCA1‐dependent DNA repair and DNA‐induced checkpoint activity.
55
 De Nicolo et al
50
 suggested that a heterozygous germline variant in the BRIP1 gene results in a truncated protein product and is associated with loss of the wild‐type BRIP1 allele in the tissues of affected breast cancer patients. Thus, the BRIP1 c.1794+1G>A/BRCA1 c.923_924del (p.Ser308Lysfs*11) double variants of patient PT1 may further increase the risk of breast cancer through the instability and functional impairment of the encoded proteins.

We further showed that BRCA1 c.5152+6T>C, classified as VUS in ClinVar (https://www.ncbi.nlm.nih.gov/clinvar/), has an in‐frame deletion in the BRCT domain resulting from exon 18 skipping. BRCT domains can form a phospho‐recognition motif that preferentially binds proteins containing phosphoserine and interact with several proteins implicated in DNA repair, including Abraxas, BRIP1, and CtIP.
56

BRCA1 c.5152+6T>C was detected in patient PT51, who was diagnosed with triple‐negative breast cancer (TNBC) IDC at 39 years of age. The proband's grandfather died of lung cancer and the grandmother died of ovarian cancer at 77 years old (Figure 2). Splicing variants in BRCT domains of BRCA1 have been reported to be associated with aberrant splicing in patients with breast/ovarian cancer.
27
, 
57
 This variant does not have frequency information in genome databases, including ExAC (http://exac.broadinstitute.org/) and gnomAD (http://gnomad.broadinstitue.org), and has not been reported in the literature. Moreover, BRCA1 c.5152+6T>C was not detected in the 393 healthy female Korean controls (Table 4). By employing a saturation‐genome‐editing technique based on CRISPR‐mediated homology‐directed repair, Findlay et al
58
 suggested that BRCA1 c.5152+6T>C is a loss‐of‐function variant. Thus, we propose that BRCA1 c.5152+6T>C be reclassified as likely pathogenic. Nevertheless, further analysis of the patient and relatives is needed to clarify the actual clinical impact of this variant. In addition, the detected pathogenic variants should have moderately established carcinogenic lifetime risk as well as appropriate counseling recommendations for the patient.

The use of customized multigene panels to confirm associations with genes other than BRCA1/2 in patients with HBOC has increased, and through this approach we have revealed additional pathogenic variants in 4.4% of cases. Although this study included fewer ovarian cancer patients than breast cancer patients due to the low consent rate, our results highlight the importance of performing multigene panel testing of patients with HBOC in the Korean population as an alternative strategy for identifying shaded P/LP variants. We also demonstrated how exon splicing analysis by conducting in silico predictions or functional studies using patient samples can be beneficial in the identification of uncharacterized intronic variants that are expected to increase HBOC risk. Finally, further analysis is warranted to determine the clinical impact and patient outcomes associated with the identification of P/LP variants in non‐BRCA1/2.

### ETHICAL CONSIDERATIONS

This study was approved by the International Review Board (IRB) of the National Cancer Center of Korea (IRB No. NCCNCS13717 and NCC2017‐0127), and written informed consent was obtained from all participating patients.

### CONFLICT OF INTEREST

None.

### Supporting information

Supplementary Materials

Click here for additional data file.



# SUPPLEMENTAL FILE 1: CAS-111-3912.pdf

# Preparing to download ...

[HHS Vulnerability Disclosure](https://www.hhs.gov/vulnerability-disclosure-policy/index.html)