# MAIN TEXT

## Analysis of high‐risk pedigrees identifies 11 candidate variants for Alzheimer's disease

### Abstract

AbstractIntroductionAnalysis of sequence data in high‐risk pedigrees is a powerful approach to detect rare predisposition variants.MethodsRare, shared candidate predisposition variants were identified from exome sequencing 19 Alzheimer's disease (AD)‐affected cousin pairs selected from high‐risk pedigrees. Variants were further prioritized by risk association in various external datasets. Candidate variants emerging from these analyses were tested for co‐segregation to additional affected relatives of the original sequenced pedigree members.ResultsAD‐affected high‐risk cousin pairs contained 564 shared rare variants. Eleven variants spanning 10 genes were prioritized in external datasets: rs201665195 (ABCA7), and rs28933981 (TTR) were previously implicated in AD pathology; rs141402160 (NOTCH3) and rs140914494 (NOTCH3) were previously reported; rs200290640 (PIDD1) and rs199752248 (PIDD1) were present in more than one cousin pair; rs61729902 (SNAP91), rs140129800 (COX6A2, AC026471), and rs191804178 (MUC16) were not present in a longevity cohort; and rs148294193 (PELI3) and rs147599881 (FCHO1) approached significance from analysis of AD‐related phenotypes. Three variants were validated via evidence of co‐segregation to additional relatives (PELI3, ABCA7, and SNAP91).DiscussionThese analyses support ABCA7 and TTR as AD risk genes, expand on previously reported NOTCH3 variant identification, and prioritize seven additional candidate variants.

### INTRODUCTION

Genetic variation significantly impacts Alzheimer's disease (AD) risk and is estimated to account for 53% of total trait variance.
1
 To date, genome‐wide association studies (GWAS) for AD have identified common variants occurring in more than 30 genes.
1
, 
2
 Despite this progress, it is estimated that more than 40% of the genetic variance of AD remains uncharacterized.
3
 While rare variation contributing to AD is difficult to detect via GWAS, it may be discoverable in pedigree‐based designs. It is well recognized that high‐risk pedigree studies are a powerful and efficient method for identification of rare predisposition variants and should be performed when these rare resources are available.
4
, 
5
, 
6
 Pedigree‐based studies are ideal for rare variant identification because rare variants can occur at a higher‐than‐population rate among related affected individuals, thereby enhancing statistical power.
7
, 
8
 Additionally, pedigree analyses limit locus heterogeneity (i.e., different gene loci or gene loci interactions causing a similar phenotype) because the rare variants are inherited by a common ancestor with a common haplotype. Although familial locus heterogeneity may exist, the success rate at identifying causal variants in large pedigrees can be much higher than in clinical settings.
9
 Pedigree‐based designs have been successfully applied to gene discovery for many phenotypes
10
, 
11
, 
12
 including AD, where notable examples are PLD3,
13

NOTCH3,
14
 and RAB10.
15
 Here, a pedigree‐based analysis using exome sequences from 19 AD‐affected cousin pairs belonging to pedigrees with a statistical excess of AD mortality was performed to identify genetic variants associated with AD risk. Candidate variants arising from this analysis were then further prioritized in publicly available, large‐scale datasets to further evaluate their risk of AD. Prioritized variants were then tested for co‐segregation in additional previously sampled relatives of the index cousin pairs.

### MATERIAL AND METHODS

The UPDB includes a genealogy of Utah, representing the founders in the mid‐1800s to their modern‐day descendants. It originated from three‐generation genealogy records
16
 and is kept current with Utah vital statistics data, which includes approximately 50,000 largely unrelated founders of European descent. More than three million individuals with at least three generations of genealogy connecting to Utah founders in the UPDB were analyzed here. The UPDB has been linked to various phenotypic data including the Utah Cancer Registry from 1966 and Utah death certificates from 1904 to 2014. This unique combination of genealogy with phenotypes has contributed to the identification, recruitment, consent, and sampling of over 30,000 individuals in thousands of high‐risk pedigrees representing many disorders. Previous analyses of the UPDB genealogy have reported that founder effects are not present among the Utah population and indicate low inbreeding levels similar to other places in the United States.
17
, 
18

Systematic review: Despite recent advances in identifying common genetic variants associated with Alzheimer's disease (AD), it is estimated that more than 40% of the genetic variance of AD remains uncharacterized. High‐risk pedigree studies offer more statistical power for identification of rare predisposition variants, and they have previously been used to identify variants in PLD3, NOTCH3, and RAB10.
Interpretation: Pedigree‐based analyses using exome sequences from 19 AD‐affected cousin pairs with extended segregation assays were performed. These analyses add support to ABCA7 and TTR as AD risk genes. Furthermore, seven other candidate variants are prioritized in addition to NOTCH3 variants that were previously reported.
Future directions: Our analyses provide support for rare variant prioritization through pedigree‐based analyses. Additional inquiries into ABCA7 and TTR as AD risk genes are warranted. Moreover, these analyses suggest rare variants in NOTCH3, PIDD1, SNAP91, COX6A2, MUC16, PELI3, and FCHO1 may contribute to AD risk.

Systematic review: Despite recent advances in identifying common genetic variants associated with Alzheimer's disease (AD), it is estimated that more than 40% of the genetic variance of AD remains uncharacterized. High‐risk pedigree studies offer more statistical power for identification of rare predisposition variants, and they have previously been used to identify variants in PLD3, NOTCH3, and RAB10.

Interpretation: Pedigree‐based analyses using exome sequences from 19 AD‐affected cousin pairs with extended segregation assays were performed. These analyses add support to ABCA7 and TTR as AD risk genes. Furthermore, seven other candidate variants are prioritized in addition to NOTCH3 variants that were previously reported.

Future directions: Our analyses provide support for rare variant prioritization through pedigree‐based analyses. Additional inquiries into ABCA7 and TTR as AD risk genes are warranted. Moreover, these analyses suggest rare variants in NOTCH3, PIDD1, SNAP91, COX6A2, MUC16, PELI3, and FCHO1 may contribute to AD risk.

Among the individuals with stored DNA from Utah pedigrees were 199 subjects whose death certificates indicated AD as a cause of death. These 199 individuals were related in 102 independent, descending pedigrees including between two and six sampled AD mortality cases. Each pedigree was tested for an excess of AD deaths among descendants as described elsewhere.
15
using all available genealogy and death certificate data. Twenty‐four high‐risk pedigrees including two to four sampled AD cases were identified. Of those pedigrees, 19 unrelated pedigrees were selected for this analysis that included at least one sampled AD‐affected approximate‐cousin pair, ranging in relationship from avuncular (n = 1; expected shared DNA = 25%) to third cousin (n = 1; expected shared DNA = 0.78%).

Two AD‐affected cousins from each of the 19 high‐risk pedigrees were exome sequenced at the Huntsman Cancer Institute's Genomics Core Facility. The Agilent SureSelect XT Human All Exon + UTR (v5) capture kit was used to prepare the DNA library from 2 μg of DNA per sample. Paired‐end reads of up to 150 base pairs were sequenced on the Illumina HiSeq 2000 sequencer. BWA‐MEM
19
, 
20
 mapped raw reads to the human genome v37 (GRCh37) reference genome. The Genome Analysis Toolkit 3.5.0 (GATK)
21
 software called variants using the Broad Institute's best practices guidelines. Variants were removed if they occurred outside the exon capture kit intended area of coverage, and the remaining variants were annotated with ANNOVAR
22
 for their predicted pathogenicity.

Pathogenicity predictions were conducted in ANNOVAR using various algorithms because deleterious predictions made by multiple prediction methods are more likely to be reliable.
23
 Table S1 in supporting information shows scores reported by SIFT,
24
, 
25
 Polyphen‐2,
26
 LRT,
23
 MutationTaster2,
27
 FATHMM,
28
 PROVEAN,
29
 VEST3,
30
 MetaSVM, MetaLR, M‐CAP,
31
 REVEL,
32
 CADD,
33
 DANN,
34
 EIGEN,
35
 GenoCanyon,
36
 and GERP++.
37
 Prioritized variants were predicted to be deleterious by at least two of those functional prediction algorithms.

Some variant prioritization (described below) was based on comparisons to the Alzheimer's Disease Genetics Consortium (ADGC) dataset. ADGC contains imputed single‐nucleotide polymorphism (SNP) array data for 28,730 subjects (11,967 males and 16,760 females), including 10,486 AD cases and 10,168 healthy controls. Of the shared rare exonic variants in the cousin pairs, 291 were sequenced or imputed in the ADGC dataset.

The Knight Alzheimer's Disease Research Center at Washington University School of Medicine (Knight ADRC) cerebrospinal fluid (CSF) dataset was used for analyses of association between candidate variants and levels of AB42, Tau, or PTau. This cohort consisted of 3963 subjects (1895 males and 1675 females), including 1479 AD cases and 1370 healthy controls. All samples were genotyped using the Illumina 610 or OmniExpress chip. Of the shared rare variants in the cousin pairs, only 12 appeared in the Knight ADRC dataset. Additional imputation on the Knight ADRC dataset would likely lead to a higher false discovery rate because rare variant imputation is much less accurate than common variant imputation.
38
 Therefore, those 12 variants were used only for additional variant prioritization on the rare variants identified in the high‐risk pedigrees.

The Wellderly study
39
 is an ongoing study that includes elderly people (age 80–105) who are cognitively healthy without medical interventions. Six hundred Wellderly individuals had their whole genomes sequenced using the Complete Genomics platform. These genomes were used to perform additional variant prioritization because functionally relevant rare genetic variants associated with AD should not be present in an elderly population that does not exhibit cognitive decline.

Polygenic risk scores for each cousin were calculated and compared to ADGC cases and controls to ensure that the excess risk for AD was not caused by an excess of common disease‐associated variants (see Note S1 in supporting information). The AD polygenic risk scores calculated from Lambert, Ibrahim‐Verbaas,
2
 and apolipoprotein E (APOE) genotype for each cousin are shown in Table S2 in supporting information. A Welch's two‐sample t‐test shows that the mean AD risk from common variants in the high‐risk cousin pairs is significantly less than the mean AD risk for ADGC cases (P‐value = 2.836 × 10−9; see Figure S1 in supporting information) and controls (P‐value = 9.486 × 10−4; see Figure S2 in supporting information). Therefore, the propensity of AD in these high‐risk pedigrees is likely caused by rare genetic variants that can be prioritized using the pipeline shown below.

Candidate variants were initially required to be present in both AD‐affected cousin pairs from at least 1 of the 19 high‐risk AD pedigrees. Further analysis in Ingenuity Variant Analysis software (QIAGEN, Inc.) ensured variants were rare by removing variants with a population minor allele frequency (MAF) greater than 0.01 in 1000 Genomes,
40
 Exome Aggregation Consortium (ExAC) non‐Finnish European,
41
 The Genome Aggregation Database (gnomAD),
41
 and NHLBI GO Exome Sequencing Project (ESP), Seattle, WA (http://evs.gs.washington.edu/EVS/) [March 2018].

Ingenuity Variant Analysis then prioritized variants based on their predicted pathogenicity. Variants considered “Pathogenic,” “Likely Pathogenic,” or “Unknown” by the American College of Medical Genetics (ACMG),
42
 or resulted in either a loss or gain of gene function by in silico functional prediction algorithms were prioritized.

Four distinct strategies were used to further prioritize candidate variants in external datasets related to AD. Each strategy was evaluated independently of the others. Variants that met any one of these criteria were prioritized as candidate variants for AD. Variants were assumed to be likely causing the excess of AD mortality in the pedigrees if they were (1) a known AD risk variant, (2) a CSF biomarker of AD with a P‐value less than .1, (3) not present in the Wellderly dataset and more prevalent in ADGC AD cases than ADGC controls, or (4) observed in more than one AD‐affected cousin pair, as described below.

A literature search of all shared rare variants was conducted in Ingenuity Variant Analysis to determine the extent to which these variants were previously implicated in AD pathology, indicating additional support from independent studies. Variants with publications supporting their association with AD were prioritized.

Three linear regressions were conducted on the Knight ADRC CSF data using PLINK
43
: one for each of the phenotypes of interest (AB42, Tau, and PTau), with age, APOE status, sex, and two principal components
2
, 
44
 as covariates. The threshold for prioritizing candidate variants was P‐value less than .1 for any of the three phenotypes. The significance threshold was relaxed because variants had a low minor allele frequency and a small sample size.

Variants positively affecting AD mortality are expected to be more prevalent in diseased elderly cohorts than healthy elderly cohorts. The Wellderly dataset was used to ensure that prioritized variants were not present in a healthy longevity cohort. Furthermore, the MAF of the variant in ADGC controls needed to be less than or equal to the MAF in ADGC AD cases to ensure that the variant was more prevalent in AD cases than controls.

Variants that were identified in more than one AD‐affected cousin pair were also prioritized as candidates. This prioritization provided additional evidence for the candidate variant as a potential AD risk factor because of its prevalence in AD subjects from multiple high‐risk pedigrees. Because all variants analyzed in this study are rare variants (MAF < 0.01), it is not likely that the same rare variant would be present in two independent pedigrees and shared by both affected cousins by random chance.

A set of 199 individuals previously sampled for Utah high‐risk disease studies whose Utah death certificate included AD as a cause of death and had at least three generations of genealogy linking to Utah founders, was assayed for 9 of the 11 variants. An assay was not available for PIDD1 rs200290640, and the assay for PIDD1 rs199752248 failed due to too many homologous regions in the genome to specify the right location. The AD cases from all affected cousin pairs were included, and all assays correctly identified the observed carriers for each of the nine candidate variants. Evidence of co‐segregation to additional affected relatives was evaluated with the RVsharing program.
45
 The RVsharing calculates the probability of an observed configuration of carriers and affection status in a pedigree having occurred by chance transmission, assuming that the variant is rare (MAF < 0.01) and has entered the pedigree only once. The RVsharing program provides a test of co‐segregation (e.g., linkage) in which the strength of evidence against the null hypothesis of no co‐segregation between the disease and a variant is expressed as an exact probability (P‐value) for a given pedigree structure and disease configuration. Because the pedigrees were pre‐screened for a statistical excess of AD, rare variants that were shared among the related individuals likely contribute to the excess in AD cases. Institutional Review Board approval was in place for all reported analyses.

### Utah Population Database (UPDB)

The UPDB includes a genealogy of Utah, representing the founders in the mid‐1800s to their modern‐day descendants. It originated from three‐generation genealogy records
16
 and is kept current with Utah vital statistics data, which includes approximately 50,000 largely unrelated founders of European descent. More than three million individuals with at least three generations of genealogy connecting to Utah founders in the UPDB were analyzed here. The UPDB has been linked to various phenotypic data including the Utah Cancer Registry from 1966 and Utah death certificates from 1904 to 2014. This unique combination of genealogy with phenotypes has contributed to the identification, recruitment, consent, and sampling of over 30,000 individuals in thousands of high‐risk pedigrees representing many disorders. Previous analyses of the UPDB genealogy have reported that founder effects are not present among the Utah population and indicate low inbreeding levels similar to other places in the United States.
17
, 
18

Systematic review: Despite recent advances in identifying common genetic variants associated with Alzheimer's disease (AD), it is estimated that more than 40% of the genetic variance of AD remains uncharacterized. High‐risk pedigree studies offer more statistical power for identification of rare predisposition variants, and they have previously been used to identify variants in PLD3, NOTCH3, and RAB10.
Interpretation: Pedigree‐based analyses using exome sequences from 19 AD‐affected cousin pairs with extended segregation assays were performed. These analyses add support to ABCA7 and TTR as AD risk genes. Furthermore, seven other candidate variants are prioritized in addition to NOTCH3 variants that were previously reported.
Future directions: Our analyses provide support for rare variant prioritization through pedigree‐based analyses. Additional inquiries into ABCA7 and TTR as AD risk genes are warranted. Moreover, these analyses suggest rare variants in NOTCH3, PIDD1, SNAP91, COX6A2, MUC16, PELI3, and FCHO1 may contribute to AD risk.

Systematic review: Despite recent advances in identifying common genetic variants associated with Alzheimer's disease (AD), it is estimated that more than 40% of the genetic variance of AD remains uncharacterized. High‐risk pedigree studies offer more statistical power for identification of rare predisposition variants, and they have previously been used to identify variants in PLD3, NOTCH3, and RAB10.

Interpretation: Pedigree‐based analyses using exome sequences from 19 AD‐affected cousin pairs with extended segregation assays were performed. These analyses add support to ABCA7 and TTR as AD risk genes. Furthermore, seven other candidate variants are prioritized in addition to NOTCH3 variants that were previously reported.

Future directions: Our analyses provide support for rare variant prioritization through pedigree‐based analyses. Additional inquiries into ABCA7 and TTR as AD risk genes are warranted. Moreover, these analyses suggest rare variants in NOTCH3, PIDD1, SNAP91, COX6A2, MUC16, PELI3, and FCHO1 may contribute to AD risk.

### RESEARCH IN CONTEXT

Systematic review: Despite recent advances in identifying common genetic variants associated with Alzheimer's disease (AD), it is estimated that more than 40% of the genetic variance of AD remains uncharacterized. High‐risk pedigree studies offer more statistical power for identification of rare predisposition variants, and they have previously been used to identify variants in PLD3, NOTCH3, and RAB10.
Interpretation: Pedigree‐based analyses using exome sequences from 19 AD‐affected cousin pairs with extended segregation assays were performed. These analyses add support to ABCA7 and TTR as AD risk genes. Furthermore, seven other candidate variants are prioritized in addition to NOTCH3 variants that were previously reported.
Future directions: Our analyses provide support for rare variant prioritization through pedigree‐based analyses. Additional inquiries into ABCA7 and TTR as AD risk genes are warranted. Moreover, these analyses suggest rare variants in NOTCH3, PIDD1, SNAP91, COX6A2, MUC16, PELI3, and FCHO1 may contribute to AD risk.

Systematic review: Despite recent advances in identifying common genetic variants associated with Alzheimer's disease (AD), it is estimated that more than 40% of the genetic variance of AD remains uncharacterized. High‐risk pedigree studies offer more statistical power for identification of rare predisposition variants, and they have previously been used to identify variants in PLD3, NOTCH3, and RAB10.

Interpretation: Pedigree‐based analyses using exome sequences from 19 AD‐affected cousin pairs with extended segregation assays were performed. These analyses add support to ABCA7 and TTR as AD risk genes. Furthermore, seven other candidate variants are prioritized in addition to NOTCH3 variants that were previously reported.

Future directions: Our analyses provide support for rare variant prioritization through pedigree‐based analyses. Additional inquiries into ABCA7 and TTR as AD risk genes are warranted. Moreover, these analyses suggest rare variants in NOTCH3, PIDD1, SNAP91, COX6A2, MUC16, PELI3, and FCHO1 may contribute to AD risk.

### Identification of sampled high‐risk AD mortality pedigrees

Among the individuals with stored DNA from Utah pedigrees were 199 subjects whose death certificates indicated AD as a cause of death. These 199 individuals were related in 102 independent, descending pedigrees including between two and six sampled AD mortality cases. Each pedigree was tested for an excess of AD deaths among descendants as described elsewhere.
15
using all available genealogy and death certificate data. Twenty‐four high‐risk pedigrees including two to four sampled AD cases were identified. Of those pedigrees, 19 unrelated pedigrees were selected for this analysis that included at least one sampled AD‐affected approximate‐cousin pair, ranging in relationship from avuncular (n = 1; expected shared DNA = 25%) to third cousin (n = 1; expected shared DNA = 0.78%).

### Exome sequencing of AD cousin pairs

Two AD‐affected cousins from each of the 19 high‐risk pedigrees were exome sequenced at the Huntsman Cancer Institute's Genomics Core Facility. The Agilent SureSelect XT Human All Exon + UTR (v5) capture kit was used to prepare the DNA library from 2 μg of DNA per sample. Paired‐end reads of up to 150 base pairs were sequenced on the Illumina HiSeq 2000 sequencer. BWA‐MEM
19
, 
20
 mapped raw reads to the human genome v37 (GRCh37) reference genome. The Genome Analysis Toolkit 3.5.0 (GATK)
21
 software called variants using the Broad Institute's best practices guidelines. Variants were removed if they occurred outside the exon capture kit intended area of coverage, and the remaining variants were annotated with ANNOVAR
22
 for their predicted pathogenicity.

Pathogenicity predictions were conducted in ANNOVAR using various algorithms because deleterious predictions made by multiple prediction methods are more likely to be reliable.
23
 Table S1 in supporting information shows scores reported by SIFT,
24
, 
25
 Polyphen‐2,
26
 LRT,
23
 MutationTaster2,
27
 FATHMM,
28
 PROVEAN,
29
 VEST3,
30
 MetaSVM, MetaLR, M‐CAP,
31
 REVEL,
32
 CADD,
33
 DANN,
34
 EIGEN,
35
 GenoCanyon,
36
 and GERP++.
37
 Prioritized variants were predicted to be deleterious by at least two of those functional prediction algorithms.

### The Alzheimer's Disease Genetics Consortium Dataset

Some variant prioritization (described below) was based on comparisons to the Alzheimer's Disease Genetics Consortium (ADGC) dataset. ADGC contains imputed single‐nucleotide polymorphism (SNP) array data for 28,730 subjects (11,967 males and 16,760 females), including 10,486 AD cases and 10,168 healthy controls. Of the shared rare exonic variants in the cousin pairs, 291 were sequenced or imputed in the ADGC dataset.

### The Knight cerebrospinal fluid dataset

The Knight Alzheimer's Disease Research Center at Washington University School of Medicine (Knight ADRC) cerebrospinal fluid (CSF) dataset was used for analyses of association between candidate variants and levels of AB42, Tau, or PTau. This cohort consisted of 3963 subjects (1895 males and 1675 females), including 1479 AD cases and 1370 healthy controls. All samples were genotyped using the Illumina 610 or OmniExpress chip. Of the shared rare variants in the cousin pairs, only 12 appeared in the Knight ADRC dataset. Additional imputation on the Knight ADRC dataset would likely lead to a higher false discovery rate because rare variant imputation is much less accurate than common variant imputation.
38
 Therefore, those 12 variants were used only for additional variant prioritization on the rare variants identified in the high‐risk pedigrees.

### Wellderly dataset

The Wellderly study
39
 is an ongoing study that includes elderly people (age 80–105) who are cognitively healthy without medical interventions. Six hundred Wellderly individuals had their whole genomes sequenced using the Complete Genomics platform. These genomes were used to perform additional variant prioritization because functionally relevant rare genetic variants associated with AD should not be present in an elderly population that does not exhibit cognitive decline.

### Assessing common genetic variants

Polygenic risk scores for each cousin were calculated and compared to ADGC cases and controls to ensure that the excess risk for AD was not caused by an excess of common disease‐associated variants (see Note S1 in supporting information). The AD polygenic risk scores calculated from Lambert, Ibrahim‐Verbaas,
2
 and apolipoprotein E (APOE) genotype for each cousin are shown in Table S2 in supporting information. A Welch's two‐sample t‐test shows that the mean AD risk from common variants in the high‐risk cousin pairs is significantly less than the mean AD risk for ADGC cases (P‐value = 2.836 × 10−9; see Figure S1 in supporting information) and controls (P‐value = 9.486 × 10−4; see Figure S2 in supporting information). Therefore, the propensity of AD in these high‐risk pedigrees is likely caused by rare genetic variants that can be prioritized using the pipeline shown below.

### Initial variant prioritization

Candidate variants were initially required to be present in both AD‐affected cousin pairs from at least 1 of the 19 high‐risk AD pedigrees. Further analysis in Ingenuity Variant Analysis software (QIAGEN, Inc.) ensured variants were rare by removing variants with a population minor allele frequency (MAF) greater than 0.01 in 1000 Genomes,
40
 Exome Aggregation Consortium (ExAC) non‐Finnish European,
41
 The Genome Aggregation Database (gnomAD),
41
 and NHLBI GO Exome Sequencing Project (ESP), Seattle, WA (http://evs.gs.washington.edu/EVS/) [March 2018].

Ingenuity Variant Analysis then prioritized variants based on their predicted pathogenicity. Variants considered “Pathogenic,” “Likely Pathogenic,” or “Unknown” by the American College of Medical Genetics (ACMG),
42
 or resulted in either a loss or gain of gene function by in silico functional prediction algorithms were prioritized.

### Prioritization in external datasets

Four distinct strategies were used to further prioritize candidate variants in external datasets related to AD. Each strategy was evaluated independently of the others. Variants that met any one of these criteria were prioritized as candidate variants for AD. Variants were assumed to be likely causing the excess of AD mortality in the pedigrees if they were (1) a known AD risk variant, (2) a CSF biomarker of AD with a P‐value less than .1, (3) not present in the Wellderly dataset and more prevalent in ADGC AD cases than ADGC controls, or (4) observed in more than one AD‐affected cousin pair, as described below.

A literature search of all shared rare variants was conducted in Ingenuity Variant Analysis to determine the extent to which these variants were previously implicated in AD pathology, indicating additional support from independent studies. Variants with publications supporting their association with AD were prioritized.

Three linear regressions were conducted on the Knight ADRC CSF data using PLINK
43
: one for each of the phenotypes of interest (AB42, Tau, and PTau), with age, APOE status, sex, and two principal components
2
, 
44
 as covariates. The threshold for prioritizing candidate variants was P‐value less than .1 for any of the three phenotypes. The significance threshold was relaxed because variants had a low minor allele frequency and a small sample size.

Variants positively affecting AD mortality are expected to be more prevalent in diseased elderly cohorts than healthy elderly cohorts. The Wellderly dataset was used to ensure that prioritized variants were not present in a healthy longevity cohort. Furthermore, the MAF of the variant in ADGC controls needed to be less than or equal to the MAF in ADGC AD cases to ensure that the variant was more prevalent in AD cases than controls.

Variants that were identified in more than one AD‐affected cousin pair were also prioritized as candidates. This prioritization provided additional evidence for the candidate variant as a potential AD risk factor because of its prevalence in AD subjects from multiple high‐risk pedigrees. Because all variants analyzed in this study are rare variants (MAF < 0.01), it is not likely that the same rare variant would be present in two independent pedigrees and shared by both affected cousins by random chance.

### Known AD risk variant

A literature search of all shared rare variants was conducted in Ingenuity Variant Analysis to determine the extent to which these variants were previously implicated in AD pathology, indicating additional support from independent studies. Variants with publications supporting their association with AD were prioritized.

### Increased AD risk in CSF dataset

Three linear regressions were conducted on the Knight ADRC CSF data using PLINK
43
: one for each of the phenotypes of interest (AB42, Tau, and PTau), with age, APOE status, sex, and two principal components
2
, 
44
 as covariates. The threshold for prioritizing candidate variants was P‐value less than .1 for any of the three phenotypes. The significance threshold was relaxed because variants had a low minor allele frequency and a small sample size.

### AD risk gradient

Variants positively affecting AD mortality are expected to be more prevalent in diseased elderly cohorts than healthy elderly cohorts. The Wellderly dataset was used to ensure that prioritized variants were not present in a healthy longevity cohort. Furthermore, the MAF of the variant in ADGC controls needed to be less than or equal to the MAF in ADGC AD cases to ensure that the variant was more prevalent in AD cases than controls.

### Multiple hit pedigrees

Variants that were identified in more than one AD‐affected cousin pair were also prioritized as candidates. This prioritization provided additional evidence for the candidate variant as a potential AD risk factor because of its prevalence in AD subjects from multiple high‐risk pedigrees. Because all variants analyzed in this study are rare variants (MAF < 0.01), it is not likely that the same rare variant would be present in two independent pedigrees and shared by both affected cousins by random chance.

### Evidence for segregation with AD risk in additional affected sampled relatives

A set of 199 individuals previously sampled for Utah high‐risk disease studies whose Utah death certificate included AD as a cause of death and had at least three generations of genealogy linking to Utah founders, was assayed for 9 of the 11 variants. An assay was not available for PIDD1 rs200290640, and the assay for PIDD1 rs199752248 failed due to too many homologous regions in the genome to specify the right location. The AD cases from all affected cousin pairs were included, and all assays correctly identified the observed carriers for each of the nine candidate variants. Evidence of co‐segregation to additional affected relatives was evaluated with the RVsharing program.
45
 The RVsharing calculates the probability of an observed configuration of carriers and affection status in a pedigree having occurred by chance transmission, assuming that the variant is rare (MAF < 0.01) and has entered the pedigree only once. The RVsharing program provides a test of co‐segregation (e.g., linkage) in which the strength of evidence against the null hypothesis of no co‐segregation between the disease and a variant is expressed as an exact probability (P‐value) for a given pedigree structure and disease configuration. Because the pedigrees were pre‐screened for a statistical excess of AD, rare variants that were shared among the related individuals likely contribute to the excess in AD cases. Institutional Review Board approval was in place for all reported analyses.

### RESULTS

Initial variant prioritization of rare variants shared by at least one cousin pair identified 400 rare variants spanning 470 genes. Ingenuity Variant Analysis subsequently prioritized 382 variants in 447 genes that were pathogenic, likely pathogenic, uncertain significance, or associated with a gain or loss of gene function. Of those shared rare variants, 117 variants in 138 genes had a biological interaction with genes implicated in AD. The numbers of variants included after each analysis are shown in Figure 1, and each variant is listed in File S1 in supporting information. That list of 117 variants was then used by four independent prioritization screens to identify 11 rare variants spanning 10 genes with the strongest support for increasing AD risk in these high‐risk pedigrees (see Table 1). Four rare variants previously associated with increased AD risk were identified using a literature search in Ingenuity Variant Analysis. The Knight ADRC CSF dataset identified two additional variants (Note: of the 117 prioritized variants, only 11 appeared in the CSF dataset) that were associated with increased risk for AB42, Tau, or Ptau in CSF. Additionally, the AD Risk Gradient identified three variants that were not present in the Wellderly dataset and were more prevalent in ADGC cases than ADGC controls. Finally, two variants were present in two independent high‐risk AD pedigrees.

Flowchart depicting variant prioritization. The number of rare variants and genes remaining at each level are shown. ADRC, Alzheimer's Disease Research Center; CSF, cerebrospinal fluid; MAF, minor allele frequency.

Prioritized variants

Notes: These variants are most likely to affect AD pathology in the high‐risk pedigrees. The accession number, affected gene, the Human Genome Variation Society (HGVS) variant annotation, impact on translation, and the final prioritization option that identified the variant are shown.

Abbreviations: AD, Alzheimer's disease; CSF, cerebrospinal fluid.

Five candidate variants assayed for segregation evidence had no additional carriers beyond the original cousin pair, including: MUC16 rs191804178, NOTCH3 rs141402160, NOTCH3 rs140914494, FCHO1 rs147599881, and COX6A2 rs140129800. An additional carrier was identified among the 199 assayed AD cases for TTR rs28933981, but the additional carrier was not related to the original AD‐affected cousin pair sharing the variant. Additional related AD‐affected carriers were observed among the 199 assayed cases for four of the nine assayed candidate variants. The other five newly identified AD case carriers were not related to any of the other carriers or to each other. Five AD‐affected carriers were observed for PELI3 rs148294193, including the original cousin pair. One of the newly identified case carriers was a cousin (and avuncular) to the original cousin pair (rare variant sharing P‐value = .006), and the other two newly identified case carriers were unrelated to all other carriers. Four AD‐affected carriers were observed for ABCA7 rs201665195 including the original cousin pair. One of the newly identified case carriers was a sibling to one of the original cousin carriers (rare variant sharing P‐value = .027), and one carrier was not related to any other carriers. Three AD‐affected carriers were observed for SNAP91 rs61729902, including the original cousin pair. The newly identified case carrier was a sibling of one of the original cousins (rare variant sharing P‐value = .027).

### Variant prioritization

Initial variant prioritization of rare variants shared by at least one cousin pair identified 400 rare variants spanning 470 genes. Ingenuity Variant Analysis subsequently prioritized 382 variants in 447 genes that were pathogenic, likely pathogenic, uncertain significance, or associated with a gain or loss of gene function. Of those shared rare variants, 117 variants in 138 genes had a biological interaction with genes implicated in AD. The numbers of variants included after each analysis are shown in Figure 1, and each variant is listed in File S1 in supporting information. That list of 117 variants was then used by four independent prioritization screens to identify 11 rare variants spanning 10 genes with the strongest support for increasing AD risk in these high‐risk pedigrees (see Table 1). Four rare variants previously associated with increased AD risk were identified using a literature search in Ingenuity Variant Analysis. The Knight ADRC CSF dataset identified two additional variants (Note: of the 117 prioritized variants, only 11 appeared in the CSF dataset) that were associated with increased risk for AB42, Tau, or Ptau in CSF. Additionally, the AD Risk Gradient identified three variants that were not present in the Wellderly dataset and were more prevalent in ADGC cases than ADGC controls. Finally, two variants were present in two independent high‐risk AD pedigrees.

Flowchart depicting variant prioritization. The number of rare variants and genes remaining at each level are shown. ADRC, Alzheimer's Disease Research Center; CSF, cerebrospinal fluid; MAF, minor allele frequency.

Prioritized variants

Notes: These variants are most likely to affect AD pathology in the high‐risk pedigrees. The accession number, affected gene, the Human Genome Variation Society (HGVS) variant annotation, impact on translation, and the final prioritization option that identified the variant are shown.

Abbreviations: AD, Alzheimer's disease; CSF, cerebrospinal fluid.

### Evidence for co‐segregation with AD risk in additional affected sampled relatives

Five candidate variants assayed for segregation evidence had no additional carriers beyond the original cousin pair, including: MUC16 rs191804178, NOTCH3 rs141402160, NOTCH3 rs140914494, FCHO1 rs147599881, and COX6A2 rs140129800. An additional carrier was identified among the 199 assayed AD cases for TTR rs28933981, but the additional carrier was not related to the original AD‐affected cousin pair sharing the variant. Additional related AD‐affected carriers were observed among the 199 assayed cases for four of the nine assayed candidate variants. The other five newly identified AD case carriers were not related to any of the other carriers or to each other. Five AD‐affected carriers were observed for PELI3 rs148294193, including the original cousin pair. One of the newly identified case carriers was a cousin (and avuncular) to the original cousin pair (rare variant sharing P‐value = .006), and the other two newly identified case carriers were unrelated to all other carriers. Four AD‐affected carriers were observed for ABCA7 rs201665195 including the original cousin pair. One of the newly identified case carriers was a sibling to one of the original cousin carriers (rare variant sharing P‐value = .027), and one carrier was not related to any other carriers. Three AD‐affected carriers were observed for SNAP91 rs61729902, including the original cousin pair. The newly identified case carrier was a sibling of one of the original cousins (rare variant sharing P‐value = .027).

### DISCUSSION

Analysis of exomes from 19 AD‐affected cousin pairs identified 400 shared rare candidate AD predisposition variants. Initial prioritization with Ingenuity Variant Analysis on likely pathogenicity and biological context reduced this list to 117 rare variants occurring in 138 genes. Further prioritization in one of four independent datasets, further prioritized 11 variants in 10 genes (Table 1). Four of these variants in three genes (ABCA7, TTR, and NOTCH3) represent replications of previous associations to AD, while the remaining eight variants are better classified as candidate variants still requiring validation, although each exhibited some level of replication through one of the four independent prioritization strategies. Finally, variants in ABCA7, SNAP91, and PELI3 showed significant evidence of segregation to other related AD cases in the high‐risk pedigrees in which they were identified.

Four identified rare variants were previously reported in the literature as associated with increasing AD risk (known AD risk variants). The first variant, NM_019112.3:c.302T > G (rs201665195), is found in the ABCA7 gene on chromosome 19. Unfortunately, this variant was not included in the ADGC validation datasets, so it could not be independently validated beyond the cousin pairs. Le Guennec et al.
46
 observed this variant, along with other rare ABCA7 variants, in AD cases and confirmed that rare, loss of function, and predicted damaging missense variants in ABCA7 are more common in patients with AD. Vardarajan et al.
47
 also identified this rare variant in two of their three late‐onset AD (LOAD) cohorts, and this variant was not seen in unaffected individuals. The cousin pairs analysis, and evidence for segregation, adds additional support to ABCA7 as a gene that impacts AD.

Variant NM_000371.3:c.416C > T (rs28933981) is located in the TTR gene on chromosome 18. Although this variant was sequenced in the ADGC validation datasets, it did not pass quality control thresholds, so it could not be used for validation. Sassi et al.
48
 observed this variant in 6 out of 332 AD cases (1.8% of cases), and this variant had a strong effect size (odds ratio = 6.19, 95% confidence interval = 1.099–63.091). It was also observed in two of their 676 cognitive‐normal control samples (0.30% of controls). TTR is known to be involved in amyloid beta (Aβ) catabolism and has no homologous proteins, suggesting that subtle changes to this protein could have strong functional implications.
48
 This analysis of high‐risk cousin pairs adds additional support to TTR as a gene that impacts AD.

Patel et al.
49
 reported missense mutations NM_000435.3:c.743G > C (rs141402160) and NM_000435:c.593C > T (rs140914494) in NOTCH3 using the pedigrees from this analysis combined with other population resources. However, these variants were not sequenced in the validation datasets. NOTCH3 plays a key role in neural development and is known to cause cerebral autosomal dominant arteriopathy with subcortical infarcts and leukoencephalopathy (CADASIL). The same region of NOTCH3 was also linked to AD in a Turkish family.
14
. NOTCH3 is known to bind PSEN1

50
, 
51
 and PSEN2.
50
 Furthermore, NOTCH3 had a P‐value of 1 × 10−4 from the gene‐based testing of AD risk in the ADGC dataset, implying that multiple variants in the gene may play a role in AD risk.

Two variants approached significance for increased AD risk in the Knight ADRC CSF dataset. A missense mutation in PELI3 (NM_001243135.1:c.115G > C; rs148294193) has a suggestive positive influence on AD risk (P‐value = .0553). PELI3 binds UBC

52
 and APP.
53

UBC is significantly downregulated in AD brains suggesting that decreased UBC function may be important in AD pathogenesis including increased neuronal death and non‐regulated APP production.
54

APP is a well‐studied AD risk gene, and many mutations in this gene are known to cause early‐onset AD. APP is cleaved into Aβ peptides, which are a major component of the amyloid plaques deposited in the brains of AD patients.
55

PELI3 encodes a scaffold protein and an intermediate signaling protein in the innate immune response pathway. Segregation of the PELI3 variant to an additional AD‐affected relative in the original high‐risk pedigree was observed.

A missense mutation NM_001161357.1:c.557G > A (rs147599881) in FCHO1 has an AD risk P‐value of .0624 in the CSF with a positive direction. FCHO1 also binds UBC

56
 and APP.
53

FCHO1 is involved in vesicle‐mediated transportation and clathrin‐mediated endocytosis.

Three variants were not present in the Wellderly dataset and were more prevalent in ADGC cases than ADGC controls. Missense variant NM_00124279.1:c.2113C > T (rs61729902) is located in SNAP91, which encodes a protein that binds UBC.
57

SNAP91 is involved in vesicle‐mediated transportation and clathrin‐mediated endocytosis. SNAP91 is a paralog of PICALM, a known top 10 AD risk gene,
58
 and they both have similarities between their functions of clathrin‐mediated endocytosis. The SNAP91 variant was observed in an additional AD case in the original high‐risk pedigree sequenced.

Another missense mutation, NM_005205.3:c.34T > G (rs140129800), was found in the COX6A2 gene, which encodes a protein that binds APP.
53

COX6A2 is a terminal enzyme in the mitochondrial respiratory chain and is part of the Kyoto Encyclopedia of Genes and Genomes (KEGG) pathway for AD.
59
, 
60
, 
61

Missense mutation NM_024690.2:c.10900C > T (rs191804178) in MUC16 affects binding of UBC,
62
 and MUC16 forms a protective mucous barrier on the apical surfaces of the epithelia. Mutations in MUC16 are associated with ovarian cancer, endometriosis, Pseudo‐Meigs syndrome, serous cystadenocarcinoma, and bronchogenic cysts. MUC16 has not previously been implicated in AD.

Three variants were observed to be shared by both members of the AD‐affected cousin pair in more than one high‐risk pedigree. Missense variant NM_145886.3:c.2044C > T (rs200290640) and splice site variant NM_145886.3:c.2042‐2A > G (rs199752248) both affect the PIDD1 gene, which is shown to bind APOE ε4.
63
 Neither of these variants were sequenced in the ADGC validation datasets. PIDD1 contains a death domain, interacts with other death domain proteins, and is suggested to be an effector of p53‐dependent apoptosis. Mutations in PIDD1 are associated with poikiloderma with neutropenia (PN) disorder that affects the skin and the immune system.

Rare variant‐sharing in AD‐affected relatives who are members of validated high‐risk pedigrees is central to prioritize candidate variants. Similar approaches capitalizing on shared genetics between related individuals in high‐risk pedigrees have been used successfully in Utah for decades to identify rare variants contributing to common diseases,
10
, 
64
, 
65
 as well as more recent adaptations.
12
, 
66
 This study design is limited by a relatively small available sample size (19 pedigrees with 38 index cases), the absence of other ethnicities besides Whites, and the fact that AD phenotyping was based solely on death certificate data. However, the approach raises statistical power by increasing the relative allele frequencies of rare variants, which allows a single pedigree in the sample set to identify rare candidate variants, a key advantage in the presence of locus heterogeneity. Replication in external datasets can be difficult to achieve, given the rare nature of such variants. However, the absence of these limitations in external dataset validations (i.e., if all variants were present in ADGC, Knight ADRC, or the Wellderly dataset) would likely lead to more candidate variants identified through this approach. The prioritization criteria were very conservative and report only the most supported variants that likely affect AD within these high‐risk pedigrees. Additional variants that may affect AD within these pedigrees may have been de‐prioritized because of the stringent nature of the analysis. All genetic variants identified at each prioritization level are reported for future research to assess the relative support of each variant. The most supported variants were also assessed for co‐segregation using a small number of additional AD‐affected relatives with the AD phenotype within the original high‐risk pedigrees. Although the UPDB population genealogy data include the possibility of undocumented relationships among subjects and inadequate phenotyping strategies, previous research shows that founder effects and inbreeding are no greater in the UPDB than in the general population,
17
, 
18
 and the UPDB has been used extensively in previous disease studies.

The presence of rare variants identified here may explain the prevalence of AD mortality in 19 of the 36 AD‐affected individuals from high‐risk pedigrees (see Table S2). The excess AD mortality observed in the remaining individuals might be due to complex interactions, heterogeneity in the pedigree, misdiagnosis of AD, non‐coding variants, or variants that were removed due to stringent prioritization criteria. Because the initial prioritization included only variants in genes known or predicted to affect AD pathology, the analysis is dependent on current understanding of AD pathology and may not encompass all disease‐causing variants. However, the purpose of this study was to identify highly supported candidate variants associated with AD risk. To that end, these analyses provide additional support to previous studies that implicate ABCA7 and TTR with AD mortality. Because many known AD risk variants were prioritized in these analyses, other prioritized variants also likely affect AD mortality, and all 11 variants spanning 9 genes should be prioritized in future analyses. These outcomes indicate that a high‐risk pedigree approach can achieve sufficient power to detect rare variants, particularly when coupled with external datasets that contain meaningful data about disease risk.

### Variants observed in more than one AD‐affected cousin pair

Three variants were observed to be shared by both members of the AD‐affected cousin pair in more than one high‐risk pedigree. Missense variant NM_145886.3:c.2044C > T (rs200290640) and splice site variant NM_145886.3:c.2042‐2A > G (rs199752248) both affect the PIDD1 gene, which is shown to bind APOE ε4.
63
 Neither of these variants were sequenced in the ADGC validation datasets. PIDD1 contains a death domain, interacts with other death domain proteins, and is suggested to be an effector of p53‐dependent apoptosis. Mutations in PIDD1 are associated with poikiloderma with neutropenia (PN) disorder that affects the skin and the immune system.

### Strengths and limitations

Rare variant‐sharing in AD‐affected relatives who are members of validated high‐risk pedigrees is central to prioritize candidate variants. Similar approaches capitalizing on shared genetics between related individuals in high‐risk pedigrees have been used successfully in Utah for decades to identify rare variants contributing to common diseases,
10
, 
64
, 
65
 as well as more recent adaptations.
12
, 
66
 This study design is limited by a relatively small available sample size (19 pedigrees with 38 index cases), the absence of other ethnicities besides Whites, and the fact that AD phenotyping was based solely on death certificate data. However, the approach raises statistical power by increasing the relative allele frequencies of rare variants, which allows a single pedigree in the sample set to identify rare candidate variants, a key advantage in the presence of locus heterogeneity. Replication in external datasets can be difficult to achieve, given the rare nature of such variants. However, the absence of these limitations in external dataset validations (i.e., if all variants were present in ADGC, Knight ADRC, or the Wellderly dataset) would likely lead to more candidate variants identified through this approach. The prioritization criteria were very conservative and report only the most supported variants that likely affect AD within these high‐risk pedigrees. Additional variants that may affect AD within these pedigrees may have been de‐prioritized because of the stringent nature of the analysis. All genetic variants identified at each prioritization level are reported for future research to assess the relative support of each variant. The most supported variants were also assessed for co‐segregation using a small number of additional AD‐affected relatives with the AD phenotype within the original high‐risk pedigrees. Although the UPDB population genealogy data include the possibility of undocumented relationships among subjects and inadequate phenotyping strategies, previous research shows that founder effects and inbreeding are no greater in the UPDB than in the general population,
17
, 
18
 and the UPDB has been used extensively in previous disease studies.

### Conclusion

The presence of rare variants identified here may explain the prevalence of AD mortality in 19 of the 36 AD‐affected individuals from high‐risk pedigrees (see Table S2). The excess AD mortality observed in the remaining individuals might be due to complex interactions, heterogeneity in the pedigree, misdiagnosis of AD, non‐coding variants, or variants that were removed due to stringent prioritization criteria. Because the initial prioritization included only variants in genes known or predicted to affect AD pathology, the analysis is dependent on current understanding of AD pathology and may not encompass all disease‐causing variants. However, the purpose of this study was to identify highly supported candidate variants associated with AD risk. To that end, these analyses provide additional support to previous studies that implicate ABCA7 and TTR with AD mortality. Because many known AD risk variants were prioritized in these analyses, other prioritized variants also likely affect AD mortality, and all 11 variants spanning 9 genes should be prioritized in future analyses. These outcomes indicate that a high‐risk pedigree approach can achieve sufficient power to detect rare variants, particularly when coupled with external datasets that contain meaningful data about disease risk.

Partial support for all data sets within the Utah Population Database (UPDB) was provided by Huntsman Cancer Institute, University of Utah and the Huntsman Cancer Institute's Cancer Center Support grant, P30 CA42014 from National Cancer Institute. LACA receives partial support from the Huntsman Cancer Institute's Cancer Center Support grant, P30 CA42014 from National Cancer Institute.

Acknowledgment is made to the donors of Alzheimer's Disease Research, a program of the BrightFocus Foundation for support of this research through grant #A2020118F (PI: Miller).

This study is part of the NHLBI Grand Opportunity Exome Sequencing Project (GO‐ESP). Funding for GO‐ESP was provided by NHLBI grants RC2 HL103010 (HeartGO), RC2 HL102923 (LungGO), and RC2 HL102924 (WHISP). The exome sequencing was performed through NHLBI grants RC2 HL102925 (BroadGO) and RC2 HL102926 (SeattleGO). HeartGO gratefully acknowledges the following groups and individuals who provided biological samples or data for this study. DNA samples and phenotypic data were obtained from the following studies supported by the NHLBI: the Atherosclerosis Risk in Communities (ARIC) study, the Coronary Artery Risk Development in Young Adults (CARDIA) study, Cardiovascular Health Study (CHS), the Framingham Heart Study (FHS), the Jackson Heart Study (JHS), and the Multi‐Ethnic Study of Atherosclerosis (MESA).

Data from Alzheimer's Disease Genetics Consortium (ADGC) were appropriately downloaded from dbGaP (accession: phs000372.v1.p1). We acknowledge the contributions of the members of the ADGC listed in Appendix: Alzheimer's Disease Genetics Consortium Collaborators.

### Supporting information

Supporting Information

Click here for additional data file.

Supporting Information

Click here for additional data file.

Supporting Information

Click here for additional data file.

