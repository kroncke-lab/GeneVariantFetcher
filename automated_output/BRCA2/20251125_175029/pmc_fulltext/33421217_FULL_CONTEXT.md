# MAIN TEXT

## Pathogenicity assessment of variants for breast cancer susceptibility genes based on BRCAness of tumor sample

### Abstract

AbstractGenes involved in the homologous recombination repair pathway—as exemplified by BRCA1, BRCA2, PALB2, ATM, and CHEK2—are frequently associated with hereditary breast and ovarian cancer syndrome. Germline mutations in the loci of these genes with loss of heterozygosity or additional somatic truncation at the WT allele lead to the development of breast cancers with characteristic clinicopathological features and prominent genomic features of homologous recombination deficiency, otherwise referred to as “BRCAness.” Although clinical genetic testing for these and other genes has increased the chances of identifying pathogenic variants, there has also been an increase in the prevalence of variants of uncertain significance, which poses a challenge to patient care because of the difficulties associated with making further clinical decisions. To overcome this challenge, we sought to develop a methodology to reclassify the pathogenicity of these unknown variants using statistical modeling of BRCAness. The model was developed with Lasso logistic regression by comparing 116 genomic attributes derived from 37 BRCA1/2 biallelic mutant and 32 homologous recombination‐quiescent breast cancer exomes. The model showed 95.8% and 86.7% accuracies in the training cohort and The Cancer Genome Atlas validation cohort, respectively. Through application of the model for variant reclassification of homologous recombination‐associated hereditary breast and ovarian cancer causal genes and further assessment with clinicopathological features, we finally identified one likely pathogenic and five likely benign variants. As such, the BRCAness model developed from the tumor exome was robust and provided a reasonable basis for variant reclassification.

### INTRODUCTION

Hereditary breast and ovarian cancer syndrome exhibits an autosomal dominant trait and is caused by a defect in susceptibility genes with tumor suppressor function. This includes the highly penetrant BRCA1, BRCA2, and PALB2, the syndromic TP53, CDH1, and NF1, and the moderately penetrant CHEK2, ATM, BARD1, and RAD51C.
1
 Recent advances in sequencing technology have enabled multigene‐panel genetic testing to identify causal genetic variants for patients with HBOC.
2
 However, these genetic tests have also resulted in the identification of numerous variants for which, in most cases, the pathogenicity and clinical significance are unknown. Such VUS pose a challenge in clinical practice because of the difficulties associated with deciding whether further medical care is warranted; for instance, surveillance, drugs, prophylactic surgery, or risk assessment for relatives.
3
 It would therefore be beneficial to develop a methodology to reclassify these VUS and correctly interpret their pathogenicity.

A subset of the causal genes associated with HBOC syndrome is known to be involved in the HR repair pathway, which repairs DNA double‐stranded breaks in the cell. Complete loss of function of a HR‐related gene—achieved by WT allele inactivation through LOH by copy loss or AST
4
—leads to HRD, which, in turn, causes characteristic chromosomal rearrangement and nucleotide substitution.
5
, 
6
, 
7
, 
8
 Such genomic features have been referred to as “BRCAness,” as the features largely resemble those in tumors with a mutation in the prototypical BRCA1 or BRCA2 gene.
9
, 
10
, 
11
 Identifying BRCAness is considered critical because it has the potential to be a biomarker, like BRCA1 or BRCA2 mutation, for predicting drug sensitivity to poly (ADP‐ribose) polymerase inhibitors or platinum.
12
 Thus far, numerous genomic features representative of BRCAness have been reported, including: BRCA1‐mutant associated copy number aberrations
13
; alterations in the numbers of intermediate‐size LOH (HRD‐LOH),
14
 subtelomeric regions with allelic imbalance (HRD‐NtAI),
15
 or large‐scale state transitions (HRD‐LST)
16
; the HRD score, as determined by the sum of HRD‐LOH, HRD‐NtAI, and HRD‐LST
17
; the COSMIC mutational signature #3
18
; and the LOH state in the BRCA1 or BRCA2 locus.
19
 Indeed, previous clinical studies have confirmed the utility of the HRD score in the prediction of platinum sensitivity.
20
, 
21
, 
22
, 
23
 Furthermore, statistical models encompassing several genomic variables detected by WGS can be used to identify BRCAness in the cell.
7
 Considering the link between HR‐related genomic changes and BRCAness, there could be the potential to use BRCAness for the reclassification of VUS identified through screens of HBOC genetic testing. However, a limited number of studies thus far have pointed to or investigated this potential utility. It remains to be elucidated whether a BRCAness biomarker can be used to interpret VUS pathogenicity, particularly variants of non‐BRCA1/2 HBOC‐causative genes.

Here, we developed a robust statistical model to estimate BRCAness in a breast cancer sample derived from 116 genomic variables obtained through WES. We validated the accuracy and versatility of the model using breast cancer data from TCGA as an independent external cohort. We then applied this model for the reclassification of VUS of HR‐associated HBOC causal genes and identified one LP and five LB variants. The current study shows that our BRCAness model can contribute to the medical care of HBOC patients by lowering the number of VUS and providing a better pathogenicity classification.

### MATERIALS AND METHODS

Samples for the training cohort were obtained by WES, whereas data for the validation cohort were obtained from TCGA (https://portal.gdc.cancer.gov/). For each dataset, we analyzed the germline and somatic genomic aberrations. Pathogenicity interpretations of germline variants were undertaken according to the ACMG–AMP guidelines
24
 using the pipeline previously described.
25
 We developed the BRCAness model by Lasso logistic regression using 116 genomic variables. See Document S1 for details regarding sample data collection and the experimental, bioinformatics, and statistical analyses.

### RESULTS

During breast cancer tumorigenesis in patients with HBOC syndrome, the WT allele of an HR‐associated gene bearing a germline mutation (HBOC‐HR gene) is inactivated as a “second hit” through LOH or AST.
6
 The tumor cells of such breast cancer typically present with a set of genomic features, frequently referred to as BRCAness, derived from HRD
11
; this is also linked with characteristic clinicopathological features.
13
 Because the presence of a pathogenic variant in one HBOC‐HR gene with LOH or AST is tightly linked with BRCAness and typical clinicopathological features, we considered that both BRCAness and clinicopathology could be used to estimate whether a VUS is pathogenic or benign (Figure 1A).

Strategy and outline of the current study. A, Tumorigenic process and pathogenicity assessment for variants of uncertain significance (VUS). Schematic representation of the tumorigenesis of homologous recombination (HR)‐deficient breast cancer is shown with a research strategy, in which “BRCAness” and several clinicopathological features of a tumor sample are used to assess the pathogenicity of a VUS. In the current study, BRCAness is defined as a set of characteristic genomic features in breast cancer samples with HR deficiency, originally described in tumors from patients with germline BRCA1 or BRCA2 mutations.
11
 B, BRCAness and its application for assessing VUS pathogenicity. The BRCAness model was built by comparing the genomic features of tumor samples with BRCA1 or BRCA2 biallelic loss (BAL) with those of HR‐quiescent tumors. We then predicted the BRCAness of the samples: these included tumors with a VUS in any of the germline BRCA1, BRCA2, or other hereditary breast and ovarian cancer‐HR repair‐associated (HBOC‐HR) genes. VUS pathogenicity was assessed based on BRCAness and clinicopathological features to identify potentially likely pathogenic (PLP) and potentially likely benign (PLB) variants. The likely pathogenic (LP) and likely benign (LB) variants were finally classified by interrogating previous reports including functional assays or case–control studies. AST, additional somatic truncation; AUC, area under the curve

Figure 1B shows the strategy of the current study. First, we developed a BRCAness model as the training step using exome data generated at the Japanese Foundation for Cancer Research (JFCR cohort). For this, we extracted “BRCAness” using Lasso logistic regression analysis by comparing the genomic features of 37 tumor samples with BRCA1/2 BAL with those of 32 tumor samples without any genomic or epigenetic alterations in HBOC‐HR genes (ie, HR‐quiescent tumors); these were used as positive and negative controls, respectively (Figures 1B and S1). Second, the model was validated with TCGA breast cancer exome data (TCGA cohort) using 18 positive and 95 negative controls (Figures 1B and S1). Third, relying on the model, BRCAness was subsequently predicted for 413 miscellaneous samples from the JFCR and TCGA cohorts. These samples were not used to develop the statistical model, had at least one miscellaneous (germline or somatic, genetic or epigenetic) alteration in an HBOC‐HR gene, and included 251 samples with VUS (Figures 1B and S1). BRCAness caused by a germline variant in an HBOC‐HR gene has been observed not only in hereditary but also sporadic disease.
26
 Therefore, we included all breast cancer samples with exome data in the current study. By comparing BRCAness probability with clinicopathological feature(s), the subsequent assessment of VUS pathogenicity from 269 unique VUS revealed five PLP and 35 PLB variants. After interrogating previously published results from functional assays or case–control studies, we finally classified one LP and five LB variants.

Of the 175 and 735 breast cancer exomes from the JFCR and TCGA cohorts, respectively, one JFCR and 314 TCGA data were excluded by stringent quality control criteria, as follows. First, we excluded 214 TCGA data that underwent WGA before exome library preparation because statistically significant differences were noted between data with and without WGA in several genomic features, such as the number of abnormal gain segments and the number of indels and HRD‐NtAI (data not shown); this might indicate that WGA introduces some artifactual genomic features into the exome. No JFCR samples were subjected to WGA before sequencing. We then excluded one JFCR and 100 TCGA data because of low tumor content (less than 0.3) or low detectability of SNVs/indels in a sample (less than 20) (Figure S1A). The remaining 174 JFCR and 421 TCGA samples were subsequently genotyped for germline or somatic, genetic or epigenetic alterations in HBOC‐HR genes (Figure S1B; see also Doc. S1). Pathogenicity of germline variants was interpreted according to the ACMG‐AMP guidelines.
24
, 
25

Thirty‐seven JFCR and 18 TCGA samples had BRCA1/2 BAL (germline pathogenic variants with LOH or AST); 32 JFCR and 95 TCGA samples were HR‐quiescent tumors (no germline or somatic, genetic or epigenetic alterations among any of the HBOC‐HR genes in the exome analyses) (Figure S1C,D). These samples were assigned as positive and negative controls at the training and validation steps during development of the BRCAness model (Figures 1B and S1D). The remaining 105 JFCR and 308 TCGA samples were miscellaneous in genotype: the samples had at least one germline or somatic, genetic or epigenetic change in an HBOC‐HR gene. These samples were combined and used for BRCAness prediction to link VUS with BRCAness and clinicopathological features (Figure S1C,D).

Figure S2 illustrates the genomic attributes derived from exome analyses of 37 BRCA1/2 BAL tumors, 32 HR‐quiescent tumors, and 105 miscellaneous tumors in the JFCR cohort, together with data for hormonal and subtype status. Consistent with previous reports,
6
, 
8
, 
27
, 
28
 clear distinctions were observed in the genomic attributes between BRCA1/2 BAL tumors and HR‐quiescent tumors; these distinctions included the numbers of SNVs, indels, and abnormal copy number segments, the frequency of nucleotide substitutions and COSMIC mutational signatures, and HRD scores. The numbers of SNVs, indels, and abnormal copy number segments were significantly higher in tumors with BRCA1/2 BAL tumors than in HR‐quiescent tumors (all P values were below .0001; data not shown). Whereas COSMIC mutational signature #3 (“BRCA” signature) was a prominent characteristic of BRCA1/2 BAL tumors, most of the HR‐quiescent tumors had signature #1 (“aging” signature) and a fraction had signatures #2 and #13 (“APOBEC” signatures [apolipoprotein B mRNA editing catalytic polypeptide‐like]) (Figures 2 and S2).
6
, 
8
 The tumor genotypes were also associated with hormonal receptor or subtype status: many BRCA1/2 BAL tumors were triple‐negative subtype, whereas HR‐quiescent tumors were luminal subtype and mostly positive for estrogen and progesterone receptors (Figure S2). Importantly, miscellaneous tumors were characterized by an intermediate set of features along a broad spectrum (Figure S2). These genomic features were reproducibly observed in the samples from TCGA cohort; however, there were different proportions of BRCA1/2 BAL, HR‐quiescent, and miscellaneous genomic changes (Figures 2 and S3).

BRCAness probabilities pertaining to breast cancer homologous recombination (HR) gene status. A, Violin plot diagram of the HR deficiency (HRD)‐LOH score, proportion of COSMIC signature #3, and number of indels in the Japanese Foundation for Cancer Research (JFCR) cohort. B, Violin plot diagram of HRD‐LOH score, proportion of COSMIC signature #3, and number of indels in The Cancer Genome Atlas (TCGA) cohort. In (A) and (B), these variables were selected in the final BRCAness model by Lasso logistic regression. C, Violin plot diagram of BRCAness probabilities in the JFCR cohort. D, Violin plot diagram of BRCAness probabilities in TCGA cohort. The threshold of BRCAness probability was 0.534 to distinguish genomic features of BRCA1/2‐biallelic loss (BAL) tumors from HR‐quiescent tumors

Clinicopathological features of BRCA1/2 BAL and HR‐quiescent tumor samples in the JFCR cohort were also highly consistent with those reported in previous studies (Figure S4A).
25

BRCA1/2 BAL tumors were characterized by juvenile onset, triple‐negative subtype, bilateral and/or multiple occurrence, solid tubular histology, and nuclear grade 3; these features were more salient in BRCA1 BAL tumors than in BRCA2 BAL tumors (Figures S2 and S4A). Conversely, older age of onset, luminal subtype, unilateral and single occurrence, papillotubular or scirrhous histology, and lower nuclear grade were more representative features of HR‐quiescent tumors (Figures S2 and S4A). Although concomitant clinical information was limited in TCGA cohort data, a similar pattern was observed in terms of age at diagnosis and molecular subtype (Figure S4B). These strong associations between tumor genotype and genomic/clinicopathological features not only validates our exome‐based genomic characterization but also provides a solid basis for the development of a statistical model reliant on tumor genotype and genomic attributes.

We selected three variables and developed a BRCAness model using Lasso logistic regression analysis as a method of machine learning (Figure 2A,B). The Lasso model was superior to the conventional multivariate logistic regression in terms of the mean accuracies of both the training and validation data (Table 1). We then applied the model to estimate the BRCAness status of 413 miscellaneous tumor samples having variable genomic/epigenetic changes in other HBOC‐HR genes (Figure 2C,D). Document S2 outlines the details regarding the development and application of the BRCAness model.

Statistics of Lasso and conventional multivariate logistic regression models to determine BRCAness of breast cancer tumor samples

Lasso logistic regression

Multivariate logistic regression

Abbreviations: HRD, homologous recombination deficiency; LOH, loss of heterozygosity; SNV, single nucleotide variant.

Figure 3 shows the relationship between Lasso BRCAness probability and the type of genomic or epigenetic change in HBOC‐HR genes, along with the associated clinicopathological parameters in 413 breast tumor samples with miscellaneous genotypes. Based on the BRCAness predictive model, samples were subdivided into BRCAnesshigh (n = 97) and BRCAnesslow (n = 316) probability tumors, with a cut‐off set to 0.534 (Figure 3A). Mann‐Whitney tests, comparing samples with and without a genomic or epigenetic change in HBOC‐HR genes, showed that probability was significantly correlated with BRCA1 or RAD51C promoter hypermethylation, germline pathogenic variant with LOH/AST, or somatic homozygous deletion (Figure 3B, Table S1), consistent with previous observations.
18
 These results validated the use of our analytical methodologies and the predicted results that these genomic or epigenetic alterations were most likely BRCAness‐causing (Table S1). Worth noting, one tumor with a germline mutation in BRCA1 c.5193 + 3insT with LOH together with a PTEN somatic homozygous deletion was categorized as BRCAnesshigh in the analysis. The BRCA1 variant was located three nucleotides away from the 3′ splice site of exon 18 in the BRCA1 gene (NM_007294) and therefore was initially classified as VUS in our informatic pipeline: a splice variant is defined only as a truncation when alterations are present within two bases of the exon‐intron boundary. Nevertheless, relying on the Myriad‐Falco database, this variant was classified as pathogenic. The prediction was also validated by the clinicopathological data supplied for the samples with miscellaneous genotypes. As seen in Figure 3C, triple‐negative subtype, nuclear grade 3, and solid tubular histology were highly enriched in BRCAnesshigh tumors, which are reminiscent of BRCA1/2 BAL tumors (Figure S4). These findings indicate tight associations between BRCAness probability and the types of genomic or epigenetic alterations that occur, along with relevant clinicopathological features.

BRCAness probability, genomic/epigenetic alterations in homologous recombination (HR) genes, and clinicopathological features of breast cancer tumors with miscellaneous genotypes, excluding BRCA1/2 biallelic loss and HR‐quiescence. A, BRCAness probability and status of genomic/epigenetic alterations in BRCA1/2 or hereditary breast and ovarian cancer‐HR repair pathway‐associated (HBOC‐HR) genes. Top panel, bar plots for BRCAness probability. Bottom panel, presence or absence of hypermethylation for BRCA1 or RAD51C promoters, followed by presence or absence of germline pathogenic variants, germline variants of uncertain significance (VUS), somatic truncations, somatic VUS with (+) or without (–) LOH or additional somatic truncation (AST) in BRCA1/2 or HBOC‐HR genes, and finally, presence or absence of homozygous deletion (copy number [CN] = 0) and CN loss (CN = 1) in BRCA1/2 or HBOC‐HR genes. The samples are sorted according to BRCAness probability and divided into BRCAnesshigh and BRCAnesslow tumors using a threshold of 0.534. B, BRCAness probabilities in tumors with BRCA1 or RAD51C promoter hypermethylation, in tumors with germline pathogenic variants with LOH or AST, and in tumors with somatic homozygous deletion in BRCA1/2 or HBOC‐HR genes. P values were computed using the Mann‐Whitney U test. C, Clinicopathological features of BRCAnesshigh and BRCAnesslow tumors. Left panel, intrinsic subtype. Middle panel, nuclear grade. Right panel, tumor histology. P values were computed using Fisher’s exact test. DCIS, ductal carcinoma in situ; H, hypermethylated; NH, not hypermethylated; MT, mutated; TN, triple‐negative

The validity of the BRCAness prediction prompted us to assess VUS pathogenicity by examining the status of the variant and the tumor clinicopathological features (Figure 4). Here, we assumed that WT allelic inactivation by LOH or AST is essential for manifesting full‐blown BRCAness. Using a tree diagram and relying on variant status, we show that 337 germline VUS (n = 337; 269 unique variants) could be classified into variants with LOH/AST (n = 77) and without LOH/AST (n = 260); those VUS lacking LOH/AST were dismissed from further assessment. Among the 77 variants with LOH/AST, 32 and 45 VUS were classified as BRCAnesshigh and BRCAnesslow tumors. These 32 BRCAnesshigh variants could be further subdivided by the absence or presence of concomitant BRCAness‐causing genomic or epigenetic alterations (Figure 4). A VUS (with LOH/AST) classified as a BRCAnesshigh tumor but lacking a concomitant BRCAness‐causing event is considered possibly pathogenic and was designated “category 1”. Comparatively, a VUS (with LOH/AST) classified as a BRCAnesshigh tumor and accompanied by a BRCAness‐causing event—that would have been assigned as “category 2”—is considered indeterminable, as the BRCAness‐causing event masks the effect of the VUS. Furthermore, we considered a VUS with LOH/AST in a BRCAnesslow tumor (“category 3”) as possibly benign, as a VUS with a homozygous state does not cause HRD (Figure 4).

Pathogenicity of variants of uncertain significance (VUS) with BRCAness probability and clinicopathological features. Breast cancer tumors with VUS with LOH or additional somatic truncations (AST) were classified into BRCAnesshigh and BRCAnesslow tumor categories according to BRCAness probability. BRCAnesshigh tumors were subdivided into tumors with and without concomitant BRCAness‐causing events: BRCA1 or RAD51C promoter hypermethylation, germline pathogenic variants with LOH/AST, or somatic homozygous deletions in BRCA1/2 or hereditary breast and ovarian cancer‐homologous recombination repair pathway‐associated genes. Germline VUS were then assessed against clinicopathological features that are characteristic of tumors with BRCA1/2 biallelic loss. One representative variant per category is shown, with tumor characteristics as examples. PLB, potentially likely benign; PLP, potentially likely pathogenic

In assessing the clinicopathological attributes, we used three features to distinguish between BRCAnesshigh and BRCAnesslow tumors: internal subtype, histology, and nuclear grade. Whereas triple‐negative subtype, solid tubular histology, and nuclear grade 3 were characteristic of BRCAnesshigh tumors, other phenotypes were linked with BRCAnesslow tumors. When a category 1 VUS consistently had clinicopathological features characteristic of BRCAnesshigh tumors, the sample was classified as PLP. A VUS of category 3 with one or more BRCAnesslow tumor characteristics (ie, not BRCAnesshigh characteristics) was interpreted as PLB (Figure 4).

Using the tree classification, we identified 19 (18 unique), 10 (10 unique), and 43 (41 unique) category 1, category 2, and category 3 variants from among the 77 VUS with LOH/AST from the original categorization. Five (two unique) variants could not be classified according to the tree diagram because they appeared in both BRCAnesshigh and BRCAnesslow tumor classes, and, instead, were assigned as “conflicting” (Table S2). Ten cases had more than one VUS (total 21 variants) accompanied by LOH/AST: three cases had six variants, two cases had five variants, and five cases had 10 variants; these were assigned to categories 1, 2, and 3, respectively (Table S2). Of note, there are three important considerations concerning these assignments: (i) the BRCAness‐causing variant cannot be inferred when two or more category 1 variants simultaneously reside in a tumor; (ii) in contrast, regardless of the number of VUS in a case, the category 2 variant is indeterminable because of the concomitant BRCAness‐causing event; and (iii) all category 3 variants are considered benign as none of them causes BRCAness. Based on these points, when there were two or more variants with LOH/AST in a tumor, the variants in categories 1 and 2 were reclassified as VUS, whereas category 3 variants were classified as PLB (Table S2). Further assessments in terms of the clinicopathological features identified six (five unique) and 36 (35 unique) VUS as PLP and PLB, respectively (Table S2). Among them, one category 1 (BLM c.T11C p.V4A) and two category 3 (ATM c.T7912G p.W2638G and CHEK2 c.G683A p.S228N) variants were recurrently detected.

We next reassessed the identified variants using previous reports comprising functional assays or case–control studies. The pathogenicity of BRCA1 c.T4951C p.S1651P as a PLP variant was supported by a functional assessment based on saturation genome editing
29
 (Table 2). A large case–control study of Japanese patients with breast cancer (7051 cases and 11 241 controls)
30
 revealed a similar prevalence of ATM c.C2879A p.P960H, BRCA2 c.A125G p.Y42C, and PALB2 c.G2509A p.E837K among cases and controls, indicating that these variants are not pathogenic. PALB2 c.C2135T p.A712V was also found to be benign in an Australian case–control study (1996 cases and 1998 controls).
31

BRCA1 c.G4535T p.S1512I and BRCA2 c.A125G p.Y42C showed comparable HR repair
32
, 
33
 and cell viability
34
, 
35
, 
36
 to that of WT proteins in multiple functional assessments, indicative of their nonpathogenicity. From these comparisons, we can conclude that our BRCAness model was useful in identifying LP and LB variants from VUS.

Proposed pathogenicity of variants of uncertain significance for breast cancer susceptibility genes

Abbreviations: LB, likely benign; LP, likely pathogenic.

### Research strategy and study design

During breast cancer tumorigenesis in patients with HBOC syndrome, the WT allele of an HR‐associated gene bearing a germline mutation (HBOC‐HR gene) is inactivated as a “second hit” through LOH or AST.
6
 The tumor cells of such breast cancer typically present with a set of genomic features, frequently referred to as BRCAness, derived from HRD
11
; this is also linked with characteristic clinicopathological features.
13
 Because the presence of a pathogenic variant in one HBOC‐HR gene with LOH or AST is tightly linked with BRCAness and typical clinicopathological features, we considered that both BRCAness and clinicopathology could be used to estimate whether a VUS is pathogenic or benign (Figure 1A).

Strategy and outline of the current study. A, Tumorigenic process and pathogenicity assessment for variants of uncertain significance (VUS). Schematic representation of the tumorigenesis of homologous recombination (HR)‐deficient breast cancer is shown with a research strategy, in which “BRCAness” and several clinicopathological features of a tumor sample are used to assess the pathogenicity of a VUS. In the current study, BRCAness is defined as a set of characteristic genomic features in breast cancer samples with HR deficiency, originally described in tumors from patients with germline BRCA1 or BRCA2 mutations.
11
 B, BRCAness and its application for assessing VUS pathogenicity. The BRCAness model was built by comparing the genomic features of tumor samples with BRCA1 or BRCA2 biallelic loss (BAL) with those of HR‐quiescent tumors. We then predicted the BRCAness of the samples: these included tumors with a VUS in any of the germline BRCA1, BRCA2, or other hereditary breast and ovarian cancer‐HR repair‐associated (HBOC‐HR) genes. VUS pathogenicity was assessed based on BRCAness and clinicopathological features to identify potentially likely pathogenic (PLP) and potentially likely benign (PLB) variants. The likely pathogenic (LP) and likely benign (LB) variants were finally classified by interrogating previous reports including functional assays or case–control studies. AST, additional somatic truncation; AUC, area under the curve

Figure 1B shows the strategy of the current study. First, we developed a BRCAness model as the training step using exome data generated at the Japanese Foundation for Cancer Research (JFCR cohort). For this, we extracted “BRCAness” using Lasso logistic regression analysis by comparing the genomic features of 37 tumor samples with BRCA1/2 BAL with those of 32 tumor samples without any genomic or epigenetic alterations in HBOC‐HR genes (ie, HR‐quiescent tumors); these were used as positive and negative controls, respectively (Figures 1B and S1). Second, the model was validated with TCGA breast cancer exome data (TCGA cohort) using 18 positive and 95 negative controls (Figures 1B and S1). Third, relying on the model, BRCAness was subsequently predicted for 413 miscellaneous samples from the JFCR and TCGA cohorts. These samples were not used to develop the statistical model, had at least one miscellaneous (germline or somatic, genetic or epigenetic) alteration in an HBOC‐HR gene, and included 251 samples with VUS (Figures 1B and S1). BRCAness caused by a germline variant in an HBOC‐HR gene has been observed not only in hereditary but also sporadic disease.
26
 Therefore, we included all breast cancer samples with exome data in the current study. By comparing BRCAness probability with clinicopathological feature(s), the subsequent assessment of VUS pathogenicity from 269 unique VUS revealed five PLP and 35 PLB variants. After interrogating previously published results from functional assays or case–control studies, we finally classified one LP and five LB variants.

### Samples and genomic changes

Of the 175 and 735 breast cancer exomes from the JFCR and TCGA cohorts, respectively, one JFCR and 314 TCGA data were excluded by stringent quality control criteria, as follows. First, we excluded 214 TCGA data that underwent WGA before exome library preparation because statistically significant differences were noted between data with and without WGA in several genomic features, such as the number of abnormal gain segments and the number of indels and HRD‐NtAI (data not shown); this might indicate that WGA introduces some artifactual genomic features into the exome. No JFCR samples were subjected to WGA before sequencing. We then excluded one JFCR and 100 TCGA data because of low tumor content (less than 0.3) or low detectability of SNVs/indels in a sample (less than 20) (Figure S1A). The remaining 174 JFCR and 421 TCGA samples were subsequently genotyped for germline or somatic, genetic or epigenetic alterations in HBOC‐HR genes (Figure S1B; see also Doc. S1). Pathogenicity of germline variants was interpreted according to the ACMG‐AMP guidelines.
24
, 
25

Thirty‐seven JFCR and 18 TCGA samples had BRCA1/2 BAL (germline pathogenic variants with LOH or AST); 32 JFCR and 95 TCGA samples were HR‐quiescent tumors (no germline or somatic, genetic or epigenetic alterations among any of the HBOC‐HR genes in the exome analyses) (Figure S1C,D). These samples were assigned as positive and negative controls at the training and validation steps during development of the BRCAness model (Figures 1B and S1D). The remaining 105 JFCR and 308 TCGA samples were miscellaneous in genotype: the samples had at least one germline or somatic, genetic or epigenetic change in an HBOC‐HR gene. These samples were combined and used for BRCAness prediction to link VUS with BRCAness and clinicopathological features (Figure S1C,D).

### Genomic and clinicopathological features in tumors with BRCA1/2 biallelic loss

Figure S2 illustrates the genomic attributes derived from exome analyses of 37 BRCA1/2 BAL tumors, 32 HR‐quiescent tumors, and 105 miscellaneous tumors in the JFCR cohort, together with data for hormonal and subtype status. Consistent with previous reports,
6
, 
8
, 
27
, 
28
 clear distinctions were observed in the genomic attributes between BRCA1/2 BAL tumors and HR‐quiescent tumors; these distinctions included the numbers of SNVs, indels, and abnormal copy number segments, the frequency of nucleotide substitutions and COSMIC mutational signatures, and HRD scores. The numbers of SNVs, indels, and abnormal copy number segments were significantly higher in tumors with BRCA1/2 BAL tumors than in HR‐quiescent tumors (all P values were below .0001; data not shown). Whereas COSMIC mutational signature #3 (“BRCA” signature) was a prominent characteristic of BRCA1/2 BAL tumors, most of the HR‐quiescent tumors had signature #1 (“aging” signature) and a fraction had signatures #2 and #13 (“APOBEC” signatures [apolipoprotein B mRNA editing catalytic polypeptide‐like]) (Figures 2 and S2).
6
, 
8
 The tumor genotypes were also associated with hormonal receptor or subtype status: many BRCA1/2 BAL tumors were triple‐negative subtype, whereas HR‐quiescent tumors were luminal subtype and mostly positive for estrogen and progesterone receptors (Figure S2). Importantly, miscellaneous tumors were characterized by an intermediate set of features along a broad spectrum (Figure S2). These genomic features were reproducibly observed in the samples from TCGA cohort; however, there were different proportions of BRCA1/2 BAL, HR‐quiescent, and miscellaneous genomic changes (Figures 2 and S3).

BRCAness probabilities pertaining to breast cancer homologous recombination (HR) gene status. A, Violin plot diagram of the HR deficiency (HRD)‐LOH score, proportion of COSMIC signature #3, and number of indels in the Japanese Foundation for Cancer Research (JFCR) cohort. B, Violin plot diagram of HRD‐LOH score, proportion of COSMIC signature #3, and number of indels in The Cancer Genome Atlas (TCGA) cohort. In (A) and (B), these variables were selected in the final BRCAness model by Lasso logistic regression. C, Violin plot diagram of BRCAness probabilities in the JFCR cohort. D, Violin plot diagram of BRCAness probabilities in TCGA cohort. The threshold of BRCAness probability was 0.534 to distinguish genomic features of BRCA1/2‐biallelic loss (BAL) tumors from HR‐quiescent tumors

Clinicopathological features of BRCA1/2 BAL and HR‐quiescent tumor samples in the JFCR cohort were also highly consistent with those reported in previous studies (Figure S4A).
25

BRCA1/2 BAL tumors were characterized by juvenile onset, triple‐negative subtype, bilateral and/or multiple occurrence, solid tubular histology, and nuclear grade 3; these features were more salient in BRCA1 BAL tumors than in BRCA2 BAL tumors (Figures S2 and S4A). Conversely, older age of onset, luminal subtype, unilateral and single occurrence, papillotubular or scirrhous histology, and lower nuclear grade were more representative features of HR‐quiescent tumors (Figures S2 and S4A). Although concomitant clinical information was limited in TCGA cohort data, a similar pattern was observed in terms of age at diagnosis and molecular subtype (Figure S4B). These strong associations between tumor genotype and genomic/clinicopathological features not only validates our exome‐based genomic characterization but also provides a solid basis for the development of a statistical model reliant on tumor genotype and genomic attributes.

### Development of BRCAness model based on tumor genotype and genomic attributes

We selected three variables and developed a BRCAness model using Lasso logistic regression analysis as a method of machine learning (Figure 2A,B). The Lasso model was superior to the conventional multivariate logistic regression in terms of the mean accuracies of both the training and validation data (Table 1). We then applied the model to estimate the BRCAness status of 413 miscellaneous tumor samples having variable genomic/epigenetic changes in other HBOC‐HR genes (Figure 2C,D). Document S2 outlines the details regarding the development and application of the BRCAness model.

Statistics of Lasso and conventional multivariate logistic regression models to determine BRCAness of breast cancer tumor samples

Lasso logistic regression

Multivariate logistic regression

Abbreviations: HRD, homologous recombination deficiency; LOH, loss of heterozygosity; SNV, single nucleotide variant.

### BRCAness in tumors with miscellaneous genomic/epigenetic alterations

Figure 3 shows the relationship between Lasso BRCAness probability and the type of genomic or epigenetic change in HBOC‐HR genes, along with the associated clinicopathological parameters in 413 breast tumor samples with miscellaneous genotypes. Based on the BRCAness predictive model, samples were subdivided into BRCAnesshigh (n = 97) and BRCAnesslow (n = 316) probability tumors, with a cut‐off set to 0.534 (Figure 3A). Mann‐Whitney tests, comparing samples with and without a genomic or epigenetic change in HBOC‐HR genes, showed that probability was significantly correlated with BRCA1 or RAD51C promoter hypermethylation, germline pathogenic variant with LOH/AST, or somatic homozygous deletion (Figure 3B, Table S1), consistent with previous observations.
18
 These results validated the use of our analytical methodologies and the predicted results that these genomic or epigenetic alterations were most likely BRCAness‐causing (Table S1). Worth noting, one tumor with a germline mutation in BRCA1 c.5193 + 3insT with LOH together with a PTEN somatic homozygous deletion was categorized as BRCAnesshigh in the analysis. The BRCA1 variant was located three nucleotides away from the 3′ splice site of exon 18 in the BRCA1 gene (NM_007294) and therefore was initially classified as VUS in our informatic pipeline: a splice variant is defined only as a truncation when alterations are present within two bases of the exon‐intron boundary. Nevertheless, relying on the Myriad‐Falco database, this variant was classified as pathogenic. The prediction was also validated by the clinicopathological data supplied for the samples with miscellaneous genotypes. As seen in Figure 3C, triple‐negative subtype, nuclear grade 3, and solid tubular histology were highly enriched in BRCAnesshigh tumors, which are reminiscent of BRCA1/2 BAL tumors (Figure S4). These findings indicate tight associations between BRCAness probability and the types of genomic or epigenetic alterations that occur, along with relevant clinicopathological features.

BRCAness probability, genomic/epigenetic alterations in homologous recombination (HR) genes, and clinicopathological features of breast cancer tumors with miscellaneous genotypes, excluding BRCA1/2 biallelic loss and HR‐quiescence. A, BRCAness probability and status of genomic/epigenetic alterations in BRCA1/2 or hereditary breast and ovarian cancer‐HR repair pathway‐associated (HBOC‐HR) genes. Top panel, bar plots for BRCAness probability. Bottom panel, presence or absence of hypermethylation for BRCA1 or RAD51C promoters, followed by presence or absence of germline pathogenic variants, germline variants of uncertain significance (VUS), somatic truncations, somatic VUS with (+) or without (–) LOH or additional somatic truncation (AST) in BRCA1/2 or HBOC‐HR genes, and finally, presence or absence of homozygous deletion (copy number [CN] = 0) and CN loss (CN = 1) in BRCA1/2 or HBOC‐HR genes. The samples are sorted according to BRCAness probability and divided into BRCAnesshigh and BRCAnesslow tumors using a threshold of 0.534. B, BRCAness probabilities in tumors with BRCA1 or RAD51C promoter hypermethylation, in tumors with germline pathogenic variants with LOH or AST, and in tumors with somatic homozygous deletion in BRCA1/2 or HBOC‐HR genes. P values were computed using the Mann‐Whitney U test. C, Clinicopathological features of BRCAnesshigh and BRCAnesslow tumors. Left panel, intrinsic subtype. Middle panel, nuclear grade. Right panel, tumor histology. P values were computed using Fisher’s exact test. DCIS, ductal carcinoma in situ; H, hypermethylated; NH, not hypermethylated; MT, mutated; TN, triple‐negative

### Reclassification of VUS with BRCAness

The validity of the BRCAness prediction prompted us to assess VUS pathogenicity by examining the status of the variant and the tumor clinicopathological features (Figure 4). Here, we assumed that WT allelic inactivation by LOH or AST is essential for manifesting full‐blown BRCAness. Using a tree diagram and relying on variant status, we show that 337 germline VUS (n = 337; 269 unique variants) could be classified into variants with LOH/AST (n = 77) and without LOH/AST (n = 260); those VUS lacking LOH/AST were dismissed from further assessment. Among the 77 variants with LOH/AST, 32 and 45 VUS were classified as BRCAnesshigh and BRCAnesslow tumors. These 32 BRCAnesshigh variants could be further subdivided by the absence or presence of concomitant BRCAness‐causing genomic or epigenetic alterations (Figure 4). A VUS (with LOH/AST) classified as a BRCAnesshigh tumor but lacking a concomitant BRCAness‐causing event is considered possibly pathogenic and was designated “category 1”. Comparatively, a VUS (with LOH/AST) classified as a BRCAnesshigh tumor and accompanied by a BRCAness‐causing event—that would have been assigned as “category 2”—is considered indeterminable, as the BRCAness‐causing event masks the effect of the VUS. Furthermore, we considered a VUS with LOH/AST in a BRCAnesslow tumor (“category 3”) as possibly benign, as a VUS with a homozygous state does not cause HRD (Figure 4).

Pathogenicity of variants of uncertain significance (VUS) with BRCAness probability and clinicopathological features. Breast cancer tumors with VUS with LOH or additional somatic truncations (AST) were classified into BRCAnesshigh and BRCAnesslow tumor categories according to BRCAness probability. BRCAnesshigh tumors were subdivided into tumors with and without concomitant BRCAness‐causing events: BRCA1 or RAD51C promoter hypermethylation, germline pathogenic variants with LOH/AST, or somatic homozygous deletions in BRCA1/2 or hereditary breast and ovarian cancer‐homologous recombination repair pathway‐associated genes. Germline VUS were then assessed against clinicopathological features that are characteristic of tumors with BRCA1/2 biallelic loss. One representative variant per category is shown, with tumor characteristics as examples. PLB, potentially likely benign; PLP, potentially likely pathogenic

In assessing the clinicopathological attributes, we used three features to distinguish between BRCAnesshigh and BRCAnesslow tumors: internal subtype, histology, and nuclear grade. Whereas triple‐negative subtype, solid tubular histology, and nuclear grade 3 were characteristic of BRCAnesshigh tumors, other phenotypes were linked with BRCAnesslow tumors. When a category 1 VUS consistently had clinicopathological features characteristic of BRCAnesshigh tumors, the sample was classified as PLP. A VUS of category 3 with one or more BRCAnesslow tumor characteristics (ie, not BRCAnesshigh characteristics) was interpreted as PLB (Figure 4).

Using the tree classification, we identified 19 (18 unique), 10 (10 unique), and 43 (41 unique) category 1, category 2, and category 3 variants from among the 77 VUS with LOH/AST from the original categorization. Five (two unique) variants could not be classified according to the tree diagram because they appeared in both BRCAnesshigh and BRCAnesslow tumor classes, and, instead, were assigned as “conflicting” (Table S2). Ten cases had more than one VUS (total 21 variants) accompanied by LOH/AST: three cases had six variants, two cases had five variants, and five cases had 10 variants; these were assigned to categories 1, 2, and 3, respectively (Table S2). Of note, there are three important considerations concerning these assignments: (i) the BRCAness‐causing variant cannot be inferred when two or more category 1 variants simultaneously reside in a tumor; (ii) in contrast, regardless of the number of VUS in a case, the category 2 variant is indeterminable because of the concomitant BRCAness‐causing event; and (iii) all category 3 variants are considered benign as none of them causes BRCAness. Based on these points, when there were two or more variants with LOH/AST in a tumor, the variants in categories 1 and 2 were reclassified as VUS, whereas category 3 variants were classified as PLB (Table S2). Further assessments in terms of the clinicopathological features identified six (five unique) and 36 (35 unique) VUS as PLP and PLB, respectively (Table S2). Among them, one category 1 (BLM c.T11C p.V4A) and two category 3 (ATM c.T7912G p.W2638G and CHEK2 c.G683A p.S228N) variants were recurrently detected.

We next reassessed the identified variants using previous reports comprising functional assays or case–control studies. The pathogenicity of BRCA1 c.T4951C p.S1651P as a PLP variant was supported by a functional assessment based on saturation genome editing
29
 (Table 2). A large case–control study of Japanese patients with breast cancer (7051 cases and 11 241 controls)
30
 revealed a similar prevalence of ATM c.C2879A p.P960H, BRCA2 c.A125G p.Y42C, and PALB2 c.G2509A p.E837K among cases and controls, indicating that these variants are not pathogenic. PALB2 c.C2135T p.A712V was also found to be benign in an Australian case–control study (1996 cases and 1998 controls).
31

BRCA1 c.G4535T p.S1512I and BRCA2 c.A125G p.Y42C showed comparable HR repair
32
, 
33
 and cell viability
34
, 
35
, 
36
 to that of WT proteins in multiple functional assessments, indicative of their nonpathogenicity. From these comparisons, we can conclude that our BRCAness model was useful in identifying LP and LB variants from VUS.

Proposed pathogenicity of variants of uncertain significance for breast cancer susceptibility genes

Abbreviations: LB, likely benign; LP, likely pathogenic.

### DISCUSSION

Numerous efforts have been made to identify genomic biomarkers that can be used to estimate the extent of HRD or BRCAness in a cell. Such biomarkers include the number of copy number variants,
13
 HRD‐LOH,
14
 HRD‐NtAI,
15
 HRD‐LST,
16
 HRD score (the sum of HRD‐LOH, HRD‐NtAI, and HRD‐LST),
17
 LOH status of BRCA1 or BRCA2 gene,
19
 COSMIC signature #3,
18
 and HRDetect.
7
 These biomarkers were identified using a range of different technologies, including multiplex ligation‐dependent probe amplification, high probe density SNP arrays (SNP 6.0 arrays), WES, and WGS. Yet, compared with these other methods, the BRCAness model developed in the current study showed superior or at least comparable statistics in terms of accuracy,
13
, 
16
 sensitivity,
7
, 
16
 specificity,
16
 and AUC.
7
, 
15
, 
16
 Among the aforementioned biomarkers, only HRDetect was developed based on a multivariate model.
7
 Whereas our BRCAness model (based on WES) comprises COSMIC signature #3, HRD‐LOH, and number of indels as the variables, HRDetect is derived from WGS and consists of COSMIC signatures #3 and #8, HRD‐LOH, the proportion of deletions with microhomology, and rearrangement signatures #3 and #8,
7
 with partially overlapping variables (COSMIC signatures #3 and HRD‐LOH). It is worth noting that rearrangement signatures are only available from WGS data and not from WES.
8
 It is also important to note that our BRCAness model reproducibly detected BRCA1/2 BAL samples in TCGA breast cancer data with very high accuracy (86.7%), and this independent validation indicates the robustness of the model.

Assessing the pathogenicity of HBOC‐causing genes is critically important when undertaking prophylactic surgery, drug treatment, surveillance, and risk assessment (eg, for at‐risk relatives). However, VUS with unknown relevance
24
 pose a challenge in clinical practice.
37
 Therefore, multiple studies, including case–control studies
38
 and functional assays,
39
 have sought to overcome the challenges associated with VUS by determining their pathogenicity. The VUS reclassification based on the hallmarks of BRCAness
26
 relies on the knowledge that many HBOC‐causing genes are biologically involved in the HR repair pathway and that their functional deficit leads to tumorigenesis. For example, of the HRD biomarkers, HRDetect and signature #3 implicate the pathogenicity of BRCA1 p.L1780P and BRCA1 p.C61G. Through stringent assessment in the current study, we propose BRCA1 p.S1651P as an LP variant, and ATM p.P960H, BRCA2 p.Y42C, PALB2 p.E837K, PALB2 p.A712V, and BRCA1 p.S1512I as LB variants. We also suggest that additional tumor sequencing data, particularly from patients with HBOC syndrome who have undergone genetic tests, will efficiently contribute to the pathogenicity classification of VUS of HBOC‐HR genes, as will more conventional case–control studies and in vitro functional assays.

There are several limitations to the current research. First, “splice site mutation” was conservatively defined and, therefore, we could have missed other critical protein truncating mutations, like BRCA1 c.5193 + 3insT, as indicated. Second, we lacked data regarding family history and the status of any additional cancers, both of which are relevant for assessing genetic risk. Third, we developed the BRCAness model using genomic features of breast cancer and, thus, it remains unclear whether the model can be applied for other HRD‐related cancers, such as ovarian, pancreatic, or prostate cancer. Finally, because we relied on VUS reclassification for BRCAness per tumor sample, the BRCAness‐causing variant cannot be inferred when two or more category 1 VUS simultaneously reside in a tumor.

In the current study, we undertook VUS reclassification for HBOC‐HR genes based on the assumption that all HBOC‐HR gene BAL tumors show the BRCAness phenotype. Whereas tumors with germline pathogenic variants with LOH or AST were enriched among BRCAnesshigh tumors, ATM BAL tumors showed low BRCAness probability and were not identified as the triple‐negative subtype (Table S1B). This finding is indeed consistent with previous observations, in which ATM‐mutated breast cancers showed low COSMIC signature #3,
18
 low HRD‐LST,
16
 and a luminal subtype.
28
, 
40
 Although ATM is a component of the HR repair pathway, ATM loss of function could exert a distinct tumorigenic program that differs from that of the other HBOC‐HR genes. These findings could imply the need for a different algorithm for the phenotypic evaluation of ATM and perhaps some of the other HBOC‐HR genes. However, given that the final classification results were supported by many previous findings derived from case–control studies and/or in vitro functional assays, we consider that the BRCAness model using tumor samples provides a reasonable basis for VUS reclassification. Further combinatorial assessment using the BRCAness model along with a case–control study and/or experimental functional analyses could provide a more detailed pathogenicity determination for VUS in breast cancer, and is warranted in future studies.

### CONFLICT OF INTEREST

Takayuki Ueno received lecture fees from Eisai, Chugai Pharmaceutical Co., and Novartis Pharma. Shinji Ohno received lecture fees from AstraZeneca, Chugai Pharmaceutical, Pfizer, Eli Lilly, Eisai, and Taiho, and research funds from Taiho. Yoshio Miki received lecture fees from AstraZeneca. Shunji Takahashi received lecture fees from MSD and Chugai Pharmaceutical and research funds from Taiho and Daiichi‐Sankyo.

### ETHICAL CONSIDERATIONS

All 171 patients received surgical treatment at the Cancer Institute Hospital of the Japanese Foundation for Cancer Research (JFCR) and provided written informed consent. Ethical approval was obtained from the internal review board of JFCR.

### Supporting information

Fig. S1‐S4

Click here for additional data file.

Table S1

Click here for additional data file.

Table S2

Click here for additional data file.

Doc. S1

Click here for additional data file.

Doc. S2

Click here for additional data file.



# SUPPLEMENTAL FILE 1: CAS-112-1310.pdf

# Preparing to download ...

[HHS Vulnerability Disclosure](https://www.hhs.gov/vulnerability-disclosure-policy/index.html)