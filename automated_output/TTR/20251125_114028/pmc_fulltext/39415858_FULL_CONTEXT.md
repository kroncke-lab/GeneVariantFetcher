# MAIN TEXT

## Bioinformatics Analysis of Differentially Expressed Genes in Carpal Tunnel Syndrome Using RNA Sequencing

### Abstract

Background:Carpal tunnel syndrome (CTS) is a common disease resulting from the median nerve entrapment at the wrist. Although CTS (prevalence=5%–10% in the general population) is the most common neuropathy, its molecular mechanisms need elucidation. We used bioinformatics to detect genes with differential expressions in CTS and introduce the molecular regulatory noncoding RNAs and signaling pathways involved.Methods:The raw files of the RNA sequencing of CTS patients and controls were obtained from GEO (accession: GSE108023), and the samples were analyzed. Differentially expressed genes were isolated using DESeq2 R. Functional analyses were conducted on the signaling pathways, biological processes, molecular functions, and cellular components of the differentially expressed genes. Additionally, interactions between the most differentially expressed genes and miRNAs and lncRNAs were investigated bioinformatically.Results:Upregulation and downregulation were observed in 790 and 922 genes, respectively. The signaling pathway analysis identified the metabolism pathways of arachidonic acid, linoleic acid, and tyrosine as the most significant pathways in CTS. Moreover, PLA2G2D and PLA2G2A with upregulated expressions and PLA2G2F, PLA2G4F, PLA2G4D, PLA2G3, and PLA2G4E with downregulated expressions were genes from the phospholipase family playing significant roles in the pathways. Further analyses demonstrated that hsa-miR-3150b-3p targeted PLA2G2A and PLA2G4F, and RP11-573D15.8-018 lncRNA had regulatory interactions with the aforementioned genes.Conclusion:Molecular studies on CTS will clarify the involved signaling pathways and provide critical data for biomedical research, drug development, and clinical applications.

### Introduction

Carpal tunnel syndrome (CTS) is a medical disorder caused by the density of the median nerve, which moves through the wrist in the carpal tunnel. The carpal tunnel is an anatomical chamber located at the base of the palm. Nine flexor tendons and median nerves pass through the carpal tunnel, surrounded on 3 sides by carpal bones. The median nerve provides sensation to the thumb, the index finger, the long finger, and half of the ring finger(1, 2).

The main symptoms of CTS are pain, numbness, and tingling in the thumb, the index finger, and the middle finger. Grip strength may be impaired, and the thenar musculature might become weak eventually. CTS affects approximately 1 in 10 persons during their lifetime and constitutes the most common stress syndrome (approximately 90% of all stress syndromes) (3–5).

About 5% of the population of the United States suffers from CTS. Caucasians have the highest risk of CTS development of all races. Women are more susceptible to CTS than men between ages 45 and 60 in a 3: 1 ratio. With only 10% of the reported CTS cases below age 30, aging is deemed a risk factor (6–8). In Iran, Moosazadeh et al (9) assessed 14 525 subjects in a meta-analysis and reported an estimated overall prevalence rate for CTS of 17.53%. Other risk factors for CTS include obesity, frequent wrist work, pregnancy, genetics, and rheumatoid arthritis. Further, experimental evidence indicates that hypothyroidism can also increase the risk of the disease (10, 11).

Genetic factors are the most significant determinants of CTS. A genome-wide association study on CTS identified 16 genomic loci significantly associated with this disease, including several loci previously associated with human height (2, 10–12). Biological factors such as genetic predisposition and anthropometric traits have significantly stronger causal relationships with CTS than occupational and environmental factors such as repetitive stress injuries, suggesting that CTS may not be simply avoided by eschewing certain activities or tasks (10).

In the present study, we drew upon bioinformatics to address the current paucity of research on the molecular mechanisms, signaling pathways, and genes involved in CTS. We also examined the potential noncoding RNAs (microRNAs, lncRNAs and circular RNAs) with regulatory role on involved genes. Further, we review literature about the expression analysis in CTS.

### Methods

The present study was conducted jointly at Bam University of Medical Sciences and Rajaie Cardiovascular, Medical and Research Center, Tehran, Iran from May 2022 to May 2023 (IR.MUBAM.REC.1401.029).

The data set was downloaded from the NCBI Gene Expression Omnibus (GEO) database. RNA-sequencing samples were obtained from 41 patients with CTS and 6 individuals without CTS from the database (GEO accession: GSE108023) (12). This research was approved by the review boards of the authors’ affiliated institutions.

Raw data files were first examined using FastQC to check the quality of the reads. Adapter sequences were removed using Trimmomatic, version 0.36, and the sequence reads were aligned to the human genome (hg38) using HISAT2, version 2.1.0. Then, the University of California Santa Cruz (UCSC) annotation file (hg38) was applied to analyze the expression patterns of the genes. The read count of each transcript was determined using HTSeq, version 0.9.1. The DESeq2 R package helped identify genes with differential expression patterns between the patient and control samples by considering 2 parameters: a log2FoldChange ≠5 and an adjusted P value <0.05.

Gene Ontology analyses, consisting of biological processes, molecular functions, and cellular components, were conducted using Enrichr (http://amp.pharm.mssm.edu/Enrichr/). Additionally, important pathways of the genes with differential expression patterns were examined using the Kyoto Encyclopedia of Genes and Genomes (KEGG). The genes exhibiting the highest upregulation and downregulation were selected for further analysis.

Following the selection of the most upregulated and downregulated genes, miRNAs with target sites on these genes were examined using TargetScan (http://www.targetscan.org/vert_80/). Next, miRNAs that were common among the highly differentially expressed genes were isolated utilizing Venny 2.1 to scrutinize their possible effects on CTS. Subsequently, interactions between lncRNAs and the aforementioned genes were explored using LncRRIsearch (http://http://rtools.cbrc.jp/LncRRIsearch), and 10 lncRNAs with the highest interactions with the genes were selected.

Google scholar (https://scholar.google.com/) and PubMed Central (https://www.ncbi.nlm.nih.gov/pmc/) cites were investigated to find all the studies on RNAs expression analysis in Carpal Tunnel Syndrome.

### Obtaining Data Sets and Samples

The data set was downloaded from the NCBI Gene Expression Omnibus (GEO) database. RNA-sequencing samples were obtained from 41 patients with CTS and 6 individuals without CTS from the database (GEO accession: GSE108023) (12). This research was approved by the review boards of the authors’ affiliated institutions.

### Bioinformatics Analysis of Raw Sample Data

Raw data files were first examined using FastQC to check the quality of the reads. Adapter sequences were removed using Trimmomatic, version 0.36, and the sequence reads were aligned to the human genome (hg38) using HISAT2, version 2.1.0. Then, the University of California Santa Cruz (UCSC) annotation file (hg38) was applied to analyze the expression patterns of the genes. The read count of each transcript was determined using HTSeq, version 0.9.1. The DESeq2 R package helped identify genes with differential expression patterns between the patient and control samples by considering 2 parameters: a log2FoldChange ≠5 and an adjusted P value <0.05.

### In Silico Functional Analysis of the Differentially Expressed Genes

Gene Ontology analyses, consisting of biological processes, molecular functions, and cellular components, were conducted using Enrichr (http://amp.pharm.mssm.edu/Enrichr/). Additionally, important pathways of the genes with differential expression patterns were examined using the Kyoto Encyclopedia of Genes and Genomes (KEGG). The genes exhibiting the highest upregulation and downregulation were selected for further analysis.

### Target Prediction of the Respective Interactions Between the Highly Differentially Expressed Genes and MicroRNAs (miRNAs) and Long Noncoding RNAs (lncRNAs)

Following the selection of the most upregulated and downregulated genes, miRNAs with target sites on these genes were examined using TargetScan (http://www.targetscan.org/vert_80/). Next, miRNAs that were common among the highly differentially expressed genes were isolated utilizing Venny 2.1 to scrutinize their possible effects on CTS. Subsequently, interactions between lncRNAs and the aforementioned genes were explored using LncRRIsearch (http://http://rtools.cbrc.jp/LncRRIsearch), and 10 lncRNAs with the highest interactions with the genes were selected.

### RNA Expression Analysis in Literature

Google scholar (https://scholar.google.com/) and PubMed Central (https://www.ncbi.nlm.nih.gov/pmc/) cites were investigated to find all the studies on RNAs expression analysis in Carpal Tunnel Syndrome.

### Results

Significant differential expression patterns were detected between the CTS group and the control group in 1712 genes, of which 790 exhibited increased expression patterns and 922 showed decreased expression patterns. The volcano plot below shows the genes whose expression patterns differed between the 2 groups (Fig. 1), and Supplementary 1 lists the names of the genes.

The volcanic plot of the genes. The red dots indicate genes with decreased expression patterns, the blue dots indicate genes with increased expression patterns, and the green dots indicate genes without significant differential expression patterns in terms of the log2FoldChange

Gene Ontology analyses, composed of biological processes, molecular functions, and cellular components, were performed on the genes whose expression patterns differed between the 2 groups. The biological process analysis showed that the genes with differential expression patterns were involved in biological processes such as inflammatory and immune responses. The molecular function analysis revealed the involvement of the differentially expressed genes in such molecular functions as protein/threonine kinase activity, magnesium ion binding, and metallopeptidase processes. The cellular component analysis demonstrated that the differentially expressed genes were present in the cytoskeleton of micro-tubules, organ membranes, peroxisomes, and nuclei (Fig. 2A–C).

A) biological processes (BP), B) molecular functions (MF), and C) cellular components (CC) of the genes with differential expression patterns

The results of the signaling pathway analysis with adjusted P values <0.05 showed that 21 genes whose expression patterns differed between the CTS and control groups were involved in the metabolism pathways of arachidonic acid, linoleic acid, and tyrosine (Table 1). As is demonstrated in Table 1, the following genes were common in the metabolism pathways of arachidonic acid, linoleic acid, and tyrosine: PLA2G2D, PLA2G2A, PLA2G2F, PLA2G4F, PLA2G4D, and PLA2G3. Moreover, PLA2G2D and PLA2G2A with upregulated expressions and PLA2G2F, PLA2G4F, PLA2G4D, PLA2G3, and PLA2G4E with downregulated expressions were genes from the phospholipase family, playing significant roles in the signaling pathways of CTS. Further, PLA2G2A with log2FoldChange =6.950397 and PLA2G4F with log2FoldChange = −15.08284 had the highest upregulation and downregulation, respectively.

Signaling pathways involved in genes with differential expression patterns between the patient and control groups

Our analysis revealed that 9 miRNAs and 53 miRNAs probably targeted the 3′-UTR of PLA2G2A and PLA2G4F, respectively (Supplementary Tables 2 & 3). Further, of all 62 predicted miRNAs, hsa-miR-3150b-3p was the only common miRNA targeting both genes (Fig. 3) with a probable regulatory role in CTS.

The common microRNAs (miRNAs) that targeted both PLA2G4F and PLA2G2A corresponding to the TargetScan database

The analysis of the interactions between lncRNAs and the PLA2G2A and PLA2G4F genes yielded 100 lncRNAs with probable interactions with each gene, respectively. LncRRIsearch predicts the local base-pairing interactions based on interaction energy. Accordingly, the threshold interaction energy was set to −16 kcal/mol, with lower scores predicting the most probable interaction. Among all lncRNAs with probable interactions with PLA2G2A and PLA2G4F, 10 lncRNAs had the highest energy score of interaction (Table 2). Interestingly, RP11-573D15.8-018 lncRNA, with probable interactions with both PLA2G2A and PLA2G4F, had the highest score.

Top 10 long noncoding RNAs (lncRNAs) that could interact with the PLA2G4F and PLA2G2A genes

Totally 10 studies were found in both investigated data bases (Google Scholar and PubMed). Five studies were about association analysis of SNPs and mutations of genes in CTS (10–14). Other five studies indicated that the expression of ITGAL, ITGAM, PECAM1, TGFBR2, PRG5, CASP8, IGFBP5, COL5A1 and WINT9a might be involved in the CTS (15–19).

No studies were found about ncRNAs expression and their regulatory roles in CTS and this study is the first study to suggest the potential regulatory roles of miRNAs and lncRNAs in CTS.

### Genes with Differential Expression Patterns in the Patient and Control Samples

Significant differential expression patterns were detected between the CTS group and the control group in 1712 genes, of which 790 exhibited increased expression patterns and 922 showed decreased expression patterns. The volcano plot below shows the genes whose expression patterns differed between the 2 groups (Fig. 1), and Supplementary 1 lists the names of the genes.

The volcanic plot of the genes. The red dots indicate genes with decreased expression patterns, the blue dots indicate genes with increased expression patterns, and the green dots indicate genes without significant differential expression patterns in terms of the log2FoldChange

### Significance of Metabolism and Immune Response Pathways in CTS

Gene Ontology analyses, composed of biological processes, molecular functions, and cellular components, were performed on the genes whose expression patterns differed between the 2 groups. The biological process analysis showed that the genes with differential expression patterns were involved in biological processes such as inflammatory and immune responses. The molecular function analysis revealed the involvement of the differentially expressed genes in such molecular functions as protein/threonine kinase activity, magnesium ion binding, and metallopeptidase processes. The cellular component analysis demonstrated that the differentially expressed genes were present in the cytoskeleton of micro-tubules, organ membranes, peroxisomes, and nuclei (Fig. 2A–C).

A) biological processes (BP), B) molecular functions (MF), and C) cellular components (CC) of the genes with differential expression patterns

The results of the signaling pathway analysis with adjusted P values <0.05 showed that 21 genes whose expression patterns differed between the CTS and control groups were involved in the metabolism pathways of arachidonic acid, linoleic acid, and tyrosine (Table 1). As is demonstrated in Table 1, the following genes were common in the metabolism pathways of arachidonic acid, linoleic acid, and tyrosine: PLA2G2D, PLA2G2A, PLA2G2F, PLA2G4F, PLA2G4D, and PLA2G3. Moreover, PLA2G2D and PLA2G2A with upregulated expressions and PLA2G2F, PLA2G4F, PLA2G4D, PLA2G3, and PLA2G4E with downregulated expressions were genes from the phospholipase family, playing significant roles in the signaling pathways of CTS. Further, PLA2G2A with log2FoldChange =6.950397 and PLA2G4F with log2FoldChange = −15.08284 had the highest upregulation and downregulation, respectively.

Signaling pathways involved in genes with differential expression patterns between the patient and control groups

### Probable Regulatory Effects of hsa-miR-3150b-3p and RP11-573D15.8-018 lncRNA on PLA2G2A and PLA2G4F

Our analysis revealed that 9 miRNAs and 53 miRNAs probably targeted the 3′-UTR of PLA2G2A and PLA2G4F, respectively (Supplementary Tables 2 & 3). Further, of all 62 predicted miRNAs, hsa-miR-3150b-3p was the only common miRNA targeting both genes (Fig. 3) with a probable regulatory role in CTS.

The common microRNAs (miRNAs) that targeted both PLA2G4F and PLA2G2A corresponding to the TargetScan database

The analysis of the interactions between lncRNAs and the PLA2G2A and PLA2G4F genes yielded 100 lncRNAs with probable interactions with each gene, respectively. LncRRIsearch predicts the local base-pairing interactions based on interaction energy. Accordingly, the threshold interaction energy was set to −16 kcal/mol, with lower scores predicting the most probable interaction. Among all lncRNAs with probable interactions with PLA2G2A and PLA2G4F, 10 lncRNAs had the highest energy score of interaction (Table 2). Interestingly, RP11-573D15.8-018 lncRNA, with probable interactions with both PLA2G2A and PLA2G4F, had the highest score.

Top 10 long noncoding RNAs (lncRNAs) that could interact with the PLA2G4F and PLA2G2A genes

### RNAs expression analysis and CTS

Totally 10 studies were found in both investigated data bases (Google Scholar and PubMed). Five studies were about association analysis of SNPs and mutations of genes in CTS (10–14). Other five studies indicated that the expression of ITGAL, ITGAM, PECAM1, TGFBR2, PRG5, CASP8, IGFBP5, COL5A1 and WINT9a might be involved in the CTS (15–19).

No studies were found about ncRNAs expression and their regulatory roles in CTS and this study is the first study to suggest the potential regulatory roles of miRNAs and lncRNAs in CTS.

### Discussion

CTS, albeit prevalent, has poorly-understood pathophysiology, particularly in terms of its genetic predisposers and molecular mechanisms (13). In the present bioinformatics study on patients with CTS, we found that the metabolism pathways of arachidonic acid, linoleic acid, α-linolenic acid, and tyrosine were influential in the development of this disease.

Arachidonic acid and its metabolites are essential for modulating the function and survival of nerve cells. Arachidonic acid derivatives are involved in the production of cyclooxygenases and lipoxygenases to generate such mediators as prostaglandins, leukotrienes, and eicosanoids, involved in various biological events, including inflammation. Additionally, arachidonic acid contributes to neuron growth via a calcium-dependent mechanism, and its metabolites are involved in the physiopathology of many central nervous system disorders and acute inflammatory disorders (14–17). A study showed that the involvement of lipoxygenases in arachidonic acid metabolism could significantly affect the incidence of rheumatoid arthritis, a systemic connective tissue disease (18).

Linoleic acid is the main essential fatty acid in the diet and must be metabolized to other substances. Linoleic acid concentrations are almost normal or slightly above normal in patients with diabetes (19, 20). According to a previous study, the concentration of the major metabolites of insulin increased significantly in the treatment of patients with diabetes, and it was concluded that the acupuncture treatment of patients with diabetes affected unique metabolic pathways, including those of α-linolenic acid, glutamine, alanine, and activated vitamin B6 (21).

Williams et al (16) examined the effects of the precursor doses of catecholamine L-tyrosine and reported that abnormal amine metabolism was often the etiology of schizophrenia. Magy et al (22) found that TTR (Tyr78Phe) mutation was associated with peripheral neuropathy, CTS, and cutaneous amyloidosis.

According to our results vis-à-vis differentially expressed genes, inflammatory and immune responses are among the major biological processes involved in CTS.

The immune system has been implicated in several neurological disorders such as peripheral neuropathies and neuropathic pains (23, 24). Information is, however, scarce regarding the immune response in patients with CTS, although neuroinflammation has been suggested as a potential contributor to its pathogenesis (25, 26). Moalem-Taylor et al (25) explored the immune profile in patients with CTS and posited an association between the syndrome and alterations in the homeostasis of memory T cells and inflammatory modulating cytokines/chemokines. Arachidonic acid and its metabolites have attracted much attention concerning inflammatory processes and diseases (27, 28). Cheng et al(29) reported that as a result of acid-induced injury (elevated H2O2 levels) in the lower esophageal sphincter, inflammation could lead to changes in arachidonic acid metabolism. Calabrese et al (30) showed that increased arachidonic acid levels in the glycerolipid pools of inflammatory cells could produce excessive cysteinyl leukotrienes in patients with asthma. Studies on animal models indicate that conjugated linoleic acid, a class of linoleic acid isomers, enhances immune functions (31–33) while ameliorating immune-mediated catabolism (34, 35) and stimulating the production of anti-inflammatory cytokines (36). Erdinest et al (37) concluded that α-linolenic acid might serve as an anti-inflammatory agent and was mediated through the signal transduction of nuclear factor κ-light-chain-enhancer of activated B cells (NF-κB). Likewise, in the present study, PLA2G2D, PLA2G2F, PLA2G4F, PLA2G4D, PLA2G3, PLA2G2A, and PLA2G4E showed significant expression changes in the mentioned pathways, resulting in CTS.

Phospholipase A2 (PLA2) catalyzes the hydrolysis of fatty acyl ester bonds at the sn-2 position, catalyzes phosphoglycerides, and releases fatty acids and lysophospholipids. One of the fatty acids that can be released from the membrane by PLA2 activity is arachidonic acid, which is a vital precursor to the biosynthesis of various eicosanoids, including prostaglandins, thromboxane, and leukotrienes (38). Moreover, the PLA2G2A gene (PLA2 IIA) is a secretory PLA2 expressed in lacrimal glands, cartilage cells, and amniotic epithelial cells (38, 39) and appears to play a variety of roles in human diseases such as colon cancer, coronary artery disease, and inflammation (40, 41). The induction of PLA2G2A is a feature of inflammatory responses, and increased PLA2G2A expression has been reported in several types of malignancies such as pancreatic and prostate cancer (39–42). The PLA2G2D gene (sPLA2-IID) is abundant in the spleen and lymph nodes, attenuates the immune responses of Th1 and Th17, and is associated with the anti-inflammatory metabolites u3 polyunsaturated fatty acids (PUFAs) (43). The PLA2G3 gene is associated with the risk of colorectal cancer and ovarian cancer metastasis (44, 45). Three members of the PLA2 family, both the soluble (PLA2G2F) and cytoplasmic (PLA2G4D and PLA2G4E) forms, comprise a vital pathogenic pathway between pityriasis rubra pilaris and psoriasis (46).

A large family of RNAs without a coding function, noncoding RNAs (ncRNAs) are composed of miRNAs, circular RNAs (circRNAs), and lncRNAs and play significant roles in the regulatory process(47–49). To our knowledge, the existing literature features no study on the role of miRNAs in CTS. Our analysis revealed that hsamiR-3150b-3p could target the PLA2G2A and PLA2G4F genes, and RP11-573D15.8-018 lncRNA could interact with both genes. Thus, these 2 miRNA and lncRNA might have regulatory roles in CTS by interfering with PLA2G2A and PLA2G4F and the relevant pathways. Nevertheless, there is no study on the role of hsamiR-3150b-3p and RP11-573D15.8-018 lncRNA in the regulation of PLA2G2A and PLA2G4F and the metabolism pathways of arachidonic acid, linoleic acid, and tyrosine in CTS, requiring further experimental analyses. Based on the findings of the present study, we posit that the down-regulation of PLA2G4F and the upregulation of PLA2G2A could be the result of miR-3150b-3p expression and could affect RP11-573D15.8-018 expression. Hence, the arachidonic acid pathway, a significant pathway in the inflammation process, could be regulated by these molecular markers and contribute to CTS progression (Fig. 4). Further experimental validation to confirm the expression patterns of the genes, microRNA and lncRNA presented in this study are proposed.

A schematic view of the suggested pathway involved in carpal tunnel syndrome (CTS) according to our findings. The downregulation of PLA2G4F and the upregulation of PLA2G2A could be the result of the altered expressions of miR-3150b-3p and RP11-573D15.8-018

The findings of the current investigation should be interpreted in light of the following limitations. Our results need confirmation via experimental analyses. Secondly, studies and bioinformatics data sets on the molecular mechanisms of CTS pathogenesis are limited. The only existed data set of high throughput sequencing for CTS was GSE108023 which we employed for the present study.

### Conclusion

The existing literature contains a paucity of information regarding the molecular mechanisms underlying CTS. In addition, for all the aforementioned investigations, the current study is the first to utilize gene expression profiling to identify the molecular pathways involved in CTS. Indeed, we succeeded in detecting genes whose expression patterns differed between patients with CTS and subjects without CTS. Nonetheless, to boost the diagnostic and therapeutic processes concerning CTS, we need more detailed studies aimed at identifying the exact molecular pathways and major regulators of this syndrome. However, our results demonstrated that hsa-miR-3150b-3p and RP11-573D15.8-018 might have regulatory role in CTS through possible interaction with PLA2G2A and PLA2G4F.

Additionally, we suggest that the expression patterns of the genes, microRNA and lncRNA introduced herein and the pathways involved be confirmed by more detailed laboratory studies to help clarify the exact mechanisms of the disease.

### Journalism Ethics considerations

Ethical issues (Including plagiarism, informed consent, misconduct, data fabrication and/or falsification, double publication and/or submission, redundancy, etc.) have been completely observed by the authors.



# SUPPLEMENTAL FILE 1: IJPH-53-1871.pdf

# Preparing to download ...

[HHS Vulnerability Disclosure](https://www.hhs.gov/vulnerability-disclosure-policy/index.html)