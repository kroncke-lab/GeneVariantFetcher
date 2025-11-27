# MAIN TEXT

## Germline variation in BRCA1/2 is highly ethnic‐specific: Evidence from over 30,000 Chinese hereditary breast and ovarian cancer patients

### Abstract

BRCA1 and BRCA2 play essential roles in maintaining the genome stability. Pathogenic germline mutations in these two genes disrupt their function, lead to genome instability and increase the risk of developing breast and ovarian cancers. BRCA mutations have been extensively screened in Caucasian populations, and the resulting information are used globally as the standard reference in clinical diagnosis, treatment and prevention of BRCA‐related cancers. Recent studies suggest that BRCA mutations can be ethnic‐specific, raising the question whether a Caucasian‐based BRCA mutation information can be used as a universal standard worldwide, or whether an ethnicity‐based BRCA mutation information system need to be developed for the corresponding ethnic populations. In this study, we used Chinese population as a model to test ethnicity‐specific BRCA mutations considering that China has one of the latest numbers of breast cancer patients therefore BRCA mutation carriers. Through comprehensive data mining, standardization and annotation, we collected 1,088 distinct BRCA variants derived from over 30,000 Chinese individuals, one of the largest BRCA data set from a non‐Caucasian population covering nearly all known BRCA variants in the Chinese population (https://dbBRCA-Chinese.fhs.umac.mo). Using this data, we performed multi‐layered analyses to determine the similarities and differences of BRCA variation between Chinese and non‐Chinese ethnic populations. The results show the substantial differences of BRCA data between Chinese and non‐Chinese ethnicities. Our study indicates that the current Caucasian population‐based BRCA data is not adequate to represent the BRCA status in non‐Caucasian populations. Therefore, ethnic‐based BRCA standards need to be established to serve for the non‐Caucasian populations.

### Introduction

Approximately 10–15% of breast cancer cases are caused by hereditary genetic mutations.1 The most penetrating mutations are those in the BRCA1 and BRCA2 (BRCA) genes,2, 3 which are essential for maintaining the genome stability. Women carrying pathogenic mutations in BRCA1 have a 72% lifetime risk of developing breast cancer, while those with BRCA2 mutations have a 69% risk.4 Mutations in BRCA also increase the risk for ovarian cancer, prostate cancer, melanoma and pancreatic cancer. Identification of BRCA mutation carriers before the development of cancer is crucial in order to protect them from developing cancer by taking preventive measures of early cancer surveillance, chemoprevention and preventive surgery.5, 6, 7, 8, 9 Extensive efforts have identified a large number of BRCA mutations, mainly in the Caucasian populations of Europe and North America. BRCA databases with well‐documented, annotated and freely accessible BRCA information have been developed and used worldwide as the references for diagnosis, treatment and prevention of BRCA associated cancers.

Studies have revealed that human BRCA1 and BRCA2 are rapidly evolving under positive selection to effectively protect genome stability.10, 11 Recent studies have further suggested that the variation in human BRCA could be ethnic‐specific in different ethnic populations. For example, BRCA variants within Latin American populations are highly heterogeneous,12 and BRCA variants in Asian populations differ substantially from those in other populations.13 Understanding the ethnic‐specificity of BRCA variation is important, as it can provide a precise genetic basis to study the relationship between human evolution and diseases. Furthermore, it will determine whether the Caucasian‐based BRCA data is adequate to serve as a universal reference to determine BRCA status in non‐Caucasian populations around the world, or whether the ethnicity‐based BRCA mutation data should be developed instead. However, the answer remains elusive owing to the lack of BRCA data from most of the non‐Caucasian populations.12, 13, 14, 15 For example, in the recently completed CIMBA study that collected BRCA data from 49 countries across six continents, the data available from any single, non‐Caucasian ethnic populations remain very limited.15 To fully prove the existence of ethnic‐specific BRCA mutations, a comprehensive data from particular ethnic populations will be required.

The Chinese population is the largest one in the world.16 Breast cancer is the most common cancer among Chinese women, with 260,000 new breast cancer cases diagnosed and 70,000 mortalities annually.17 With nearly two decades of BRCA studies in China, in particular in the recent years owing largely to the adoption of next‐generation sequencing technologies, BRCA data exclusively derived from the Chinese population are increasingly reported.17, 18, 19 Thus, the Chinese population can serve as an ideal model to test the presence of ethnic‐specific BRCA mutations. Through a comprehensive data mining, standardization and annotation, we collected nearly all BRCA variant data currently available for the Chinese population and developed the data into a public database dbBRCA‐Chinese (https://dbBRCA-Chinese.fhs.umac.mo). Using this rich data set, we studied the similarities and differences in BRCA variation between Chinese and non‐Chinese populations. The data from our study provides a convincing evidence for the existence of ethnic‐specific BRCA mutation. Here, we report detailed information from the study.

### Materials and Methods

We searched for resources reporting BRCA variation data from individuals of Chinese ethnicity, including publications in PubMed, the China National Knowledge Infrastructure (CNKI) database (http://oversea.cnki.net/kns55/default.aspx) and WanFang (http://www.wanfangdata.com/COJ/intr.asp#China Online), as well as Chinese‐derived BRCA variants in existing BRCA databases.19 For the collected BRCA variants, we performed extensive standardization and reannotation, following HGVS20 and ACMG guidelines.21 The reference sequences used for BRCA1 analysis were: cDNA NM_007294.3, protein NP_009225.1, genome hg19, BIC cDNA: U14680.1 and BIC protein: AAA73985.1; those used for BRCA2 were: cDNA NM_000059.3, protein NP_000050.2, genome hg19, BIC cDNA: U43746.1 and BIC protein: AAB07223.1, respectively22 (https://doi.org/10.1093/nar/gkw1070). The transcript‐oriented position of each variant was converted to its respective genome position in hg19 using the Position Converter tool in Mutalyzer,23 and the consistency with HGVS nomenclature was confirmed using the Name Checker DNA tool.24 The variants were annotated using ANNOVAR with eight reference databases namely: RefGene, dbSNP (version 150), 1000genome, ESP6500, ExAC, ClinVar, InterVar and DBNSFP.24 Following BRCA databases were used for the comparative analysis: BIC 25 (https://research.nhgri.nih.gov/bic/, accessed February 20, 2018), ClinVar26(http://www.ncbi.nlm.nih.gov/clinvar/, accessed February 20, 2018), BRCA Exchange (http://brcaexchange.org, accessed February 20, 2018), ENIGMA27 (downloaded from BED database, http://brcaexchange.org, accessed February 20, 2018), BMD (http://www.arup.utah.edu/database/BRCA/Home/BRCA1_landing.php,http://www.arup.utah.edu/database/BRCA/Home/BRCA2_landing.php)(https://www.aruplab.com/topics/breast-cancer/brcadatabase, accessed February 20, 2018), LOVD28 (http://www.lovd.nl/3.0/home, accessed February 20, 2018) and CIMBA.15 The Chinese variants present in these databases were classified as known variants; those absent were classified as novel variants and deposited in ClinVar database29 (accession number nstd165). Five categories based on ACMG guidelines were used: pathogenic, likely pathogenic, uncertain significance, likely benign, and benign for the variant classification of known variants using ClinVar database.30, 36 For those variants not matched with the variant list using the above‐mentioned resources, their classifications were predicted from InterVar database using ANNOVAR annotation tool.30 For these variants not being able to classified, they were included as “Unclassified” group. BRCA variants in Latin population were extracted from Villarreal‐Garza et al.,12
BRCA variants in Asian populations were extracted from,13 and BRCA variants in the Indian population were extracted from,13, 31 respectively.

Fisher's exact test was used to analyze the differences between the Chinese and non‐Chinese variant data. p < 0.05 was considered as significant difference.

### Variant data collection and analysis

We searched for resources reporting BRCA variation data from individuals of Chinese ethnicity, including publications in PubMed, the China National Knowledge Infrastructure (CNKI) database (http://oversea.cnki.net/kns55/default.aspx) and WanFang (http://www.wanfangdata.com/COJ/intr.asp#China Online), as well as Chinese‐derived BRCA variants in existing BRCA databases.19 For the collected BRCA variants, we performed extensive standardization and reannotation, following HGVS20 and ACMG guidelines.21 The reference sequences used for BRCA1 analysis were: cDNA NM_007294.3, protein NP_009225.1, genome hg19, BIC cDNA: U14680.1 and BIC protein: AAA73985.1; those used for BRCA2 were: cDNA NM_000059.3, protein NP_000050.2, genome hg19, BIC cDNA: U43746.1 and BIC protein: AAB07223.1, respectively22 (https://doi.org/10.1093/nar/gkw1070). The transcript‐oriented position of each variant was converted to its respective genome position in hg19 using the Position Converter tool in Mutalyzer,23 and the consistency with HGVS nomenclature was confirmed using the Name Checker DNA tool.24 The variants were annotated using ANNOVAR with eight reference databases namely: RefGene, dbSNP (version 150), 1000genome, ESP6500, ExAC, ClinVar, InterVar and DBNSFP.24 Following BRCA databases were used for the comparative analysis: BIC 25 (https://research.nhgri.nih.gov/bic/, accessed February 20, 2018), ClinVar26(http://www.ncbi.nlm.nih.gov/clinvar/, accessed February 20, 2018), BRCA Exchange (http://brcaexchange.org, accessed February 20, 2018), ENIGMA27 (downloaded from BED database, http://brcaexchange.org, accessed February 20, 2018), BMD (http://www.arup.utah.edu/database/BRCA/Home/BRCA1_landing.php,http://www.arup.utah.edu/database/BRCA/Home/BRCA2_landing.php)(https://www.aruplab.com/topics/breast-cancer/brcadatabase, accessed February 20, 2018), LOVD28 (http://www.lovd.nl/3.0/home, accessed February 20, 2018) and CIMBA.15 The Chinese variants present in these databases were classified as known variants; those absent were classified as novel variants and deposited in ClinVar database29 (accession number nstd165). Five categories based on ACMG guidelines were used: pathogenic, likely pathogenic, uncertain significance, likely benign, and benign for the variant classification of known variants using ClinVar database.30, 36 For those variants not matched with the variant list using the above‐mentioned resources, their classifications were predicted from InterVar database using ANNOVAR annotation tool.30 For these variants not being able to classified, they were included as “Unclassified” group. BRCA variants in Latin population were extracted from Villarreal‐Garza et al.,12
BRCA variants in Asian populations were extracted from,13 and BRCA variants in the Indian population were extracted from,13, 31 respectively.

### Statistical data analysis

Fisher's exact test was used to analyze the differences between the Chinese and non‐Chinese variant data. p < 0.05 was considered as significant difference.

### Results

The majority of Chinese BRCA variation data was gathered from mainland China (90.2%), and the remainder were from Hong Kong (6.4%), Taiwan (1.6%), Singapore (0.9) and Malaysia (0.4%) (Table S1A, S1B, Supporting Information). In total, 31,689 cases with Chinese ethnicity were tested for BRCA mutations and their results were reported between 1999 and 2017; of these, 69.3% were reported between 2016 and 2017. Nearly all BRCA variants were from breast and/or ovarian cancer patients under different clinical criteria, except for four variants that were derived from 1,043 healthy control individuals.32

We collected and summarized the clinical information reported from each of the reference studies (Fig. S1 and Table S2, Supporting Information). The data showed several unique features as follows: 1). BMI is considered as one of the risk factors for breast cancer. However, according to our data, 83.7% of Chinese patients were within the normal range of 18.5–22.9 BMI and 10.5% were seen to have even lower than 18.5 BMI, indicating that obesity had no significant role in increasing risk for breast and ovarian cancer in these Chinese patients; 2). Family history is considered as a high‐risk factor for hereditary breast cancer. However, 72.3% of patients from our study did not report any family history of cancer, suggesting that family history was not an essential factor in this disease cohort; 3). Stage II breast cancer cases accounted for 66.2% whereas stage I cases only 6.3%, indicating the lack of earlier diagnosis.

A total of 3,791 BRCA mutation carriers were identified (12%, 2,123 BRCA1 and 1,688 BRCA2), of which 1,978 carriers (52.2%, 990 BRCA1 and 988 BRCA2) were within the clinically reported mutation categories of pathogenic or likely pathogenic. By standardizing and re‐annotating all the variants following Human Genome Variation Society (HGVS)‐recommended nomenclature, we identified a total of 1,088 distinct BRCA variants (557 in BRCA1 and 531 in BRCA2) of which 519 (47.7%, 278 in BRCA1 and 241 in BRCA2) were recurrent, and 50% were either detected or validated by Sanger sequencing (Tables S4–S5, Supporting Information). Except the 26 variants specifically from Uygur ethnic population, all variants were from Han Chinese population. We developed the dbBRCA‐Chinese database as an open source to host the entire set of BRCA variants and their annotation information (https://dbBRCA-Chinese.fhs.umac.mo).

The 1,088 BRCA variants were distributed with different frequencies among the 3,791 BRCA variant carriers. The age distribution data show that 52.3% of the cancer developed at early age (<40 years) in these BRCA variant carriers (Table S1C, Supporting Information). About half of the variants were detected only once in single individuals (46.5% in BRCA1, 54.6% in BRCA2), and the rest variants were distributed between 2 to 100 individuals with increased frequencies (Table S1D, Supporting Information).

We classified the BRCA variants into five classes—pathogenic, likely pathogenic, uncertain significance, likely benign, and benign—following the American College of Medical Genetics and Genomics (ACMG) standards and guidelines. Pathogenic and likely pathogenic variants, which are clinically reportable, accounted for 46% of BRCA1 variants and 52% of BRCA2 variants. Such higher rates do not necessarily imply high BRCA variation rate in Chinese breast and ovarian cancer patient population but reflect the fact that the cancer patients included in many of the studies were selected from high‐risk patients of strong family history or early age of cancer development. Importantly, 13% and 8% of BRCA1 and BRCA2 were classified as of uncertain significance, and 30% of both BRCA1 and BRCA2 variants remained as unclassified variants (Fig. 1
a). This fact indicates that much effort needs to be made in order to determine the function of these BRCA variants existing in Chinese population.

BRCA data. (a) Clinical classification of Chinese BRCA data as pathogenic, likely pathogenic, uncertain significance, likely benign, benign and unclassified. The pathogenic and likely pathogenic variants accounted for 49.5%. (b) Relationship between population sizes and their contribution to current BRCA data. The proportions of different human ethnic populations were from the database (http://www.worldometers.info/world‐population/), the ethnic origins of BRCA data were from different BRCA databases. It shows that the current BRCA data is not proportional to the human ethnic populations. (c) Comparison of BRCA data between Chinese and non‐Chinese populations. A total of 557 BRCA1 and 531 BRCA2 Chinese variants were compared to 6,344 BRCA1 and 8,886 BRCA2 non‐Chinese variants compiled from all existing BRCA databases. The results show that 38% of Chinese BRCA variants were present only in the Chinese population.

We compiled ethnic origins of currently existing BRCA variation data to know the status of BRCA study across ethnic human populations. By combining the BRCA variants from the Breast Cancer Information Core (BIC), ClinVar, BRCA Exchange Database (BED), BRCA1 and BRCA2 Mutation Database (BMD), Leiden Open Variation Database (LOVD) and ENIGMA (Evidence‐based Network for the Interpretation of Germline Mutant Alleles), we generated a single BRCA variation data set containing 6,343 distinct BRCA1 and 8,884 distinct BRCA2 variants (S6A). Classification of the ethnic origins of these variants showed that 62% were from Caucasian populations and 15% from the Ashkenazi Jewish population. The remaining 23% were originated from non‐Chinese Asian (13%), Latino (5%), African (3%), and Chinese (2%) populations (Fig. 1
b and Tables S6B, S6C, Supporting Information). Of the 5,925 CIMBA BRCA variants data with defined ethnicities, 80.2% were from Caucasian, 2.4% from Ashkenazi Jews, 10.8% from Asian, 3.8% from Hispanic and 2.8% from African American.15 The analysis shows that the current BRCA variant data contains very limited information from non‐Caucasian populations.

We performed the multi‐layer analyses to investigate the similarities and differences in BRCA variation between Chinese and non‐Chinese populations.

In order to determine the similarities and differences of BRCA variation between Chinese and non‐Chinese populations, we made a comprehensive comparison between Chinese and any available non‐Chinese BRCA data as represented below. The rationale for comparing with each database are:GnomAD: It contains extensive normal population variation data collected from human population by the largest exome and whole‐genome sequencing projects. The comparison aimed to determine the similarities and differences of the normal BRCA variants present between Chinese and non‐Chinese populations;BIC, BED, BMD, ClinVar, ENIGMA and LOVD: These are the major BRCA databases, with the data mostly derived from Caucasians as indicated by our analysis (Fig. 1
b). Comparison with these databases aimed to determine the similarities and differences of BRCA mutation between Chinese and Caucasian (mostly) populations;CIMBA data: It contains BRCA data collected from 49 countries across six continents. The comparison aimed to determine the similarities and differences between Chinese and worldwide non‐Chinese populations including more non‐Caucasian data;Latin America and the Caribbean data: The data were from Latin American population of Argentina, Bahamas, Brazil, Chile, Colombia, Costa Rica, Cuba, Mexico, Peru, Puerto Rico, Uruguay, Venezuela and the Hispanic population in the United States.12 The comparison aimed to determine the similarities and differences between Chinese and Latin America populations;Non‐Chinese Asian populations: BRCA data are available from Bangladeshi, Filipino, Iranian, Israeli, Japanese, Korean, Lebanese, Malay, Oman, Pakistani, Sri Lankan, Thai and Turkish populations.13 These populations were genetically and geographically closer to the Chinese population than other non‐Chinese populations. The comparison aimed to determine the similarities and differences between Chinese and the non‐Chinese Asian populations;Indian population: India has the 2nd largest population in the world, with highly diversified genetic background. Several large‐scale BRCA studies were reported recently with substantial BRCA data collected from the Indian patients.13, 31 The comparison aimed to determine the similarities and differences between Chinese and Indian populations, the two largest populations in the world.

GnomAD: It contains extensive normal population variation data collected from human population by the largest exome and whole‐genome sequencing projects. The comparison aimed to determine the similarities and differences of the normal BRCA variants present between Chinese and non‐Chinese populations;

BIC, BED, BMD, ClinVar, ENIGMA and LOVD: These are the major BRCA databases, with the data mostly derived from Caucasians as indicated by our analysis (Fig. 1
b). Comparison with these databases aimed to determine the similarities and differences of BRCA mutation between Chinese and Caucasian (mostly) populations;

CIMBA data: It contains BRCA data collected from 49 countries across six continents. The comparison aimed to determine the similarities and differences between Chinese and worldwide non‐Chinese populations including more non‐Caucasian data;

Latin America and the Caribbean data: The data were from Latin American population of Argentina, Bahamas, Brazil, Chile, Colombia, Costa Rica, Cuba, Mexico, Peru, Puerto Rico, Uruguay, Venezuela and the Hispanic population in the United States.12 The comparison aimed to determine the similarities and differences between Chinese and Latin America populations;

Non‐Chinese Asian populations: BRCA data are available from Bangladeshi, Filipino, Iranian, Israeli, Japanese, Korean, Lebanese, Malay, Oman, Pakistani, Sri Lankan, Thai and Turkish populations.13 These populations were genetically and geographically closer to the Chinese population than other non‐Chinese populations. The comparison aimed to determine the similarities and differences between Chinese and the non‐Chinese Asian populations;

Indian population: India has the 2nd largest population in the world, with highly diversified genetic background. Several large‐scale BRCA studies were reported recently with substantial BRCA data collected from the Indian patients.13, 31 The comparison aimed to determine the similarities and differences between Chinese and Indian populations, the two largest populations in the world.

Matching the Chinese BRCA data with those in GnomAD from the largest exome and whole‐genome sequence data collection33 showed that only 76 of the 557 (13.6%) Chinese BRCA1 had matches in 2,476 (3%) BRCA1 variants in GnomAD and 97 of the 531 (18.2%) Chinese BRCA2 had matches in 3,674 variants (2.6%) in GnomAD. The results indicate that the vast majority of the Chinese BRCA variants were not present in the population data provided by current exome and whole‐genome sequencing studies (Table 1A and Table S7A, Supporting Information). For those with the matches, their abundance as judged by the East Asia population frequencies were mostly at lower levels [62/76 (81.5%) BRCA1 variants and 73/97 (72.3%) BRCA2 variants <0.001], highlighting their pathogenic potential.

Comparison of BRCA variants between Chinese and other populations

Total refers to the numbers in each reference database.

Proportion = Shared variants / total variants (557 in BRCA1 or 531 in BRCA2) in Chinese population.

Variants can be overlapped among different databases and populations.

Indian BRCA variants in Refs. 13, 22 were combined for the comparison.

Comparing the Chinese BRCA data with non‐Chinese BRCA data in these major BRCA databases shows that 38% of BRCA variants were present only in the Chinese population [186 (33.4%) of 557 BRCA1 variants and 226 (42.6%) of 531 BRCA2 variants] (Fig. 1
c, Table 1B). Of all databases used for the comparison, the ClinVar database had the highest matching rates of 60.1% and 53.4% in BRCA1 and BRCA2, respectively due to its large data collection.

Comparing the Chinese BRCA data with the recent CIMBA data enriching non‐Caucasian data shows that 17.6% of the 1,088 Chinese BRCA variants had matches [BRCA1: 106/557 (19%) and BRCA2: 86/531 (16.2%). There are 15 Chinese BRCA1 variants and 16 Chinese BRCA2 variants included in both our Chinese data and the CIMBA data. These shared Chinese variants were removed in order to know the similarity and differences between Chinese and non‐Chinese populations included in the CIMBA data. After the removal, the overall matched rate decreased to 14.9% [BRCA1: 91/557 (16.3%) and BRCA2: 71/531 (13.4%)] (Table 1C and Table S7B, Supporting Information).

Comparing the Chinese BRCA data with those from Latin America and the Caribbean found that 97.8% of total variants [544 (97.7%) BRCA1 variants and 520 (97.9%) BRCA2] were specific only to the Chinese population (Table 1D and Table S7C, Supporting Information).

Comparing the Chinese BRCA data with those in non‐Chinese Asian populations found that 78.6% of total variants [436 (78.2%) of BRCA1 and 420 (79.2%) of BRCA2] were present only in the Chinese population (Table 1E and Table S7D, Supporting Information).

Comparison shows that only 23 (4.1%) BRCA1 and two (0.4%) BRCA2 variants were shared between the Chinese and Indian populations (Table 1F and Table S7E, Supporting Information).

Through these extensive comparisons, we were able to determine that around 40% of the BRCA data in Chinese was explicit and absent from the current BRCA data derived from non‐Chinese population.

We compared the variant distribution across BRCA1 and BRCA2 exons between Chinese and non‐Chinese populations using the data from BIC database as a testing model. For both BRCA1 and BRCA2, the distribution of variants differed significantly in multiple exons between the two data sets (Fig. S2, Supporting Information). For BRCA1, the proportions of variants in exons 2, 11D, 16, 20 and 24 were higher in the BIC data than in the Chinese data, whereas the proportions in exons 11B and 11C in the Chinese data were significantly higher than in the BIC data; in BRCA2, the proportions in exons 11A, 25 and 27 were higher in the BIC data than in the Chinese data, whereas the proportions in exons 2, 11F, 14, 21 and 22 were higher in the Chinese data than in the BIC data, respectively. The results show the presence of differences of exon distribution between Chinese and non‐Chinese BRCA variants.

We compared single‐base changes and other variant types between Chinese and BIC data. The results show the significant differences present in multiple types of base changes between the two data sets, including G > A, C > T and delT in BRCA1, and A > C, delA, delG and delC in BRCA2. In BRCA1, delT had higher frequency in the Chinese population than in the BIC data (6.8% versus 3.5%, p < 0.0034); in BRCA2, delA, delG and delC were more frequent in the Chinese data than in the BIC data (11.5% versus 4.1%, 6.7% versus 1.6% and 5.1% versus 2.4%, respectively, p < 0.000, 0.000, 0.004, accordingly) (Table 2A). Significant differences were also present in the missense, nonsense, stop gain, splice variants and intronic variant types in both BRCA1 and BRCA2, and in frameshifts in BRCA2 (Table 2B).

Comparison of BRCA variation types between Chinese and BIC data

To comprimise naming differences between Chinese data and BIC data for comparison, the names in Chinese data were converted as: Frameshift deletion and frameshift insertion were combined as frameshift; Nonsynonymous SNV and missense were combined as Missense, Nonframeshift deletion was converted as In Frame Deletion, Splice site to Splice. Statistical comparison was performed using Fisher's exact test.

We compared the clinical categories: pathogenic, likely pathogenic, uncertain significance, likely benign, and benign between Chinese BRCA data and BIC data. The results showed significant differences in multiple categories between the two data sets. For example, 12.7% of BRCA1 variants in Chinese were variants of uncertain significance, which was much higher than the value of 0.57% in the BIC data (p < 0.0000); the proportion of BRCA2 pathogenic variants in the Chinese population was also significantly higher than that in the BIC data (p < 0.0005), and the proportions of unclassified variants in both BRCA1 and BRCA2 were much higher in the BIC data than in the Chinese variants (Table 3).

Comparison of Clinical classification between Chinese and BIC BRCA variants

Firstly, we checked if the Chinese variant data contained any BRCA founder mutations known in other populations, including BRCA1 c.66_67delAG (185delAG), c.5263_5264insC (5382insC), and BRCA2 c.5946delT (6174delT) in Ashkenazi Jews;34
BRCA1 c.‐58C > G (C61G), c.4153delA (c.4035delA), and c.5263_5264insC (5382insC) in Poles;35
BRCA1 c.303 T > G, c.1623dupG, c.4122_4123delTG and c.5324 T > G in Africans;14
BRCA1 ex9‐12del in Mexicans,36
BRCA1 390C > A in Koreans and Japanese; and BRCA1 c.470_471delCT, BRCA2 c.7480C > T in Koreans.37 We observed that many of these founder mutations were either absent or present at low prevalence, hence they could not be considered as founder mutations in the Chinese population. Secondly, we searched for potential BRCA founder mutation candidates in the Chinese population by referring to 1) the abundance as calculated by the total number of variant carriers divided by the total number of individuals tested although haplotype data will be required to finally determine the true founder mutations; and 2) additional criteria for removing the variants unlikely to be founder mutations in order to focus on the variants as potential candidates: more than 100 tested individuals (a founder mutation should have a reasonable prevalence in a given population. 100 was set as a minimal cut‐off for the population size); at least two detected variant carriers (as a founder mutation, it cannot be only present in a single individual. Therefore, 2 cases were set as the minimal number of mutation carriers. In this way, all variants detected only in single individuals will be eliminated); carrier frequency > 1% (a precondition for founder mutation as pathogenic one is its lower population frequency in population. We set 1% of mutation carrier as the minimal cut‐off to eliminate these with high population frequency, which are mostly normal polymorphism); and variants in the categories of pathogenic, likely pathogenic, uncertain significance, or unclassified (Founder mutations must be pathogenic. Restricted the candidates to these classes will narrow down the founder mutation candidates by eliminating the benign and likely benign variants as they do not increase cancer risk). Using these conditions, we tested whether the data could support the Chinese BRCA founder mutations proposed by previous studies including BRCA1 c.981_982delAT (1100delAT),38
BRCA1 1081delT(1081delG),39
BRCA1c.5154G > A and BRCA1c.5468‐1del8;40 and BRCA2 c.3109C > T, BRCA2 c.7436_7805del370 and BRCA2c.9097_9098insA.38 With an exception for BRCA1c.5154G > A variant, our data do not support the above‐mentioned variants as the founder mutations. Next, we searched for high‐frequency variants meeting the same criteria as above and identified a total of 16 pathogenic, two likely pathogenic, 22 unclassified variants, and five of uncertain significance in BRCA1; and ten pathogenic, one likely pathogenic, and 11 unclassified variants in BRCA2 (Table 4 and Table S8, Supporting Information), respectively. The higher prevalence and clinical pathogenicity of these variants supported them as potential candidates for BRCA1 and BRCA2 founder mutations in the Chinese population. Despite of the higher prevalence, the unclassified variants or those of uncertain significance cannot be regarded as potential founder mutations unless their pathogenicity is determined.

High frequenct BRCA variants in Chinese population

Case tested refers to the total cases included in each study.

The most significant variants found were the following:
BRCA1 c.5154G > A. This variant had the highest prevalence of 5.6% (15 out of 266 detected by five studies). It is a stop‐gain pathogenic mutation, present in the BIC, BMD, LOVD, ClinVar, BED databases and was reported as a Chinese founder mutation by a previous study (31).
BRCA1 c.4258C > T. This variant had a prevalence of 2.5% (3 out of 118), is pathogenic, and is present in the BIC, BMD, LOVD, ClinVar and BED databases.
BRCA1 c.3296delC. This variant had a prevalence of 2.4% (3 out of 124), is pathogenic, and is present in the BIC, BMD, ClinVar, BED and LOVD databases.
BRCA1 c.5533_5540delATTGGGCA/delTACCAGTG. This variant had a prevalence of 2.5% (3 out of 125), is pathogenic, and is absent from other BRCA databases.
BRCA2 c.7655_7658delTTAA. This variant had a prevalence of 3.7% (4 out of 107, reported by four studies), is a pathogenic frameshift deletion and is present in the BIC, ClinVar and BED databases.
BRCA2 c.2636_2637delCT. This variant had a prevalence of 2.2% (4 out of 180), is pathogenic, and is present in the BIC, BMD, ClinVar and BED databases.
BRCA2 c.2339C > G. This variant had a prevalence of 2% (2 out of 99), is pathogenic, and is present in the BMD, ClinVar, BED and LOVD databases.

BRCA1 c.5154G > A. This variant had the highest prevalence of 5.6% (15 out of 266 detected by five studies). It is a stop‐gain pathogenic mutation, present in the BIC, BMD, LOVD, ClinVar, BED databases and was reported as a Chinese founder mutation by a previous study (31).

BRCA1 c.4258C > T. This variant had a prevalence of 2.5% (3 out of 118), is pathogenic, and is present in the BIC, BMD, LOVD, ClinVar and BED databases.

BRCA1 c.3296delC. This variant had a prevalence of 2.4% (3 out of 124), is pathogenic, and is present in the BIC, BMD, ClinVar, BED and LOVD databases.

BRCA1 c.5533_5540delATTGGGCA/delTACCAGTG. This variant had a prevalence of 2.5% (3 out of 125), is pathogenic, and is absent from other BRCA databases.

BRCA2 c.7655_7658delTTAA. This variant had a prevalence of 3.7% (4 out of 107, reported by four studies), is a pathogenic frameshift deletion and is present in the BIC, ClinVar and BED databases.

BRCA2 c.2636_2637delCT. This variant had a prevalence of 2.2% (4 out of 180), is pathogenic, and is present in the BIC, BMD, ClinVar and BED databases.

BRCA2 c.2339C > G. This variant had a prevalence of 2% (2 out of 99), is pathogenic, and is present in the BMD, ClinVar, BED and LOVD databases.

Although these high‐frequent BRCA mutations suggests the presence of certain potential founder mutations in Chinese population, it is also obvious from the data that there are unlikely to be high‐prevalence founder mutations in the Chinese population as these in the Ashkenazi Jewish population. However, the Chinese population is composed of highly heterogeneous ethnic groups with different genetic features, and even the dominant Han ethnic group is not homogeneous. Therefore, there remains a possibility for the presence of certain high‐prevalence founder mutations in certain specific ethnic groups, and in certain populations located at specific geographical locations.

It is of interest to know whether BRCA ethnic‐specificity exists within the Chinese population, given the fact that it has 56 ethnic groups with divergent genetic backgrounds.41 We tested this possibility by using Uygur group as a model. Uygur group is the largest minority group in Xinjiang, northwestern China, with its unique genetic features.41 A series of BRCA studies have been carried out in Uygur group, with the identification of 70 BRCA variants. Of these, 20 BRCA1 variants and 6 BRCA2 variant were present only in the Uygur group. Of the 26 Uygur‐specific BRCA1 variants, one was likely pathogenic, one was uncertain significance, one was likely benign, and 23 remain unclassified (Table S9, Supporting Information). The results indicate the presence of Uygur‐specific BRCA variants within the Chinese population.

### Collection of Chinese BRCA data

The majority of Chinese BRCA variation data was gathered from mainland China (90.2%), and the remainder were from Hong Kong (6.4%), Taiwan (1.6%), Singapore (0.9) and Malaysia (0.4%) (Table S1A, S1B, Supporting Information). In total, 31,689 cases with Chinese ethnicity were tested for BRCA mutations and their results were reported between 1999 and 2017; of these, 69.3% were reported between 2016 and 2017. Nearly all BRCA variants were from breast and/or ovarian cancer patients under different clinical criteria, except for four variants that were derived from 1,043 healthy control individuals.32

We collected and summarized the clinical information reported from each of the reference studies (Fig. S1 and Table S2, Supporting Information). The data showed several unique features as follows: 1). BMI is considered as one of the risk factors for breast cancer. However, according to our data, 83.7% of Chinese patients were within the normal range of 18.5–22.9 BMI and 10.5% were seen to have even lower than 18.5 BMI, indicating that obesity had no significant role in increasing risk for breast and ovarian cancer in these Chinese patients; 2). Family history is considered as a high‐risk factor for hereditary breast cancer. However, 72.3% of patients from our study did not report any family history of cancer, suggesting that family history was not an essential factor in this disease cohort; 3). Stage II breast cancer cases accounted for 66.2% whereas stage I cases only 6.3%, indicating the lack of earlier diagnosis.

A total of 3,791 BRCA mutation carriers were identified (12%, 2,123 BRCA1 and 1,688 BRCA2), of which 1,978 carriers (52.2%, 990 BRCA1 and 988 BRCA2) were within the clinically reported mutation categories of pathogenic or likely pathogenic. By standardizing and re‐annotating all the variants following Human Genome Variation Society (HGVS)‐recommended nomenclature, we identified a total of 1,088 distinct BRCA variants (557 in BRCA1 and 531 in BRCA2) of which 519 (47.7%, 278 in BRCA1 and 241 in BRCA2) were recurrent, and 50% were either detected or validated by Sanger sequencing (Tables S4–S5, Supporting Information). Except the 26 variants specifically from Uygur ethnic population, all variants were from Han Chinese population. We developed the dbBRCA‐Chinese database as an open source to host the entire set of BRCA variants and their annotation information (https://dbBRCA-Chinese.fhs.umac.mo).

### General features of Chinese BRCA data

The 1,088 BRCA variants were distributed with different frequencies among the 3,791 BRCA variant carriers. The age distribution data show that 52.3% of the cancer developed at early age (<40 years) in these BRCA variant carriers (Table S1C, Supporting Information). About half of the variants were detected only once in single individuals (46.5% in BRCA1, 54.6% in BRCA2), and the rest variants were distributed between 2 to 100 individuals with increased frequencies (Table S1D, Supporting Information).

We classified the BRCA variants into five classes—pathogenic, likely pathogenic, uncertain significance, likely benign, and benign—following the American College of Medical Genetics and Genomics (ACMG) standards and guidelines. Pathogenic and likely pathogenic variants, which are clinically reportable, accounted for 46% of BRCA1 variants and 52% of BRCA2 variants. Such higher rates do not necessarily imply high BRCA variation rate in Chinese breast and ovarian cancer patient population but reflect the fact that the cancer patients included in many of the studies were selected from high‐risk patients of strong family history or early age of cancer development. Importantly, 13% and 8% of BRCA1 and BRCA2 were classified as of uncertain significance, and 30% of both BRCA1 and BRCA2 variants remained as unclassified variants (Fig. 1
a). This fact indicates that much effort needs to be made in order to determine the function of these BRCA variants existing in Chinese population.

BRCA data. (a) Clinical classification of Chinese BRCA data as pathogenic, likely pathogenic, uncertain significance, likely benign, benign and unclassified. The pathogenic and likely pathogenic variants accounted for 49.5%. (b) Relationship between population sizes and their contribution to current BRCA data. The proportions of different human ethnic populations were from the database (http://www.worldometers.info/world‐population/), the ethnic origins of BRCA data were from different BRCA databases. It shows that the current BRCA data is not proportional to the human ethnic populations. (c) Comparison of BRCA data between Chinese and non‐Chinese populations. A total of 557 BRCA1 and 531 BRCA2 Chinese variants were compared to 6,344 BRCA1 and 8,886 BRCA2 non‐Chinese variants compiled from all existing BRCA databases. The results show that 38% of Chinese BRCA variants were present only in the Chinese population.

### Age and abundance of variants

The 1,088 BRCA variants were distributed with different frequencies among the 3,791 BRCA variant carriers. The age distribution data show that 52.3% of the cancer developed at early age (<40 years) in these BRCA variant carriers (Table S1C, Supporting Information). About half of the variants were detected only once in single individuals (46.5% in BRCA1, 54.6% in BRCA2), and the rest variants were distributed between 2 to 100 individuals with increased frequencies (Table S1D, Supporting Information).

### Pathogenic and nonpathogenic variants

We classified the BRCA variants into five classes—pathogenic, likely pathogenic, uncertain significance, likely benign, and benign—following the American College of Medical Genetics and Genomics (ACMG) standards and guidelines. Pathogenic and likely pathogenic variants, which are clinically reportable, accounted for 46% of BRCA1 variants and 52% of BRCA2 variants. Such higher rates do not necessarily imply high BRCA variation rate in Chinese breast and ovarian cancer patient population but reflect the fact that the cancer patients included in many of the studies were selected from high‐risk patients of strong family history or early age of cancer development. Importantly, 13% and 8% of BRCA1 and BRCA2 were classified as of uncertain significance, and 30% of both BRCA1 and BRCA2 variants remained as unclassified variants (Fig. 1
a). This fact indicates that much effort needs to be made in order to determine the function of these BRCA variants existing in Chinese population.

BRCA data. (a) Clinical classification of Chinese BRCA data as pathogenic, likely pathogenic, uncertain significance, likely benign, benign and unclassified. The pathogenic and likely pathogenic variants accounted for 49.5%. (b) Relationship between population sizes and their contribution to current BRCA data. The proportions of different human ethnic populations were from the database (http://www.worldometers.info/world‐population/), the ethnic origins of BRCA data were from different BRCA databases. It shows that the current BRCA data is not proportional to the human ethnic populations. (c) Comparison of BRCA data between Chinese and non‐Chinese populations. A total of 557 BRCA1 and 531 BRCA2 Chinese variants were compared to 6,344 BRCA1 and 8,886 BRCA2 non‐Chinese variants compiled from all existing BRCA databases. The results show that 38% of Chinese BRCA variants were present only in the Chinese population.

### Ethnic origins of current BRCA data

We compiled ethnic origins of currently existing BRCA variation data to know the status of BRCA study across ethnic human populations. By combining the BRCA variants from the Breast Cancer Information Core (BIC), ClinVar, BRCA Exchange Database (BED), BRCA1 and BRCA2 Mutation Database (BMD), Leiden Open Variation Database (LOVD) and ENIGMA (Evidence‐based Network for the Interpretation of Germline Mutant Alleles), we generated a single BRCA variation data set containing 6,343 distinct BRCA1 and 8,884 distinct BRCA2 variants (S6A). Classification of the ethnic origins of these variants showed that 62% were from Caucasian populations and 15% from the Ashkenazi Jewish population. The remaining 23% were originated from non‐Chinese Asian (13%), Latino (5%), African (3%), and Chinese (2%) populations (Fig. 1
b and Tables S6B, S6C, Supporting Information). Of the 5,925 CIMBA BRCA variants data with defined ethnicities, 80.2% were from Caucasian, 2.4% from Ashkenazi Jews, 10.8% from Asian, 3.8% from Hispanic and 2.8% from African American.15 The analysis shows that the current BRCA variant data contains very limited information from non‐Caucasian populations.

We performed the multi‐layer analyses to investigate the similarities and differences in BRCA variation between Chinese and non‐Chinese populations.

### Comparison with existing BRCA data

In order to determine the similarities and differences of BRCA variation between Chinese and non‐Chinese populations, we made a comprehensive comparison between Chinese and any available non‐Chinese BRCA data as represented below. The rationale for comparing with each database are:GnomAD: It contains extensive normal population variation data collected from human population by the largest exome and whole‐genome sequencing projects. The comparison aimed to determine the similarities and differences of the normal BRCA variants present between Chinese and non‐Chinese populations;BIC, BED, BMD, ClinVar, ENIGMA and LOVD: These are the major BRCA databases, with the data mostly derived from Caucasians as indicated by our analysis (Fig. 1
b). Comparison with these databases aimed to determine the similarities and differences of BRCA mutation between Chinese and Caucasian (mostly) populations;CIMBA data: It contains BRCA data collected from 49 countries across six continents. The comparison aimed to determine the similarities and differences between Chinese and worldwide non‐Chinese populations including more non‐Caucasian data;Latin America and the Caribbean data: The data were from Latin American population of Argentina, Bahamas, Brazil, Chile, Colombia, Costa Rica, Cuba, Mexico, Peru, Puerto Rico, Uruguay, Venezuela and the Hispanic population in the United States.12 The comparison aimed to determine the similarities and differences between Chinese and Latin America populations;Non‐Chinese Asian populations: BRCA data are available from Bangladeshi, Filipino, Iranian, Israeli, Japanese, Korean, Lebanese, Malay, Oman, Pakistani, Sri Lankan, Thai and Turkish populations.13 These populations were genetically and geographically closer to the Chinese population than other non‐Chinese populations. The comparison aimed to determine the similarities and differences between Chinese and the non‐Chinese Asian populations;Indian population: India has the 2nd largest population in the world, with highly diversified genetic background. Several large‐scale BRCA studies were reported recently with substantial BRCA data collected from the Indian patients.13, 31 The comparison aimed to determine the similarities and differences between Chinese and Indian populations, the two largest populations in the world.

GnomAD: It contains extensive normal population variation data collected from human population by the largest exome and whole‐genome sequencing projects. The comparison aimed to determine the similarities and differences of the normal BRCA variants present between Chinese and non‐Chinese populations;

BIC, BED, BMD, ClinVar, ENIGMA and LOVD: These are the major BRCA databases, with the data mostly derived from Caucasians as indicated by our analysis (Fig. 1
b). Comparison with these databases aimed to determine the similarities and differences of BRCA mutation between Chinese and Caucasian (mostly) populations;

CIMBA data: It contains BRCA data collected from 49 countries across six continents. The comparison aimed to determine the similarities and differences between Chinese and worldwide non‐Chinese populations including more non‐Caucasian data;

Latin America and the Caribbean data: The data were from Latin American population of Argentina, Bahamas, Brazil, Chile, Colombia, Costa Rica, Cuba, Mexico, Peru, Puerto Rico, Uruguay, Venezuela and the Hispanic population in the United States.12 The comparison aimed to determine the similarities and differences between Chinese and Latin America populations;

Non‐Chinese Asian populations: BRCA data are available from Bangladeshi, Filipino, Iranian, Israeli, Japanese, Korean, Lebanese, Malay, Oman, Pakistani, Sri Lankan, Thai and Turkish populations.13 These populations were genetically and geographically closer to the Chinese population than other non‐Chinese populations. The comparison aimed to determine the similarities and differences between Chinese and the non‐Chinese Asian populations;

Indian population: India has the 2nd largest population in the world, with highly diversified genetic background. Several large‐scale BRCA studies were reported recently with substantial BRCA data collected from the Indian patients.13, 31 The comparison aimed to determine the similarities and differences between Chinese and Indian populations, the two largest populations in the world.

Matching the Chinese BRCA data with those in GnomAD from the largest exome and whole‐genome sequence data collection33 showed that only 76 of the 557 (13.6%) Chinese BRCA1 had matches in 2,476 (3%) BRCA1 variants in GnomAD and 97 of the 531 (18.2%) Chinese BRCA2 had matches in 3,674 variants (2.6%) in GnomAD. The results indicate that the vast majority of the Chinese BRCA variants were not present in the population data provided by current exome and whole‐genome sequencing studies (Table 1A and Table S7A, Supporting Information). For those with the matches, their abundance as judged by the East Asia population frequencies were mostly at lower levels [62/76 (81.5%) BRCA1 variants and 73/97 (72.3%) BRCA2 variants <0.001], highlighting their pathogenic potential.

Comparison of BRCA variants between Chinese and other populations

Total refers to the numbers in each reference database.

Proportion = Shared variants / total variants (557 in BRCA1 or 531 in BRCA2) in Chinese population.

Variants can be overlapped among different databases and populations.

Indian BRCA variants in Refs. 13, 22 were combined for the comparison.

Comparing the Chinese BRCA data with non‐Chinese BRCA data in these major BRCA databases shows that 38% of BRCA variants were present only in the Chinese population [186 (33.4%) of 557 BRCA1 variants and 226 (42.6%) of 531 BRCA2 variants] (Fig. 1
c, Table 1B). Of all databases used for the comparison, the ClinVar database had the highest matching rates of 60.1% and 53.4% in BRCA1 and BRCA2, respectively due to its large data collection.

Comparing the Chinese BRCA data with the recent CIMBA data enriching non‐Caucasian data shows that 17.6% of the 1,088 Chinese BRCA variants had matches [BRCA1: 106/557 (19%) and BRCA2: 86/531 (16.2%). There are 15 Chinese BRCA1 variants and 16 Chinese BRCA2 variants included in both our Chinese data and the CIMBA data. These shared Chinese variants were removed in order to know the similarity and differences between Chinese and non‐Chinese populations included in the CIMBA data. After the removal, the overall matched rate decreased to 14.9% [BRCA1: 91/557 (16.3%) and BRCA2: 71/531 (13.4%)] (Table 1C and Table S7B, Supporting Information).

Comparing the Chinese BRCA data with those from Latin America and the Caribbean found that 97.8% of total variants [544 (97.7%) BRCA1 variants and 520 (97.9%) BRCA2] were specific only to the Chinese population (Table 1D and Table S7C, Supporting Information).

Comparing the Chinese BRCA data with those in non‐Chinese Asian populations found that 78.6% of total variants [436 (78.2%) of BRCA1 and 420 (79.2%) of BRCA2] were present only in the Chinese population (Table 1E and Table S7D, Supporting Information).

Comparison shows that only 23 (4.1%) BRCA1 and two (0.4%) BRCA2 variants were shared between the Chinese and Indian populations (Table 1F and Table S7E, Supporting Information).

Through these extensive comparisons, we were able to determine that around 40% of the BRCA data in Chinese was explicit and absent from the current BRCA data derived from non‐Chinese population.

### GnomAD

Matching the Chinese BRCA data with those in GnomAD from the largest exome and whole‐genome sequence data collection33 showed that only 76 of the 557 (13.6%) Chinese BRCA1 had matches in 2,476 (3%) BRCA1 variants in GnomAD and 97 of the 531 (18.2%) Chinese BRCA2 had matches in 3,674 variants (2.6%) in GnomAD. The results indicate that the vast majority of the Chinese BRCA variants were not present in the population data provided by current exome and whole‐genome sequencing studies (Table 1A and Table S7A, Supporting Information). For those with the matches, their abundance as judged by the East Asia population frequencies were mostly at lower levels [62/76 (81.5%) BRCA1 variants and 73/97 (72.3%) BRCA2 variants <0.001], highlighting their pathogenic potential.

Comparison of BRCA variants between Chinese and other populations

Total refers to the numbers in each reference database.

Proportion = Shared variants / total variants (557 in BRCA1 or 531 in BRCA2) in Chinese population.

Variants can be overlapped among different databases and populations.

Indian BRCA variants in Refs. 13, 22 were combined for the comparison.

### BIC, ClinVar, BED, BMD, LOVD and ENIGMA

Comparing the Chinese BRCA data with non‐Chinese BRCA data in these major BRCA databases shows that 38% of BRCA variants were present only in the Chinese population [186 (33.4%) of 557 BRCA1 variants and 226 (42.6%) of 531 BRCA2 variants] (Fig. 1
c, Table 1B). Of all databases used for the comparison, the ClinVar database had the highest matching rates of 60.1% and 53.4% in BRCA1 and BRCA2, respectively due to its large data collection.

### CIMBA

Comparing the Chinese BRCA data with the recent CIMBA data enriching non‐Caucasian data shows that 17.6% of the 1,088 Chinese BRCA variants had matches [BRCA1: 106/557 (19%) and BRCA2: 86/531 (16.2%). There are 15 Chinese BRCA1 variants and 16 Chinese BRCA2 variants included in both our Chinese data and the CIMBA data. These shared Chinese variants were removed in order to know the similarity and differences between Chinese and non‐Chinese populations included in the CIMBA data. After the removal, the overall matched rate decreased to 14.9% [BRCA1: 91/557 (16.3%) and BRCA2: 71/531 (13.4%)] (Table 1C and Table S7B, Supporting Information).

### Latin America and the Caribbean data

Comparing the Chinese BRCA data with those from Latin America and the Caribbean found that 97.8% of total variants [544 (97.7%) BRCA1 variants and 520 (97.9%) BRCA2] were specific only to the Chinese population (Table 1D and Table S7C, Supporting Information).

### Non‐Chinese Asian populations

Comparing the Chinese BRCA data with those in non‐Chinese Asian populations found that 78.6% of total variants [436 (78.2%) of BRCA1 and 420 (79.2%) of BRCA2] were present only in the Chinese population (Table 1E and Table S7D, Supporting Information).

### Indian population

Comparison shows that only 23 (4.1%) BRCA1 and two (0.4%) BRCA2 variants were shared between the Chinese and Indian populations (Table 1F and Table S7E, Supporting Information).

Through these extensive comparisons, we were able to determine that around 40% of the BRCA data in Chinese was explicit and absent from the current BRCA data derived from non‐Chinese population.

### Comparison in exon distribution

We compared the variant distribution across BRCA1 and BRCA2 exons between Chinese and non‐Chinese populations using the data from BIC database as a testing model. For both BRCA1 and BRCA2, the distribution of variants differed significantly in multiple exons between the two data sets (Fig. S2, Supporting Information). For BRCA1, the proportions of variants in exons 2, 11D, 16, 20 and 24 were higher in the BIC data than in the Chinese data, whereas the proportions in exons 11B and 11C in the Chinese data were significantly higher than in the BIC data; in BRCA2, the proportions in exons 11A, 25 and 27 were higher in the BIC data than in the Chinese data, whereas the proportions in exons 2, 11F, 14, 21 and 22 were higher in the Chinese data than in the BIC data, respectively. The results show the presence of differences of exon distribution between Chinese and non‐Chinese BRCA variants.

### Comparison in base changes and variant types

We compared single‐base changes and other variant types between Chinese and BIC data. The results show the significant differences present in multiple types of base changes between the two data sets, including G > A, C > T and delT in BRCA1, and A > C, delA, delG and delC in BRCA2. In BRCA1, delT had higher frequency in the Chinese population than in the BIC data (6.8% versus 3.5%, p < 0.0034); in BRCA2, delA, delG and delC were more frequent in the Chinese data than in the BIC data (11.5% versus 4.1%, 6.7% versus 1.6% and 5.1% versus 2.4%, respectively, p < 0.000, 0.000, 0.004, accordingly) (Table 2A). Significant differences were also present in the missense, nonsense, stop gain, splice variants and intronic variant types in both BRCA1 and BRCA2, and in frameshifts in BRCA2 (Table 2B).

Comparison of BRCA variation types between Chinese and BIC data

To comprimise naming differences between Chinese data and BIC data for comparison, the names in Chinese data were converted as: Frameshift deletion and frameshift insertion were combined as frameshift; Nonsynonymous SNV and missense were combined as Missense, Nonframeshift deletion was converted as In Frame Deletion, Splice site to Splice. Statistical comparison was performed using Fisher's exact test.

### Comparison in clinical categories

We compared the clinical categories: pathogenic, likely pathogenic, uncertain significance, likely benign, and benign between Chinese BRCA data and BIC data. The results showed significant differences in multiple categories between the two data sets. For example, 12.7% of BRCA1 variants in Chinese were variants of uncertain significance, which was much higher than the value of 0.57% in the BIC data (p < 0.0000); the proportion of BRCA2 pathogenic variants in the Chinese population was also significantly higher than that in the BIC data (p < 0.0005), and the proportions of unclassified variants in both BRCA1 and BRCA2 were much higher in the BIC data than in the Chinese variants (Table 3).

Comparison of Clinical classification between Chinese and BIC BRCA variants

### Comparison in founder mutations

Firstly, we checked if the Chinese variant data contained any BRCA founder mutations known in other populations, including BRCA1 c.66_67delAG (185delAG), c.5263_5264insC (5382insC), and BRCA2 c.5946delT (6174delT) in Ashkenazi Jews;34
BRCA1 c.‐58C > G (C61G), c.4153delA (c.4035delA), and c.5263_5264insC (5382insC) in Poles;35
BRCA1 c.303 T > G, c.1623dupG, c.4122_4123delTG and c.5324 T > G in Africans;14
BRCA1 ex9‐12del in Mexicans,36
BRCA1 390C > A in Koreans and Japanese; and BRCA1 c.470_471delCT, BRCA2 c.7480C > T in Koreans.37 We observed that many of these founder mutations were either absent or present at low prevalence, hence they could not be considered as founder mutations in the Chinese population. Secondly, we searched for potential BRCA founder mutation candidates in the Chinese population by referring to 1) the abundance as calculated by the total number of variant carriers divided by the total number of individuals tested although haplotype data will be required to finally determine the true founder mutations; and 2) additional criteria for removing the variants unlikely to be founder mutations in order to focus on the variants as potential candidates: more than 100 tested individuals (a founder mutation should have a reasonable prevalence in a given population. 100 was set as a minimal cut‐off for the population size); at least two detected variant carriers (as a founder mutation, it cannot be only present in a single individual. Therefore, 2 cases were set as the minimal number of mutation carriers. In this way, all variants detected only in single individuals will be eliminated); carrier frequency > 1% (a precondition for founder mutation as pathogenic one is its lower population frequency in population. We set 1% of mutation carrier as the minimal cut‐off to eliminate these with high population frequency, which are mostly normal polymorphism); and variants in the categories of pathogenic, likely pathogenic, uncertain significance, or unclassified (Founder mutations must be pathogenic. Restricted the candidates to these classes will narrow down the founder mutation candidates by eliminating the benign and likely benign variants as they do not increase cancer risk). Using these conditions, we tested whether the data could support the Chinese BRCA founder mutations proposed by previous studies including BRCA1 c.981_982delAT (1100delAT),38
BRCA1 1081delT(1081delG),39
BRCA1c.5154G > A and BRCA1c.5468‐1del8;40 and BRCA2 c.3109C > T, BRCA2 c.7436_7805del370 and BRCA2c.9097_9098insA.38 With an exception for BRCA1c.5154G > A variant, our data do not support the above‐mentioned variants as the founder mutations. Next, we searched for high‐frequency variants meeting the same criteria as above and identified a total of 16 pathogenic, two likely pathogenic, 22 unclassified variants, and five of uncertain significance in BRCA1; and ten pathogenic, one likely pathogenic, and 11 unclassified variants in BRCA2 (Table 4 and Table S8, Supporting Information), respectively. The higher prevalence and clinical pathogenicity of these variants supported them as potential candidates for BRCA1 and BRCA2 founder mutations in the Chinese population. Despite of the higher prevalence, the unclassified variants or those of uncertain significance cannot be regarded as potential founder mutations unless their pathogenicity is determined.

High frequenct BRCA variants in Chinese population

Case tested refers to the total cases included in each study.

The most significant variants found were the following:
BRCA1 c.5154G > A. This variant had the highest prevalence of 5.6% (15 out of 266 detected by five studies). It is a stop‐gain pathogenic mutation, present in the BIC, BMD, LOVD, ClinVar, BED databases and was reported as a Chinese founder mutation by a previous study (31).
BRCA1 c.4258C > T. This variant had a prevalence of 2.5% (3 out of 118), is pathogenic, and is present in the BIC, BMD, LOVD, ClinVar and BED databases.
BRCA1 c.3296delC. This variant had a prevalence of 2.4% (3 out of 124), is pathogenic, and is present in the BIC, BMD, ClinVar, BED and LOVD databases.
BRCA1 c.5533_5540delATTGGGCA/delTACCAGTG. This variant had a prevalence of 2.5% (3 out of 125), is pathogenic, and is absent from other BRCA databases.
BRCA2 c.7655_7658delTTAA. This variant had a prevalence of 3.7% (4 out of 107, reported by four studies), is a pathogenic frameshift deletion and is present in the BIC, ClinVar and BED databases.
BRCA2 c.2636_2637delCT. This variant had a prevalence of 2.2% (4 out of 180), is pathogenic, and is present in the BIC, BMD, ClinVar and BED databases.
BRCA2 c.2339C > G. This variant had a prevalence of 2% (2 out of 99), is pathogenic, and is present in the BMD, ClinVar, BED and LOVD databases.

BRCA1 c.5154G > A. This variant had the highest prevalence of 5.6% (15 out of 266 detected by five studies). It is a stop‐gain pathogenic mutation, present in the BIC, BMD, LOVD, ClinVar, BED databases and was reported as a Chinese founder mutation by a previous study (31).

BRCA1 c.4258C > T. This variant had a prevalence of 2.5% (3 out of 118), is pathogenic, and is present in the BIC, BMD, LOVD, ClinVar and BED databases.

BRCA1 c.3296delC. This variant had a prevalence of 2.4% (3 out of 124), is pathogenic, and is present in the BIC, BMD, ClinVar, BED and LOVD databases.

BRCA1 c.5533_5540delATTGGGCA/delTACCAGTG. This variant had a prevalence of 2.5% (3 out of 125), is pathogenic, and is absent from other BRCA databases.

BRCA2 c.7655_7658delTTAA. This variant had a prevalence of 3.7% (4 out of 107, reported by four studies), is a pathogenic frameshift deletion and is present in the BIC, ClinVar and BED databases.

BRCA2 c.2636_2637delCT. This variant had a prevalence of 2.2% (4 out of 180), is pathogenic, and is present in the BIC, BMD, ClinVar and BED databases.

BRCA2 c.2339C > G. This variant had a prevalence of 2% (2 out of 99), is pathogenic, and is present in the BMD, ClinVar, BED and LOVD databases.

Although these high‐frequent BRCA mutations suggests the presence of certain potential founder mutations in Chinese population, it is also obvious from the data that there are unlikely to be high‐prevalence founder mutations in the Chinese population as these in the Ashkenazi Jewish population. However, the Chinese population is composed of highly heterogeneous ethnic groups with different genetic features, and even the dominant Han ethnic group is not homogeneous. Therefore, there remains a possibility for the presence of certain high‐prevalence founder mutations in certain specific ethnic groups, and in certain populations located at specific geographical locations.

### Comparison within the Chinese population

It is of interest to know whether BRCA ethnic‐specificity exists within the Chinese population, given the fact that it has 56 ethnic groups with divergent genetic backgrounds.41 We tested this possibility by using Uygur group as a model. Uygur group is the largest minority group in Xinjiang, northwestern China, with its unique genetic features.41 A series of BRCA studies have been carried out in Uygur group, with the identification of 70 BRCA variants. Of these, 20 BRCA1 variants and 6 BRCA2 variant were present only in the Uygur group. Of the 26 Uygur‐specific BRCA1 variants, one was likely pathogenic, one was uncertain significance, one was likely benign, and 23 remain unclassified (Table S9, Supporting Information). The results indicate the presence of Uygur‐specific BRCA variants within the Chinese population.

### Discussion

By using the rich Chinese BRCA variation data as a representative of non‐Caucasian populations, our study provides solid evidence to conclude the presence of ethnic‐specific BRCA mutation. This likely reflects the human evolutionary history of genetic diversity and environmental adaptation.11 The variants shared between different ethnic populations were likely originated before their diversification, whereas the ethnic‐specific variants were likely generated after their diversification. Since BRCA reference data plays key roles in identifying the mutation carriers, lack of ethnic‐specific data in the present references implies that they have inadequate power in locating the mutation carriers with non‐Caucasian ethnic background. This is vividly exemplified by the presence of only 16 BRCA variants derived from mainland Chinese among the 3,791 BRCA variants in the BIC database.19 In order to identify the mutation carriers with various ethnic background, ethnic‐specific BRCA references need to be developed. Combined usage of both ethnic‐specific and existing BRCA reference databases should provide comprehensive identification of BRCA mutation carriers in different ethnic populations, a critical step towards precision medicine. Developing ethnic‐specific BRCA references will certainly be a challenge both scientifically and financially, but this task needs to be completed sooner or later for the sake of prevention of BRCA‐related cancers in the non‐Caucasian populations. The issue of ethnic‐specific germline mutation could also exist in other cancer predisposition genes. Experiences from developing ethnic‐specific BRCA references should provide a valuable example to address the same issue in these genes.

### Supporting information

Table S1A Original publications of BRCA1/2 variation in Chinese population

Table S1B. Origins of Chinese taken BRCA test between 1999 and 2017

Table S1C. Age distribution of breast/ovarian cancer with BRCA variants*

Table S1D. Abundance of BRCA variants in Chinese patients

Click here for additional data file.

Table S2 Clinical information for the patients reported by the studies*

Click here for additional data file.

Table S3 Exon positions in BRCA1 and BRCA2

Click here for additional data file.

Table S4 BRCA1 variants identified in Chinese population

Click here for additional data file.

Table S5 BRCA2 variants identified in Chinese population

Click here for additional data file.

Table S6A Total BRCA1/2 variants in current BRCA databases (BED, BIC, BMD, clinvar, LOVD, ENIGMA)

Click here for additional data file.

Table S7E Chinese BRCA variants shared with Indian population

Click here for additional data file.

Table S8 High frequent BRCA variants in Chinese population*

Click here for additional data file.

Table S9A Summary of BRCA1 variants in Uygur of Chinese population

Table S9B. BRCA variants reported in Xinjiang population

Click here for additional data file.

Fig S1 Status of BRCA studies in Chinese populations. (A) Regions where BRCA testing was conducted. The regions highlighted in yellow are the ones where the tests were conducted between 1999 and 2017; the populations of the Beijing and Shanghai regions are among the most extensively tested, but rapid progress was made in 2016–2017 regarding covering more inland areas (green area). (B) Number of individuals who underwent BRCA testing. BRCA testing in China was initiated in 1999. The numbers remained steady until 2015, then increased rapidly owing to the use of next‐generation sequencing technologies.

Click here for additional data file.

Fig S2 Comparison of exon distribution between Chinese and non‐Chinese BRCA variants (BIC). The rate of variation in each exon was calculated by dividing the number of variants in the exon by the total number of variants in BRCA1 or BRCA2. *p < 0.05 by Fisher's exact test as significance.

Click here for additional data file.

Appendix S1: Supporting Information

Click here for additional data file.



# SUPPLEMENTAL FILE 1: IJC-145-962.pdf

# Preparing to download ...

[HHS Vulnerability Disclosure](https://www.hhs.gov/vulnerability-disclosure-policy/index.html)