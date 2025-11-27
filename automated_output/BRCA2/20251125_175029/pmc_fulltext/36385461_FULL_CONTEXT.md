# MAIN TEXT

## Ethnic‐specificity, evolution origin and deleteriousness of Asian BRCA variation revealed by over 7500 BRCA variants derived from Asian population

### Abstract

AbstractPathogenic variation in BRCA1 and BRCA2 (BRCA) causes high risk of breast and ovarian cancer, and BRCA variation data are important markers for BRCA‐related clinical cancer applications. However, comprehensive BRCA variation data are lacking from the Asian population despite its large population size, heterogenous genetic background and diversified living environment across the Asia continent. We performed a systematic study on BRCA variation in Asian population including extensive data mining, standardization, annotation and characterization. We identified 7587 BRCA variants from 685 592 Asian individuals in 40 Asia countries and regions, including 1762 clinically actionable pathogenic variants and 4915 functionally unknown variants (https://genemutation.fhs.um.edu.mo/Asian-BRCA/). We observed the highly ethnic‐specific nature of Asian BRCA variants between Asian and non‐Asian populations and within Asian populations, highlighting that the current European descendant population‐based BRCA data is inadequate to reflect BRCA variation in the Asian population. We also provided archeological evidence for the evolutionary origin and arising time of Asian BRCA variation. We further provided structural‐based evidence for the deleterious variants enriched within the functionally unknown Asian BRCA variants. The data from our study provide a current view of BRCA variation in the Asian population and a rich resource to guide clinical applications of BRCA‐related cancer for the Asian population.

### INTRODUCTION

BRCA1 and BRCA2 (BRCA hereafter) are two of the well‐characterized tumor suppressor genes. BRCA repairs double‐strand DNA breaks through homologous recombination to maintain genome stability.
1
 The pathogenic variation in BRCA damages the function of BRCA, causing genome instability and leading to an increased risk of cancer mostly in breast and ovarian.
2
, 
3
 Since the establishment of the relationship between pathogenic BRCA variation and cancer risk, extensive efforts have been made in characterizing BRCA variation, resulting in the collection, characterization and classification of massive human BRCA variation data
4
, 
5
 that nearly 70 000 human BRCA variants have been identified (https://brcaexchange.org/factsheet, accessed on September 21, 2022). Such high‐frequent genetic variation is rarely seen in other cancer‐related human genes. Positive selection in human BRCA is considered an explanation.
8
 The rich knowledge of BRCA variation has greatly enhanced our understanding of the genetic basis of cancer, and BRCA variation data are used as valuable markers for clinical diagnosis, treatment and prevention of BRCA‐related cancer.
6
, 
7

It is well determined that human BRCA variation is highly ethnic‐specific in human populations.
9
, 
10
, 
11
, 
12
, 
13
, 
14
, 
15
, 
16
, 
17
, 
18
, 
19
, 
20
, 
21
 Like in many human genetic studies,
22
 however, the majority of current BRCA variation data were derived from the European descendant populations.
23
, 
24
, 
25
 With highly heterogeneous genetic background, Asia has a population size of 4.7 billion or 60% of the world population (https://www.worldometers.info/world-population/asia-population/), living in over 50 countries and regions across the geographically highly diversified Asia continent. It would be expected that BRCA variation in the Asian population may have distinguishing features from non‐Asian populations. While attempts have been made to address the issue,
12
, 
15
, 
26

BRCA variation in the Asian population remains poorly characterized and fragmented, in isolated populations. Lack of data prevents a deep understanding of the genetic basis of BRCA‐related cancer in the Asian population and heavily relies on the use of European descendant‐derived BRCA variation data to guide clinical cancer applications for the Asian population.

In the present study, we performed a systematic analysis for BRCA variation in the Asian population. Through extensive data mining, standardization, annotation and characterization, we identified a total of 7587 BRCA variants from 685 592 cancer and non‐cancer individuals in 40 Asia countries and regions (https://genemutation.fhs.um.edu.mo/Asian-BRCA/), which is the largest collection of BRCA variant data at the continental level besides the data from European descendant populations. Our analysis of the rich data reveals the highly ethnic‐specific nature of BRCA variation between Asian and non‐Asian populations and within Asian populations. We also determined the evolutionary origin and arising time for Asian BRCA variation. We further provided protein structure‐based evidence for the enrichment of deleterious variants in the abundant unclassified Asian BRCA variants.

### MATERIALS AND METHODS

We collected the Asian BRCA variants from both literature and genomics databases. Pre‐conditions were set to ensure that the variants originated from the Asian population: (a) The original study must be performed at an institution located in an Asian country (countries), with clearly indicated region, city and hospital within the country; (b) The sample sources used in the reports must be indicated with a geographic location within the country; (c) The sequence data from genome databases must be derived from Asian samples with sample ID and (d) The study performed by non‐Asia laboratories must indicate their sample sources of Asian origin with a geographic location within the country. For the collected variants, we performed extensive standardization following HGVS nomenclature
27
 and ACMG guidelines.
28
 The following reference sequences were used to annotate BRCA variants: BRCA1: genome assembly: hg38: 43044295 to 43 125 364, HGVSg: NC_000017.11, HGVSc: NM_007294.3, HGVSp: NP_ 009225.1, BIC cDNA: U14680.1, BIC protein: AAA73985.1; BRCA2: genome assembly: hg38: 32315508 to 32 400 268, HGVSg: NC_000013.11, HGVSc: NM_000059.3, HGVSp: NP_000050.2, BIC cDNA: U43746.1, BIC protein: AAB07223.1. BIC‐oriented position was converted to HGVSc defined position, with minus 119 for BRCA1 and minus 228 for BRCA2. HGVSc was converted to HGVSg position in hg38 using the Position Converter in Mutalyzer
29
 and validated using the Name Checker in Mutalyzer. LiftOver
30
 was used for genome position conversion from genome assembly hg19 to hg38 and for extracting the reference sequences at the corresponding genome position. Variants were annotated using ANNOVAR
31
 referring to multiple databases including RefGene, dbSNP (version 150), gnomAD exome collection (v2.1.1)
32
 and gnomAD genome collection (v3.0) and dbNSFP (v4.2a).
33
 Minor Allele Frequency for variants was from gnomAD exome collection and gnomAD genome collection. Nonsynonymous frameshift variant was annotated with dbNSFP containing 19 protein prediction and 9 evolution analysis tools. Clinical significance was annotated with ClinVar (accessed on March 2, 2022)
24
; The variants without matches in ClinVar were classified as unclassified. The novelty of variants was determined by comparing them to BRCA Exchange (accessed on March 8, 2022). The compliance with HGVS nomenclature for all variants was validated by Mutalyzer (https://mutalyzer.nl).

An open‐access database “dbBRCA‐Asian” was constructed to host the Asian BRCA variation data and their annotation information. The database was divided into a front‐end and a back end. The back end was implemented with LAMP stack including Linux operating system (CentOS 7, https://www.centos.org/), Apache HTTP Server (Apache 2.4.52, https://httpd.apache.org/), MySQL relational database management system (5.6.50, https://www.mysql.com/) and PHP (PHP 7.3, https://www.php.net/) programming language, for managing user's search and retrieval requests. The database implemented a user interface using web front‐end languages (HTML, CSS and JavaScript). HTML presented the contents for users and received the data entered by the user to the back‐end, CSS managed the layout typography of the page and adjusted the style of the elements in HTML, and JavaScript managed the front‐end communication.

The Asian population was divided into five regions of East Asia (China mainland, Hong Kong, Japan, Korea, Mongolia, Macau, North Korea, Taiwan), South‐East Asia (Brunei, Cambodia, East Timor, Indonesia, Laos, Malaysia, Myanmar, Philippines, Singapore, Thailand, Vietnam), South Asia (Bangladesh, Bhutan, India, Maldives, Nepal, Pakistan, Sri Lanka), West Asia (Arabic, Armenia, Azerbaijan, Bahrain, Georgia, Iran, Iraq, Israel, Jordan, Kuwait, Lebanon, Oman, Palestine, Qatar, Saudi Arabia, Syria, Turkey, United Arab Emirates, Yemen) and Central Asia (Afghanistan, Kazakhstan, Kyrgyzstan, Russia, Tajikistan, Turkmenistan, Uzbekistan). Saturation of BRCA variant detection was estimated by comparing the pathogenic/likely pathogenic variant distribution in BRCA1 BRCT1/2 domain (aa position: 1650‐1724/ 1758‐1842) and BRCA2 DNA binding domain (aa position: 250‐500) between the Asian variant data and these reported in the two recent studies.
34
, 
35
 Venn was used to visualize the comparative results using the R package “VennDiagram.”
36

BRCA sequences from 4817 ancient human individuals dated between 80 000 to 100 BP (before present), 8 Neanderthal genomes dated between 38 310 and 80 000 BP and a Denisovan genome dated 78 000 BP were extracted from literature and the Allen Ancient DNA Resource (AADR) database (https://reich.hms.harvard.edu/allen‐ancient‐dna‐resource‐aadr‐downloadable‐genotypes‐present‐day‐and‐ancient‐dna‐data, accessed on April 22, 2022, V50.0). We also searched Google Scholar, PubMed and Web of Sciences with keywords “Ancient human,” “Ancient DNA,” “Ancient genome,” “Ancient sequence,” “Ancient variant” and “Ancient mutation” to identify potential BRCA variants in ancient humans. The mapping process followed the procedures (Li et al, 2022).

We collected BRCA founder variants from the Asian population reported in the literature using the keywords of “BRCA1,” “BRCA2,” “founder mutation,” “variant,” “haplotype,” “Asian” or each Asian country or region. Only the variants with haplotype evidence were considered founder variants.

The RPMDS performance followed the procedures.
37
 Briefly, BRCA1 RING domain (PDB ID: 1JM7, 1‐103 residues),
38
 BRCA1 BRCT domain (PDB ID: 1JNX, 1649‐1859 residues)
39
 and BRCA2 BRC4 domain (PDB ID: 1N0W, 1519‐1551 residues)
40
 were used as the templates to build the mutant structure for each selected missense variant following the procedure described.
41
, 
42
, 
43
 The mutant structure was used as the starting configuration for MD simulations using GROMACS.
44
 Each simulation held approximately 99 000 atoms for the BRCA1 RING domain and BRCT domain, and 125 000 atoms for the BRCA2 BRC4 domain. Three independent simulations for the wildtype BRCA1 RING domain were performed for 100 ns.
45
 Hydrogen bonds were constrained at equilibrium length by using the LINC algorithm.
46
 The final 100 frames of the trajectories were used for variant classification. A total of 21 wildtypes in 10 ns‐frames and 33 pathogenic variants in BRCA1 RING domain, 11 benign and pathogenic variants in BRCA1 BRCT domain and 9 benign and pathogenic variants in BRCA2 BRC4 domain were used to determine the cut‐off values to classify the unknown variants into deleterious or tolerated variants. BRCA1 RING domain used wild‐type bases for benign due to the lack of benign variants. The structure deviation of BRCA1 RING domain pathogenic variants was 38.5% ± 5.8%, with a range of 36.4% and 40.3%; the structural deviation for BRCA1 BRCT pathogenic variants had a mean of 25.9% ± 8.9%, with a range of 18.6% and 33.3%; the structural deviation for BRCA2 BRC4 pathogenic variants had the mean of 23.5% ± 4.9%, with a range of 20.1% and 28.4%. The cut‐off values were set at 33.3, 30.3 and 22.9 for BRCA1 RING, BRCA1 BRCT and BRCA2 BRC4 domains, respectively. Wilcoxon test was used to test the difference between pathogenic variant and benign variant distribution in BRCA1 RING domain, BRCA1 BRCT1/2 domain and BRCA2 DNA binding domain.
41

### Analysis of Asian 
BRCA
 variant data

We collected the Asian BRCA variants from both literature and genomics databases. Pre‐conditions were set to ensure that the variants originated from the Asian population: (a) The original study must be performed at an institution located in an Asian country (countries), with clearly indicated region, city and hospital within the country; (b) The sample sources used in the reports must be indicated with a geographic location within the country; (c) The sequence data from genome databases must be derived from Asian samples with sample ID and (d) The study performed by non‐Asia laboratories must indicate their sample sources of Asian origin with a geographic location within the country. For the collected variants, we performed extensive standardization following HGVS nomenclature
27
 and ACMG guidelines.
28
 The following reference sequences were used to annotate BRCA variants: BRCA1: genome assembly: hg38: 43044295 to 43 125 364, HGVSg: NC_000017.11, HGVSc: NM_007294.3, HGVSp: NP_ 009225.1, BIC cDNA: U14680.1, BIC protein: AAA73985.1; BRCA2: genome assembly: hg38: 32315508 to 32 400 268, HGVSg: NC_000013.11, HGVSc: NM_000059.3, HGVSp: NP_000050.2, BIC cDNA: U43746.1, BIC protein: AAB07223.1. BIC‐oriented position was converted to HGVSc defined position, with minus 119 for BRCA1 and minus 228 for BRCA2. HGVSc was converted to HGVSg position in hg38 using the Position Converter in Mutalyzer
29
 and validated using the Name Checker in Mutalyzer. LiftOver
30
 was used for genome position conversion from genome assembly hg19 to hg38 and for extracting the reference sequences at the corresponding genome position. Variants were annotated using ANNOVAR
31
 referring to multiple databases including RefGene, dbSNP (version 150), gnomAD exome collection (v2.1.1)
32
 and gnomAD genome collection (v3.0) and dbNSFP (v4.2a).
33
 Minor Allele Frequency for variants was from gnomAD exome collection and gnomAD genome collection. Nonsynonymous frameshift variant was annotated with dbNSFP containing 19 protein prediction and 9 evolution analysis tools. Clinical significance was annotated with ClinVar (accessed on March 2, 2022)
24
; The variants without matches in ClinVar were classified as unclassified. The novelty of variants was determined by comparing them to BRCA Exchange (accessed on March 8, 2022). The compliance with HGVS nomenclature for all variants was validated by Mutalyzer (https://mutalyzer.nl).

### Database construction

An open‐access database “dbBRCA‐Asian” was constructed to host the Asian BRCA variation data and their annotation information. The database was divided into a front‐end and a back end. The back end was implemented with LAMP stack including Linux operating system (CentOS 7, https://www.centos.org/), Apache HTTP Server (Apache 2.4.52, https://httpd.apache.org/), MySQL relational database management system (5.6.50, https://www.mysql.com/) and PHP (PHP 7.3, https://www.php.net/) programming language, for managing user's search and retrieval requests. The database implemented a user interface using web front‐end languages (HTML, CSS and JavaScript). HTML presented the contents for users and received the data entered by the user to the back‐end, CSS managed the layout typography of the page and adjusted the style of the elements in HTML, and JavaScript managed the front‐end communication.

### Comparative analysis

The Asian population was divided into five regions of East Asia (China mainland, Hong Kong, Japan, Korea, Mongolia, Macau, North Korea, Taiwan), South‐East Asia (Brunei, Cambodia, East Timor, Indonesia, Laos, Malaysia, Myanmar, Philippines, Singapore, Thailand, Vietnam), South Asia (Bangladesh, Bhutan, India, Maldives, Nepal, Pakistan, Sri Lanka), West Asia (Arabic, Armenia, Azerbaijan, Bahrain, Georgia, Iran, Iraq, Israel, Jordan, Kuwait, Lebanon, Oman, Palestine, Qatar, Saudi Arabia, Syria, Turkey, United Arab Emirates, Yemen) and Central Asia (Afghanistan, Kazakhstan, Kyrgyzstan, Russia, Tajikistan, Turkmenistan, Uzbekistan). Saturation of BRCA variant detection was estimated by comparing the pathogenic/likely pathogenic variant distribution in BRCA1 BRCT1/2 domain (aa position: 1650‐1724/ 1758‐1842) and BRCA2 DNA binding domain (aa position: 250‐500) between the Asian variant data and these reported in the two recent studies.
34
, 
35
 Venn was used to visualize the comparative results using the R package “VennDiagram.”
36

### Evolution analysis

BRCA sequences from 4817 ancient human individuals dated between 80 000 to 100 BP (before present), 8 Neanderthal genomes dated between 38 310 and 80 000 BP and a Denisovan genome dated 78 000 BP were extracted from literature and the Allen Ancient DNA Resource (AADR) database (https://reich.hms.harvard.edu/allen‐ancient‐dna‐resource‐aadr‐downloadable‐genotypes‐present‐day‐and‐ancient‐dna‐data, accessed on April 22, 2022, V50.0). We also searched Google Scholar, PubMed and Web of Sciences with keywords “Ancient human,” “Ancient DNA,” “Ancient genome,” “Ancient sequence,” “Ancient variant” and “Ancient mutation” to identify potential BRCA variants in ancient humans. The mapping process followed the procedures (Li et al, 2022).

### Asian 
BRCA
 founder variants

We collected BRCA founder variants from the Asian population reported in the literature using the keywords of “BRCA1,” “BRCA2,” “founder mutation,” “variant,” “haplotype,” “Asian” or each Asian country or region. Only the variants with haplotype evidence were considered founder variants.

### RPMDS analysis

The RPMDS performance followed the procedures.
37
 Briefly, BRCA1 RING domain (PDB ID: 1JM7, 1‐103 residues),
38
 BRCA1 BRCT domain (PDB ID: 1JNX, 1649‐1859 residues)
39
 and BRCA2 BRC4 domain (PDB ID: 1N0W, 1519‐1551 residues)
40
 were used as the templates to build the mutant structure for each selected missense variant following the procedure described.
41
, 
42
, 
43
 The mutant structure was used as the starting configuration for MD simulations using GROMACS.
44
 Each simulation held approximately 99 000 atoms for the BRCA1 RING domain and BRCT domain, and 125 000 atoms for the BRCA2 BRC4 domain. Three independent simulations for the wildtype BRCA1 RING domain were performed for 100 ns.
45
 Hydrogen bonds were constrained at equilibrium length by using the LINC algorithm.
46
 The final 100 frames of the trajectories were used for variant classification. A total of 21 wildtypes in 10 ns‐frames and 33 pathogenic variants in BRCA1 RING domain, 11 benign and pathogenic variants in BRCA1 BRCT domain and 9 benign and pathogenic variants in BRCA2 BRC4 domain were used to determine the cut‐off values to classify the unknown variants into deleterious or tolerated variants. BRCA1 RING domain used wild‐type bases for benign due to the lack of benign variants. The structure deviation of BRCA1 RING domain pathogenic variants was 38.5% ± 5.8%, with a range of 36.4% and 40.3%; the structural deviation for BRCA1 BRCT pathogenic variants had a mean of 25.9% ± 8.9%, with a range of 18.6% and 33.3%; the structural deviation for BRCA2 BRC4 pathogenic variants had the mean of 23.5% ± 4.9%, with a range of 20.1% and 28.4%. The cut‐off values were set at 33.3, 30.3 and 22.9 for BRCA1 RING, BRCA1 BRCT and BRCA2 BRC4 domains, respectively. Wilcoxon test was used to test the difference between pathogenic variant and benign variant distribution in BRCA1 RING domain, BRCA1 BRCT1/2 domain and BRCA2 DNA binding domain.
41

### RESULTS

We collected BRCA variants originated from Asia population by mining through 571 publications published by May 2022 reporting Asian BRCA variation data and 18 human genomic databases containing Asian BRCA genomic sequence data from any of the 49 Asia countries (Russia's Asia part) and 4 regions of China (mainland, Hong Kong, Macau and Taiwan) (Table S1). After standardization, annotation and classification of the collected BRCA variation data, we identified a total of 7587 BRCA variants including 3245 (42.8%) in BRCA1 and 4342 (57.2%) in BRCA2, from 685 592 individuals including 324 144 (47.2%) cancer cases and 361 448 (52.8%) non‐cancer individuals, in 40 Asian countries and all 4 regions of China (Figure 1A, Tables S2A, S2B and S3A‐S3C). BRCA1 has 24 exons, with an mRNA size of 7224 bases and a protein size of 1863 amino acid residues; BRCA2 has 27 exons, with an mRNA size of 11 386 bases and protein size of 3418 amino acid residues. On average, 1 in 2.2‐bases in mRNA and half of the 1863 amino acid residues in BRCA1, and 1 in 2.6‐bases in mRNA and every of the 3418 amino acid residues in BRCA2 were altered by these variants (Figure 1B). Over three quarters (5815, 76.6%) of the 7587 variants were exonic (Figure 1C), over a third (2591, 34.1%) were nonsynonymous SNVs (Figure 1D), over a fifth (1762, 23.3%) were pathogenic and likely pathogenic variants, and nearly two‐third (4915, 64.8%) were unknown variants including the unclassified, classified as uncertain significance and conflicting interpretation (Figure 1E). Breast cancer, hereditary breast and ovarian cancer syndrome (HBOC) and ovarian cancer were the major cancer types (92.8%) used by the original studies (48.2%, 31.0% and 13.6% accordingly). However, the information is not equivalent to the origin of cancer types for individual BRCA variants identified, as many original studies did not provide one‐to‐one information between the original cancer cases and the identified BRCA variants.

Distribution of BRCA variants in Asian population. (A) BRCA variants identified in Asia countries and regions. It shows the quantity of BRCA variants from different Asia countries. Color intensities represent the number of variants identified in each Asia country and region. Only the Asia part of Russia's territory is shown. See Table S3 for details. (B) Variant distribution in BRCA1 and BRCA2 coding regions. The outer circle shows the exons of BRCA; the second circle shows the density of BRCA variants with each dot representing the frequency of variant at corresponding cDNA positions; the third circle shows the frequency of BRCA variants at coding‐changing amino‐acid residues position with a bar chart, the pathogenic variant is marked in red; the fourth circle shows the high‐prevalent (carrier number >100, present in both cancer and noncancer populations) coding‐changing pathogenic BRCA variants in the Asian population with blue bar; the inner circle shows the functional domains in BRCA, three Ashkenazi Jew founder variants (BRCA1 c.68_69del, BRCA1 c.5266dupC and BRCA2 c.5946del were marked). See Table S7C for details. (C) Location of Asian BRCA variants. (D) Variation type of Asian BRCA variants. (E) Clinical significance of Asian BRCA variants

Of the countries with variant data, Chinese and Japanese populations contributed the largest quantities of 4808 (64.6%) and 1457 (19.6%) of the total variants, respectively. Other Asia countries contributed modest BRCA variant data, and eight Asia countries had no BRCA variant data reported so far (Table S3A). We estimated the saturation level of the detection in Chinese and Japanese populations by referring to the BRCA pathogenic variant data with these from two recent studies on Chinese
34
 and Japanese populations
35
 not included in our study (but the data were included after the comparison). We observed that 61.5% of the 358 Chinese pathogenic BRCA variants and 41% of the 315 Japanese pathogenic BRCA variants in these studies were included in our dataset (Tables S4A‐S4C). The results highlighted the need to continuously collect BRCA variation data to reach higher coverage even in the Chinese and Japanese populations.

We developed an open‐access database, dbBRCA‐Asian, to share the Asian BRCA variation data with the community (https://genemutation.fhs.um.edu.mo/Asian‐BRCA/). The database contains information for the variants identified from the study, including the origin, annotation, functional classification and prevalence in the cancer cohort and non‐cancer population.

The Asian population has a highly diversified genetic background in reflecting its evolutionary history of adaptation to different natural environments in the Asian continent. We tested the regional differences within the Asian population by dividing the Asia continent into East, Southeast, South, West and Central geographic regions (see Section 2 for the covered countries/region in each region) and compared the variation data between the five regions. The results showed that BRCA variation was highly different among the regions that most of the variants were notably regional present only in a single region or shared only with a few other regions. Of the 37 variants commonly shared among all five regions, 35 were benign and only 2 were pathogenic including a founder variant BRCA1 c.5266dup in Ashkenazi Jews,
17
 which originated 1800 years ago in Scandinavia population and integrated into Ashkenazi Jews 400 to 500 years ago,
47
 and BRCA2 c.2808_2811del common in breast cancer in Western European and African (Figure 2A,B, Table S5).
48
 Data from Indian and Pakistani were typical examples that of the 796 BRCA variants identified in Indian, only 90 (45 in BRCA1, 46 in BRCA2) were shared with Pakistani's 158 BRCA variants (Figure 2C,D and Tables S6A‐S6C).

Comparison of BRCA variants within Asian populations. (A) BRCA1 variants among the populations in five Asian regions. East Asia (EA), South‐East Asia (SEA), South Asia (SA), West Asia (WA) and Central Asia (CA). See data for the BRCA1 variants common in five regions in Table S5. (B) BRCA2 variants among the populations in five Asian regions. See data for the BRCA2 variants common in five regions in Table S5. (C) BRCA1 variants between Indian and Pakistani. See actual data in Table S6. (D) BRCA2 variants between Indian and Pakistani. See actual data in Table S6

Pathogenic BRCA variant is a valuable marker for clinical diagnosis and the use of synthetic lethal‐based PARPi therapy.
7
 A total of 1762 pathogenic variants were identified from the 7587 (23.3%) Asian BRCA variants (Tables S7A and S7B). c.5496_5506delinsA was the most common Asian BRCA1 pathogenic variant (11.3%, 1981 carriers in the tested 15 880 cancer and 1581 non‐cancer individuals) and c.7480C>T was the most common Asian BRCA2 pathogenic variant (4%, 6815 carriers in the tested 104 289 cancer and 67 467 non‐cancer individuals). Table 1 lists the most common pathogenic BRCA variants in the Asian population. We tested the ethnic specificity of the 1762 Asian pathogenic BRCA variants by comparing them with the reported ones with clear ethnic origin information, including the 431 pathogenic BRCA variants collected from the Latino/Hispanic American population,
16
 and the 466 pathogenic BRCA variants identified from Middle Eastern, North African and South European countries.
49
 The results showed that only 192 of the 431 variants were shared accounting for 10.9% of the 1762 Asian pathogenic BRCA variants or 44.5% of the 431 Latin/Hispanic pathogenic BRCA variants (Figure 3A and Tables S8A‐S8C); and only 264 of the 466 variants were shared accounting for 15% of the 1762 Asian pathogenic BRCA variants and 56.7% of the 466 pathogenic BRCA variants from Middle Eastern, North African and South European countries (Figure 3B and Tables S8D‐S8F). The differences were also reflected by the significant difference of pathogenic variants in BRCA1 BRCT1 domain (P value = .002) and BRCT2 domain (P = .002) (Figure 3C and Table S8). As an additional coverage test, we also compared the 208 pathogenic BRCA variants only in the Middle Eastern population from the same study
49
 with the Asian variant data. We observed that 153 (73.6%) of the variants and 135 (64.9%) were included in the 1762 Asian variants, demonstrating high inclusion of Asian data in covering the BRCA pathogenic variants in the Middle Eastern population.

Common pathogenic BRCA variants in Asian population
a

Carrier number >100, screened in both cancer and noncancer populations.

Number of carrier/total cases × 100.

Founder variants.

Ethnic specificity of Asian pathogenic BRCA variants. (A) BRCA1 and BRCA2 pathogenic variants between Asian and Latin/Hispanic American populations. (B) BRCA1 and BRCA2 pathogenic variants between Asian and Middle Eastern, North African and South European populations. See Table S8 for the detailed data. (C) Differences of pathogenic variant density between the Asian and Latin/Hispanic American populations in BRCA1 BRCT domains and BRCA2 DNA binding domain. For example, the Latin/Hispanic American population had a lower occurrence in BRCT1 but a higher occurrence in BRCT2 than the Asian population. Red line: Asian population; Blue line: Latin/Hispanic American population. (D) Distribution of BRCA1 founder variant c.68_69del in Asian countries. It shows that the c.68_69del, a founder variant in Ashkenazi Jews and the Middle Eastern population, was widely present in Asian populations. (E) Highly variable c.67‐69 locus in Asian population. Of the eight variants present in the Asian population, two were not presented in BRCA Exchange database but only in Asian including the c.67_68delinsAG in 13 of 579 cancer patients of Singapore and c.69del in 1 of 61 cancer patients of Saudi Arabia. Green: RING domain

Many BRCA founder variants are pathogenic and present in specific ethnic populations.
11
, 
17
, 
18
, 
19
, 
20
, 
21
 We analyzed the distribution of the BRCA founder variants reported in Asian populations. Of the 49 Asian BRCA founder variants reported with haplotype evidence, 45 (90.8%) (31 in BRCA1 and 14 in BRCA2) were included in the Asian pathogenic BRCA variant list (Tables 2 and S9). As reported by the original studies, a part of the founder variants was only present in one or a few ethnic Asian population(s). For example, BRCA1 c.5154G>A in Chinese, BRCA1 c.5339T>C in Korean, BRCA2 c.6685G>T in Palestine/Jordans and BRCA2 c.5574_5577del in Japanese. However, a part of founder variants was present in more Asia populations than originally reported. For example, BRCA1 c.981_982del and BRCA2 c.7480C>T were highly prevalent in Asia (Table 1); BRCA1 c.2269del reported in Pakistani was also present in Arabic, China (mainland, Hong Kong), India, Iran, Israel, Jordan, Korea populations; BRCA1 c.5470_5477del reported in Chinese was also present in Korea and Japan; BRCA2 c.4037_4038del reported in Filipino was also present in Armenia, China (mainland, Hong Kong), India, Iran, Korea, Malaysia. BRCA1 c.68_69del, the founder variant in Ashkenazi Jews and certain Arabic populations, was widely present in Asian populations including Azerbaijan, China, India, Iran, Iraq, Japan, Malaysia, Mongolia, Nepal, Pakistan, Philippines, Qatar, Russia, Saudi Arabia, Sri Lanka, Singapore and Thailand (Figure 3D). For some highly prevalent variants, their actual distribution can be very different in different populations. For example, the locus around c.68‐69 was highly variable, surrounded by multiple variants including c.67G>A p.(Glu23Lys), c.67_68delinsAG p.(Glu23Arg), c.68dup p.(Cys24Valfs*17), c.68A>G p.(Glu23Gly), c.68del p.(Glu23Glyfs*8), c.68_69del p.(Glu23Valfs*17), c.69G>C p.(Glu23Asp), c.69del p.(Glu23Aspfs*8). Of these variants, c.67_68delinsAG p.(Glu23Arg) was only present in 13 of 579 Singapore cancer patients and c.69del p.(Glu23Aspfs*8) was only present in 1 of 61 Saudi Arabia cancer patients (Figure 3E and Table S2).

Distribution of founder variants in Asian population

Frequency = Carrier number/total cases × 100.

Our previous study indicated that human pathogenic BRCA variants mostly arose within the last 10 000 years.
50
 To identify the ancestral origins and arising time of Asian BRCA pathogenic variants, we searched their presence in the 4817 ancient human genome data dated up to 80 000 years before present (BP), in 8 Neanderthal genomes dated between 38 310 to 80 000 BP and a Denisovan genome dated 78 000 BP. We identified 153 Asian pathogenic BRCA variants (85 in BRCA1 and 68 in BRCA2) in the ancient human genomes (Figure 4A and Table S10A) dated mostly within the last 5000 BP (Figure 4B). In BRCA1, c.4485‐1G>A was the oldest identified in a 36 260‐year‐old individual in Voronezh Oblast, Russia and c.191G>A was the youngest identified in a 308‐year‐old individual in Fujian, China; In BRCA2, c.7617+1G>A was the oldest in an >10 000 BP individual in South America and c.1688G>A was the youngest in a 305‐year‐old individual in Vanuatu. As the control, we also performed the same analyses for the benign variants and identified 350 (161 in BRCA1 and 189 in BRCA2) in ancient individuals between 37 470 and 145 BP, with similar arising timing as the pathogenic variants (Figure S1 and Table S10B). While we found no evidence for the presence of human pathogenic BRCA variants in the published 9 Neanderthal and Denisovan genomes, we did identify 5 benign BRCA1 and 16 benign BRCA2 variants in 3 Neanderthals dated between 44 500 to 80 000 BP and 16 BRCA1 benign and 6 BRCA2 benign variants in the Denisovan. These contributed 4.7% of the 912 Asian benign BRCA variants (Table S10C).

Ancient fossils containing Asian pathogenic BRCA variants and timing. (A) Ancient human genomes containing Asian pathogenic BRCA variants. It showed the locations of the ancient fossil genomes and their dated times before presence (BP). Red dot line: the boundary between Europe and Asian continents; circle dot: the ancient fossil genome containing åBRCA1 pathogenic variants; triangle dot: the ancient fossil genome containing BRCA2 pathogenic variant; overlapped dot: the ancient fossil genome containing both BRCA1 and BRCA2 pathogenic variants. Color in icon: year before present. See the detailed data in Table S10. (B) Arising times of Asian pathogenic BRCA variants. It shows that the number of variants increased in the past 5000 years following the increased ancient samples available during this period. Blue: BRCA1 variants; yellow: BRCA2 variants; dot line: availability of ancient samples

Of the 7587 Asian BRCA variants identified, 4915 (64.8%) were functionally unknown as unclassified, uncertain significance or with conflicting interpretations (Figure 1E). Determination of the functional impact of the unknown variants in non‐European populations is challenging due to the lack of reference information. Previously, we developed the Ramachandran Plot‐Molecular Dynamics Simulation (RPMDS) method for the functional classification of missense variants. The method predicts the deleteriousness of missense variants based on the impact of the variant on protein structure stability.
39
, 
43
 Using this method, we analyzed the unknown missense variants in the structurally known BRCA1 RING and BRCT domains and the BRCA2 BRC4 domain. Of the 361 missense variants analyzed, 161 (44.6%) showed significant disturbance of the structural stability in the corresponding domain (Figure 5A and Table S11).

Classification of unknown BRCA missense variants by RP‐MDS method. (A) Summary of the deleterious variants classified from the unknown variants. See the actual data in Table S11. (B) c.5156T>C p(Val1719Ala) classified as a deleterious variant in the BRCA1 BRCT domain. The top left image showed the wildtype structure with VAL1719, the bottom left image showed the structure altered by Ala1719; the right images show the zoom‐in details corresponding to the structure in Val1719 and the structure in Ala1719

Figure 5B showed an example of the disturbed structure in the BRCA1 BRCT domain by a missense variant. BRCT domain (PDB ID: 1JNX) has two head‐to‐tail repeats with a similar structure. Each domain comprises of four parallel β‐sheet (residues 1650‐1655; 1675‐1677; 1685‐1689 and 1712‐1719) and three α‐helix (α‐helix‐1: residues 1659‐1669; α‐helix‐2: residues 1701‐1707; α‐helix‐3:1716‐1725). The α‐β sheet structure is maintained through the hydrophobic residues (Phe1734, Val1741, Leu1764, Met1783 and Leu1839). The wildtype Val1719 interacts with 7 residues (Leu1664, Val1665, Phe1668, Val1688, Tyr1716, Ser1722 and Ile1723), stabilizes the α‐helix‐3 by interacting with Tyr1716, Ser1722 and Ile1723 and stabilized the α‐helix‐1 to α‐helix‐3 by interacting with Leu1664 and Val1665. However, the Val1719Ala by the unclassified missense variant c.5156T>C changed BRCT structure by altering the interactions of Ala1719 to Phe1668 and Tyr1716, destabilized the structural integrity of α‐helix‐3 and the localized region. Therefore, c.5156T>C p.(Val1719Ala) was classified as a deleterious variant.

### Asian 
BRCA
 variation data

We collected BRCA variants originated from Asia population by mining through 571 publications published by May 2022 reporting Asian BRCA variation data and 18 human genomic databases containing Asian BRCA genomic sequence data from any of the 49 Asia countries (Russia's Asia part) and 4 regions of China (mainland, Hong Kong, Macau and Taiwan) (Table S1). After standardization, annotation and classification of the collected BRCA variation data, we identified a total of 7587 BRCA variants including 3245 (42.8%) in BRCA1 and 4342 (57.2%) in BRCA2, from 685 592 individuals including 324 144 (47.2%) cancer cases and 361 448 (52.8%) non‐cancer individuals, in 40 Asian countries and all 4 regions of China (Figure 1A, Tables S2A, S2B and S3A‐S3C). BRCA1 has 24 exons, with an mRNA size of 7224 bases and a protein size of 1863 amino acid residues; BRCA2 has 27 exons, with an mRNA size of 11 386 bases and protein size of 3418 amino acid residues. On average, 1 in 2.2‐bases in mRNA and half of the 1863 amino acid residues in BRCA1, and 1 in 2.6‐bases in mRNA and every of the 3418 amino acid residues in BRCA2 were altered by these variants (Figure 1B). Over three quarters (5815, 76.6%) of the 7587 variants were exonic (Figure 1C), over a third (2591, 34.1%) were nonsynonymous SNVs (Figure 1D), over a fifth (1762, 23.3%) were pathogenic and likely pathogenic variants, and nearly two‐third (4915, 64.8%) were unknown variants including the unclassified, classified as uncertain significance and conflicting interpretation (Figure 1E). Breast cancer, hereditary breast and ovarian cancer syndrome (HBOC) and ovarian cancer were the major cancer types (92.8%) used by the original studies (48.2%, 31.0% and 13.6% accordingly). However, the information is not equivalent to the origin of cancer types for individual BRCA variants identified, as many original studies did not provide one‐to‐one information between the original cancer cases and the identified BRCA variants.

Distribution of BRCA variants in Asian population. (A) BRCA variants identified in Asia countries and regions. It shows the quantity of BRCA variants from different Asia countries. Color intensities represent the number of variants identified in each Asia country and region. Only the Asia part of Russia's territory is shown. See Table S3 for details. (B) Variant distribution in BRCA1 and BRCA2 coding regions. The outer circle shows the exons of BRCA; the second circle shows the density of BRCA variants with each dot representing the frequency of variant at corresponding cDNA positions; the third circle shows the frequency of BRCA variants at coding‐changing amino‐acid residues position with a bar chart, the pathogenic variant is marked in red; the fourth circle shows the high‐prevalent (carrier number >100, present in both cancer and noncancer populations) coding‐changing pathogenic BRCA variants in the Asian population with blue bar; the inner circle shows the functional domains in BRCA, three Ashkenazi Jew founder variants (BRCA1 c.68_69del, BRCA1 c.5266dupC and BRCA2 c.5946del were marked). See Table S7C for details. (C) Location of Asian BRCA variants. (D) Variation type of Asian BRCA variants. (E) Clinical significance of Asian BRCA variants

Of the countries with variant data, Chinese and Japanese populations contributed the largest quantities of 4808 (64.6%) and 1457 (19.6%) of the total variants, respectively. Other Asia countries contributed modest BRCA variant data, and eight Asia countries had no BRCA variant data reported so far (Table S3A). We estimated the saturation level of the detection in Chinese and Japanese populations by referring to the BRCA pathogenic variant data with these from two recent studies on Chinese
34
 and Japanese populations
35
 not included in our study (but the data were included after the comparison). We observed that 61.5% of the 358 Chinese pathogenic BRCA variants and 41% of the 315 Japanese pathogenic BRCA variants in these studies were included in our dataset (Tables S4A‐S4C). The results highlighted the need to continuously collect BRCA variation data to reach higher coverage even in the Chinese and Japanese populations.

We developed an open‐access database, dbBRCA‐Asian, to share the Asian BRCA variation data with the community (https://genemutation.fhs.um.edu.mo/Asian‐BRCA/). The database contains information for the variants identified from the study, including the origin, annotation, functional classification and prevalence in the cancer cohort and non‐cancer population.

### Ethnic specificity of Asian 
BRCA
 variants

The Asian population has a highly diversified genetic background in reflecting its evolutionary history of adaptation to different natural environments in the Asian continent. We tested the regional differences within the Asian population by dividing the Asia continent into East, Southeast, South, West and Central geographic regions (see Section 2 for the covered countries/region in each region) and compared the variation data between the five regions. The results showed that BRCA variation was highly different among the regions that most of the variants were notably regional present only in a single region or shared only with a few other regions. Of the 37 variants commonly shared among all five regions, 35 were benign and only 2 were pathogenic including a founder variant BRCA1 c.5266dup in Ashkenazi Jews,
17
 which originated 1800 years ago in Scandinavia population and integrated into Ashkenazi Jews 400 to 500 years ago,
47
 and BRCA2 c.2808_2811del common in breast cancer in Western European and African (Figure 2A,B, Table S5).
48
 Data from Indian and Pakistani were typical examples that of the 796 BRCA variants identified in Indian, only 90 (45 in BRCA1, 46 in BRCA2) were shared with Pakistani's 158 BRCA variants (Figure 2C,D and Tables S6A‐S6C).

Comparison of BRCA variants within Asian populations. (A) BRCA1 variants among the populations in five Asian regions. East Asia (EA), South‐East Asia (SEA), South Asia (SA), West Asia (WA) and Central Asia (CA). See data for the BRCA1 variants common in five regions in Table S5. (B) BRCA2 variants among the populations in five Asian regions. See data for the BRCA2 variants common in five regions in Table S5. (C) BRCA1 variants between Indian and Pakistani. See actual data in Table S6. (D) BRCA2 variants between Indian and Pakistani. See actual data in Table S6

Pathogenic BRCA variant is a valuable marker for clinical diagnosis and the use of synthetic lethal‐based PARPi therapy.
7
 A total of 1762 pathogenic variants were identified from the 7587 (23.3%) Asian BRCA variants (Tables S7A and S7B). c.5496_5506delinsA was the most common Asian BRCA1 pathogenic variant (11.3%, 1981 carriers in the tested 15 880 cancer and 1581 non‐cancer individuals) and c.7480C>T was the most common Asian BRCA2 pathogenic variant (4%, 6815 carriers in the tested 104 289 cancer and 67 467 non‐cancer individuals). Table 1 lists the most common pathogenic BRCA variants in the Asian population. We tested the ethnic specificity of the 1762 Asian pathogenic BRCA variants by comparing them with the reported ones with clear ethnic origin information, including the 431 pathogenic BRCA variants collected from the Latino/Hispanic American population,
16
 and the 466 pathogenic BRCA variants identified from Middle Eastern, North African and South European countries.
49
 The results showed that only 192 of the 431 variants were shared accounting for 10.9% of the 1762 Asian pathogenic BRCA variants or 44.5% of the 431 Latin/Hispanic pathogenic BRCA variants (Figure 3A and Tables S8A‐S8C); and only 264 of the 466 variants were shared accounting for 15% of the 1762 Asian pathogenic BRCA variants and 56.7% of the 466 pathogenic BRCA variants from Middle Eastern, North African and South European countries (Figure 3B and Tables S8D‐S8F). The differences were also reflected by the significant difference of pathogenic variants in BRCA1 BRCT1 domain (P value = .002) and BRCT2 domain (P = .002) (Figure 3C and Table S8). As an additional coverage test, we also compared the 208 pathogenic BRCA variants only in the Middle Eastern population from the same study
49
 with the Asian variant data. We observed that 153 (73.6%) of the variants and 135 (64.9%) were included in the 1762 Asian variants, demonstrating high inclusion of Asian data in covering the BRCA pathogenic variants in the Middle Eastern population.

Common pathogenic BRCA variants in Asian population
a

Carrier number >100, screened in both cancer and noncancer populations.

Number of carrier/total cases × 100.

Founder variants.

Ethnic specificity of Asian pathogenic BRCA variants. (A) BRCA1 and BRCA2 pathogenic variants between Asian and Latin/Hispanic American populations. (B) BRCA1 and BRCA2 pathogenic variants between Asian and Middle Eastern, North African and South European populations. See Table S8 for the detailed data. (C) Differences of pathogenic variant density between the Asian and Latin/Hispanic American populations in BRCA1 BRCT domains and BRCA2 DNA binding domain. For example, the Latin/Hispanic American population had a lower occurrence in BRCT1 but a higher occurrence in BRCT2 than the Asian population. Red line: Asian population; Blue line: Latin/Hispanic American population. (D) Distribution of BRCA1 founder variant c.68_69del in Asian countries. It shows that the c.68_69del, a founder variant in Ashkenazi Jews and the Middle Eastern population, was widely present in Asian populations. (E) Highly variable c.67‐69 locus in Asian population. Of the eight variants present in the Asian population, two were not presented in BRCA Exchange database but only in Asian including the c.67_68delinsAG in 13 of 579 cancer patients of Singapore and c.69del in 1 of 61 cancer patients of Saudi Arabia. Green: RING domain

Many BRCA founder variants are pathogenic and present in specific ethnic populations.
11
, 
17
, 
18
, 
19
, 
20
, 
21
 We analyzed the distribution of the BRCA founder variants reported in Asian populations. Of the 49 Asian BRCA founder variants reported with haplotype evidence, 45 (90.8%) (31 in BRCA1 and 14 in BRCA2) were included in the Asian pathogenic BRCA variant list (Tables 2 and S9). As reported by the original studies, a part of the founder variants was only present in one or a few ethnic Asian population(s). For example, BRCA1 c.5154G>A in Chinese, BRCA1 c.5339T>C in Korean, BRCA2 c.6685G>T in Palestine/Jordans and BRCA2 c.5574_5577del in Japanese. However, a part of founder variants was present in more Asia populations than originally reported. For example, BRCA1 c.981_982del and BRCA2 c.7480C>T were highly prevalent in Asia (Table 1); BRCA1 c.2269del reported in Pakistani was also present in Arabic, China (mainland, Hong Kong), India, Iran, Israel, Jordan, Korea populations; BRCA1 c.5470_5477del reported in Chinese was also present in Korea and Japan; BRCA2 c.4037_4038del reported in Filipino was also present in Armenia, China (mainland, Hong Kong), India, Iran, Korea, Malaysia. BRCA1 c.68_69del, the founder variant in Ashkenazi Jews and certain Arabic populations, was widely present in Asian populations including Azerbaijan, China, India, Iran, Iraq, Japan, Malaysia, Mongolia, Nepal, Pakistan, Philippines, Qatar, Russia, Saudi Arabia, Sri Lanka, Singapore and Thailand (Figure 3D). For some highly prevalent variants, their actual distribution can be very different in different populations. For example, the locus around c.68‐69 was highly variable, surrounded by multiple variants including c.67G>A p.(Glu23Lys), c.67_68delinsAG p.(Glu23Arg), c.68dup p.(Cys24Valfs*17), c.68A>G p.(Glu23Gly), c.68del p.(Glu23Glyfs*8), c.68_69del p.(Glu23Valfs*17), c.69G>C p.(Glu23Asp), c.69del p.(Glu23Aspfs*8). Of these variants, c.67_68delinsAG p.(Glu23Arg) was only present in 13 of 579 Singapore cancer patients and c.69del p.(Glu23Aspfs*8) was only present in 1 of 61 Saudi Arabia cancer patients (Figure 3E and Table S2).

Distribution of founder variants in Asian population

Frequency = Carrier number/total cases × 100.

### Ancestral origin and arising time of Asian pathogenic 
BRCA
 variants

Our previous study indicated that human pathogenic BRCA variants mostly arose within the last 10 000 years.
50
 To identify the ancestral origins and arising time of Asian BRCA pathogenic variants, we searched their presence in the 4817 ancient human genome data dated up to 80 000 years before present (BP), in 8 Neanderthal genomes dated between 38 310 to 80 000 BP and a Denisovan genome dated 78 000 BP. We identified 153 Asian pathogenic BRCA variants (85 in BRCA1 and 68 in BRCA2) in the ancient human genomes (Figure 4A and Table S10A) dated mostly within the last 5000 BP (Figure 4B). In BRCA1, c.4485‐1G>A was the oldest identified in a 36 260‐year‐old individual in Voronezh Oblast, Russia and c.191G>A was the youngest identified in a 308‐year‐old individual in Fujian, China; In BRCA2, c.7617+1G>A was the oldest in an >10 000 BP individual in South America and c.1688G>A was the youngest in a 305‐year‐old individual in Vanuatu. As the control, we also performed the same analyses for the benign variants and identified 350 (161 in BRCA1 and 189 in BRCA2) in ancient individuals between 37 470 and 145 BP, with similar arising timing as the pathogenic variants (Figure S1 and Table S10B). While we found no evidence for the presence of human pathogenic BRCA variants in the published 9 Neanderthal and Denisovan genomes, we did identify 5 benign BRCA1 and 16 benign BRCA2 variants in 3 Neanderthals dated between 44 500 to 80 000 BP and 16 BRCA1 benign and 6 BRCA2 benign variants in the Denisovan. These contributed 4.7% of the 912 Asian benign BRCA variants (Table S10C).

Ancient fossils containing Asian pathogenic BRCA variants and timing. (A) Ancient human genomes containing Asian pathogenic BRCA variants. It showed the locations of the ancient fossil genomes and their dated times before presence (BP). Red dot line: the boundary between Europe and Asian continents; circle dot: the ancient fossil genome containing åBRCA1 pathogenic variants; triangle dot: the ancient fossil genome containing BRCA2 pathogenic variant; overlapped dot: the ancient fossil genome containing both BRCA1 and BRCA2 pathogenic variants. Color in icon: year before present. See the detailed data in Table S10. (B) Arising times of Asian pathogenic BRCA variants. It shows that the number of variants increased in the past 5000 years following the increased ancient samples available during this period. Blue: BRCA1 variants; yellow: BRCA2 variants; dot line: availability of ancient samples

### Deleteriousness of the unknown Asian 
BRCA
 variants

Of the 7587 Asian BRCA variants identified, 4915 (64.8%) were functionally unknown as unclassified, uncertain significance or with conflicting interpretations (Figure 1E). Determination of the functional impact of the unknown variants in non‐European populations is challenging due to the lack of reference information. Previously, we developed the Ramachandran Plot‐Molecular Dynamics Simulation (RPMDS) method for the functional classification of missense variants. The method predicts the deleteriousness of missense variants based on the impact of the variant on protein structure stability.
39
, 
43
 Using this method, we analyzed the unknown missense variants in the structurally known BRCA1 RING and BRCT domains and the BRCA2 BRC4 domain. Of the 361 missense variants analyzed, 161 (44.6%) showed significant disturbance of the structural stability in the corresponding domain (Figure 5A and Table S11).

Classification of unknown BRCA missense variants by RP‐MDS method. (A) Summary of the deleterious variants classified from the unknown variants. See the actual data in Table S11. (B) c.5156T>C p(Val1719Ala) classified as a deleterious variant in the BRCA1 BRCT domain. The top left image showed the wildtype structure with VAL1719, the bottom left image showed the structure altered by Ala1719; the right images show the zoom‐in details corresponding to the structure in Val1719 and the structure in Ala1719

Figure 5B showed an example of the disturbed structure in the BRCA1 BRCT domain by a missense variant. BRCT domain (PDB ID: 1JNX) has two head‐to‐tail repeats with a similar structure. Each domain comprises of four parallel β‐sheet (residues 1650‐1655; 1675‐1677; 1685‐1689 and 1712‐1719) and three α‐helix (α‐helix‐1: residues 1659‐1669; α‐helix‐2: residues 1701‐1707; α‐helix‐3:1716‐1725). The α‐β sheet structure is maintained through the hydrophobic residues (Phe1734, Val1741, Leu1764, Met1783 and Leu1839). The wildtype Val1719 interacts with 7 residues (Leu1664, Val1665, Phe1668, Val1688, Tyr1716, Ser1722 and Ile1723), stabilizes the α‐helix‐3 by interacting with Tyr1716, Ser1722 and Ile1723 and stabilized the α‐helix‐1 to α‐helix‐3 by interacting with Leu1664 and Val1665. However, the Val1719Ala by the unclassified missense variant c.5156T>C changed BRCT structure by altering the interactions of Ala1719 to Phe1668 and Tyr1716, destabilized the structural integrity of α‐helix‐3 and the localized region. Therefore, c.5156T>C p.(Val1719Ala) was classified as a deleterious variant.

### DISCUSSION

BRCA variation has become one of the most used markers in clinical cancer applications. However, the potential of BRCA variation remains to be fully explored to benefit public health,
51
 of which include the lack of BRCA variation information from non‐European populations.
22
 Our study attempts to find the largely missed BRCA variation information in the Asian population to facilitate the cancer risk study and clinical applications in the population.

The data from our study reveals the highly ethnic‐specific nature of Asian BRCA variation. Multiple factors can attribute to this feature. The positive selection uniquely imposed on human BRCA leads to the constantly arising of new BRCA variants,
8
 the adaptation of different ethnic populations to their diversified natural environment in the Asia continent promotes highly ethnic‐specific BRCA variation in Asian populations. Furthermore, bottleneck, genetic drift and founder variants can also play roles in different ethnic populations living in specific geographic environments. The highly ethnic‐specific BRCA variation in Asian population questions the use of the European descendant population‐derived BRCA variants as the sole reference for clinical applications in Asian population, as they do not include many Asian‐specific BRCA variants. A similar situation may also be present in other ethnic populations such as African and Latin American populations, in which BRCA variation is known to be highly ethnic‐specific.
11
, 
17
, 
18
, 
19
, 
20
, 
21

Identification of the evolutionary origin of BRCA variants in modern humans can reveal the evolution selection of BRCA variation. Human pathogenic BRCA variants mostly arose within the last 10 000 years after the last ice age,
34
 and nearly all Asian human BRCA founder variants identified so far arose in the past few thousand years, for example, BRCA1 c.5470‐5477del arose 2090 years ago,
34

BRCA1 c.68_69del, the founder variant in Ashkenazi Jews, arose 1200 years ago.
52
 Our archeological study identified a group of the Asian BRCA pathogenic variants in ancient individuals. However, caution must be taken in interpreting the data. First, not being observed in existing ancient genomes does not imply non‐existence because of the rarity of ancient genomes. Second, the dated age of the matched sample refers to the variant present in the individual carrier, but that variant could already arise earlier in the ancestor of the carrier. In addition, most published ancient genomes were from the samples collected in recent 10 000 years. This may also cause bias for the actual arising date of the variants. Although the modern human genome contains a few percentages of the DNA from Neanderthal and Denisovan,
53
 we found no evidence for the presence of Asian pathogenic BRCA variants in the 9 Neanderthal or Denisovan except several benign BRCA variants. However, this could be related to the limited number of Neanderthal and Denisovan cases. Analysis of more Neanderthal and Denisovan data when available should help to further address this issue.

The presence of the same variant in different ethnic populations does not imply that the variant arose from the same ancestor unless they shared the same haplotype. For example, the BRCA1 c.68_69del is a founder variant in Ashkenazi Jews and Iraq as confirmed by haplotype analysis.
52
 Our data showed that c.68_69del was also widely present in Asian populations (Figure 3C). There can be two possible sources for the sharing c.68_69del in these populations, as proposed by previous studies
52
, 
54
: it could be inherited from the common ancestor of the different populations as supported by haplotype evidence, or it could arise in different populations sporadically by co‐incidence as the locus is highly variable (Figure 3D). For instance, haplotype evidence revealed that c.68_69del in Ashkenazi Jews and Iraq populations shared the same haplotype dated 1200 years ago, whereas the haplotypes for c.68_69del in Malaysia, East India and other ethnic cases were different.
52
 As there was no haplotype information for c.68_69del in most of the Asian populations sharing this variant, we speculate that the c.68_69del in these populations likely arose by co‐incidence due to the unstable c.67‐69 locus.

It is a challenge to annotate Asian unknown BRCA variants due to the lack of reference information. Since human BRCA variation was not originated through cross‐species conservation, evolution conservation‐based in silico methods are not suitable to annotate human BRCA variants. We developed the RPMDS method to predict the deleteriousness of the unknown missense variants by referring to their impact on protein structural stability.
41
 Although structure‐based deleterious information cannot be directly translated as pathogenicity for clinical application, it does provide structure‐based functional evidence to support their pathogenicity.

Data from previous and current studies showed a high prevalence of pathogenic BRCA variants in non‐cancer human populations between 0.2% and 0.5%, or 1 pathogenic variant carrier per a few hundreds of healthy individuals, such as the 0.18% in Malaysian,
55
 0.26% in Japanese,
56
 0.38% in Han Chinese,
57
0.38% in Mexicans
58
 and 0.53% in American.
59
 And pathogenic BRCA2 variation was the highest among 169 DNA damage repair (DDR) genes in global populations.
60
, 
61
, 
62
 This raises an interesting biological question: why evolution selection allows the apparently deleterious variants to be present at such a high level in the human population? A possible explanation is that there could be beneficial effects of the pathogenic BRCA variants besides causing high cancer risk, considering multiple essential biological functions of BRCA involved including reproduction, immunity, neural development and so on.
63
, 
64
 The situation could be similar to the beneficial variation in human lactase that allows humans to consume animal milk as a protein source,
65
 and the variation in human G6PD in preventing humans from malaria infection.
66

There are certain limitations of the study. The use of country as the unit for variant collection has a limitation as certain countries like Singapore has multiple ethnic populations. Many Asia countries contributed a modest quantity of BRCA variation data. This limited in‐depth analysis of specific questions such as the relationship between BRCA variation and consanguinity tradition in the Arabic population.
67
 No BRCA data have been reported in eight Asian countries, for which we do not have clues for the features of BRCA variation in these populations. Limited ethnic/race information is provided in current BRCA databases including BRCA Exchange and ClinVar databases. This largely limits comprehensive analysis of the ethnic‐specific issue for BRCA variation between different populations. BRCA structure has not been fully determined except for a few functional domains. Therefore, the missense variants outsides these domains cannot be classified by using the structure‐based methods. While our study generated a centralized Asia BRCA variation data set, the scale of 7587 BRCA variants still lags behind the nearly 70 000 BRCA variants achieved in non‐Asian populations, mostly the European descent populations (https://brcaexchange.org/factsheet).
25
 it is also necessary to indicate that the different categories based on the combined data may not reflect the differences accurately (Figure 1C‐E), because different sequencing methods were used in original BRCA sequencing studies, including whole genome sequencing, exome sequencing, BRCA‐targeted sequencing, BRCA hotspot mutation‐targeted sequencing and so on. Data from each sequencing method can be biased to identify the variants in different regions. For example, whole genome sequencing covers both coding and non‐coding regions, exome sequencing mainly covers the coding and splicing sites, hotspot mutation‐targeted sequencing detects only a few spots of the coding region and so on.

Nevertheless, the variant data from our study provide a current view of BRCA variation and a rich resource for further study and clinical applications of BRCA‐related cancer in Asian population.

### AUTHOR CONTRIBUTIONS

Zixin Qin: Methodology; Formal Analysis; Data Curation; Visualization; Writing ‐ Review & Editing; Jiaheng Li: Methodology; Formal Analysis; Data Curation; Visualization; Writing ‐ Review & Editing; Benjamin Tam: Methodology; Software; Formal Analysis; Data Curation; Visualization; Writing ‐ Review & Editing; Siddharth Sinha: Methodology; Software; Formal Analysis; Data Curation; Visualization; Writing ‐ Review & Editing; Bojin Zhao: Software; Formal Analysis; Data Curation; Visualization; Writing ‐ Review & Editing; Shanmuga Priya Bhaskaran: Software; Formal Analysis; Data Curation; Visualization; Writing ‐ Review & Editing; Teng Huang: Software; Formal Analysis; Visualization; Writing ‐ Review & Editing; Xiaobing Wu: Formal Analysis; Writing ‐ Review & Editing; Jia Sheng Chian: Formal Analysis; Writing ‐ Review & Editing; Maoni Guo: Formal Analysis; Visualization; Writing ‐ Review & Editing; Si Hoi Kou: Formal Analysis; Writing ‐ Review & Editing; Huijun Lei: Formal Analysis; Writing ‐ Review & Editing; Li Zhang: Formal Analysis; Writing ‐ Review & Editing; Xiaoyu Wang: Formal Analysis; Writing ‐ Review & Editing; Philip Naderev P. Lagniton: Formal Analysis; Writing ‐ Review & Editing; Fengxia Xiao: Formal Analysis; Writing ‐ Review & Editing; Xinyang Jiang: Formal Analysis; Writing ‐ Review & Editing; San Ming Wang: Conceptualization; Formal Analysis; Writing ‐ Original Draft; Writing ‐ Review & Editing; Supervision and funding; with input from all authors. The work reported in the article has been performed by the authors, unless clearly specified in the text.

### FUNDING INFORMATION

Our study was supported by grants from the Macau Science and Technology Development Fund (085/2017/A2, 0077/2019/AMJ, 0032/2022/A1), the University of Macau (SRG2017‐00097‐FHS, MYRG2019‐00018‐FHS, MYRG2020‐00094‐FHS), the Faculty of Health Sciences, University of Macau (FHSIG/SW/0007/2020P and a startup fund) to SMW. BT is the recipient of the University of Macau Postdoctoral Fellowship Class A of the Macao Talent Program.

### CONFLICT OF INTEREST

The authors declare no conflict of interest in the study.

### Supporting information

Figure S1. Ancient individuals with Asian benign BRCA variants and their arising timing.

Click here for additional data file.

Table S1. Sources of BRCA variants.

Click here for additional data file.

Table S2. List of Asian BRCA variants.

Click here for additional data file.

Table S3. List of Asian BRCA variants by country and region.

Click here for additional data file.

Table S4. Saturation test.

Click here for additional data file.

Table S5. Comparison of BRCA variants among five Asia regions.

Click here for additional data file.

Table S6. Comparison of BRCA variants between India and Pakistan populations.

Click here for additional data file.

Table S7. List of Asian pathogenic BRCA variants.

Click here for additional data file.

Table S8. Comparison of pathogenic BRCA variants between Asian and Latin/Hispanic American populations.

Click here for additional data file.

Table S9. List of Asian founder variants.

Click here for additional data file.

Table S10. Ancestral origin and arising time of Asian pathogenic BRCA variants.

Click here for additional data file.

Table S11. Deleterious BRCA variants predicted by RPMDS.

Click here for additional data file.



# SUPPLEMENTAL FILE 1: IJC-152-1159.pdf

# Preparing to download ...

[HHS Vulnerability Disclosure](https://www.hhs.gov/vulnerability-disclosure-policy/index.html)