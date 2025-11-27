# MAIN TEXT

## Familial trigeminal neuralgia – a systematic clinical study with a genomic screen of the neuronal electrogenisome

### Abstract

ObjectiveThis cross-sectional study examined, for the first time, a large cohort of patients with trigeminal neuralgia, to ascertain the occurrence of familial cases, providing a systematic description of clinical features of familial disease. Since there is evidence linking hyperexcitability of trigeminal ganglion neurons to trigeminal neuralgia, we also carried out an exploratory genetic analysis of the neuronal electrogenisome in these patients.MethodsWe recorded familial occurrence by systematically interviewing all patients with a definite diagnosis of classical or idiopathic trigeminal neuralgia. We found 12 occurrences of trigeminal neuralgia with positive family history out of 88 enrolled patients. Whole-exome sequencing was carried out in 11 patients. We concentrated on the genetic variants within a 173-gene panel, comprising channel genes encoding sodium, potassium, calcium, chloride, transient receptor potential channels, and gap junction channels. Gene expression profiles were based on published RNA sequencing datasets of rodent/human trigeminal ganglia tissues, with a focus on genes related to neuronal excitability.ResultsIn patients with familial trigeminal neuralgia, pain was more often located in the right, second division. All patients reported triggers. Four patients experienced concomitant continuous pain. Whole-exome sequencing analysis within the trigeminal ganglion electrogenisome identified 41 rare variants in ion channels, consisting of variants in sodium channels (6), potassium channels (10), chloride channels (5), calcium channels (7), transient receptor potential channels (12), and gap junction channels (1). In one patient, a previously profiled gain-of-function mutation in SCN10A (Nav1.8 p.Ala1304Thr), previously reported in painful neuropathy, was found; this variant was not present in unaffected siblings.ConclusionsOur results suggest that familial occurrence of trigeminal neuralgia is more common than previously considered. Although our results demonstrate variants in genes encoding voltage-gated ion channels and transient receptor potential channels within these patients, further study will be needed to determine their roles in the pathogenesis of trigeminal neuralgia.

### Background

Trigeminal neuralgia (TN) is a unique neuropathic facial pain condition, characterized by unilateral paroxysmal pain, evoked by trigger factors (1). A subgroup of patients suffers from TN with concomitant continuous pain, described as dull, burning or aching (atypical TN). Recent diagnostic criteria distinguish TN as “classical”; that is, related to neurovascular compression producing morphological changes on the trigeminal root, “secondary” to a major neurological disease, or “idiopathic” with unknown aetiology (1,2). Genetic factors may play a role in the pathophysiology of TN: this hypothesis is suggested by the younger age of TN patients without neurovascular conflict (3).

Although TN usually occurs in a sporadic fashion, rare familial occurrences have been described (4–9). The occurrence of familial TN strongly supports a genetic, presumably multifactorial and polygenic, origin of TN. In this cross-sectional study, we set out to determine, for the first time, the occurrence of familial cases in a large cohort of patients with classical and idiopathic TN, and to provide a systematic description of clinical features in these familial cases. Since hyperexcitability of trigeminal ganglion neurons may contribute to the pathophysiology of TN (10–13), we carried out whole exome sequencing, filtered to focus on the electrogenisome of trigeminal ganglion neurons, in patients with familial occurrence. Given the postulated multifactorial and polygenic origin of TN, this exploratory genetic analysis provides a stepping-stone toward future studies on larger numbers of genes and with larger numbers of patients, unrelated controls and healthy family members.

### Methods

The study was approved by the Ethic Committee of Sapienza University (reference 649/17). Written informed consent was obtained from each participant.

In this cross-sectional study, we prospectively screened consecutive patients attending the Center for Neuropathic Pain at Sapienza University, from January 2017 to January 2019. The inclusion criterion was a definite diagnosis of TN according to ICHD-3 beta criteria and the 2018 update, including 13.1.1.1 Classical TN and 13.1.1.3 Idiopathic TN. Exclusion criteria were a diagnosis of secondary TN, related to a major neurological disease, a diagnosis of orofacial pain other than TN, and cognitive disturbances or communication barriers. The diagnosis of TN was confirmed by two clinicians.

Each patient underwent a clinical examination with sensory profiling using bedside tools (14,15), trigeminal reflex testing (16,17), and a dedicated 3T MRI, with specifications optimized for identification of the TN etiology (2,18).

Familial occurrence was systematically investigated using a dedicated questionnaire, including demographic data and clinical characteristics, encompassing details of trigger factors and possible concomitant continuous pain (Supplemental Data). When a familial occurrence was reported, more information about the diagnosis of TN was acquired (i.e. the center where the diagnosis was carried out, the clinical notes, the prescribed pharmacological or surgical treatment) and verified by at least two clinicians. We included patients in the familial group only when a certified diagnosis provided by a dedicated center for the management of facial pain was available. As a further criterion, we considered the response to specific drugs; that is, carbamazepine, oxcarbazepine or phenytoin. Familial occurrence was corroborated by the pedigree analysis showing a familial occurrence in at least one first- or second-degree relative. In any cases where the level of the available information was deemed as insufficient to support the evidence of a TN diagnosis, the record was not included for the purposes of our analysis. Trigger zone overlap profiling in patients with familial and sporadic occurrence of TN was provided using dedicated software (19).

Patients 1 to 10 and 12 underwent blood sampling for the genetic analysis.

Genomic DNA was extracted from whole blood of 11 cases and two unaffected siblings (a brother and a sister) of TN4 (TN1-10 and TN12) with informed consent using NucleoSpin Tissue kit (#740952, Machery-Nagel, Duren, Germany) according to the manufacturer’s instructions. Exome library preparation and sequencing were performed by the Yale Center for Genome Analysis. Briefly, 1 µg of genomic DNA was sheared into fragments with the Covaris E220 system (Covaris, Woburn, MA, United States), followed by end-repair, A-tailing, ligation of multiplexing adaptors, and pre-capture ligation-mediated PCR amplification. The quantification and insert size distribution of product were determined using Caliper LabChip GX system (PerkinElmer Inc., Santa Clara, CA, USA). Exome library was captured using IDT xGen Exome Research Panel V1.0 (Integrated DNA Technologies, Coralville, IA, USA), and quantified by qRT-PCR using a commercially available kit (KAPA Biosystems, Wilmington MA, USA). Libraries were then loaded onto S4 flowcell and sequenced on the Illumina NovaSeq6000 platform using 101 bp paired-end sequencing reads according to Illumina protocols (Illumina Inc., San Diego, CA, USA).

The obtained reads were aligned to the human genome reference (UCSC Genome Browser, hg19) using the BWA-MEM aligner. Picard MarkDuplicates was used to mark PCR duplicates, and then the GATK indel realignment and base quality score recalibration tools were applied to generate the final BAM files. GATK HaplotypeCaller was utilized to generate GVCF files for each sample, then they were joint-called to create the initial VCF file of variant calls. The VCF files were then processed for filtering and annotation using Biomedical Genomics Workbench (Qiagen, Hilden, Germany) and in-house R scripts. A coverage depth cutoff of 10-fold and alteration frequency cutoff of 30% were applied. All variants detected in our patients were checked against the databases of gnomAD (The Genome Aggregation Database, https://gnomad.broadinstitute.org/), ExAC (The Exome Aggregation Consortium, http://exac.broadinstitute.org), 1000 genomes (https://www.ncbi.nlm.nih.gov/variation/tools/1000genomes/), dbSNP (https://www.ncbi.nlm.nih.gov/snp/), ESP (http://evs.gs.washington.edu/EVS/), UK10K (https://www.uk10k.org/data_access.html), and Yale whole-exome database. Five in-silico tools were utilized to predict the pathogenicity of these variants: SIFT, PolyPhen-2, Mutation Assessor (MA), FATHMM (Functional Analysis through Hidden Markov Models), and CONDEL (Consensus deleteriousness score) (20–24). We also annotated all genes using the Human Pain Genes Database (https://humanpaingenetics.org/hpgdb/), which contains >330 pain-associated genes. Gene expression profiles were obtained from published RNA sequencing datasets of rodent/human trigeminal ganglia tissues. These datasets consist of mouse (GSE100175) (25), rat (GSE96765) (26), and human (27). Gene expression level was estimated using the mean value of normalized FPKM (Fragments Per Kilobase of gene per Million mapped reads) separately in each dataset. Because the sequencing depth of mouse/rat/human samples varied dramatically, the absolute values should not be compared between different species.

Genomic DNA was amplified using High Fidelity, AccuPrime Taq DNA Polymerase according to manufacturer’s protocol (ThermoFisher Cat#12346-086). Thermal cycling was initiated at 94℃ for 2 minutes followed by 35 cycles of 15 s at 94℃, annealing for 30 s at 55℃, and an extension for 60 s at 68℃. PCR amplicons were sequenced at the Keck DNA Sequencing facility at Yale University. The sequencing data was analyzed using Biomedical Genomics Workbench (Qiagen, Hilden, Germany).

The Mann-Whitney test was used to compare the onset age of TN between groups; that is, familial and sporadic TN, given that this demographic data did not have a normal distribution. For comparisons of categorical data, such as the clinical characteristics, we used Fisher’s exact test or Chi-squared tests.

### Standard protocol approvals, registrations, and patient consents

The study was approved by the Ethic Committee of Sapienza University (reference 649/17). Written informed consent was obtained from each participant.

### Clinical assessment

In this cross-sectional study, we prospectively screened consecutive patients attending the Center for Neuropathic Pain at Sapienza University, from January 2017 to January 2019. The inclusion criterion was a definite diagnosis of TN according to ICHD-3 beta criteria and the 2018 update, including 13.1.1.1 Classical TN and 13.1.1.3 Idiopathic TN. Exclusion criteria were a diagnosis of secondary TN, related to a major neurological disease, a diagnosis of orofacial pain other than TN, and cognitive disturbances or communication barriers. The diagnosis of TN was confirmed by two clinicians.

Each patient underwent a clinical examination with sensory profiling using bedside tools (14,15), trigeminal reflex testing (16,17), and a dedicated 3T MRI, with specifications optimized for identification of the TN etiology (2,18).

Familial occurrence was systematically investigated using a dedicated questionnaire, including demographic data and clinical characteristics, encompassing details of trigger factors and possible concomitant continuous pain (Supplemental Data). When a familial occurrence was reported, more information about the diagnosis of TN was acquired (i.e. the center where the diagnosis was carried out, the clinical notes, the prescribed pharmacological or surgical treatment) and verified by at least two clinicians. We included patients in the familial group only when a certified diagnosis provided by a dedicated center for the management of facial pain was available. As a further criterion, we considered the response to specific drugs; that is, carbamazepine, oxcarbazepine or phenytoin. Familial occurrence was corroborated by the pedigree analysis showing a familial occurrence in at least one first- or second-degree relative. In any cases where the level of the available information was deemed as insufficient to support the evidence of a TN diagnosis, the record was not included for the purposes of our analysis. Trigger zone overlap profiling in patients with familial and sporadic occurrence of TN was provided using dedicated software (19).

Patients 1 to 10 and 12 underwent blood sampling for the genetic analysis.

### Whole-exome sequencing (WES)

Genomic DNA was extracted from whole blood of 11 cases and two unaffected siblings (a brother and a sister) of TN4 (TN1-10 and TN12) with informed consent using NucleoSpin Tissue kit (#740952, Machery-Nagel, Duren, Germany) according to the manufacturer’s instructions. Exome library preparation and sequencing were performed by the Yale Center for Genome Analysis. Briefly, 1 µg of genomic DNA was sheared into fragments with the Covaris E220 system (Covaris, Woburn, MA, United States), followed by end-repair, A-tailing, ligation of multiplexing adaptors, and pre-capture ligation-mediated PCR amplification. The quantification and insert size distribution of product were determined using Caliper LabChip GX system (PerkinElmer Inc., Santa Clara, CA, USA). Exome library was captured using IDT xGen Exome Research Panel V1.0 (Integrated DNA Technologies, Coralville, IA, USA), and quantified by qRT-PCR using a commercially available kit (KAPA Biosystems, Wilmington MA, USA). Libraries were then loaded onto S4 flowcell and sequenced on the Illumina NovaSeq6000 platform using 101 bp paired-end sequencing reads according to Illumina protocols (Illumina Inc., San Diego, CA, USA).

### WES data analysis

The obtained reads were aligned to the human genome reference (UCSC Genome Browser, hg19) using the BWA-MEM aligner. Picard MarkDuplicates was used to mark PCR duplicates, and then the GATK indel realignment and base quality score recalibration tools were applied to generate the final BAM files. GATK HaplotypeCaller was utilized to generate GVCF files for each sample, then they were joint-called to create the initial VCF file of variant calls. The VCF files were then processed for filtering and annotation using Biomedical Genomics Workbench (Qiagen, Hilden, Germany) and in-house R scripts. A coverage depth cutoff of 10-fold and alteration frequency cutoff of 30% were applied. All variants detected in our patients were checked against the databases of gnomAD (The Genome Aggregation Database, https://gnomad.broadinstitute.org/), ExAC (The Exome Aggregation Consortium, http://exac.broadinstitute.org), 1000 genomes (https://www.ncbi.nlm.nih.gov/variation/tools/1000genomes/), dbSNP (https://www.ncbi.nlm.nih.gov/snp/), ESP (http://evs.gs.washington.edu/EVS/), UK10K (https://www.uk10k.org/data_access.html), and Yale whole-exome database. Five in-silico tools were utilized to predict the pathogenicity of these variants: SIFT, PolyPhen-2, Mutation Assessor (MA), FATHMM (Functional Analysis through Hidden Markov Models), and CONDEL (Consensus deleteriousness score) (20–24). We also annotated all genes using the Human Pain Genes Database (https://humanpaingenetics.org/hpgdb/), which contains >330 pain-associated genes. Gene expression profiles were obtained from published RNA sequencing datasets of rodent/human trigeminal ganglia tissues. These datasets consist of mouse (GSE100175) (25), rat (GSE96765) (26), and human (27). Gene expression level was estimated using the mean value of normalized FPKM (Fragments Per Kilobase of gene per Million mapped reads) separately in each dataset. Because the sequencing depth of mouse/rat/human samples varied dramatically, the absolute values should not be compared between different species.

### Variant validation

Genomic DNA was amplified using High Fidelity, AccuPrime Taq DNA Polymerase according to manufacturer’s protocol (ThermoFisher Cat#12346-086). Thermal cycling was initiated at 94℃ for 2 minutes followed by 35 cycles of 15 s at 94℃, annealing for 30 s at 55℃, and an extension for 60 s at 68℃. PCR amplicons were sequenced at the Keck DNA Sequencing facility at Yale University. The sequencing data was analyzed using Biomedical Genomics Workbench (Qiagen, Hilden, Germany).

### Statistical analysis

The Mann-Whitney test was used to compare the onset age of TN between groups; that is, familial and sporadic TN, given that this demographic data did not have a normal distribution. For comparisons of categorical data, such as the clinical characteristics, we used Fisher’s exact test or Chi-squared tests.

### Findings

We have screened 112 patients with a definite diagnosis of TN. We excluded 22 patients: 17 with TN secondary to multiple sclerosis, three patients with a benign tumor of the posterior fossa, one patient with TN secondary to a superior cerebellar artery aneurysm, and one patient with megadolicobasilar artery. Thus, we included in the analysis 88 patients with classical or idiopathic TN (62 females, 26 men; age 65 ± 12.4). Fourteen cases of possible familial occurrence of TN were found, but two patients were excluded from the analysis given that a certified diagnosis provided by a dedicated center for the management of facial pain could not be obtained. In the end, 12 occurrences of possibly hereditary TN (seven female, five males, 10 families; age 64 ± 12) were included in the familial group. We defined these cases as possibly hereditary on the basis of the diagnosis of TN in at least one first- or second-degree relative, according to the pedigree (Supplemental Data, Figure e-1). Eleven index patients and two unaffected relatives underwent genetic analysis.

Clinical characteristics of case series are summarized in Table 1.
Table 1.Case series of familial TN.CaseClassical/ IdiopathicAffected relativesAgeGenderAge at onsetAffected divisionSideConcomitant continuous pain (Y/N), division, sideRemission periods (Y/N)Drug/responseTN1ClassicalTwo brothers61Female60V2LeftNNOXC 600 mg/pain reliefTN2ClassicalGrandmother45Female39V2-V3RightNNOXC 1200 mg/pain reliefTN3IdiopathicSister65Male43V1-V2-V3RightY, V2, RightYRefractoryTN4ClassicalFather73Male63V1-V2RightY, V1-V2, RightYCBZ 1000 mg/pain reliefTN5UnknownDaughter72Male64V1-V2RightNNRefractoryTN6IdiopathicFather43Female22V1-V2-V3RightY, V1-V2-V3, RightNRefractoryTN7ClassicalMother61Male46V2-V3LeftY, V3, LeftYOXC 1800 mg/pain reliefTN8ClassicalAunt76Female50V2-V3LeftY, V2-V3, LeftYOXC 900 mg/pain reliefTN9IdiopathicFather, father’s brother78Female64V2-V3RightNYDrop out for AEsTN10IdiopathicFather77Female71V1-V2-V3RightNYDrop out for AEsTN11ClassicalSister53Female43V3RightNYRefractoryTN12ClassicalFather64Male59V1-V2LeftNNOXC 1200 mg/pain relief

Case series of familial TN.

Case 1 was a 61-year-old woman whose two brothers suffered from TN. Paroxysmal pain was evoked by talking, chewing and the light touch of the nasal wing. MRI showed a neurovascular conflict producing dislocation and atrophy of the trigeminal root.

Case 2 was a 45-year-old woman whose grandmother suffered from TN. Paroxysmal pain was evoked by eating, drinking, tooth brushing and light touch of the cheek and lower eyelid. MRI showed dislocation and atrophy of the trigeminal root at the site of the neurovascular conflict.

Case 3 (brother of case 1) was a 65-year-old man suffering from atypical TN. Paroxysms were evoked by talking, chewing, tooth brushing, washing the face and tactile stimulation of trigger zones including conjunctival fornix, nasal wing, upper lip, chin, alveolar gingiva and scalp. MRI showed a contact between trigeminal root and superior cerebellar artery, without morphological changes on the trigeminal root. The patient was refractory to pharmacological treatment and underwent microvascular decompression in the posterior fossa with reduction in number of attacks per day.

Case 4 was a 73-year-old man whose father suffered from TN. TN started with purely paroxysmal pain and, after an interval of nine years without any symptoms, the patient developed an atypical form of TN. Paroxysms were evoked by chewing and light touch of supraorbital region, external side of the eye and upper lip. MRI identified a dislocation of the trigeminal root with two conflicting vessels: the superior cerebellar artery located supero-medially to the root and a vein located superiorly.

Cases 5 and 6 were father and daughter respectively. The father suffered from purely paroxysmal pain evoked by talking, chewing and light touch of lower eyelid, nasal wing and upper lip. Standard MRI excluded secondary TN related to major neurological diseases but a dedicated study of the site of neurovascular conflict was not available. The age of onset of TN in the daughter was 22, an extremely uncommon feature in classical TN. She suffered from atypical TN. Trigger maneuvers included chewing and tooth brushing; trigger zones were located in the upper and lower lip. MRI did not identify a neurovascular conflict. The patient was refractory to first line treatment and underwent stereotactic radiosurgery with complete pain relief.

Case 7 was a 61-year-old man, whose mother was affected. Paroxysmal pain was evoked by talking, chewing, washing his face and gently touching chin, cheek and upper lip. After 15 years, concomitant continuous pain developed in the third trigeminal division. MRI showed a neurovascular conflict producing dislocation and atrophy of the trigeminal root.

Case 8 was a 76-year-old woman whose aunt was affected by TN. She suffered from atypical TN. Paroxysms were evoked by chewing, washing her face, pronouncing labial letters and gently touching the scalp, nasal wing, upper and lower lip. MRI showed a bilateral neurovascular conflict with morphological changes on the affected side.

Case 9 was a 78-year-old woman whose father and uncle (father’s brother) were affected by TN. As trigger maneuvers, patient reported talking, swallowing and washing her face. Trigger zones were reported on the upper and lower lip and chin. MRI did not identify a neurovascular conflict. The patient could not assume full dosage of first line drugs due to intolerable side effects and a surgical procedure was recommended; however, follow up data were not available.

Case 10 was a 77-year-old woman whose father was affected by TN. Medical history was notable for idiopathic hemifacial spasm on the right side. Pain was triggered by talking and eating but no extraoral triggers were detected. Dedicated MRI did not identify a neurovascular conflict. Pharmacological treatment with oxcarbazepine produced side effects on the central nervous system that required a dosage reduction to an unsatisfactory level.

Case 11 was a 53-year-old woman whose sister was affected by TN. Pain was triggered by talking, eating and gently touching the face. Dedicated MRI identified a neurovascular conflict producing dislocation and atrophy on the trigeminal root. The patient was refractory to pharmacological standard treatment.

Case 12 was a 64-year-old man whose father was affected by TN. Pain was triggered by eating and tooth brushing; trigger zones included cheek and supraorbital region. Dedicated MRI identified a neurovascular conflict with a venous vessel, producing dislocation and atrophy of the trigeminal root.

Although the age at onset in patients with familial occurrences (median 55; IQR: 43÷64) was earlier in comparison with that observed in sporadic cases (median 60; IQR: 49÷64), this difference did not reach statistical significance. In one patient with familial TN, paroxysmal pain started at a very young age, 22 years old. Apart from the onset at young age in two patients with familial TN (39 and 22 years old), there was a substantial overlap between the two groups in the onset age distribution.

A comparison of demographic and clinical characteristics between sporadic and familial occurrences of TN is provided in Table 2.
Table 2.Demographics and clinical characteristics in familial and sporadic groups.PatientsGender (F/M)Affected side (L/R)Affected division (% of patients)ConcomitantTriggerRemissionTreatmentV1V2V3p (effect size)continuous pain (Y/N)factors (Y/N)periods (Y/N)response (Y/N)Familial cases n = 127/54/86 (50)11 (92)8 (67)0.0825/712/07/54/8Sporadic cases n = 7655/2120/5620 (26)56 (74)44 (58)<0.01 (0.49)231/4574/243/339/67
p
0.3310.7310.1710.2810.751>0.501>0.501>0.501>0.0511Fisher’s exact test. Reported p-values refer to the comparison between the sporadic and the familial groups.2Chi-square test. Reported p-values refer to the comparison, relevant to the most affected division, within the same group of patients.

Demographics and clinical characteristics in familial and sporadic groups.

Fisher’s exact test. Reported p-values refer to the comparison between the sporadic and the familial groups.

Chi-square test. Reported p-values refer to the comparison, relevant to the most affected division, within the same group of patients.

In patients with familial TN, as well as in patients with sporadic TN, the second trigeminal division was the most affected (11/12 and 56/76 respectively), followed by the third (8/12 and 44/76) and the first (6/12 and 20/76). In both groups, the right side was more frequently involved than the left side (8/12 and 56/76).

All patients with familial TN reported trigger factors, similarly to sporadic TN, where trigger factors were reported in 74 cases out of 76. Overlap profiling of trigger zones in patients with and without a family history of TN is provided in Figure 1.
Figure 1.Trigger zones overlap profiling in patients with sporadic (a) and familial (b) TN. The number of superimpositions ranged from 2 (dark cyan) to 15 (dark orange), in sporadic forms, and between 2 (dark cyan) and 7 (dark orange) in familial forms.

Trigger zones overlap profiling in patients with sporadic (a) and familial (b) TN. The number of superimpositions ranged from 2 (dark cyan) to 15 (dark orange), in sporadic forms, and between 2 (dark cyan) and 7 (dark orange) in familial forms.

A comparable percentage of patients in both groups (5/12 and 31/76) reported concomitant continuous pain between the paroxysmal attacks.

The percentage of remission periods in patients with familial and sporadic occurrences was 7/12 (58%) and 47/76 (57%), respectively.

In both classical and idiopathic TN, the same first-choice pharmacological treatment (CBZ or OXC) was administered (17).

Four out of 12 patients with familial occurrence and 9/76 patients with sporadic cases were refractory to standard pharmacological treatment. In 2/12 patients with familial TN and 15/76 patients with sporadic cases, side-effects produced interruption of treatment or required a dosage reduction to an unsatisfactory level.

In 4/12 and 21/76 cases with familial and sporadic TN respectively, dedicated MRI did not identify a neurovascular conflict producing morphological changes on the trigeminal root: these cases were classified as idiopathic.

Due to the purported role of hyperexcitability along the trigeminal pathway in TN (10–13), we focused our analysis on 173 ion channel genes that are known to regulate neuronal excitability in peripheral sensory neurons such as trigeminal ganglion neurons (Supplemental Data, Table e-1). We identified 41 rare variants (minor allele frequency < 1%) in 36 genes from 11 familial TN cases (Table 3). The data did not show a common variant, or variants, in the same gene in more than three unrelated probands. Future studies, testing these variants in heathy family members and unrelated controls, are needed. However, we identified potentially pathogenic variants, and applied additional filters to narrow down the number of variants for future functional assessment.
Table 3.All variants detected in our ion channel gene panel from 11 cases with familial TN.GeneCaseVariantProteinPaingnomADMouse TGRat TGHuman TGSIFTPP2MAFAMCON
SCN2A
TN774T > CIle25Thr−1/251012ND282.760.750.010.5973.445–4.570.68
SCN3A
TN846C > TLeu16Phe−1/2512580.76134.770.300.180.013–0.315–3.890.47
SCN5A
TN81019G > AArg340Gln+21/2798380.5069.800.520.260.5520.490–3.810.52
SCN7A
TN73610G > AGly1204Arg−29/17643221.4910333.3632.830.040.4842.555–4.160.61
SCN9A
TN43269T > GMet1090Arg+25/27642042.335164.2912.100.170.0080.345–1.680.48
SCN10A
TN43910G > AAla1304Thr+13/28279076.546030.398.310.000.7823.600–4.400.69
TRPA1
TN5/61178G > AArg393Gln+6/24684017.921693.4910.850.470.0011.610.980.45
TRPC6
TN2864G > TLeu288Phe−/7.5154.640.260.001.0002.775–2.070.64
TRPM2
TN5/61327A > GSer443Gly−34/2827582.25308.984.130.310.0302.505–0.160.59
TRPM3
TN104934C > TAla1645Val−280/2828242.561236.186.350.220.0000.0000.570.43
TRPM4
TN81009C > TArg337Cys+36/2803225.72261.084.150.000.6441.6500.060.49
TRPM7
TN42791G > AAla931Thr−/2.881341.877.750.030.9943.125–5.070.66
TRPM8
TN189G > AArg30Gln+9/2826888.952794.571.610.320.6891.2400.210.49
TRPS1
TN123263C > TAla1088Val−/0.1988.403.570.090.0020.000–5.070.52
TRPS1
TN8138G > TLys46Asn−12/2492240.1988.403.570.010.0110.000–5.130.53
TRPV4
TN12847T > ATyr283Asn+29/2826000.3011.431.010.000.8431.940–0.460.51
TRPV5
TN368T > CLeu23Pro−343/282814ND3.000.130.040.0071.085–1.570.53
TRPV6
TN121203C > GSer401Arg−2/2436120.06183.240.080.510.0031.320–1.950.54
KCNA5
TN4251A > CGlu84Ala−29/2357140.3292.940.310.700.0010.000–4.350.51
KCNC3
TN1/3642_650dupDuplication−27/1903941.522463.330.70/////
KCNC3
TN121706C > TPro569Leu−16/2318261.522463.330.700.610.008–0.405–4.510.48
KCND2
TN1/31467G > ASplice site+3/2504420.15ND0.04/////
KCNH2
TN11489G > AGly497Ser−/2.331287.308.800.560.0560.000–5.400.53
KCNH7
TN21120G > AVal374Met−805/2822062.40155.041.580.060.9932.430–5.390.63
KCNJ6
TN41142T > CLeu381Pro+34/2489480.4875.150.090.220.9081.100–2.490.53
KCNK17
TN2806C > TThr269Met−2/282860NDND1.660.080.0060.0001.850.39
KCNS2
TN5850A > GThr284Ala−420/2828580.266.570.370.490.0091.040–4.900.57
KCNV1
TN81418C > TAla473Val−16/2828143.0294.912.060.010.3140.695–4.490.55
CACNA1A
TN77364C > APro2455His+/1.601544.425.670.010.1350.345–3.810.51
CACNA1D
TN16032G > AArg2011Gln−4/2824480.25459.040.940.270.0021.6100.540.47
CACNA1G
TN72117G > AArg706Gln−21/2800060.0482.350.100.070.9201.355–4.280.56
CACNA1I
TN64603G > TCys1905Ser−2/2492300.50147.270.170.150.4440.700–4.890.56
CACNA1I
TN105714G > CVal1535Leu−116/2806920.50147.270.170.750.2671.610–3.990.56
CACNA1S
TN25407G > AAsp1803Asn−/0.151.870.020.000.9892.2750.270.54
CACNB1
TN101337C > TPro446Leu−6/1982444.46706.9910.750.230.0111.100–1.050.53
CLCN1
TN486A > CHis29Pro−708/282226ND116.970.040.040.4971.700–1.980.54
CLCN1
TN9501C > GPhe167Leu−308/282890ND116.970.040.560.0550.690–2.920.51
CLCN1
TN102666A > GAsn889Ser−2/248918ND116.970.040.270.0010.490–1.880.49
CLCN2
TN12704G > AArg235Gln−310/282620NDND9.250.030.7622.095–3.490.57
CLIC5
TN4799C > THis267Tyr+/0.03202.851.260.170.2421.440–3.260.55
GJB5
TN265G > AArg22His−3/251452ND2.730.060.000.9992.390–5.630.63+/−: gene listed/not listed in the Human Pain Gene Database; Mouse/rat/human TG: mean value of normalized FPKM (Fragments Per Kilobase of gene per Million mapped reads) of trigeminal ganglia; ND: not detected; Score: SIFT (cutoff < 0.05), PP2 (PolyPhen-2, cutoff > 0.5), MA (Mutation Accessor, cutoff > 0.65), FA (FATHMM, cutoff < −1.5), CON (CONDEL, cutoff > 0.5). /: no data.

All variants detected in our ion channel gene panel from 11 cases with familial TN.

+/−: gene listed/not listed in the Human Pain Gene Database; Mouse/rat/human TG: mean value of normalized FPKM (Fragments Per Kilobase of gene per Million mapped reads) of trigeminal ganglia; ND: not detected; Score: SIFT (cutoff < 0.05), PP2 (PolyPhen-2, cutoff > 0.5), MA (Mutation Accessor, cutoff > 0.65), FA (FATHMM, cutoff < −1.5), CON (CONDEL, cutoff > 0.5). /: no data.

It is notable that TN1 and TN3 are sister/brother, and they share a 9-bp duplication (c.642_650dup) in the potassium channel gene KCNC3 and a synonymous splice site variant (c.1467G > A) in KCND2. This splice site variant is located at the last nucleotide of a coding exon, and this alteration is predicted to reduce the splicing efficiency using multiple online tools (28–30). Although TN1 carries a variant in TRPM8 (Arg30Gln), a gene which is known to contribute to pain (Human Pain Genes Database), it is not shared by her sibling.

We found rare variants in four known pain-associated genes in TN4: SCN9A (encodes voltage-gated sodium channel Nav1.7), SCN10A (Nav1.8), CLIC5 (encodes Chloride Intracellular Channel 5), and KCNJ6 (encodes ATP-sensitive inward rectifier potassium channel 2) (Table 3). Among these genes, SCN9A and SCN10A have been implicated in human pain disorders (31,32). Pedigree analysis of the unaffected siblings of TN4 showed that only the p.AAla1304Thr variant in SCN10A was unique to the affected proband; this SCN10A variant is known to increase sensory neuron excitability and was previously identified as pathogenic in a patient with painful peripheral neuropathy (33).

TN5 and TN6 are a father/daughter pair and they share variants in genes encoding transient receptor potential (TRP) channels, TRPA1 (Arg393Gln) and TRPM2 (Ser443Gly). Rare variants of TRP genes were identified in all other TN cases except for TN7 and TN9 (Table 3). All of these variants are single-nucleotide substitutions: TRPC6 (TN2: p.Leu288Phe), TRPM3 (TN10: p.Ala1645Val), TRPM4 (TN8: p.Arg337Cys), TRPM7 (TN4: p.Ala931Thr), TRPM8 (TN1: p.Arg30Gln), TRPS1 (TN8: p.Lys46Asn; TN12: p.Ala1088Val), TRPV4 (TN12: pTyr283Asn), and TRPV5 (TN3: p.Leu23Pro), and TRPV6 (TN12: Ser401Arg). Only one variant was not scored as damaging by any of the five in silico analyses, and three variants are classified as pathogenic by all five algorithms, while the rest were scored as potentially pathogenic by one or more of these tests (Table 3).

Rare variants of one gene, CLCN1, were seen in three patients. Mutations of this gene have been linked to autosomal dominant/recessive myotonia congenita (34), but while expressed in trigeminal ganglia, its expression level is low. No muscle weakness or myotonia was recorded in these patients. However, given that different mutations in genes can cause different diseases, we cannot exclude a priori that these rare variants in the CLCN1 gene may contribute to TN.

### Case series of familial TN

Clinical characteristics of case series are summarized in Table 1.
Table 1.Case series of familial TN.CaseClassical/ IdiopathicAffected relativesAgeGenderAge at onsetAffected divisionSideConcomitant continuous pain (Y/N), division, sideRemission periods (Y/N)Drug/responseTN1ClassicalTwo brothers61Female60V2LeftNNOXC 600 mg/pain reliefTN2ClassicalGrandmother45Female39V2-V3RightNNOXC 1200 mg/pain reliefTN3IdiopathicSister65Male43V1-V2-V3RightY, V2, RightYRefractoryTN4ClassicalFather73Male63V1-V2RightY, V1-V2, RightYCBZ 1000 mg/pain reliefTN5UnknownDaughter72Male64V1-V2RightNNRefractoryTN6IdiopathicFather43Female22V1-V2-V3RightY, V1-V2-V3, RightNRefractoryTN7ClassicalMother61Male46V2-V3LeftY, V3, LeftYOXC 1800 mg/pain reliefTN8ClassicalAunt76Female50V2-V3LeftY, V2-V3, LeftYOXC 900 mg/pain reliefTN9IdiopathicFather, father’s brother78Female64V2-V3RightNYDrop out for AEsTN10IdiopathicFather77Female71V1-V2-V3RightNYDrop out for AEsTN11ClassicalSister53Female43V3RightNYRefractoryTN12ClassicalFather64Male59V1-V2LeftNNOXC 1200 mg/pain relief

Case series of familial TN.

Case 1 was a 61-year-old woman whose two brothers suffered from TN. Paroxysmal pain was evoked by talking, chewing and the light touch of the nasal wing. MRI showed a neurovascular conflict producing dislocation and atrophy of the trigeminal root.

Case 2 was a 45-year-old woman whose grandmother suffered from TN. Paroxysmal pain was evoked by eating, drinking, tooth brushing and light touch of the cheek and lower eyelid. MRI showed dislocation and atrophy of the trigeminal root at the site of the neurovascular conflict.

Case 3 (brother of case 1) was a 65-year-old man suffering from atypical TN. Paroxysms were evoked by talking, chewing, tooth brushing, washing the face and tactile stimulation of trigger zones including conjunctival fornix, nasal wing, upper lip, chin, alveolar gingiva and scalp. MRI showed a contact between trigeminal root and superior cerebellar artery, without morphological changes on the trigeminal root. The patient was refractory to pharmacological treatment and underwent microvascular decompression in the posterior fossa with reduction in number of attacks per day.

Case 4 was a 73-year-old man whose father suffered from TN. TN started with purely paroxysmal pain and, after an interval of nine years without any symptoms, the patient developed an atypical form of TN. Paroxysms were evoked by chewing and light touch of supraorbital region, external side of the eye and upper lip. MRI identified a dislocation of the trigeminal root with two conflicting vessels: the superior cerebellar artery located supero-medially to the root and a vein located superiorly.

Cases 5 and 6 were father and daughter respectively. The father suffered from purely paroxysmal pain evoked by talking, chewing and light touch of lower eyelid, nasal wing and upper lip. Standard MRI excluded secondary TN related to major neurological diseases but a dedicated study of the site of neurovascular conflict was not available. The age of onset of TN in the daughter was 22, an extremely uncommon feature in classical TN. She suffered from atypical TN. Trigger maneuvers included chewing and tooth brushing; trigger zones were located in the upper and lower lip. MRI did not identify a neurovascular conflict. The patient was refractory to first line treatment and underwent stereotactic radiosurgery with complete pain relief.

Case 7 was a 61-year-old man, whose mother was affected. Paroxysmal pain was evoked by talking, chewing, washing his face and gently touching chin, cheek and upper lip. After 15 years, concomitant continuous pain developed in the third trigeminal division. MRI showed a neurovascular conflict producing dislocation and atrophy of the trigeminal root.

Case 8 was a 76-year-old woman whose aunt was affected by TN. She suffered from atypical TN. Paroxysms were evoked by chewing, washing her face, pronouncing labial letters and gently touching the scalp, nasal wing, upper and lower lip. MRI showed a bilateral neurovascular conflict with morphological changes on the affected side.

Case 9 was a 78-year-old woman whose father and uncle (father’s brother) were affected by TN. As trigger maneuvers, patient reported talking, swallowing and washing her face. Trigger zones were reported on the upper and lower lip and chin. MRI did not identify a neurovascular conflict. The patient could not assume full dosage of first line drugs due to intolerable side effects and a surgical procedure was recommended; however, follow up data were not available.

Case 10 was a 77-year-old woman whose father was affected by TN. Medical history was notable for idiopathic hemifacial spasm on the right side. Pain was triggered by talking and eating but no extraoral triggers were detected. Dedicated MRI did not identify a neurovascular conflict. Pharmacological treatment with oxcarbazepine produced side effects on the central nervous system that required a dosage reduction to an unsatisfactory level.

Case 11 was a 53-year-old woman whose sister was affected by TN. Pain was triggered by talking, eating and gently touching the face. Dedicated MRI identified a neurovascular conflict producing dislocation and atrophy on the trigeminal root. The patient was refractory to pharmacological standard treatment.

Case 12 was a 64-year-old man whose father was affected by TN. Pain was triggered by eating and tooth brushing; trigger zones included cheek and supraorbital region. Dedicated MRI identified a neurovascular conflict with a venous vessel, producing dislocation and atrophy of the trigeminal root.

### Clinical characteristics, etiology and drug response in patients with familial and sporadic cases

Although the age at onset in patients with familial occurrences (median 55; IQR: 43÷64) was earlier in comparison with that observed in sporadic cases (median 60; IQR: 49÷64), this difference did not reach statistical significance. In one patient with familial TN, paroxysmal pain started at a very young age, 22 years old. Apart from the onset at young age in two patients with familial TN (39 and 22 years old), there was a substantial overlap between the two groups in the onset age distribution.

A comparison of demographic and clinical characteristics between sporadic and familial occurrences of TN is provided in Table 2.
Table 2.Demographics and clinical characteristics in familial and sporadic groups.PatientsGender (F/M)Affected side (L/R)Affected division (% of patients)ConcomitantTriggerRemissionTreatmentV1V2V3p (effect size)continuous pain (Y/N)factors (Y/N)periods (Y/N)response (Y/N)Familial cases n = 127/54/86 (50)11 (92)8 (67)0.0825/712/07/54/8Sporadic cases n = 7655/2120/5620 (26)56 (74)44 (58)<0.01 (0.49)231/4574/243/339/67
p
0.3310.7310.1710.2810.751>0.501>0.501>0.501>0.0511Fisher’s exact test. Reported p-values refer to the comparison between the sporadic and the familial groups.2Chi-square test. Reported p-values refer to the comparison, relevant to the most affected division, within the same group of patients.

Demographics and clinical characteristics in familial and sporadic groups.

Fisher’s exact test. Reported p-values refer to the comparison between the sporadic and the familial groups.

Chi-square test. Reported p-values refer to the comparison, relevant to the most affected division, within the same group of patients.

In patients with familial TN, as well as in patients with sporadic TN, the second trigeminal division was the most affected (11/12 and 56/76 respectively), followed by the third (8/12 and 44/76) and the first (6/12 and 20/76). In both groups, the right side was more frequently involved than the left side (8/12 and 56/76).

All patients with familial TN reported trigger factors, similarly to sporadic TN, where trigger factors were reported in 74 cases out of 76. Overlap profiling of trigger zones in patients with and without a family history of TN is provided in Figure 1.
Figure 1.Trigger zones overlap profiling in patients with sporadic (a) and familial (b) TN. The number of superimpositions ranged from 2 (dark cyan) to 15 (dark orange), in sporadic forms, and between 2 (dark cyan) and 7 (dark orange) in familial forms.

Trigger zones overlap profiling in patients with sporadic (a) and familial (b) TN. The number of superimpositions ranged from 2 (dark cyan) to 15 (dark orange), in sporadic forms, and between 2 (dark cyan) and 7 (dark orange) in familial forms.

A comparable percentage of patients in both groups (5/12 and 31/76) reported concomitant continuous pain between the paroxysmal attacks.

The percentage of remission periods in patients with familial and sporadic occurrences was 7/12 (58%) and 47/76 (57%), respectively.

In both classical and idiopathic TN, the same first-choice pharmacological treatment (CBZ or OXC) was administered (17).

Four out of 12 patients with familial occurrence and 9/76 patients with sporadic cases were refractory to standard pharmacological treatment. In 2/12 patients with familial TN and 15/76 patients with sporadic cases, side-effects produced interruption of treatment or required a dosage reduction to an unsatisfactory level.

In 4/12 and 21/76 cases with familial and sporadic TN respectively, dedicated MRI did not identify a neurovascular conflict producing morphological changes on the trigeminal root: these cases were classified as idiopathic.

### Genetic analysis

Due to the purported role of hyperexcitability along the trigeminal pathway in TN (10–13), we focused our analysis on 173 ion channel genes that are known to regulate neuronal excitability in peripheral sensory neurons such as trigeminal ganglion neurons (Supplemental Data, Table e-1). We identified 41 rare variants (minor allele frequency < 1%) in 36 genes from 11 familial TN cases (Table 3). The data did not show a common variant, or variants, in the same gene in more than three unrelated probands. Future studies, testing these variants in heathy family members and unrelated controls, are needed. However, we identified potentially pathogenic variants, and applied additional filters to narrow down the number of variants for future functional assessment.
Table 3.All variants detected in our ion channel gene panel from 11 cases with familial TN.GeneCaseVariantProteinPaingnomADMouse TGRat TGHuman TGSIFTPP2MAFAMCON
SCN2A
TN774T > CIle25Thr−1/251012ND282.760.750.010.5973.445–4.570.68
SCN3A
TN846C > TLeu16Phe−1/2512580.76134.770.300.180.013–0.315–3.890.47
SCN5A
TN81019G > AArg340Gln+21/2798380.5069.800.520.260.5520.490–3.810.52
SCN7A
TN73610G > AGly1204Arg−29/17643221.4910333.3632.830.040.4842.555–4.160.61
SCN9A
TN43269T > GMet1090Arg+25/27642042.335164.2912.100.170.0080.345–1.680.48
SCN10A
TN43910G > AAla1304Thr+13/28279076.546030.398.310.000.7823.600–4.400.69
TRPA1
TN5/61178G > AArg393Gln+6/24684017.921693.4910.850.470.0011.610.980.45
TRPC6
TN2864G > TLeu288Phe−/7.5154.640.260.001.0002.775–2.070.64
TRPM2
TN5/61327A > GSer443Gly−34/2827582.25308.984.130.310.0302.505–0.160.59
TRPM3
TN104934C > TAla1645Val−280/2828242.561236.186.350.220.0000.0000.570.43
TRPM4
TN81009C > TArg337Cys+36/2803225.72261.084.150.000.6441.6500.060.49
TRPM7
TN42791G > AAla931Thr−/2.881341.877.750.030.9943.125–5.070.66
TRPM8
TN189G > AArg30Gln+9/2826888.952794.571.610.320.6891.2400.210.49
TRPS1
TN123263C > TAla1088Val−/0.1988.403.570.090.0020.000–5.070.52
TRPS1
TN8138G > TLys46Asn−12/2492240.1988.403.570.010.0110.000–5.130.53
TRPV4
TN12847T > ATyr283Asn+29/2826000.3011.431.010.000.8431.940–0.460.51
TRPV5
TN368T > CLeu23Pro−343/282814ND3.000.130.040.0071.085–1.570.53
TRPV6
TN121203C > GSer401Arg−2/2436120.06183.240.080.510.0031.320–1.950.54
KCNA5
TN4251A > CGlu84Ala−29/2357140.3292.940.310.700.0010.000–4.350.51
KCNC3
TN1/3642_650dupDuplication−27/1903941.522463.330.70/////
KCNC3
TN121706C > TPro569Leu−16/2318261.522463.330.700.610.008–0.405–4.510.48
KCND2
TN1/31467G > ASplice site+3/2504420.15ND0.04/////
KCNH2
TN11489G > AGly497Ser−/2.331287.308.800.560.0560.000–5.400.53
KCNH7
TN21120G > AVal374Met−805/2822062.40155.041.580.060.9932.430–5.390.63
KCNJ6
TN41142T > CLeu381Pro+34/2489480.4875.150.090.220.9081.100–2.490.53
KCNK17
TN2806C > TThr269Met−2/282860NDND1.660.080.0060.0001.850.39
KCNS2
TN5850A > GThr284Ala−420/2828580.266.570.370.490.0091.040–4.900.57
KCNV1
TN81418C > TAla473Val−16/2828143.0294.912.060.010.3140.695–4.490.55
CACNA1A
TN77364C > APro2455His+/1.601544.425.670.010.1350.345–3.810.51
CACNA1D
TN16032G > AArg2011Gln−4/2824480.25459.040.940.270.0021.6100.540.47
CACNA1G
TN72117G > AArg706Gln−21/2800060.0482.350.100.070.9201.355–4.280.56
CACNA1I
TN64603G > TCys1905Ser−2/2492300.50147.270.170.150.4440.700–4.890.56
CACNA1I
TN105714G > CVal1535Leu−116/2806920.50147.270.170.750.2671.610–3.990.56
CACNA1S
TN25407G > AAsp1803Asn−/0.151.870.020.000.9892.2750.270.54
CACNB1
TN101337C > TPro446Leu−6/1982444.46706.9910.750.230.0111.100–1.050.53
CLCN1
TN486A > CHis29Pro−708/282226ND116.970.040.040.4971.700–1.980.54
CLCN1
TN9501C > GPhe167Leu−308/282890ND116.970.040.560.0550.690–2.920.51
CLCN1
TN102666A > GAsn889Ser−2/248918ND116.970.040.270.0010.490–1.880.49
CLCN2
TN12704G > AArg235Gln−310/282620NDND9.250.030.7622.095–3.490.57
CLIC5
TN4799C > THis267Tyr+/0.03202.851.260.170.2421.440–3.260.55
GJB5
TN265G > AArg22His−3/251452ND2.730.060.000.9992.390–5.630.63+/−: gene listed/not listed in the Human Pain Gene Database; Mouse/rat/human TG: mean value of normalized FPKM (Fragments Per Kilobase of gene per Million mapped reads) of trigeminal ganglia; ND: not detected; Score: SIFT (cutoff < 0.05), PP2 (PolyPhen-2, cutoff > 0.5), MA (Mutation Accessor, cutoff > 0.65), FA (FATHMM, cutoff < −1.5), CON (CONDEL, cutoff > 0.5). /: no data.

All variants detected in our ion channel gene panel from 11 cases with familial TN.

+/−: gene listed/not listed in the Human Pain Gene Database; Mouse/rat/human TG: mean value of normalized FPKM (Fragments Per Kilobase of gene per Million mapped reads) of trigeminal ganglia; ND: not detected; Score: SIFT (cutoff < 0.05), PP2 (PolyPhen-2, cutoff > 0.5), MA (Mutation Accessor, cutoff > 0.65), FA (FATHMM, cutoff < −1.5), CON (CONDEL, cutoff > 0.5). /: no data.

It is notable that TN1 and TN3 are sister/brother, and they share a 9-bp duplication (c.642_650dup) in the potassium channel gene KCNC3 and a synonymous splice site variant (c.1467G > A) in KCND2. This splice site variant is located at the last nucleotide of a coding exon, and this alteration is predicted to reduce the splicing efficiency using multiple online tools (28–30). Although TN1 carries a variant in TRPM8 (Arg30Gln), a gene which is known to contribute to pain (Human Pain Genes Database), it is not shared by her sibling.

We found rare variants in four known pain-associated genes in TN4: SCN9A (encodes voltage-gated sodium channel Nav1.7), SCN10A (Nav1.8), CLIC5 (encodes Chloride Intracellular Channel 5), and KCNJ6 (encodes ATP-sensitive inward rectifier potassium channel 2) (Table 3). Among these genes, SCN9A and SCN10A have been implicated in human pain disorders (31,32). Pedigree analysis of the unaffected siblings of TN4 showed that only the p.AAla1304Thr variant in SCN10A was unique to the affected proband; this SCN10A variant is known to increase sensory neuron excitability and was previously identified as pathogenic in a patient with painful peripheral neuropathy (33).

TN5 and TN6 are a father/daughter pair and they share variants in genes encoding transient receptor potential (TRP) channels, TRPA1 (Arg393Gln) and TRPM2 (Ser443Gly). Rare variants of TRP genes were identified in all other TN cases except for TN7 and TN9 (Table 3). All of these variants are single-nucleotide substitutions: TRPC6 (TN2: p.Leu288Phe), TRPM3 (TN10: p.Ala1645Val), TRPM4 (TN8: p.Arg337Cys), TRPM7 (TN4: p.Ala931Thr), TRPM8 (TN1: p.Arg30Gln), TRPS1 (TN8: p.Lys46Asn; TN12: p.Ala1088Val), TRPV4 (TN12: pTyr283Asn), and TRPV5 (TN3: p.Leu23Pro), and TRPV6 (TN12: Ser401Arg). Only one variant was not scored as damaging by any of the five in silico analyses, and three variants are classified as pathogenic by all five algorithms, while the rest were scored as potentially pathogenic by one or more of these tests (Table 3).

Rare variants of one gene, CLCN1, were seen in three patients. Mutations of this gene have been linked to autosomal dominant/recessive myotonia congenita (34), but while expressed in trigeminal ganglia, its expression level is low. No muscle weakness or myotonia was recorded in these patients. However, given that different mutations in genes can cause different diseases, we cannot exclude a priori that these rare variants in the CLCN1 gene may contribute to TN.

### WES data analysis and variant validation

Due to the purported role of hyperexcitability along the trigeminal pathway in TN (10–13), we focused our analysis on 173 ion channel genes that are known to regulate neuronal excitability in peripheral sensory neurons such as trigeminal ganglion neurons (Supplemental Data, Table e-1). We identified 41 rare variants (minor allele frequency < 1%) in 36 genes from 11 familial TN cases (Table 3). The data did not show a common variant, or variants, in the same gene in more than three unrelated probands. Future studies, testing these variants in heathy family members and unrelated controls, are needed. However, we identified potentially pathogenic variants, and applied additional filters to narrow down the number of variants for future functional assessment.
Table 3.All variants detected in our ion channel gene panel from 11 cases with familial TN.GeneCaseVariantProteinPaingnomADMouse TGRat TGHuman TGSIFTPP2MAFAMCON
SCN2A
TN774T > CIle25Thr−1/251012ND282.760.750.010.5973.445–4.570.68
SCN3A
TN846C > TLeu16Phe−1/2512580.76134.770.300.180.013–0.315–3.890.47
SCN5A
TN81019G > AArg340Gln+21/2798380.5069.800.520.260.5520.490–3.810.52
SCN7A
TN73610G > AGly1204Arg−29/17643221.4910333.3632.830.040.4842.555–4.160.61
SCN9A
TN43269T > GMet1090Arg+25/27642042.335164.2912.100.170.0080.345–1.680.48
SCN10A
TN43910G > AAla1304Thr+13/28279076.546030.398.310.000.7823.600–4.400.69
TRPA1
TN5/61178G > AArg393Gln+6/24684017.921693.4910.850.470.0011.610.980.45
TRPC6
TN2864G > TLeu288Phe−/7.5154.640.260.001.0002.775–2.070.64
TRPM2
TN5/61327A > GSer443Gly−34/2827582.25308.984.130.310.0302.505–0.160.59
TRPM3
TN104934C > TAla1645Val−280/2828242.561236.186.350.220.0000.0000.570.43
TRPM4
TN81009C > TArg337Cys+36/2803225.72261.084.150.000.6441.6500.060.49
TRPM7
TN42791G > AAla931Thr−/2.881341.877.750.030.9943.125–5.070.66
TRPM8
TN189G > AArg30Gln+9/2826888.952794.571.610.320.6891.2400.210.49
TRPS1
TN123263C > TAla1088Val−/0.1988.403.570.090.0020.000–5.070.52
TRPS1
TN8138G > TLys46Asn−12/2492240.1988.403.570.010.0110.000–5.130.53
TRPV4
TN12847T > ATyr283Asn+29/2826000.3011.431.010.000.8431.940–0.460.51
TRPV5
TN368T > CLeu23Pro−343/282814ND3.000.130.040.0071.085–1.570.53
TRPV6
TN121203C > GSer401Arg−2/2436120.06183.240.080.510.0031.320–1.950.54
KCNA5
TN4251A > CGlu84Ala−29/2357140.3292.940.310.700.0010.000–4.350.51
KCNC3
TN1/3642_650dupDuplication−27/1903941.522463.330.70/////
KCNC3
TN121706C > TPro569Leu−16/2318261.522463.330.700.610.008–0.405–4.510.48
KCND2
TN1/31467G > ASplice site+3/2504420.15ND0.04/////
KCNH2
TN11489G > AGly497Ser−/2.331287.308.800.560.0560.000–5.400.53
KCNH7
TN21120G > AVal374Met−805/2822062.40155.041.580.060.9932.430–5.390.63
KCNJ6
TN41142T > CLeu381Pro+34/2489480.4875.150.090.220.9081.100–2.490.53
KCNK17
TN2806C > TThr269Met−2/282860NDND1.660.080.0060.0001.850.39
KCNS2
TN5850A > GThr284Ala−420/2828580.266.570.370.490.0091.040–4.900.57
KCNV1
TN81418C > TAla473Val−16/2828143.0294.912.060.010.3140.695–4.490.55
CACNA1A
TN77364C > APro2455His+/1.601544.425.670.010.1350.345–3.810.51
CACNA1D
TN16032G > AArg2011Gln−4/2824480.25459.040.940.270.0021.6100.540.47
CACNA1G
TN72117G > AArg706Gln−21/2800060.0482.350.100.070.9201.355–4.280.56
CACNA1I
TN64603G > TCys1905Ser−2/2492300.50147.270.170.150.4440.700–4.890.56
CACNA1I
TN105714G > CVal1535Leu−116/2806920.50147.270.170.750.2671.610–3.990.56
CACNA1S
TN25407G > AAsp1803Asn−/0.151.870.020.000.9892.2750.270.54
CACNB1
TN101337C > TPro446Leu−6/1982444.46706.9910.750.230.0111.100–1.050.53
CLCN1
TN486A > CHis29Pro−708/282226ND116.970.040.040.4971.700–1.980.54
CLCN1
TN9501C > GPhe167Leu−308/282890ND116.970.040.560.0550.690–2.920.51
CLCN1
TN102666A > GAsn889Ser−2/248918ND116.970.040.270.0010.490–1.880.49
CLCN2
TN12704G > AArg235Gln−310/282620NDND9.250.030.7622.095–3.490.57
CLIC5
TN4799C > THis267Tyr+/0.03202.851.260.170.2421.440–3.260.55
GJB5
TN265G > AArg22His−3/251452ND2.730.060.000.9992.390–5.630.63+/−: gene listed/not listed in the Human Pain Gene Database; Mouse/rat/human TG: mean value of normalized FPKM (Fragments Per Kilobase of gene per Million mapped reads) of trigeminal ganglia; ND: not detected; Score: SIFT (cutoff < 0.05), PP2 (PolyPhen-2, cutoff > 0.5), MA (Mutation Accessor, cutoff > 0.65), FA (FATHMM, cutoff < −1.5), CON (CONDEL, cutoff > 0.5). /: no data.

All variants detected in our ion channel gene panel from 11 cases with familial TN.

+/−: gene listed/not listed in the Human Pain Gene Database; Mouse/rat/human TG: mean value of normalized FPKM (Fragments Per Kilobase of gene per Million mapped reads) of trigeminal ganglia; ND: not detected; Score: SIFT (cutoff < 0.05), PP2 (PolyPhen-2, cutoff > 0.5), MA (Mutation Accessor, cutoff > 0.65), FA (FATHMM, cutoff < −1.5), CON (CONDEL, cutoff > 0.5). /: no data.

It is notable that TN1 and TN3 are sister/brother, and they share a 9-bp duplication (c.642_650dup) in the potassium channel gene KCNC3 and a synonymous splice site variant (c.1467G > A) in KCND2. This splice site variant is located at the last nucleotide of a coding exon, and this alteration is predicted to reduce the splicing efficiency using multiple online tools (28–30). Although TN1 carries a variant in TRPM8 (Arg30Gln), a gene which is known to contribute to pain (Human Pain Genes Database), it is not shared by her sibling.

We found rare variants in four known pain-associated genes in TN4: SCN9A (encodes voltage-gated sodium channel Nav1.7), SCN10A (Nav1.8), CLIC5 (encodes Chloride Intracellular Channel 5), and KCNJ6 (encodes ATP-sensitive inward rectifier potassium channel 2) (Table 3). Among these genes, SCN9A and SCN10A have been implicated in human pain disorders (31,32). Pedigree analysis of the unaffected siblings of TN4 showed that only the p.AAla1304Thr variant in SCN10A was unique to the affected proband; this SCN10A variant is known to increase sensory neuron excitability and was previously identified as pathogenic in a patient with painful peripheral neuropathy (33).

TN5 and TN6 are a father/daughter pair and they share variants in genes encoding transient receptor potential (TRP) channels, TRPA1 (Arg393Gln) and TRPM2 (Ser443Gly). Rare variants of TRP genes were identified in all other TN cases except for TN7 and TN9 (Table 3). All of these variants are single-nucleotide substitutions: TRPC6 (TN2: p.Leu288Phe), TRPM3 (TN10: p.Ala1645Val), TRPM4 (TN8: p.Arg337Cys), TRPM7 (TN4: p.Ala931Thr), TRPM8 (TN1: p.Arg30Gln), TRPS1 (TN8: p.Lys46Asn; TN12: p.Ala1088Val), TRPV4 (TN12: pTyr283Asn), and TRPV5 (TN3: p.Leu23Pro), and TRPV6 (TN12: Ser401Arg). Only one variant was not scored as damaging by any of the five in silico analyses, and three variants are classified as pathogenic by all five algorithms, while the rest were scored as potentially pathogenic by one or more of these tests (Table 3).

Rare variants of one gene, CLCN1, were seen in three patients. Mutations of this gene have been linked to autosomal dominant/recessive myotonia congenita (34), but while expressed in trigeminal ganglia, its expression level is low. No muscle weakness or myotonia was recorded in these patients. However, given that different mutations in genes can cause different diseases, we cannot exclude a priori that these rare variants in the CLCN1 gene may contribute to TN.

### Discussion

In this cross-sectional study in a large cohort of patients with classical and idiopathic TN, we identified for the first time a significant subset of patients with familial occurrences; 11% of patients with a definite diagnosis of TN reported at least one affected relative. In view of prior evidence that hyperexcitability of trigeminal ganglion neurons contributes to the pathophysiology of TN (10–13), we focused our genetic assessment on the electrogenisome. In these patients with familial TN, our exploratory analysis showed rare variants of genes encoding voltage-gated channels and TRP channels. These findings support the need for further study of possible genetic contributions to disease pathogenesis in patients with familial TN.

Although no statistically significant difference was found in age at onset between groups in this exploratory study, two patients with familial TN showed an onset of disease at young age (before 40 years old), an extremely rare occurrence in non-secondary TN (16). The other clinical characteristics including affected division, trigger factors, remission phases and the development of concomitant continuous pain did not differ significantly between patients with familial or sporadic occurrences of TN (Table 2).

The percentage of patients refractory to first line treatments (CBZ or OXC) in the sporadic TN group was in line with previous data (35), reporting up to 9% of refractory patients. Conversely, the percentage of refractory patients in the familial TN group was higher in comparison with sporadic cases and with previous data in a sample of 200 patients with classical TN (35). This finding, approaching statistical significance (Table 2), suggests that the possible mechanisms underlying familial TN might reduce the effectiveness of voltage-gated sodium channel blockers, such as CBZ and OXC.

Genetic analysis identified rare variants in multiple genes which encode members of the TRP channel family. Two of these channels, TRPA1 and TRPM8, have been previously linked to pain (36), while the other members of this family have not been known to contribute to pain pathology. All of the TRP genes identified in our analysis are known to be expressed in peripheral sensory neurons, although some have been reported to be expressed at relatively low levels (25–27). TN5 and TN6 are a father/daughter pair and they share mutations in TRPA1 and TRPM2. The familial relationship adds weight to the possibility that these rare variants might contribute to the pathophysiology of TN in these patients. Thus, while the pathogenic role of the variants in these TRP channels is still unknown, their expression in sensory neurons suggests a possible link to pain in TN.

In one patient (case 4), we identified a previously profiled gain-of-function mutation in SCN10A (Nav1.8, p.Ala1304Thr), which was previously reported in a patient with painful peripheral neuropathy (33). Voltage-clamp analysis revealed that this mutation shifts channel activation 6 mV in a hyperpolarizing direction. Current-clamp analysis showed that expression of Nav1.8-A1304T in DRG neurons reduces current threshold, increases firing frequency in response to suprathreshold stimuli, and depolarizes resting potential, thus causing neuronal hyperexcitability of DRG neurons (33). At a clinical level, these changes might be expected to lower the threshold for evoked pain and increase pain intensity. Although this patient carried variants in other ion channels (Table 3), genomic analysis of two unaffected siblings showed that only the Ala1304Thr variant is private for the proband. Considering the multifactorial/polygenetic origin of TN, the presence or absence of one variant in a family member has limited consequence on the pathogenesis of the variant. On the contrary, the assessment of the presence of this variant in the affected father would have further supported the conclusion that this mutation directly contributes to TN in this patient. Future functional studies will be required to define the possible pathogenicity of the variants described in this paper, from genes of voltage gated channels and of the TRP channel family. As with peripheral neuropathies where mutations of ion channels appear to contribute to the pathophysiology of a disorder with a time-dependent (adult-onset) course (37), definitive explication of a role of these variants in the pathophysiology of TN will require demonstration of a “multi-hit model” that can explain both adult onset and the focal nature of the disorder.

Notably, although in the majority of patients with familial TN we found a neurovascular conflict producing significant morphological changes on the trigeminal root (38), in 4/12 familial cases MRI investigation did not show neurovascular compression.

A suggested pathophysiological mechanism of classical TN is the focal demyelination of primary afferents near the entry of the trigeminal root into the pons, a “locus minoris resistentiae” where Schwann cells are substituted by oligodendroglia in providing the myelin sheath (10,39–40). There is substantial evidence that demyelination can render axons hyperexcitable and increase their susceptibility to ectopic excitation, ephaptic transmission and generation of high frequency discharges (41,42). Neurovascular contact was found in a previous TN patient with a Nav1.6 mutation (43); this gain-of-function mutation would be expected to exacerbate any hyperexcitability due to vascular compression and/or demyelination. A gain-of-function mutation of Nav1.8 was found in one of the patients in the present series. While not normally present at most nodes in normal tissue, Nav1.8 is present at nodes and heminodes in some rodent models of demyelination (44). Nav1.6 is the predominant sodium channel at nodes of Ranvier but, if Nav1.8 channels are expressed in demyelinated axons in this patient with TN, they would be expected to contribute to ectopic impulse generation and/or cross-talk.

The findings of the present study do not establish a causal relationship, but are consistent with the hypothesis that ion channel variants might contribute to the pathogenesis in at least some cases of TN. This mechanism, if confirmed, raises the question of whether the risk:reward calculus for surgery may be different in familial cases of TN compared with non-familial cases, and further raises the question of whether optimal treatment of TN patients carrying mutations of ion channels might differ from that in patients who do not carry channel mutations; genetic factors might have an impact on response to pharmacological treatment, especially in patients without severe neurovascular conflict. Further studies will be needed to answer these questions.

Focusing on the electrogenisome, we identified rare variants in 36 genes in the 11 cases of familial TN. However, future studies, assessing a wider array of genes and testing these variants in wider samples of participants, including healthy family members and unrelated controls are needed. Given the plausible multifactorial and polygenic origin of TN, our genetic analysis provides exploratory information, paving the way to future studies which can examine a larger number of genes in a larger number of patients and controls.

The MRI investigations in the affected relatives used heterogeneous protocols and did not allow definitive assessment of neurovascular conflict. In particular, we did not have the opportunity to revise MRI images for assessment of the presence of possible trigeminal root morphological changes at the site of conflict. Therefore, we could not verify whether our patients and the affected relatives suffered from the same type of TN (classical or idiopathic). Nevertheless, given the hypothesis of a multifactorial pathogenesis of TN, the presence of a conflict with morphological abnormalities might not play a crucial role in these families. The concurrent presence of familial risk factors may increase the susceptibility of the trigeminal root or cell body to ectopic excitation and high frequency discharge generation.

### Limitations

Focusing on the electrogenisome, we identified rare variants in 36 genes in the 11 cases of familial TN. However, future studies, assessing a wider array of genes and testing these variants in wider samples of participants, including healthy family members and unrelated controls are needed. Given the plausible multifactorial and polygenic origin of TN, our genetic analysis provides exploratory information, paving the way to future studies which can examine a larger number of genes in a larger number of patients and controls.

The MRI investigations in the affected relatives used heterogeneous protocols and did not allow definitive assessment of neurovascular conflict. In particular, we did not have the opportunity to revise MRI images for assessment of the presence of possible trigeminal root morphological changes at the site of conflict. Therefore, we could not verify whether our patients and the affected relatives suffered from the same type of TN (classical or idiopathic). Nevertheless, given the hypothesis of a multifactorial pathogenesis of TN, the presence of a conflict with morphological abnormalities might not play a crucial role in these families. The concurrent presence of familial risk factors may increase the susceptibility of the trigeminal root or cell body to ectopic excitation and high frequency discharge generation.

### Conclusions

This cross-sectional study in a large cohort of patients with TN identified for the first time a significant subset of patients with familial occurrence of TN (11%). Although the difference in age at onset between groups did not reach statistical significance, in one case of familial TN, onset occurred at the notable age of 22. In one patient, a gain-of-function mutation in SCN10A (Nav1.8, p.Ala1304Thr), previously reported in a patient with an episodic pain syndrome (33), was detected. We also identified rare variants in other voltage-gated channels and TRP channels that are part of the neuronal electrogenisome. We conclude that familial TN is not uncommon and propose the hypothesis that ion channel variants might contribute to pathogenesis in at least some cases of TN.

### Article highlights

11% of patients with trigeminal neuralgia reported familial occurrences, thus suggesting that TN occurs as a familial disorder in an unexpectedly large number of cases.WES, filtered to focus on the neuronal electrogenisome, revealed multiple variants in ion channels.The results raise the question of whether optimal treatment for patients with familial TN might be different than for patients without a family history. More work is needed in this area.

11% of patients with trigeminal neuralgia reported familial occurrences, thus suggesting that TN occurs as a familial disorder in an unexpectedly large number of cases.

WES, filtered to focus on the neuronal electrogenisome, revealed multiple variants in ion channels.

The results raise the question of whether optimal treatment for patients with familial TN might be different than for patients without a family history. More work is needed in this area.

### Supplemental Material

Click here for additional data file.

Supplemental material, CEP897623 Supplemental Material1 for Familial trigeminal neuralgia – a systematic clinical study with a genomic screen of the neuronal electrogenisome by Giulia Di Stefano, Jun-Hui Yuan, Giorgio Cruccu, Stephen G Waxman, Sulayman D Dib-Hajj and Andrea Truini in Cephalalgia

Click here for additional data file.

Supplemental material, CEP897623 Supplemental Material2 for Familial trigeminal neuralgia – a systematic clinical study with a genomic screen of the neuronal electrogenisome by Giulia Di Stefano, Jun-Hui Yuan, Giorgio Cruccu, Stephen G Waxman, Sulayman D Dib-Hajj and Andrea Truini in Cephalalgia

Click here for additional data file.

Supplemental material, CEP897623 Supplemental Material3 for Familial trigeminal neuralgia – a systematic clinical study with a genomic screen of the neuronal electrogenisome by Giulia Di Stefano, Jun-Hui Yuan, Giorgio Cruccu, Stephen G Waxman, Sulayman D Dib-Hajj and Andrea Truini in Cephalalgia



# SUPPLEMENTAL FILE 1: 10.1177_0333102419897623.pdf

# Preparing to download ...

[HHS Vulnerability Disclosure](https://www.hhs.gov/vulnerability-disclosure-policy/index.html)