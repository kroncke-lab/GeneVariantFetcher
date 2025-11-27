# MAIN TEXT

## Whole Exome Analysis to Select Targeted Therapies for Patients with Metastatic Breast Cancer – A Feasibility Study

### Abstract

IntroductionThe purpose of this feasibility study was to select targeted therapies according to “ESMO Scale for Clinical Actionability of molecular Targets (ESCAT)”. Data interpretation was further
supported by a browser-based Treatment Decision Support platform (MH Guide, Molecular Health, Heidelberg, Germany).PatientsWe applied next generation sequencing based whole exome sequencing of tumor tissue and peripheral blood of patients with metastatic breast cancer (n = 44) to detect somatic as well as
germline mutations.Results
In 32 metastatic breast cancer patients, data interpretation was feasible. We identified 25 genomic alterations with ESCAT Level of Evidence I or II in 18/32 metastatic breast cancer
patients, which were available for evaluation: three copy number gains in
HER2
, two g
BRCA1
, two g
BRCA2
, six
PIK3CA,
one
ESR1
, three
PTEN
, one
AKT1
and two
HER2
mutations. In addition, five samples displayed Microsatellite instability high-H.
ConclusionsResulting treatment options were discussed in a tumor board and could be recommended in a small but relevant proportion of patients with metastatic breast cancer (7/18). Thus, this study
is a valuable preliminary work for the establishment of a molecular tumor board within the German initiative “Center for Personalized Medicine” which aims to shorten time for analyses and
optimize selection of targeted therapies.

### Abbreviations

BReast CAncer-Gene

Copy number gain

Copy number loss

Copy number variation

Dihydropyrimidine dehydrogenase

ESMO Scale for Clinical Actionability of molecular Targets

Estrogen receptor

European Society for Medical Oncology

Genomic alteration

Germline

Human epidermal growth factor receptor 2

Hormone receptor

Level of evidence

Metastatic breast cancer

Microsatellite instability

Phosphoinositide-3-kinase

Poly-ADP-Ribose-Polymerase

Primary tumor

Retinoblastoma 1

Triple negative breast cancer

Whole-exome sequencing

### Introduction

Metastatic breast cancer (MBC) is still an incurable disease
1
. Routine therapy options are limited or accompanied by relevant side effects. Therefore, targeted therapy to optimize treatment adherence and outcome has emerged as the
preferred approach in past years.

Alongside the decoding of the human genome “precision medicine” moved into the focus and became a reachable goal. The molecular genetic profile of a tumor, a metastatic lesion or even
germline status of a patient can give a valuable information in this context
2
. Whole exome sequencing (WES) by identifying genomic alterations (GAs) can give crucial insights into the cellular features which cause and drive cancer development and
progression, respectively, and may inform about potential therapeutic targets
3
. GAs detected by next generation sequencing (NGS) include mutations (single nucleotide variants such as missense, nonsense, splice-site mutations) as well as copy number
variations (CNV). Identifying cancer-causing and driving genes activated due to somatic CNV partly also led to specific therapeutic approaches
3
.

The first targeted therapy for
BReast CAncer
-Gene (
BRCA
)-mutated ovarian cancer was approved in 2014 being the Poly-ADP-Ribose-Polymerase (PARP)-inhibitor olaparib
4
. Meanwhile, further PARP-inhibitors were developed to treat ovarian cancer (niraparib
5
and rucaparib
6
) or approval was gained for g
BRCA
-mutated breast cancer (olaparib and talazoparib)
7
8
, as well as for
BRCA
-mutated metastatic castration-resistant prostate
9
and pancreatic
10
cancer (olaparib). Therefore, in individual cases targeted therapies approved in other tumor types may also be an option for off-label use. Other GAs provide hints for
ineffectiveness or toxicity of certain therapies. There is for example good evidence for mutations in retinoblastoma 1 (Rb1) causing resistance of CDK4/6 inhibitors in MBC patients
11
12
. Moreover since 2020 it is recommended to test patients for lack of dihydropyrimidine dehydrogenase (DPD) before starting a therapy with fluorouracil or related medicines
such as capecitabine
13
. The enzyme DPD is needed to break down fluorouracil, whose enrichment in the blood can lead to severe and life-threatening side effects
13
.

However, obtaining information of molecular genetic properties to choose a targeted therapy is still technically challenging and clinical variant interpretation is not trivial. Further the
actual challenge of personalized medicine is to combine the clinical information of a specific patient and the molecular properties of his/her tumor with the existing biomedical knowledge to
offer a personalized treatment.

To get a comprehensive picture of a patient’s GAs we present a feasibility study that explores both the somatic (tumor/metastatic tissues) and germline (white blood cells) mutational status.
Therefore, we performed WES of patients with MBC and recommended individualized targeted therapy accordingly.

### Patients, Materials and Methods

We identified patients with MBC in our weekly tumor conference. From 2017 to 2020 44 patients with MBC signed written informed consent to participate in our feasibility study “Fighting
therapy resistance in patients with solid tumors” (approved by the Ethics Committee of the Medical Faculty of the Heinrich Heine University Düsseldorf; Ref-No: 5673). Inclusion criteria
were: metastatic patient with solid tumor; no standard therapy option; Karnowsky-Index > 70%; expected life span > 6 months. Accordingly, our patient collective was heavily pretreated
with several lines of chemo- and/ or endocrine therapies.

Hormone receptor (HR) status and Human epidermal growth factor receptor 2 (HER2) amplification was determined according to German clinical routine (immunohistochemistry and/or in situ
hybridization analysis). 26/44 (59%) had a HR+ (HER2 non amplified), 7/44 (16%) a HER2 positive (HER2+; HR+ or HR-) and 10/44 (23%) had a triple negative breast cancer (TNBC). They
experienced bone and visceral (n = 16), visceral only (n = 9), cerebral (n = 9), bone only (n = 1) or other metastatic lesions. Clinical patient data are shown in
Table S1
.

Due to insufficient DNA quality and quantity or failed sequencing procedure only 32 out of the recruited 44 patients were available for WES and data interpretation (
Fig. S1
).

Peripheral blood was obtained by routine vain puncture (two 10 ml EDTA container). Suitable tumor tissue for analyses was identified according to patients’ history (preferable most recent
metastatic lesion). The tissue was obtained in clinical routine diagnostics. From archived FFPE tissues 5 µm thick sections were prepared and slides were analyzed by experienced pathologists
to determine tumor content (at least 20%) by staining with hematoxylin/ eosin and immunohistochemistry. We analyzed tumor tissue from breast, liver, lymph node, brain, skin, lung, or
bone.

DNA was purified from FFPE tissue with the GeneRead DNA FFPE Kit (Qiagen, Venlo, the Netherlands). DNA of matched blood and tumor tissue samples were sequenced using the Agilent Sure Select
XT V7 or the Illumina IDT Exome Analysis Kit. The libraries were sequenced on a HiSeq 3000 system (Illumina).

Anonymized NGS data were transferred via secure VPN channels to a server hosted by Molecular Health GmbH, Heidelberg, Germany. Data were analyzed with the CE-marked in vitro diagnostic
software (IVD) MH GUIDE (MH Guide). The MH Guide Variant detection pipeline (MH Guide VDP) uses input sequencing data in (raw) FASTQ format, MH Guide requires for the analysis of paired
exomes (WES) of somatic samples a minimum average real coverage of > 200× (the average real coverage is defined as on-target coverage after removal of duplicate read pairs), in which at
least 80% of the target region has > 100× average real coverage in the tumor sample. SNVs and Indels that passed the quality filters fulfilled the following quality parameters: PHRED
score > 28.5, coverage > 20×, allele frequency > 5%, population frequency < 1%. The raw sequencing data (reads) are aligned through the MH Guide VDP using standard (GRCh37, HG19)
or proprietary population specific human reference genomes (MH PHREGs) based on data from the 1000 genomes project for alignment, using LoFreq (PMID: 23066108) for variant calling for SNVs
and atlas and freebayes (
https://arxiv.org/abs/1207.3907
) for Indel calling. The MH Guide VDP provides all detected variants in VCF format
for transcript and protein mapping by MH Guide. In paired analyses, all variants detected in the control sample are considered to be germline variants.

For MSI detection MH Guide uses the tool MANTIS (Microsatellite Analysis for Normal-Tumor Instability) to detect MSI biomarkers from FASTQ input [MANTIS tool available at:
https://github.com/OSU-SRLab/MANTIS
]
14
. If the stepwise difference is ≥ 0.3, an MSI-H biomarker is automatically added to the variants list. This threshold of 0.3 is based on validation of 40 TCGA (Cancer
Genome Atlas) cases from three cancer entities.

MH Guide system screens all GAs identified against the reference information in the proprietary knowledge platform, Dataome. The core of this platform is a manually curated database with
evidence-based biomarker information on peer-reviewed published evidence – the so-called clinical variant interpretations (CVIs).

Information captured during the curation process of the CVIs include I) the variant – i.e. the type of genomic aberration (e.g. SNV, Insertion, Deletion etc.); II) the drug or treatment
used in the underlying published peer-reviewed evidence (preclinical studies and clinical trials) for which the data source is mainly PubMed; III) the effect of the variant on treatment –
i.e. response, resistance or safety; IV) the quantity of effect – e.g. strong, medium, weak; V) the observation context (i.e. the disease/disease stage or model system); VI) a link to the
underlying evidence and a grading of its reliability. Based on the information provided by the MH Guide platform a GA was classified as “effective” (potential target with therapeutic
options), “ineffective” (GA with evidence of less effectiveness or ineffectiveness when certain drugs are applied) and “safety” (GA raising concerns about potential toxicity when
administering some medication).

Evidence was further underpinned using the European Society for Medical Oncology (ESMO) classification “
ESMO Scale for Clinical Actionability of molecular Targets
(ESCAT)”
15
. Clinical decision on therapy options were based on these actual ESCAT evidence tier I and II
16
. In short:

ESCAT evidence tier

I: “alteration-drug match is associated with improved outcome in clinical trials”

II: “alteration-drug match is associated with antitumor activity, but magnitude of benefit is unknown”

III: “alteration-drug match suspected to improve outcome based on clinical trial data in other tumor type(s) or with similar molecular alteration”

IV: “pre-clinical evidence of actionability”

X: “lack of evidence for actionability”

### Patients

We identified patients with MBC in our weekly tumor conference. From 2017 to 2020 44 patients with MBC signed written informed consent to participate in our feasibility study “Fighting
therapy resistance in patients with solid tumors” (approved by the Ethics Committee of the Medical Faculty of the Heinrich Heine University Düsseldorf; Ref-No: 5673). Inclusion criteria
were: metastatic patient with solid tumor; no standard therapy option; Karnowsky-Index > 70%; expected life span > 6 months. Accordingly, our patient collective was heavily pretreated
with several lines of chemo- and/ or endocrine therapies.

Hormone receptor (HR) status and Human epidermal growth factor receptor 2 (HER2) amplification was determined according to German clinical routine (immunohistochemistry and/or in situ
hybridization analysis). 26/44 (59%) had a HR+ (HER2 non amplified), 7/44 (16%) a HER2 positive (HER2+; HR+ or HR-) and 10/44 (23%) had a triple negative breast cancer (TNBC). They
experienced bone and visceral (n = 16), visceral only (n = 9), cerebral (n = 9), bone only (n = 1) or other metastatic lesions. Clinical patient data are shown in
Table S1
.

Due to insufficient DNA quality and quantity or failed sequencing procedure only 32 out of the recruited 44 patients were available for WES and data interpretation (
Fig. S1
).

### Blood and tumor tissue

Peripheral blood was obtained by routine vain puncture (two 10 ml EDTA container). Suitable tumor tissue for analyses was identified according to patients’ history (preferable most recent
metastatic lesion). The tissue was obtained in clinical routine diagnostics. From archived FFPE tissues 5 µm thick sections were prepared and slides were analyzed by experienced pathologists
to determine tumor content (at least 20%) by staining with hematoxylin/ eosin and immunohistochemistry. We analyzed tumor tissue from breast, liver, lymph node, brain, skin, lung, or
bone.

### DNA extraction and whole-exome sequencing

DNA was purified from FFPE tissue with the GeneRead DNA FFPE Kit (Qiagen, Venlo, the Netherlands). DNA of matched blood and tumor tissue samples were sequenced using the Agilent Sure Select
XT V7 or the Illumina IDT Exome Analysis Kit. The libraries were sequenced on a HiSeq 3000 system (Illumina).

### Data Analysis

Anonymized NGS data were transferred via secure VPN channels to a server hosted by Molecular Health GmbH, Heidelberg, Germany. Data were analyzed with the CE-marked in vitro diagnostic
software (IVD) MH GUIDE (MH Guide). The MH Guide Variant detection pipeline (MH Guide VDP) uses input sequencing data in (raw) FASTQ format, MH Guide requires for the analysis of paired
exomes (WES) of somatic samples a minimum average real coverage of > 200× (the average real coverage is defined as on-target coverage after removal of duplicate read pairs), in which at
least 80% of the target region has > 100× average real coverage in the tumor sample. SNVs and Indels that passed the quality filters fulfilled the following quality parameters: PHRED
score > 28.5, coverage > 20×, allele frequency > 5%, population frequency < 1%. The raw sequencing data (reads) are aligned through the MH Guide VDP using standard (GRCh37, HG19)
or proprietary population specific human reference genomes (MH PHREGs) based on data from the 1000 genomes project for alignment, using LoFreq (PMID: 23066108) for variant calling for SNVs
and atlas and freebayes (
https://arxiv.org/abs/1207.3907
) for Indel calling. The MH Guide VDP provides all detected variants in VCF format
for transcript and protein mapping by MH Guide. In paired analyses, all variants detected in the control sample are considered to be germline variants.

For MSI detection MH Guide uses the tool MANTIS (Microsatellite Analysis for Normal-Tumor Instability) to detect MSI biomarkers from FASTQ input [MANTIS tool available at:
https://github.com/OSU-SRLab/MANTIS
]
14
. If the stepwise difference is ≥ 0.3, an MSI-H biomarker is automatically added to the variants list. This threshold of 0.3 is based on validation of 40 TCGA (Cancer
Genome Atlas) cases from three cancer entities.

### Data interpretation

MH Guide system screens all GAs identified against the reference information in the proprietary knowledge platform, Dataome. The core of this platform is a manually curated database with
evidence-based biomarker information on peer-reviewed published evidence – the so-called clinical variant interpretations (CVIs).

Information captured during the curation process of the CVIs include I) the variant – i.e. the type of genomic aberration (e.g. SNV, Insertion, Deletion etc.); II) the drug or treatment
used in the underlying published peer-reviewed evidence (preclinical studies and clinical trials) for which the data source is mainly PubMed; III) the effect of the variant on treatment –
i.e. response, resistance or safety; IV) the quantity of effect – e.g. strong, medium, weak; V) the observation context (i.e. the disease/disease stage or model system); VI) a link to the
underlying evidence and a grading of its reliability. Based on the information provided by the MH Guide platform a GA was classified as “effective” (potential target with therapeutic
options), “ineffective” (GA with evidence of less effectiveness or ineffectiveness when certain drugs are applied) and “safety” (GA raising concerns about potential toxicity when
administering some medication).

Evidence was further underpinned using the European Society for Medical Oncology (ESMO) classification “
ESMO Scale for Clinical Actionability of molecular Targets
(ESCAT)”
15
. Clinical decision on therapy options were based on these actual ESCAT evidence tier I and II
16
. In short:

ESCAT evidence tier

I: “alteration-drug match is associated with improved outcome in clinical trials”

II: “alteration-drug match is associated with antitumor activity, but magnitude of benefit is unknown”

III: “alteration-drug match suspected to improve outcome based on clinical trial data in other tumor type(s) or with similar molecular alteration”

IV: “pre-clinical evidence of actionability”

X: “lack of evidence for actionability”

### Results

Out of the total number of 44 screened MBC patients 32 were available for data evaluation: 21 of these harbored hormone receptor positive (HR+) disease, five presented with HER2+ MBC and
six with TNBC. In total we detected 481 GAs in 77 different genes. This comprises 253 mutations, 223 CNVs and five samples showed high microsatellite instability (MSI-H). Most alterations
were found in the
XPC
gene (n = 22), followed by the genes
FCGR3A
(n = 19) and
MTHFR
(n = 17) as well as
CYP2C19
(n = 17), which were exclusively mutations. We
detected a minimum of three and a maximum of 32 GAs per patient (median 16.5, mean 15.7).

According to ESMO recommendations 25 GAs from 18 patients belong to ESCAT LoE I or II
16
(see
Fig. 1
). Most of the patients (12/18; 67%) had one actionable GA, five had two therapeutic options and one patient had three GAs with two possible therapy recommendations.
In
Table 1
the detected mutations and the clinical consequences which were taken are listed in detail. The detected alterations of Level IA were: three copy number gains
(CNG) in
HER2,
two g
BRCA1
, two g
BRCA2
and six
PIK3CA
mutations. MSI-H as a molecular target of level 1C was identified in five samples. One
ESR1
mutation,
three
PTEN
mutations, one
AKT1
mutation and two
HER2
mutations were grouped in Level II A+B.
Fig. 1
shows the detected GAs ESCAT LoE I and II according to tumor subtype.

Twenty-five genomic alterations with potential therapeutical consequences (Level of evidence I-II according to ESMO/ESCAT) detected in 18 metastatic breast cancer
patients by next generation sequencing based whole exome sequencing. Five patients had two and one patient three different alterations (CNG = copy number gain, mut = mutation).

Only one out of five patients with HER2+ MBC (HER2 status according to clinical routine assessment) harbored a CNG in
HER2,
based on NGS CNV analysis. Contrarily, CNG in
HER2
,
based on NGS CNV analysis, were reported in two of the patients with clinically HR+/HER2− tumor (see
Table 1
). Further, in one patient diagnosed with HR+/HER2− tumors we detected two somatic
HER2
mutations (
ERBB2 p.V777L
and
ERBB2 p.V842I
).

In the HR+/HER2− subgroup one patient was carrier of a
gBRCA1
mutation (
BRCA1 p.R1645fs
) and two patients were identified with a
gBRCA2
mutation (
BRCA2 p.S2835*,
BRCA2 p.S1271*
), respectively. Interestingly the patient with
gBRCA1
mutation additionally had a somatic
BRCA2
mutation (
BRCA2 p.S3250*
, ESCAT LoE IIIA) and one
patient with a
gBRCA2
mutation (
BRCA2 p.S2835*
) also had an
PTEN
mutation (
PTEN p.T319fs
) (please refer to
Table 1
and
Table 2
). Further, one
gBRCA1
mutation (
BRCA1 p.Q1756fs
) was found in a patient with HER2+ MBC. No
gBRCA
mutations were seen in the subgroup of TNBC
patients. All detected
gBRCA
mutations were also confirmed by routine genetic analysis due to young age or family history according the criteria for genetic counselling in Germany
19
.

All spotted
PIK3CA
mutations (four times
PIK3CA p.H1047R
and twice
PIK3CA p.E545K
) were seen in the subgroup of HR+ patients (29% of HR+ patients) as listed in
Table 1
.

MSI-H tumors were detected in all subgroups (three times HR+, once HER2+ and TNBC respectively). The
ESR1 p.Y537C
and the
AKT1 p.E17K
mutation and two of the three
PTEN
mutations (
p.Y225fs
and
p.T319fs
) were found in HR+ MBC patients. One
PTEN
nonsense mutation (
p.R130*
) was detected in a TNBC patient (
Fig. 1
and
Table 1
).

We detected 20 GAs grouped in LoE III or IV according to ESCAT. In Level IIIA we detected two somatic
BRCA2
mutations (
p.S3250*; p.S368*
) in the subgroup of HR+/HER2− patients
and one somatic
BRCA2
mutation (
p.E1879fs
) in a HER2+ tumor. In LoE IV we found one
ARID1A
, one
ATM
, five
CDH1
, two
NF1
and seven
TP53
mutations, as well as one
MYC
CNG (for detail see
Table 2
).

The following known onco-targets
20
with copy number loss were detected in our MBC patients:
RB1
(n = 5)
, PTEN
(n = 4)
, CDKN2B
(n = 2)
, NF1
(n = 5)
, SMAD4
(n = 1)
,
BRCA1
(n = 9),
BRCA2
(n = 4). CNG were identified for
HER2
(n = 3, as already described)
, EGFR
(n = 1)
, MYC
(n = 1)
, PIK3CA
(n = 7)
, FGFR1
(n = 5),
FGFR2
(n = 3)
, KRAS
(n = 5)
, CCND1
(n = 4)
, MET
(n = 3)
,
and
CDK6
(n = 4).

In addition to the actionable targets, we have also identified 92 GAs in 25 different genes that indicate ineffectiveness of a particular drug in patients with MBC (mutations n = 22, MSI-H
n = 5, CNV n = 65). From a clinical point of view
11
21
22
the most relevant GAs were:
ESR1
(mutation, n = 1; as already described),
RB1
(mutation, n = 2),
CCNE1
(CNG, n = 5),
FGFR1
(CNG, n = 5)
and PIK3CA
(mutation, n = 6; as stated before). In three patients with a
PIK3CA
mutation, we saw evidence of endocrine resistance (whereby various endocrine substances were
used in multiple lines of therapy). In the other three cases, endocrine therapy was combined with either a CDK4/6 inhibitor or fulvestrant or HER2-targeted therapy when the subtype changed
during the course of the disease, and therefore probably showed a good effect of the respective therapy. Other GAs of clinical interest due to possible ineffectiveness
11
23
24
25
26
detected in our collective are:
NF1
(CNL, n = 5),
AR
(CNG, n = 4),
CDK6
(CNG, n = 4),
MET
(CNG, n = 3),
HER2
(mutation, n = 2; as
specified above),
LRP1B
(CNL, n = 2) and
FGFR2
(CNG, n = 1) (see
Fig. 2
).

Genomic alterations with potential ineffectiveness concerning certain therapies according to MH Guide detected by next generation sequencing based whole exome
sequencing in our metastatic breast cancer patients (multiple genomic alterations per patient possible; CNG = copy number gain, CNL = copy number loss, mut = mutation).

The
RB1
mutations (
RB1p.L199fs
and
RB1p.E184
) and three of the
CCNE1
CNG were seen in the HR+/HER2− subgroup which may cause failure of CDK4/6 treatment. Further
ESR1
mutation and
FGFR1
CNG can be associated with endocrine resistance and was also altered exclusively in this subgroup. However, PIK3CA mutations (as mentioned above), which
may indicate trastuzumab resistance, were also detected exclusively in the HR+/HER2− subgroup, therefore without clinical consequences in our patient collective (
Fig. 2
shows the GAs with potential ineffectiveness according to MBC subtype).

Toxicity concerns were raised for 214 GAs in 22 different genes, with mutations in
XPC
(n = 22),
FCGR3A
(n = 19),
MTHFR
(n = 17),
CYP2 C19
(n = 17) and
CYP2
D6
(n = 16) being the most altered ones (see
Fig. 3
). A
DPD
mutation was identified in four patients (two HR+ and two TNBC), which is highly important from a clinical perspective as it increases the likelihood
of relevant adverse events during therapy with fluorouracil. Further we detected three
EGFR
mutations, which may be a driver of tumorigenesis in breast cancer and/or may indicate
therapy resistance since these GAs have been found to develop under drug pressure
27
(
Fig. 3
shows the GAs with potential toxicity according to MBC subtype).

Genomic alterations with potential toxicity concerning certain therapies according to MH Guide detected by next generation sequencing based whole exome sequencing in
our metastatic breast cancer patients (multiple genomic alterations per patient possible).

### GAs detected by WES in patients with MBC

Out of the total number of 44 screened MBC patients 32 were available for data evaluation: 21 of these harbored hormone receptor positive (HR+) disease, five presented with HER2+ MBC and
six with TNBC. In total we detected 481 GAs in 77 different genes. This comprises 253 mutations, 223 CNVs and five samples showed high microsatellite instability (MSI-H). Most alterations
were found in the
XPC
gene (n = 22), followed by the genes
FCGR3A
(n = 19) and
MTHFR
(n = 17) as well as
CYP2C19
(n = 17), which were exclusively mutations. We
detected a minimum of three and a maximum of 32 GAs per patient (median 16.5, mean 15.7).

### ESCAT Level of evidence (LoE) I or II in MBC patients

According to ESMO recommendations 25 GAs from 18 patients belong to ESCAT LoE I or II
16
(see
Fig. 1
). Most of the patients (12/18; 67%) had one actionable GA, five had two therapeutic options and one patient had three GAs with two possible therapy recommendations.
In
Table 1
the detected mutations and the clinical consequences which were taken are listed in detail. The detected alterations of Level IA were: three copy number gains
(CNG) in
HER2,
two g
BRCA1
, two g
BRCA2
and six
PIK3CA
mutations. MSI-H as a molecular target of level 1C was identified in five samples. One
ESR1
mutation,
three
PTEN
mutations, one
AKT1
mutation and two
HER2
mutations were grouped in Level II A+B.
Fig. 1
shows the detected GAs ESCAT LoE I and II according to tumor subtype.

Twenty-five genomic alterations with potential therapeutical consequences (Level of evidence I-II according to ESMO/ESCAT) detected in 18 metastatic breast cancer
patients by next generation sequencing based whole exome sequencing. Five patients had two and one patient three different alterations (CNG = copy number gain, mut = mutation).

Only one out of five patients with HER2+ MBC (HER2 status according to clinical routine assessment) harbored a CNG in
HER2,
based on NGS CNV analysis. Contrarily, CNG in
HER2
,
based on NGS CNV analysis, were reported in two of the patients with clinically HR+/HER2− tumor (see
Table 1
). Further, in one patient diagnosed with HR+/HER2− tumors we detected two somatic
HER2
mutations (
ERBB2 p.V777L
and
ERBB2 p.V842I
).

In the HR+/HER2− subgroup one patient was carrier of a
gBRCA1
mutation (
BRCA1 p.R1645fs
) and two patients were identified with a
gBRCA2
mutation (
BRCA2 p.S2835*,
BRCA2 p.S1271*
), respectively. Interestingly the patient with
gBRCA1
mutation additionally had a somatic
BRCA2
mutation (
BRCA2 p.S3250*
, ESCAT LoE IIIA) and one
patient with a
gBRCA2
mutation (
BRCA2 p.S2835*
) also had an
PTEN
mutation (
PTEN p.T319fs
) (please refer to
Table 1
and
Table 2
). Further, one
gBRCA1
mutation (
BRCA1 p.Q1756fs
) was found in a patient with HER2+ MBC. No
gBRCA
mutations were seen in the subgroup of TNBC
patients. All detected
gBRCA
mutations were also confirmed by routine genetic analysis due to young age or family history according the criteria for genetic counselling in Germany
19
.

All spotted
PIK3CA
mutations (four times
PIK3CA p.H1047R
and twice
PIK3CA p.E545K
) were seen in the subgroup of HR+ patients (29% of HR+ patients) as listed in
Table 1
.

MSI-H tumors were detected in all subgroups (three times HR+, once HER2+ and TNBC respectively). The
ESR1 p.Y537C
and the
AKT1 p.E17K
mutation and two of the three
PTEN
mutations (
p.Y225fs
and
p.T319fs
) were found in HR+ MBC patients. One
PTEN
nonsense mutation (
p.R130*
) was detected in a TNBC patient (
Fig. 1
and
Table 1
).

### ESCAT LoE III or IV in MBC

We detected 20 GAs grouped in LoE III or IV according to ESCAT. In Level IIIA we detected two somatic
BRCA2
mutations (
p.S3250*; p.S368*
) in the subgroup of HR+/HER2− patients
and one somatic
BRCA2
mutation (
p.E1879fs
) in a HER2+ tumor. In LoE IV we found one
ARID1A
, one
ATM
, five
CDH1
, two
NF1
and seven
TP53
mutations, as well as one
MYC
CNG (for detail see
Table 2
).

### CNV in known onco-targets

The following known onco-targets
20
with copy number loss were detected in our MBC patients:
RB1
(n = 5)
, PTEN
(n = 4)
, CDKN2B
(n = 2)
, NF1
(n = 5)
, SMAD4
(n = 1)
,
BRCA1
(n = 9),
BRCA2
(n = 4). CNG were identified for
HER2
(n = 3, as already described)
, EGFR
(n = 1)
, MYC
(n = 1)
, PIK3CA
(n = 7)
, FGFR1
(n = 5),
FGFR2
(n = 3)
, KRAS
(n = 5)
, CCND1
(n = 4)
, MET
(n = 3)
,
and
CDK6
(n = 4).

### GAs indicating ineffectiveness or safety concerns according to MH Guide

In addition to the actionable targets, we have also identified 92 GAs in 25 different genes that indicate ineffectiveness of a particular drug in patients with MBC (mutations n = 22, MSI-H
n = 5, CNV n = 65). From a clinical point of view
11
21
22
the most relevant GAs were:
ESR1
(mutation, n = 1; as already described),
RB1
(mutation, n = 2),
CCNE1
(CNG, n = 5),
FGFR1
(CNG, n = 5)
and PIK3CA
(mutation, n = 6; as stated before). In three patients with a
PIK3CA
mutation, we saw evidence of endocrine resistance (whereby various endocrine substances were
used in multiple lines of therapy). In the other three cases, endocrine therapy was combined with either a CDK4/6 inhibitor or fulvestrant or HER2-targeted therapy when the subtype changed
during the course of the disease, and therefore probably showed a good effect of the respective therapy. Other GAs of clinical interest due to possible ineffectiveness
11
23
24
25
26
detected in our collective are:
NF1
(CNL, n = 5),
AR
(CNG, n = 4),
CDK6
(CNG, n = 4),
MET
(CNG, n = 3),
HER2
(mutation, n = 2; as
specified above),
LRP1B
(CNL, n = 2) and
FGFR2
(CNG, n = 1) (see
Fig. 2
).

Genomic alterations with potential ineffectiveness concerning certain therapies according to MH Guide detected by next generation sequencing based whole exome
sequencing in our metastatic breast cancer patients (multiple genomic alterations per patient possible; CNG = copy number gain, CNL = copy number loss, mut = mutation).

The
RB1
mutations (
RB1p.L199fs
and
RB1p.E184
) and three of the
CCNE1
CNG were seen in the HR+/HER2− subgroup which may cause failure of CDK4/6 treatment. Further
ESR1
mutation and
FGFR1
CNG can be associated with endocrine resistance and was also altered exclusively in this subgroup. However, PIK3CA mutations (as mentioned above), which
may indicate trastuzumab resistance, were also detected exclusively in the HR+/HER2− subgroup, therefore without clinical consequences in our patient collective (
Fig. 2
shows the GAs with potential ineffectiveness according to MBC subtype).

Toxicity concerns were raised for 214 GAs in 22 different genes, with mutations in
XPC
(n = 22),
FCGR3A
(n = 19),
MTHFR
(n = 17),
CYP2 C19
(n = 17) and
CYP2
D6
(n = 16) being the most altered ones (see
Fig. 3
). A
DPD
mutation was identified in four patients (two HR+ and two TNBC), which is highly important from a clinical perspective as it increases the likelihood
of relevant adverse events during therapy with fluorouracil. Further we detected three
EGFR
mutations, which may be a driver of tumorigenesis in breast cancer and/or may indicate
therapy resistance since these GAs have been found to develop under drug pressure
27
(
Fig. 3
shows the GAs with potential toxicity according to MBC subtype).

Genomic alterations with potential toxicity concerning certain therapies according to MH Guide detected by next generation sequencing based whole exome sequencing in
our metastatic breast cancer patients (multiple genomic alterations per patient possible).

### Discussion

Here we present a feasibility study performing NGS-based WES in patients with MBC. We have established a working procedure for patient recruitment, sample collection, collaboration between
different institutions and discussion of the results in a tumor board. Due to the availability of new targeted therapies whose indications are based on GAs, there is an increasing need for
genetic analyses at an early timepoint of metastatic/advanced breast cancer. Because of tumor heterogeneity and changes during carcinogenesis or the metastatic process, several analyses during
the course of disease may even be useful and necessary (e.g. regarding a somatic
PIK3CA
mutation)
28
29
. In 2017 the FDA approved two comprehensive mid-size panels for genetic testing in cancer (Integrated Mutation Profiling of Actionable Cancer Targets (MSK-IMPACT) and
FoundationOne CDx) which addressed the unmet need for precision oncology
30
. Nevertheless, it is a question of time, costs, and technical equipment to be able to perform genetic tests and, above all, to make them available for clinical use. For
instance the relatively high dropout rate (25%) of our patient population due to sample errors (too little tumor tissue or insufficient DNA quality) is in line with data published by another
research group
31
. For clinical decision-making it is further necessary to know the patient’s medical history in detail, which is essential and must be included in the discussion of
treatment decisions in a molecular tumor board.

We identified 481 GAs in patients with MBC. These include CNV, mutations and MSI-H. According to the 2019 ESMO criteria for MBC, 25 GAs were classified as ESCAT LoE I–II. In 7/18 patients
with MBC (39%), we were able to recommend targeted therapy accordingly
32
with one patient having two options and one patient having two GAs leading to the same therapy recommendation (see
Table 1
). Precisely we suggested four times pembrolizumab if MSI-H was detected, once alpelisib for a
PIK3CA
mutation, once capivasertib plus fulvestrant if there was
a
PTEN
loss-of-function mutation, once capivasertib plus fulvestrant for an
AKT
mutation and after reevaluation of the “likely pathogenic”
HER2
mutations we could
recommend neratinib in one patient diagnosed as HER2- MBC. In six patients, the possible targeted therapy was not suitable due to the individual medical history. This means that concomitant
diseases did not allow the administration of the matching drug, or the concomitant medication did not fulfil the approval. Once, the high probability of side effects and the patient’s advanced
age spoke against the optional targeted therapy. Unfortunately, the general health of six patients at the time the genetic test result was available did not allow the implementation of
potential targeted therapy. In three patients no additional therapeutic recommendations were made after WES, mainly because appropriate targeted therapies were already in use based on findings
of routine diagnostic (HER2-targeted therapies or the PARP-inhibitor olaparib). However, this also means that we recognize these findings in clinical routine and do not miss these GAs (please
refer to
Table 1
for details).

The differences in HER2 status detected by NGS (
HER2
CNG) compared to immunohistochemistry or FISH could be explained by tumor heterogeneity
33
34
or by changes in HER2 status during the course of the disease
35
, which are common phenomena in MBC. In addition, intermetastatic heterogeneity is a described phenomenon that arises under therapy pressure, as different cell clones may
respond differently or escape certain therapies. Personalized therapy aims to treat a tumor according to its genetic characteristics, which are homogeneous in the vast majority of cancer cases
at diagnosis, rather than applying a standard therapy to all patients with a specific cancer type.
36
.
HER2
mutations occur in approximately 2% of patients with breast cancer and could expand treatment options to targeted HER2 therapy in patients with HER2− PT. In
our collective, we detected two likely pathogenic
HER2
mutations (6%), which is significantly more compared to the current literature
37
. This may be explained mainly due to our small number of patients. Hempel et al. detected an
HER2
amplification by NGS in nine out of 41 patients with advanced
breast cancer (out of which seven had a HER2+ tumor by immunohistochemistry and/or in situ hybridization) and
HER2
mutation in two patients. They conclude that a threshold must be
defined on basis of CNV to allow a better interpretation of NGS based amplification analysis. Further, performing NGS analyses they could recommend promising treatment options according ESCAT
Level I in 58.5% of their patients (
HER2
mutation n = 9,
PIK3CA
mutation n = 14, MSI n = 1)
38
.

Consequently, according to our results and those of others, targeted therapy based on molecular findings is appropriate for a significant proportion of patients
39
. In this context, the selection of suitable patients and definition of targets are crucial
31
38
40
41
.

Due to the late stage of disease or irrelevance for the current treatment option none of the inefficacy or safety concerns based on GAs changed our clinical therapy decision. However, in
earlier stage of disease the detected GAs may give valuable hints concerning possible ineffectiveness or relevant toxicity of certain therapies. In this context, Reinhardt et al. screened 1270
patients with early breast cancer for the presence of
PIK3CA
mutations and found a significant influence on the effect of adjuvant aromatase inhibitors but no influence on the efficacy
of adjuvant tamoxifen
18
. Similar to ESCAT a validated classification of these GAs is highly needed for clinical routine. According to the experts, due to its high clinical relevance, the
ESR1
mutation, for example, although indicating endocrine resistance, was included in ESCAT.

Focus of our feasibility study was to identify possibly druggable molecular targets. Besides, using WES we further detected several CNV which until now have no suitable therapeutic agent or
at least there is no sufficient clinical data. Their relevance stems from the fact that the expression level of a gene strongly correlates with its copy number
17
20
. Somatic CNV typically arise during carcinogenesis and, if resulting in the deletion of tumor suppressor genes or the amplification of oncogenes, they are usually
pathogenic
20
.

Although we were able to implement the NGS based diagnostic workflow for patients with MBC in our center, several limitations must be considered. Different tumor cell clones harboring
different somatic mutations can divide a cancer into several subgroups (heterogeneity). In this context, the prognosis and response to treatment can be very unique for each clone
20
. In addition, the genomic profile of biopsy tissues provides a picture that is limited to only a single point in space and time. This may lead to an under-representation
of intratumor heterogeneity
42
, limiting the predictive value of a single tissue biopsy as in our study. Other limitations of our work are, first, the heterogeneous patient population recruited in
routine clinical practice without randomization, follow-up, or survival data. Secondly, missing or very late tissue samples of varying quality and highly variable time points in the disease
course and thirdly, the long processing time as well as the software for data evaluation, which could overlook GAs.

In line with ESMO recommendations, our data show that large gene panels and routine use of NGS result in few clinically significant responders. Nevertheless, the patient and physician may
decide together to perform a large gene panel if the patient is informed that the likelihood of benefit is rather low. According to ESMO, the use of off-label drugs matched to GAs is only
recommended when an access program and decision-making process are available
39
.

Consequently, the German initiative “Center for Personalized Medicine” established structured access to a molecular tumor board and interdisciplinary case discussion.

### Conclusion

In this feasibility study we demonstrate that WES using NGS for patients with MBC is technically possible and feasible. However, in oncology practice there are no recommendations from
scientific societies about its use in daily routine
39
. Due to the low detection rate of truly actionable GAs leading to therapeutic consequences genetic testing is recommended only for a selected patient collective. In this
context the molecular tumor board is an essential instrument to choose appropriate patients and discuss the results in an interdisciplinary setting. The identification of drugable molecular
alterations among many possible targets and especially among the various alterations of the same target, is crucial. ESCAT provides a useful continuously updated tool for classifying GAs
according to their clinical relevance which assists clinicians in delivering accurate and individualized indications for the patients
32
. The browser-based Treatment Decision Support platform MH Guide can further help interpret the abundance of genetic data.

### Supplement

Fig. S1
: Study cohort and flowchart (MBC = metastatic breast cancer, WES = whole exome sequencing).

Table S1
: Patient collective.

