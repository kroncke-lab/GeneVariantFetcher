# MAIN TEXT

## Computational analysis of long QT syndrome type 2 and the therapeutic effects of KCNQ1 antibodies

### Abstract

ObjectiveLong QT interval syndrome (LQTS) is a highly dangerous cardiac disease that can lead to sudden cardiac death; however, its underlying mechanism remains largely unknown. This study is conceived to investigate the impact of two general genotypes of LQTS type 2, and also the therapeutic effects of an emerging immunology-based treatment named KCNQ1 antibody.MethodsA multiscale virtual heart is developed, which contains multiple biological levels ranging from ion channels to a three-dimensional cardiac structure with realistic geometry. Critical biomarkers at different biological levels are monitored to investigate the remodeling of cardiac electrophysiology induced by mutations.ResultsSimulations revealed multiple important mechanisms that are hard to capture via conventional clinical techniques, including the augmented dispersion of repolarization, the increased vulnerability to arrhythmias, the impaired adaptability in tissue to high heart rates, and so on. An emerging KCNQ1 antibody-based therapy could rescue the prolonged QT interval but did not reduce the vulnerable window.ConclusionsTiny molecular alterations can lead to cardiac electrophysiological remodeling at multiple biological levels, which in turn contributes to higher susceptibility to lethal arrhythmias in long QT syndrome type 2 patients. The KCNQ1 antibody-based therapy has proarrhythmic risks notwithstanding its QT-rescuing effects.

### Introduction

Long QT syndrome (LQTS) may lead to ventricular arrhythmias, particularly torsade de pointes ventricular tachycardia, syncope, or sudden death, which pose a significant threat to health. In China, the epidemiological and genetic data on patients are much less well-defined.
1
 Among several prevailing subtypes, the long QT syndrome type 2 (LQT2) is frequently reported in existing epidemiological reports,2,3 and it is inadequately responsive to medications, costly to manage with devices, and entails various surgical complications. LQT2 is primarily associated with mutations in the KCNH2 gene, which encodes a potassium channel known as the human ether-à-go-go-related gene (hERG) channel.
4
 The hERG channel plays a crucial role in regulating the repolarizing current IKr, and loss-of-function mutations to KCNH2 will result in a decrease in IKr in cardiomyocytes, thus delaying ventricular repolarization.
5
 On an electrocardiogram (ECG), this delay of repolarization manifests as a prolonged QT interval, which suggests that the heart takes more time to prepare itself for the next heartbeat. However, the underlying mechanism for arrhythmogenesis in patients with LQT2 remains to be elucidated due to the low prevalence of LQTS in general population and the great diversity across different subtypes and individuals.

The currently available treatments for LQTS are β-blocker drugs, left cardiac sympathetic denervation (LCSD), and implementable cardioverter-defibrillator (ICD). β-Blocker drugs are the first-line therapy for treating LQTS.
6
 The advantages of β-blocker therapy include its effectiveness in controlling heart rates and well tolerance in many patients. However, the β-blocker breakthrough events are not rare.7,8 Particularly, Priori et al.
9
 showed that LQT2 and LQT3 patients experience a significantly higher rate of cardiac events during β-blocker treatment than LQT1 patients. James et al.
10
 concluded that, despite an overall 71% risk reduction of cardiac events associated with exercise triggers, the β-blocker treatment could not fully protect the non-exercise events in LQT2 patients. LCSD is a surgical procedure that involves removing or cutting the nerves that supply the heart with sympathetic innervation. This procedure is typically reserved for patients who are unresponsive to β-blocker therapy,
11
 but it could also cause complications such as Horner's syndrome, which can cause drooping of the eyelid and constriction of the pupil.
12
 Finally, implantation of an ICD is another invasive procedure that involves placing a small device under the skin of the chest that monitors the heart's rhythm and delivers an electrical shock if a life-threatening arrhythmia is detected.
13
 Despite its contributions to the LQTS population, particularly those who are at the highest risk, as a surgical intervention, it may cause potential complications such as infection, bleeding, and device malfunction.
14
 Recently, Boutjdir and Lazzerini
15
 and Maguy et al.
16
 provided a novel antibody-based antiarrhythmic approach for the treatment of LQT2. They demonstrated that KCNQ1 antibodies, as a specific IKs agonist, could reverse the prolonged repolarization due to the loss of IKr, and abolish arrhythmic activities.
16
 Different from conventional treatments, antibody-based therapy has the advantage that it is targeted on the ion channel level per se to correct the ECG phenotype
15
; therefore, it provides a novel form of treatment and has a promising translational outcome for LQT2. However, it is still a proof-of-concept study, with all validations performed at the cellular level. The reasonable therapeutic dose of KCNQ1 antibody along with its actual effects on arrhythmogenesis biomarkers at tissue and organ levels remains to be elucidated.

Aiming at the above problems, in this study, we developed a multiscale computational cardiac model that ranges from ion channel level to three-dimensional (3D) organ level. To replicate the manifestations of the different LQT2 genotypes, we investigated both heterozygous and homozygous cases by extending the experimentally observed functional alterations at the subcellular level to the tissue and organ levels and analyzed their differences using metrics on different scales. Finally, we evaluated the immunotherapeutic effects of the KCNQ1 antibody on LQT2 by incorporating its modulations on ion channels into the virtual heart.

### Methods

The human ventricular myocyte model developed by O’hara et al.
17
 (O’Hara-Rudy dynamics model or ORd model) was utilized to simulate the electrophysiological behavior of human cardiomyocytes. The ORd model is a comprehensive cellular model developed based on human experimental data. To overcome its non-physiologically slow conduction velocity,
18
 the initial INa in the ORd model was replaced with the biophysically detailed model by ten Tusscher et al.
19
 The equation demonstrating the electrophysiology on a cellular level is as follows:(1)∂Vm∂t=−(Iion+Istim)Cmwhere Vm represents the membrane potential, Iion and Istim denote the total ionic current and the stimulating current separately, and Cm indicates the membrane capacitance. The ORd model was stimulated using 50 supra-threshold stimuli (−52 pA/pF, 1 ms) with a frequency of 1 Hz to reach its steady state, which corresponds to a normal heart rate of 60 beats/minute. The parameters of the last action potential (AP), such as action potential duration (APD), overshoot, and so on, were measured.

The effective refractory period (ERP) was assessed using a standard S1–S2 protocol. In particular, 50 S1 stimuli were delivered at 1 Hz prior to a premature stimulus (S2). The S2 had identical amplitude and duration as the S1 and was incapable of inducing a new AP if administered during the refractory period after the last S1. The procedure was performed in a step-by-step manner by progressively decreasing the S1–S2 interval and ERP was determined as the minimum diastolic interval at which the S2 stimulus was able to elicit an AP with an amplitude of 80% of the overshoot of the previous AP.

For wild-type (WT) cells, the maximal conductance of the rapidly delayed rectifier potassium channel remained unchanged to reflect a normal condition of its conducting current (IKr). As for the mutant cells, existing reports suggest that most pathogenic KCNH2 variants present the dominant negative effect with IKr <50% of WT IKr,20–24 while homozygous cases are extremely severe, with most cases ending up with embryonic lethality.20,25–28 To reflect the dominant negative effect, we considered a 75% loss of IKr for the heterozygous mutant type in this study, which was in consistent with two previous reports.20,21 In summary, for heterozygous and homozygous mutant types, the maximum conductance of IKr was reduced to 25% and 0% of the control value, respectively, to reflect the partial and complete loss of IKr.

According to the study by Maguy et al.,
16
 the KCNQ1 antibody increased IKs in a concentration-dependent manner, and such effect was achieved by modulating the channel activation properties (i.e. shifting the activation curve) without affecting the expression of IKs channel proteins. Based on this observation, a Hill function was fitted to mathematically describe the relationship between the voltage shift of the activation curve and the concentration of antibody:(2)Y=1201.0+(X/118.0930)−1.7499where Y represents the voltage shift of the activation curve and X denotes the concentration of the KCNQ1 antibody.

The corresponding concentration-dependent curve is plotted in Figure 1. This relationship was coupled into the ion channel model of IKs to simulate the therapeutic effects of different antibody concentrations on LQT2.

The fitted curve of the relationship between the voltage shift and the concentration of KCNQ1 antibody.

The propagation of excitation waves was described using the monodomain equation:(3)∂Vm∂t=∇⋅D(∇Vm)−Iion+IstimCmwhere D is the diffusion coefficient tensor for describing the intercellular electrical coupling via gap junctions. For 1D simulations, we constructed a transmural 1D strand consisting of 100 nodes with a spacing of 0.15 mm. The strand's length of 15 mm aligns with the typical range of human transmural ventricle width (4–14 mm) found in previous studies. The ratio of ENDO : MCELL : EPI (which represents endocardial, midmyocardial, and epicardial myocytes) was established at 25 : 35 : 40, following a previous study.
29
 For 1D simulation, the diffusion tensor degrades to a diffusion coefficient D. In this study, the value of D was set to 0.154 mm²/ms to achieve a conduction velocity of 71.9 cm/s, which closely approximates the recorded conduction velocity of 70 cm/s in human myocardium.

Despite the computational efficiency of the 1D model, it could not reflect the potential influence of fiber orientations on the propagation of excitation waves. Therefore, a 3D heterogeneous human left ventricular wedge model with fiber orientations extracted from diffusion tensor-magnetic resonance imaging was used.
30
 Orthotropic conduction was considered for this geometry.
31
 Particularly, given a local Cartesian coordinate system with reference to a specific myocardial fiber, the diffusion tensor can be represented as:(4)D~=(D∥000D⊥1000D⊥2)where 
D∥
, 
D⊥1
, and 
D⊥2
 are a group of diffusion coefficients describing the diffusion along the fiber, within the tissue sheet but perpendicular to the fiber, and normal to the sheet, respectively. The values were determined by adjusting conduction velocities to a ratio of 6 : 3 : 1, with the velocity along the fiber achieving 70 cm/s.
30

Due to that fibers in ventricles have their own orientations, it is required to project the local diffusion tensor 
D~
 to a global Cartesian coordinate system. This can be conducted by performing the inverse diagonalization:(5)D=AD~A−1where D is the diffusion tensor with regard to the global coordinate system and A is the matrix of eigenvectors f, s, n representing the three components of fiber orientations under the global coordinate system. Combining equations (3) to (5), the D is finally calculated as:(6)D=D∥ffT+D⊥1ssT+D⊥2nnT

The vulnerable window (VW) is a particular phase around the refractory tail of the preceding excitation wave. During the VW, a premature stimulus will evoke a phenomenon known as a unidirectional conduction block. In this study, we measured the VW by applying a standard S1–S2 protocol to a 1D transmural strand. Particularly, a series of S1 stimuli with an intensity of −52 pA/pF and a duration of 3 ms was applied to the endocardial end of the strand to initiate regular conditioning waves. The pacing cycle length of S1 stimuli was set to 1000 ms to mimic a sinus rhythm.
32
 During the final cycle of S1, a premature stimulus S2 with identical stimulating intensity and duration was applied to a local area on the strand at a predetermined time point to mimic the extra- and ectopic beat. If the S2 was applied too early, a bidirectional conduction block would occur, or if the S2 was applied too late, a bidirectional conduction would occur. But, if the S2 was applied within a certain period, it would induce the unidirectional conduction block. Therefore, the period when the unidirectional conduction block occurs was recorded as VW for that position. The aforementioned process was iterated across the entire strand to ascertain the comprehensive distribution of VWs.

The pseudo-ECG was generated using equations (7) and (8):(7)Φ(x′,y′,z′)=a2σi4σe∫(−∇Vm)⋅[∇1r]dΩ(8)r=[(x−x′)2+(y−y′)2+(z−z′)2]1/2where Φ is a unipolar potential generated by the tissue, r is the distance between a source point and the virtual electrode, σ
i
 and σe stand for intracellular and extracellular conductivities, respectively, and ∫ is the domain of integration. The QT interval was measured as the period from the earliest onset of the depolarization wave to the end of the T wave. The measured QT was corrected for heart rate using the Bazette equation
33
:(9)QTC=QT/RRwhere RR is the interval between two consecutive R waves measured in milliseconds.

In addition to the QT interval, two biomarkers regarding the T-wave are also measured. The first one is the duration from the peak to the end of the T-wave, denoted as Tpe. It is measured in milliseconds and is corrected for heart rate using the Bazette (TpeB) equations:(10)TpeB=Tpe/RRThe second one is the amplitude of the T-wave, which indicates the maximum voltage value of the T-wave on pseudo-ECG.

To ensure the robustness of our findings, we conducted parallel experiments using the human ventricular myocyte model developed by Ten Tusscher et al.
19
 (the TNNP06 model). Similarly, the IKr was reduced to 25% and 0% of the control value to model the heterozygous and homozygous conditions. These effects were extended to a 3D level, and the simulation results can be found in the Supplemental Material.

### Single-cell model

The human ventricular myocyte model developed by O’hara et al.
17
 (O’Hara-Rudy dynamics model or ORd model) was utilized to simulate the electrophysiological behavior of human cardiomyocytes. The ORd model is a comprehensive cellular model developed based on human experimental data. To overcome its non-physiologically slow conduction velocity,
18
 the initial INa in the ORd model was replaced with the biophysically detailed model by ten Tusscher et al.
19
 The equation demonstrating the electrophysiology on a cellular level is as follows:(1)∂Vm∂t=−(Iion+Istim)Cmwhere Vm represents the membrane potential, Iion and Istim denote the total ionic current and the stimulating current separately, and Cm indicates the membrane capacitance. The ORd model was stimulated using 50 supra-threshold stimuli (−52 pA/pF, 1 ms) with a frequency of 1 Hz to reach its steady state, which corresponds to a normal heart rate of 60 beats/minute. The parameters of the last action potential (AP), such as action potential duration (APD), overshoot, and so on, were measured.

The effective refractory period (ERP) was assessed using a standard S1–S2 protocol. In particular, 50 S1 stimuli were delivered at 1 Hz prior to a premature stimulus (S2). The S2 had identical amplitude and duration as the S1 and was incapable of inducing a new AP if administered during the refractory period after the last S1. The procedure was performed in a step-by-step manner by progressively decreasing the S1–S2 interval and ERP was determined as the minimum diastolic interval at which the S2 stimulus was able to elicit an AP with an amplitude of 80% of the overshoot of the previous AP.

For wild-type (WT) cells, the maximal conductance of the rapidly delayed rectifier potassium channel remained unchanged to reflect a normal condition of its conducting current (IKr). As for the mutant cells, existing reports suggest that most pathogenic KCNH2 variants present the dominant negative effect with IKr <50% of WT IKr,20–24 while homozygous cases are extremely severe, with most cases ending up with embryonic lethality.20,25–28 To reflect the dominant negative effect, we considered a 75% loss of IKr for the heterozygous mutant type in this study, which was in consistent with two previous reports.20,21 In summary, for heterozygous and homozygous mutant types, the maximum conductance of IKr was reduced to 25% and 0% of the control value, respectively, to reflect the partial and complete loss of IKr.

### Modeling the effects of KCNQ1 antibodies

According to the study by Maguy et al.,
16
 the KCNQ1 antibody increased IKs in a concentration-dependent manner, and such effect was achieved by modulating the channel activation properties (i.e. shifting the activation curve) without affecting the expression of IKs channel proteins. Based on this observation, a Hill function was fitted to mathematically describe the relationship between the voltage shift of the activation curve and the concentration of antibody:(2)Y=1201.0+(X/118.0930)−1.7499where Y represents the voltage shift of the activation curve and X denotes the concentration of the KCNQ1 antibody.

The corresponding concentration-dependent curve is plotted in Figure 1. This relationship was coupled into the ion channel model of IKs to simulate the therapeutic effects of different antibody concentrations on LQT2.

The fitted curve of the relationship between the voltage shift and the concentration of KCNQ1 antibody.

### One-dimensional (1D) strand

The propagation of excitation waves was described using the monodomain equation:(3)∂Vm∂t=∇⋅D(∇Vm)−Iion+IstimCmwhere D is the diffusion coefficient tensor for describing the intercellular electrical coupling via gap junctions. For 1D simulations, we constructed a transmural 1D strand consisting of 100 nodes with a spacing of 0.15 mm. The strand's length of 15 mm aligns with the typical range of human transmural ventricle width (4–14 mm) found in previous studies. The ratio of ENDO : MCELL : EPI (which represents endocardial, midmyocardial, and epicardial myocytes) was established at 25 : 35 : 40, following a previous study.
29
 For 1D simulation, the diffusion tensor degrades to a diffusion coefficient D. In this study, the value of D was set to 0.154 mm²/ms to achieve a conduction velocity of 71.9 cm/s, which closely approximates the recorded conduction velocity of 70 cm/s in human myocardium.

### Three-dimensional (3D) heterogeneous left ventricle wedge model

Despite the computational efficiency of the 1D model, it could not reflect the potential influence of fiber orientations on the propagation of excitation waves. Therefore, a 3D heterogeneous human left ventricular wedge model with fiber orientations extracted from diffusion tensor-magnetic resonance imaging was used.
30
 Orthotropic conduction was considered for this geometry.
31
 Particularly, given a local Cartesian coordinate system with reference to a specific myocardial fiber, the diffusion tensor can be represented as:(4)D~=(D∥000D⊥1000D⊥2)where 
D∥
, 
D⊥1
, and 
D⊥2
 are a group of diffusion coefficients describing the diffusion along the fiber, within the tissue sheet but perpendicular to the fiber, and normal to the sheet, respectively. The values were determined by adjusting conduction velocities to a ratio of 6 : 3 : 1, with the velocity along the fiber achieving 70 cm/s.
30

Due to that fibers in ventricles have their own orientations, it is required to project the local diffusion tensor 
D~
 to a global Cartesian coordinate system. This can be conducted by performing the inverse diagonalization:(5)D=AD~A−1where D is the diffusion tensor with regard to the global coordinate system and A is the matrix of eigenvectors f, s, n representing the three components of fiber orientations under the global coordinate system. Combining equations (3) to (5), the D is finally calculated as:(6)D=D∥ffT+D⊥1ssT+D⊥2nnT

### Measurement of a vulnerable window to the unidirectional conduction block

The vulnerable window (VW) is a particular phase around the refractory tail of the preceding excitation wave. During the VW, a premature stimulus will evoke a phenomenon known as a unidirectional conduction block. In this study, we measured the VW by applying a standard S1–S2 protocol to a 1D transmural strand. Particularly, a series of S1 stimuli with an intensity of −52 pA/pF and a duration of 3 ms was applied to the endocardial end of the strand to initiate regular conditioning waves. The pacing cycle length of S1 stimuli was set to 1000 ms to mimic a sinus rhythm.
32
 During the final cycle of S1, a premature stimulus S2 with identical stimulating intensity and duration was applied to a local area on the strand at a predetermined time point to mimic the extra- and ectopic beat. If the S2 was applied too early, a bidirectional conduction block would occur, or if the S2 was applied too late, a bidirectional conduction would occur. But, if the S2 was applied within a certain period, it would induce the unidirectional conduction block. Therefore, the period when the unidirectional conduction block occurs was recorded as VW for that position. The aforementioned process was iterated across the entire strand to ascertain the comprehensive distribution of VWs.

### Calculation of pseudo-ECG

The pseudo-ECG was generated using equations (7) and (8):(7)Φ(x′,y′,z′)=a2σi4σe∫(−∇Vm)⋅[∇1r]dΩ(8)r=[(x−x′)2+(y−y′)2+(z−z′)2]1/2where Φ is a unipolar potential generated by the tissue, r is the distance between a source point and the virtual electrode, σ
i
 and σe stand for intracellular and extracellular conductivities, respectively, and ∫ is the domain of integration. The QT interval was measured as the period from the earliest onset of the depolarization wave to the end of the T wave. The measured QT was corrected for heart rate using the Bazette equation
33
:(9)QTC=QT/RRwhere RR is the interval between two consecutive R waves measured in milliseconds.

In addition to the QT interval, two biomarkers regarding the T-wave are also measured. The first one is the duration from the peak to the end of the T-wave, denoted as Tpe. It is measured in milliseconds and is corrected for heart rate using the Bazette (TpeB) equations:(10)TpeB=Tpe/RRThe second one is the amplitude of the T-wave, which indicates the maximum voltage value of the T-wave on pseudo-ECG.

### Model-dependent test

To ensure the robustness of our findings, we conducted parallel experiments using the human ventricular myocyte model developed by Ten Tusscher et al.
19
 (the TNNP06 model). Similarly, the IKr was reduced to 25% and 0% of the control value to model the heterozygous and homozygous conditions. These effects were extended to a 3D level, and the simulation results can be found in the Supplemental Material.

### Results

This study modeled the results of WT and LQT2 mutants. We first conducted cellular-level simulations to obtain the membrane potential profiles of the three transmural cell types under pathological conditions. As shown in Figure 2, the model is able to capture the electrophysiological remodeling within heterozygous and homozygous types at the cellular level. First, it can be observed that all three cell types failed to repolarize under homozygous conditions (Figure 2A to C), indicating a rather severe or even lethal condition of homozygous cases. This observation is in line with previous clinical reports showing that homozygous cases are rare and most of them end up with embryonic lethality.26–28 For the heterozygous condition, repolarization was delayed in ENDO and EPI cells, as indicated by a significant prolongation of the repolarization time on the membrane potential traces (Figure 2D and F). The MCELL again failed to repolarize (Figure 2E), indicating the uneven distribution of IKr channels across different types of cells. Besides, the degree of APD prolongation was not the same (Figure 2G and H). For ENDO cells, APD was prolonged by 237.6 ms, while for EPI cells, APD was prolonged by 209.1 ms. This uneven prolongation was also seen in ERP.

Simulation results of action potentials. (A–C) Steady-state action potentials of ENDO, MCELL, and EPI cells under homozygous (red) conditions. (D–F) Steady-state action potentials of ENDO, MCELL, and EPI cells under wild-type (black) and heterozygous (orange) conditions. (G) APD90 of EPI, MCELL, and ENDO cells under wild-type (black) and heterozygous (orange) conditions. (H) Effective refractory period (ERP) of EPI, MCELL, and ENDO under wild-type (black) and heterozygous (orange) conditions. The asterisks (‘*’) in the bottom two panels indicate failure of repolarization.

We then built a 1D transmural strand model by coupling cell models and paced it with a cycle length of 1000 ms which corresponds to a heart rate of 60 beats/min. The model was paced several times to reach its steady state, and subsequently, we computed the pseudo-ECG based on the last cycle. Additionally, we measured the ECG, as well as the corresponding corrected QT using Bazett's formula (QTc) and the QTc, for both WT and heterozygous models over a BCL range from 600 to 1500 ms. The simulation results are shown in Figure 3.

Simulated pseudo-electrocardiogram (ECG) results and the restitution properties of ECG. (A and B) Spatiotemporal maps of excitation wave propagation under wild-type and heterozygous conditions. (C and D) The corresponding pseudo-ECGs under wild-type and heterozygous conditions.

Due to the failure of repolarization under homozygous conditions, we only tested the heterozygous and WT conditions. It can be observed in Figure 3A and B that the reduction in IKr resulted in a longer wavelength, suggesting a delayed depolarization in the mutant tissue. Note that there was no early afterdepolarization in the MCELL segment, this was due to the source-sink effect by cell coupling.

The pseudo-ECG in Figure 3C shows that the transmural strand model generated a proper ECG morphology, particularly the upward T-wave that indicated the right order of repolarization. For the WT condition, the model predicted a QT interval of 316 ms, while for the heterozygous condition, the QT interval was significantly prolonged to 566 ms. This manifestation reflected a longer repolarization process in the mutant strand than in the WT. In addition, the morphology of the T-wave was also altered. Specifically, the T-wave in the heterozygous condition displayed a smaller amplitude compared to that under the WT condition with a decrease ratio of approximately 5.2%. Besides, the peak-to-end interval of the T-wave (Tpe) was prolonged by nearly 22% (from 46 to 56 ms). It is well known that the increased Tpe reflected an augmented transmural dispersion of repolarization.34–36 These findings are consistent with the observations in clinical practice,22,24,37–39 suggesting the reliability of the developed pathological model. Detailed measurements are listed in Table 1.

Characteristics of the ECG morphology.

ECG: electrocardiogram; QTc: corrected QT using Bazett's formula; Tpe: Tpeak-to-Tend interval; TpeB: corrected Tpe using Bazett's formula.

To investigate how the QT prolongation responds to the changing heart rates, we measured the QT interval for WT and heterozygous conditions over a BCL ranging from 600 to 1500 ms, corresponding to heart rates from 100 to 40 bpm. Simulation results are shown in Figure 4. It can be found that, at any heart rate, the QT interval in heterozygous was significantly longer than in WT. In addition, the QT lengthening (i.e. ΔQT = QTHETE − QTWT) became larger (from 207 to 282 ms), suggesting that the condition became more severe and susceptible with heart rate slowed. In terms of the corrected QT interval (i.e. QTc), both WT and heterozygous cases showed decreased QTc with the heart rate going slower; however, the QTc ratio of heterozygous to WT, that is, QTcHETE/QTcWT, showed a consistent increase. Taken together, our simulations highlight that the QT prolongation in LQT is more pronounced during bradycardia, which is consistent with clinical observations.40–42

The restitution properties of QT interval. (A and B) Simulated pseudo-ECG results for wild-type and heterozygous conditions over a BCL range of 600–1500 ms. Arrows indicate the gradual prolongation of QT intervals as BCL increases. (C and D) Changes of QT, QTc, ΔQT, and ΔQTc as BCL increases. Black and orange traces indicate wild-type and heterozygous conditions, respectively. (E) The QTc ratio of heterozygous long QT (LQT) to wild-type, that is, QTcHETE/QTcWT.

The unidirectional conduction block is a critical phenomenon that predisposes to reentry arrhythmias. Therefore, the time period when a unidirectional conduction block can be evoked, termed VW, is an important biomarker to reflect the temporal vulnerability of the heart to reentry arrhythmias. The VW was measured using the S1–S2 protocol (see the Methods section), and the S2-induced unidirectional conduction block within the VW is illustrated in Figure 5.

Visualization of the vulnerable window (VW) and its distribution. (A and B) Panels from left to right show in order of the bidirectional conduction block (S2 applied too early), the unidirectional conduction block (within the VW), and the bidirectional conduction (S2 applied too late). (C) Distribution of VWs under wild-type (black) and heterozygous (orange) conditions. (D) Comparison of the average VW width under different conditions. The average VW widths for the wild-type and heterozygous groups were 5.22 and 31.10 ms, respectively.

To achieve an overall and comprehensive evaluation of the vulnerability, we tested the VW globally across the strand instead of only measuring that of a local site. The overall distribution of different groups, together with the averaged VW width, are plotted in Figure 5C and D. It can be observed that the occurrence of VWs in heterozygous conditions was delayed by approximately 30 ms compared to that in WT. Furthermore, the average VW in the heterozygous group was nearly five times greater compared to the WT group, indicating a significantly increased susceptibility to reentry arrhythmias. It is noteworthy that this increase is observed across all locations.

The above simulations based on idealized models have partially revealed the arrhythmogenesis mechanisms underlying LQT2. However, the conclusion is limited due to the absence of realistic ventricular anatomical structures. Therefore, we further conducted simulations showing the propagation of excitation propagation in a realistic human left ventricular wedge. The simulation results are shown in Figure 6. It can be observed that the three groups did not exhibit any difference during the depolarization phase, which was evidenced by the almost identical propagation pattern before 200 ms. The unaffected depolarization is consistent with the fact that IKr is only involved in the repolarization phase and its loss-of-function mutation brings no influence on the ventricular depolarization. As expected, discrepancies began to emerge between the two groups when entering the repolarization phase. Particularly, the ventricle of WT showed signs of repolarization in its middle part at around 300 ms, but this was not observable in the heterozygous group. By 425 ms, the WT ventricle had completed most of its repolarization process; in contrast, the heterozygote type was only partially repolarized (see the snapshots in Figure 6A and B for comparisons). The delayed repolarization was also reflected in the time length that back to the resting states. According to the simulation results, the WT ventricle returned to its resting state at around 425 ms, while the time was 775 ms under heterozygous conditions.

Propagation of excitation waves in the three-dimensional (3D) human left ventricular wedge model: (A) wild-type and (B) heterozygous conditions.

After the investigation of arrhythmia in LQT2, we then evaluated the therapeutic effects of the KCNQ1 antibody on it. Two concentrations of antibodies, 60 and 300 µg/mL, were assessed in terms of cellular and tissue-level biomarkers. As shown in the top panels of Figure 7, both concentrations of KCNQ1 antibody were able to accelerate repolarizations in all three cell types. However, the particular effects of the KCNQ1 antibody depended on its concentration and the cell types to which it was applied. At a concentration of 60 µg/mL, the KCNQ1 antibody suppressed the early afterdepolarization activities in MCELL, but its effects were not obvious in the other two types of cells, with only a slight reduction in APD. We then used an antibody concentration of 300 µg/mL. Although the APD was not fully recovered to its control level under this concentration, it was reduced substantially by approximately 32.4%.

Simulation results of the action of KCNQ1 antibody on cellular action potentials and pseudo-electrocardiograms (ECGs). (A–C) Steady-state action potentials of (A) ENDO, (B) MCELL, and (C) EPI under different conditions, including the wild-type (gray), heterozygous (orange), heterozygous with 60 µg/mL KCNQ1 antibody (green), and heterozygous with 300 µg/mL KCNQ1 antibody (blue). (D) Pseudo-ECGs under wild-type (gray), heterozygous (orange), and heterozygous with 60 µg/mL KCNQ1 antibody (green) conditions. (E) Pseudo-ECGs under wild-type (gray), heterozygous (orange), and heterozygous with 300 µg/mL KCNQ1 antibody (blue) conditions.

To provide a more clinically relevant evaluation of KCNQ1 antibody therapy, we tested the two concentrations using the transmural 1D strand and generated the pseudo-ECGs for the heterozygous group (Figure 7D and E). Overall, the use of the antibody resulted in an abbreviation of the QT interval in a concentration-dependent manner. In line with the cellular level observations, the KCNQ1 antibody at 60 µg/mL exerted only slight effects on the QT interval, abbreviating it from 566 to 526 ms. When the 300 µg/mL antibody was applied, the prolonged QT was significantly rescued, which was reduced by nearly 75% to 413 ms.

The rescuing effects of the KCNQ1 antibody were also straightforwardly observed in simulations of the excitation propagation in the 3D ventricular wedge. Figure 8 illustrates the repolarization phases of different groups, where the three rows represent the WT condition, the heterozygous condition, and the heterozygous wedge treated with 300 µg/mL KCNQ1 antibody, respectively. It can be observed that the KCNQ1 antibody rescued the delayed repolarization, indicating its capability to compensate for the abnormal repolarization process in the heart of LQT patients.

Propagation of excitation waves in the ventricular wedge with the application of KCNQ1 antibodies. (A) Wild-type condition. (B) Heterozygous condition. (C) Heterozygous group with 300 µg/mL KCNQ1 antibody.

In the previous section, we have shown that the mutation increased the susceptibility to arrhythmias (i.e. the increased VW under the mutation condition) by augmenting the transmural dispersion of repolarization. Therefore, we further measured the distribution of VW in the KCNQ1 antibody group to investigate the effect of antibodies on the transmural dispersion of repolarization. Simulation results are shown in Figure 9.

Vulnerable windows (VWs) under different conditions. (A) Distributions of VWs under different conditions. (B) The comparison of average VW width under different conditions.

Interestingly, we found that the width of VW was further increased after the application of the KCNQ1 antibody. As shown in Figure 9A, the VW was increased for most positions on the strand. Quantitative measurements in terms of the average VW width showed a 35% increase compared to the heterozygous group (from 31.10 to 42.02 ms) and became almost eight times the value in the WT condition. Since the increased VW is an important indicator of the tissue's temporal susceptibility to arrhythmias, the KCNQ1 antibody may be proarrhythmic in this aspect.

To sum up, our simulations suggest that the KCNQ1 antibody is indeed able to rescue the prolonged QT interval; however, it may also increase the temporal vulnerability to arrhythmias by further enlarging the VW.

### Simulation results of action potentials

This study modeled the results of WT and LQT2 mutants. We first conducted cellular-level simulations to obtain the membrane potential profiles of the three transmural cell types under pathological conditions. As shown in Figure 2, the model is able to capture the electrophysiological remodeling within heterozygous and homozygous types at the cellular level. First, it can be observed that all three cell types failed to repolarize under homozygous conditions (Figure 2A to C), indicating a rather severe or even lethal condition of homozygous cases. This observation is in line with previous clinical reports showing that homozygous cases are rare and most of them end up with embryonic lethality.26–28 For the heterozygous condition, repolarization was delayed in ENDO and EPI cells, as indicated by a significant prolongation of the repolarization time on the membrane potential traces (Figure 2D and F). The MCELL again failed to repolarize (Figure 2E), indicating the uneven distribution of IKr channels across different types of cells. Besides, the degree of APD prolongation was not the same (Figure 2G and H). For ENDO cells, APD was prolonged by 237.6 ms, while for EPI cells, APD was prolonged by 209.1 ms. This uneven prolongation was also seen in ERP.

Simulation results of action potentials. (A–C) Steady-state action potentials of ENDO, MCELL, and EPI cells under homozygous (red) conditions. (D–F) Steady-state action potentials of ENDO, MCELL, and EPI cells under wild-type (black) and heterozygous (orange) conditions. (G) APD90 of EPI, MCELL, and ENDO cells under wild-type (black) and heterozygous (orange) conditions. (H) Effective refractory period (ERP) of EPI, MCELL, and ENDO under wild-type (black) and heterozygous (orange) conditions. The asterisks (‘*’) in the bottom two panels indicate failure of repolarization.

### Simulation results of pseudo-ECGs under different conditions

We then built a 1D transmural strand model by coupling cell models and paced it with a cycle length of 1000 ms which corresponds to a heart rate of 60 beats/min. The model was paced several times to reach its steady state, and subsequently, we computed the pseudo-ECG based on the last cycle. Additionally, we measured the ECG, as well as the corresponding corrected QT using Bazett's formula (QTc) and the QTc, for both WT and heterozygous models over a BCL range from 600 to 1500 ms. The simulation results are shown in Figure 3.

Simulated pseudo-electrocardiogram (ECG) results and the restitution properties of ECG. (A and B) Spatiotemporal maps of excitation wave propagation under wild-type and heterozygous conditions. (C and D) The corresponding pseudo-ECGs under wild-type and heterozygous conditions.

Due to the failure of repolarization under homozygous conditions, we only tested the heterozygous and WT conditions. It can be observed in Figure 3A and B that the reduction in IKr resulted in a longer wavelength, suggesting a delayed depolarization in the mutant tissue. Note that there was no early afterdepolarization in the MCELL segment, this was due to the source-sink effect by cell coupling.

The pseudo-ECG in Figure 3C shows that the transmural strand model generated a proper ECG morphology, particularly the upward T-wave that indicated the right order of repolarization. For the WT condition, the model predicted a QT interval of 316 ms, while for the heterozygous condition, the QT interval was significantly prolonged to 566 ms. This manifestation reflected a longer repolarization process in the mutant strand than in the WT. In addition, the morphology of the T-wave was also altered. Specifically, the T-wave in the heterozygous condition displayed a smaller amplitude compared to that under the WT condition with a decrease ratio of approximately 5.2%. Besides, the peak-to-end interval of the T-wave (Tpe) was prolonged by nearly 22% (from 46 to 56 ms). It is well known that the increased Tpe reflected an augmented transmural dispersion of repolarization.34–36 These findings are consistent with the observations in clinical practice,22,24,37–39 suggesting the reliability of the developed pathological model. Detailed measurements are listed in Table 1.

Characteristics of the ECG morphology.

ECG: electrocardiogram; QTc: corrected QT using Bazett's formula; Tpe: Tpeak-to-Tend interval; TpeB: corrected Tpe using Bazett's formula.

### The restitution properties of QT interval

To investigate how the QT prolongation responds to the changing heart rates, we measured the QT interval for WT and heterozygous conditions over a BCL ranging from 600 to 1500 ms, corresponding to heart rates from 100 to 40 bpm. Simulation results are shown in Figure 4. It can be found that, at any heart rate, the QT interval in heterozygous was significantly longer than in WT. In addition, the QT lengthening (i.e. ΔQT = QTHETE − QTWT) became larger (from 207 to 282 ms), suggesting that the condition became more severe and susceptible with heart rate slowed. In terms of the corrected QT interval (i.e. QTc), both WT and heterozygous cases showed decreased QTc with the heart rate going slower; however, the QTc ratio of heterozygous to WT, that is, QTcHETE/QTcWT, showed a consistent increase. Taken together, our simulations highlight that the QT prolongation in LQT is more pronounced during bradycardia, which is consistent with clinical observations.40–42

The restitution properties of QT interval. (A and B) Simulated pseudo-ECG results for wild-type and heterozygous conditions over a BCL range of 600–1500 ms. Arrows indicate the gradual prolongation of QT intervals as BCL increases. (C and D) Changes of QT, QTc, ΔQT, and ΔQTc as BCL increases. Black and orange traces indicate wild-type and heterozygous conditions, respectively. (E) The QTc ratio of heterozygous long QT (LQT) to wild-type, that is, QTcHETE/QTcWT.

### Temporal vulnerability to the unidirectional conduction block in the transmural 1D strand

The unidirectional conduction block is a critical phenomenon that predisposes to reentry arrhythmias. Therefore, the time period when a unidirectional conduction block can be evoked, termed VW, is an important biomarker to reflect the temporal vulnerability of the heart to reentry arrhythmias. The VW was measured using the S1–S2 protocol (see the Methods section), and the S2-induced unidirectional conduction block within the VW is illustrated in Figure 5.

Visualization of the vulnerable window (VW) and its distribution. (A and B) Panels from left to right show in order of the bidirectional conduction block (S2 applied too early), the unidirectional conduction block (within the VW), and the bidirectional conduction (S2 applied too late). (C) Distribution of VWs under wild-type (black) and heterozygous (orange) conditions. (D) Comparison of the average VW width under different conditions. The average VW widths for the wild-type and heterozygous groups were 5.22 and 31.10 ms, respectively.

To achieve an overall and comprehensive evaluation of the vulnerability, we tested the VW globally across the strand instead of only measuring that of a local site. The overall distribution of different groups, together with the averaged VW width, are plotted in Figure 5C and D. It can be observed that the occurrence of VWs in heterozygous conditions was delayed by approximately 30 ms compared to that in WT. Furthermore, the average VW in the heterozygous group was nearly five times greater compared to the WT group, indicating a significantly increased susceptibility to reentry arrhythmias. It is noteworthy that this increase is observed across all locations.

### Simulation results using a realistic 3D human ventricular wedge model

The above simulations based on idealized models have partially revealed the arrhythmogenesis mechanisms underlying LQT2. However, the conclusion is limited due to the absence of realistic ventricular anatomical structures. Therefore, we further conducted simulations showing the propagation of excitation propagation in a realistic human left ventricular wedge. The simulation results are shown in Figure 6. It can be observed that the three groups did not exhibit any difference during the depolarization phase, which was evidenced by the almost identical propagation pattern before 200 ms. The unaffected depolarization is consistent with the fact that IKr is only involved in the repolarization phase and its loss-of-function mutation brings no influence on the ventricular depolarization. As expected, discrepancies began to emerge between the two groups when entering the repolarization phase. Particularly, the ventricle of WT showed signs of repolarization in its middle part at around 300 ms, but this was not observable in the heterozygous group. By 425 ms, the WT ventricle had completed most of its repolarization process; in contrast, the heterozygote type was only partially repolarized (see the snapshots in Figure 6A and B for comparisons). The delayed repolarization was also reflected in the time length that back to the resting states. According to the simulation results, the WT ventricle returned to its resting state at around 425 ms, while the time was 775 ms under heterozygous conditions.

Propagation of excitation waves in the three-dimensional (3D) human left ventricular wedge model: (A) wild-type and (B) heterozygous conditions.

### Simulation results of the therapeutic effects of KCNQ1 antibody

After the investigation of arrhythmia in LQT2, we then evaluated the therapeutic effects of the KCNQ1 antibody on it. Two concentrations of antibodies, 60 and 300 µg/mL, were assessed in terms of cellular and tissue-level biomarkers. As shown in the top panels of Figure 7, both concentrations of KCNQ1 antibody were able to accelerate repolarizations in all three cell types. However, the particular effects of the KCNQ1 antibody depended on its concentration and the cell types to which it was applied. At a concentration of 60 µg/mL, the KCNQ1 antibody suppressed the early afterdepolarization activities in MCELL, but its effects were not obvious in the other two types of cells, with only a slight reduction in APD. We then used an antibody concentration of 300 µg/mL. Although the APD was not fully recovered to its control level under this concentration, it was reduced substantially by approximately 32.4%.

Simulation results of the action of KCNQ1 antibody on cellular action potentials and pseudo-electrocardiograms (ECGs). (A–C) Steady-state action potentials of (A) ENDO, (B) MCELL, and (C) EPI under different conditions, including the wild-type (gray), heterozygous (orange), heterozygous with 60 µg/mL KCNQ1 antibody (green), and heterozygous with 300 µg/mL KCNQ1 antibody (blue). (D) Pseudo-ECGs under wild-type (gray), heterozygous (orange), and heterozygous with 60 µg/mL KCNQ1 antibody (green) conditions. (E) Pseudo-ECGs under wild-type (gray), heterozygous (orange), and heterozygous with 300 µg/mL KCNQ1 antibody (blue) conditions.

To provide a more clinically relevant evaluation of KCNQ1 antibody therapy, we tested the two concentrations using the transmural 1D strand and generated the pseudo-ECGs for the heterozygous group (Figure 7D and E). Overall, the use of the antibody resulted in an abbreviation of the QT interval in a concentration-dependent manner. In line with the cellular level observations, the KCNQ1 antibody at 60 µg/mL exerted only slight effects on the QT interval, abbreviating it from 566 to 526 ms. When the 300 µg/mL antibody was applied, the prolonged QT was significantly rescued, which was reduced by nearly 75% to 413 ms.

The rescuing effects of the KCNQ1 antibody were also straightforwardly observed in simulations of the excitation propagation in the 3D ventricular wedge. Figure 8 illustrates the repolarization phases of different groups, where the three rows represent the WT condition, the heterozygous condition, and the heterozygous wedge treated with 300 µg/mL KCNQ1 antibody, respectively. It can be observed that the KCNQ1 antibody rescued the delayed repolarization, indicating its capability to compensate for the abnormal repolarization process in the heart of LQT patients.

Propagation of excitation waves in the ventricular wedge with the application of KCNQ1 antibodies. (A) Wild-type condition. (B) Heterozygous condition. (C) Heterozygous group with 300 µg/mL KCNQ1 antibody.

### The effects of KCNQ1 antibody on transmural dispersion of repolarization

In the previous section, we have shown that the mutation increased the susceptibility to arrhythmias (i.e. the increased VW under the mutation condition) by augmenting the transmural dispersion of repolarization. Therefore, we further measured the distribution of VW in the KCNQ1 antibody group to investigate the effect of antibodies on the transmural dispersion of repolarization. Simulation results are shown in Figure 9.

Vulnerable windows (VWs) under different conditions. (A) Distributions of VWs under different conditions. (B) The comparison of average VW width under different conditions.

Interestingly, we found that the width of VW was further increased after the application of the KCNQ1 antibody. As shown in Figure 9A, the VW was increased for most positions on the strand. Quantitative measurements in terms of the average VW width showed a 35% increase compared to the heterozygous group (from 31.10 to 42.02 ms) and became almost eight times the value in the WT condition. Since the increased VW is an important indicator of the tissue's temporal susceptibility to arrhythmias, the KCNQ1 antibody may be proarrhythmic in this aspect.

To sum up, our simulations suggest that the KCNQ1 antibody is indeed able to rescue the prolonged QT interval; however, it may also increase the temporal vulnerability to arrhythmias by further enlarging the VW.

### Discussion

Compared with the complicated channelopathies underlying various acquired arrhythmias,43–47 the LQTS is relatively straightforward in terms of its alterations to the ion channel functions. In most cases, there is only one type of channel involved in LQTS, and the sole factor in LQT2 is the loss of IKr. Despite the relatively simple cause of LQT2, the interplay among myriad ion channels together with the various phenotypes of different cell types complicate the underlying mechanisms responsible for arrhythmias. By constructing a multiscale computational cardiac model, we were able to simulate how the microscopic ionic alteration affects the electrophysiological activities at tissue and organ levels, and also to quantitative the risks to arrhythmias by measuring the biomarkers of different levels. The main findings are as follows: (i) at the cellular level, despite an increase of ERP in all three cell types (which was usually considered protective against arrhythmias), the increment was uneven, and therefore led to an augmented transmural dispersion of repolarization. Besides, the diminished repolarization reserve due to the loss of IKr predisposed cells to early afterdepolarizations; (ii) the increased VW suggested that patients with LQT2 had higher risks of developing unidirectional conduction blocks and reentry arrhythmias. The reentrant excitation waves could further evolve into Torsade de Pointes (TdP) if wandering around the ventricle; (iii) the restitution property was remodeled in mutant tissue, which was manifested as the more pronounced QT prolongation in slow heart rates. Importantly, this observation agrees well with the clinical findings that the LQT2 patients experience higher risks under the bradycardia condition; (iv) in terms of ECG, in addition to the prolonged QT interval as expected, the loss of IKr also increased the Tpe. Increased Tpe was widely acknowledged as a surrogate biomarker for increased dispersion of repolarization; (v) observations that the mutation delayed the repolarization phase were consistent among different scales, regardless of the geometric model used for simulation, which suggested that other factors such as fiber orientations, anatomical structures, or ventricular heterogeneity would not affect the electrophysiological remodeling effects by the loss-of-function mutations in LQT2.

Our simulation results, together with the experimental observations, demonstrate that KCNQ1 antibodies are indeed able to abbreviate the cellular APD and ERP, therefore restoring the prolonged QT interval. These findings indicate that the increased IKs by KCNQ1 antibodies accelerate the repolarization process of the AP, which compensates for IKr deficiency and offsets the APD prolongation in mutation groups. Further simulations using the 3D heterogeneous ventricular wedge show consistent observation. Despite these positive effects, we found that the KCNQ1 antibody may not be beneficial regarding the dispersion of repolarization. Particularly, the KCNQ1 antibody did not reduce the width of the VW; instead, the VW was further enlarged (see Figure 9B). Widened VW means an increased risk to the unidirectional conduction block and the reentry arrhythmia, in other words, the KCNQ1 antibody actually exerts proarrhythmic effects in this aspect. Manifested on the ECG, the augmented VW results in an increased TpeB, and the latter is a common clinical marker hinting at an increased risk of developing arrhythmias.48–50 Taking together, whether the restored QTc represents a full recovery of arrhythmogenesis warrants further investigations.

In this study, we have used two different cell models, that is, TNNP06 and ORd, from two individual research groups to verify the robustness of our findings (refer to the Supplemental material for the simulation results based on TNNP06). Although most of the results were consistent, these two models reflected different therapeutic concentrations of KCNQ1 antibody. Specifically, applying the antibody of 60 µg/mL was enough to recover both the APD and the QT interval to its control level (see Supplemental Figure S5); however, its effects were not obvious according to the ORd model, where only a slight reduction in APD, and the restoration of QT interval was not satisfactory either (Figure 7). In fact, the concentration was elevated to 300 µg/mL to eventually achieve a therapeutic effect, according to the simulation results from the ORd model.

After a comprehensive investigation, we found that the above discrepancy arose from the relatively small amplitude of IKs in the ORd model, where the IKr dominated the repolarization force. This was evidenced by the fact that blocking IKr led to the failure of repolarization in all cell types of the ORd model. In contrast, 100% IKr blocking caused limited impacts on APD in the TNNP06 model, while IKs were relatively large and played a major role in the repolarization phase. The above observation explained why the effects of enhancing IKs (by KCNQ1 antibody) were more obvious in the TNNP06 model. However, the assumption in TNNP06 was not able to reflect the severe condition of homozygous cases in previous reports.20,25–28 Therefore, we finally chose the ORd model as the basic cell model for the whole simulation study. Nevertheless, we have conducted parallel experiments using both the two models and provided all the simulation results from cell to organ levels. Further investigations regarding the therapeutic concentration are warranted.

### Mechanisms of arrhythmogenesis in LQT2

Compared with the complicated channelopathies underlying various acquired arrhythmias,43–47 the LQTS is relatively straightforward in terms of its alterations to the ion channel functions. In most cases, there is only one type of channel involved in LQTS, and the sole factor in LQT2 is the loss of IKr. Despite the relatively simple cause of LQT2, the interplay among myriad ion channels together with the various phenotypes of different cell types complicate the underlying mechanisms responsible for arrhythmias. By constructing a multiscale computational cardiac model, we were able to simulate how the microscopic ionic alteration affects the electrophysiological activities at tissue and organ levels, and also to quantitative the risks to arrhythmias by measuring the biomarkers of different levels. The main findings are as follows: (i) at the cellular level, despite an increase of ERP in all three cell types (which was usually considered protective against arrhythmias), the increment was uneven, and therefore led to an augmented transmural dispersion of repolarization. Besides, the diminished repolarization reserve due to the loss of IKr predisposed cells to early afterdepolarizations; (ii) the increased VW suggested that patients with LQT2 had higher risks of developing unidirectional conduction blocks and reentry arrhythmias. The reentrant excitation waves could further evolve into Torsade de Pointes (TdP) if wandering around the ventricle; (iii) the restitution property was remodeled in mutant tissue, which was manifested as the more pronounced QT prolongation in slow heart rates. Importantly, this observation agrees well with the clinical findings that the LQT2 patients experience higher risks under the bradycardia condition; (iv) in terms of ECG, in addition to the prolonged QT interval as expected, the loss of IKr also increased the Tpe. Increased Tpe was widely acknowledged as a surrogate biomarker for increased dispersion of repolarization; (v) observations that the mutation delayed the repolarization phase were consistent among different scales, regardless of the geometric model used for simulation, which suggested that other factors such as fiber orientations, anatomical structures, or ventricular heterogeneity would not affect the electrophysiological remodeling effects by the loss-of-function mutations in LQT2.

### KCNQ1 antibody rescued the prolonged QT interval but not the augmented dispersion of repolarization

Our simulation results, together with the experimental observations, demonstrate that KCNQ1 antibodies are indeed able to abbreviate the cellular APD and ERP, therefore restoring the prolonged QT interval. These findings indicate that the increased IKs by KCNQ1 antibodies accelerate the repolarization process of the AP, which compensates for IKr deficiency and offsets the APD prolongation in mutation groups. Further simulations using the 3D heterogeneous ventricular wedge show consistent observation. Despite these positive effects, we found that the KCNQ1 antibody may not be beneficial regarding the dispersion of repolarization. Particularly, the KCNQ1 antibody did not reduce the width of the VW; instead, the VW was further enlarged (see Figure 9B). Widened VW means an increased risk to the unidirectional conduction block and the reentry arrhythmia, in other words, the KCNQ1 antibody actually exerts proarrhythmic effects in this aspect. Manifested on the ECG, the augmented VW results in an increased TpeB, and the latter is a common clinical marker hinting at an increased risk of developing arrhythmias.48–50 Taking together, whether the restored QTc represents a full recovery of arrhythmogenesis warrants further investigations.

### Potential limitations of this study

In this study, we have used two different cell models, that is, TNNP06 and ORd, from two individual research groups to verify the robustness of our findings (refer to the Supplemental material for the simulation results based on TNNP06). Although most of the results were consistent, these two models reflected different therapeutic concentrations of KCNQ1 antibody. Specifically, applying the antibody of 60 µg/mL was enough to recover both the APD and the QT interval to its control level (see Supplemental Figure S5); however, its effects were not obvious according to the ORd model, where only a slight reduction in APD, and the restoration of QT interval was not satisfactory either (Figure 7). In fact, the concentration was elevated to 300 µg/mL to eventually achieve a therapeutic effect, according to the simulation results from the ORd model.

After a comprehensive investigation, we found that the above discrepancy arose from the relatively small amplitude of IKs in the ORd model, where the IKr dominated the repolarization force. This was evidenced by the fact that blocking IKr led to the failure of repolarization in all cell types of the ORd model. In contrast, 100% IKr blocking caused limited impacts on APD in the TNNP06 model, while IKs were relatively large and played a major role in the repolarization phase. The above observation explained why the effects of enhancing IKs (by KCNQ1 antibody) were more obvious in the TNNP06 model. However, the assumption in TNNP06 was not able to reflect the severe condition of homozygous cases in previous reports.20,25–28 Therefore, we finally chose the ORd model as the basic cell model for the whole simulation study. Nevertheless, we have conducted parallel experiments using both the two models and provided all the simulation results from cell to organ levels. Further investigations regarding the therapeutic concentration are warranted.

### Conclusion

In this study, we developed a multiscale virtual heart for LQT2 and investigated its arrhythmogenesis mechanisms by monitoring the changes in biomarkers at different biological levels. Simulations revealed the pathological consequence due to the tiny alteration at the molecular level, and how these cardiac electrical remodeling may contribute to higher susceptibility to lethal arrhythmias in LQT2 patients. Besides, we also evaluated the therapeutic effects of an emerging treatment technique, KCNQ1 antibody, along with its potential proarrhythmic risks. Overall, this study establishes a causal link between genetic mutations and organ-level ventricular arrhythmias as well as clinical ECG manifestations and provides an in silico way to assess treatment options for patients with LQT2.

### Supplemental Material

Supplemental material, sj-pdf-1-dhj-10.1177_20552076241277032 for Computational analysis of long QT syndrome type 2 and the therapeutic effects of KCNQ1 antibodies by Zhujun Pan, Qi Fu, Huasen Jiang, Zhiqiang Wei and Shugang Zhang in DIGITAL HEALTH



# SUPPLEMENTAL FILE 1: 10.1177_20552076241277032.pdf

# Preparing to download ...

[HHS Vulnerability Disclosure](https://www.hhs.gov/vulnerability-disclosure-policy/index.html)