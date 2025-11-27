# MAIN TEXT

## Study of the parent-of-origin effect in monogenic diseases with variable age of onset. Application on ATTRv

### Abstract

In genetic diseases with variable age of onset, an accurate estimation of the survival function for the mutation carriers and also modifying factors effects estimations are important for the management of asymptomatic gene carriers across life. Among the modifying factors, the gender of the parent transmitting the mutation (i.e. the parent-of-origin effect) has been shown to have a significant effect on survival curve estimation on transthyretin familial amyloid polyneuropathy (ATTRv) families. However, as most genotypes are unknown, the parent-of-origin must be calculated through a probability estimated from the pedigree. We propose in this article to extend the method providing mutation carrier survival estimates in order to estimate the parent-of-origin effect. The method is both validated on simulated data and applied to familly samples with ATTRv.

### Introduction

In variable age of onset diseases caused by a deleterious gene mutation, knowledge of the survival function for carrier individual is important both to understand the underlying mechanism of the disease (like identification of potential factors that modulate this age of onset) and for prevention strategies. In the literature on this subject, the age-specific cumulative distribution function (CDF), named also penetrance function, is preferentially used. In this paper, like in [1], we will use the classical survival function (which is simply the complementary of the CDF) to assess the probability of not being affected by the disease according to the age for mutation carrier individuals. Note that our survival function hence corresponds to the cause-specific survival (disease diagnosis) and not to the overall survival.

Recently, a semi-parametric method based on a Cox model was developed to estimate survival function from familial data ascertained through affected individual. This method is very efficient and handle the unknown genotypes through the sum-product algorithm in Bayesian network. The method is described in detail in [1]. This method has been applied on families affected by a mendelian disease: the ATTRv, that is the most frequent familial amyloidosis, with autosomal dominant transmission. A fatal outcome occurs after an average duration of 10-13 years [2, 3]. This severe diseases shows important differences in age of onset and thus on survival curve, according to different covariates as countries, gender, mutations [4–6].

Furthermore, differences of age of onset according to gender of the transmitting parent have been mentioned in Portuguese families [7]. In Sweden families, the age of onset is variable among pedigrees and is significantly higher with paternal inheritance compared with that of maternal inheritance [8]. These first observations lead to the study of the heterogeneity in the survival curves according to gender of the transmitting parent and such heterogeneity have been noted on Portugueses ATTRv families [7, 9]. In particular, Hellman et al. [4] have find that the risk of disease in the carriers was significantly higher when the mutation was inherited from the mother than from the father.

In 2009, Bonaïti et al [10] investigate the parent-of-origin effect in a sample of French and Portuguese families that have already been described in [11]. The covariate for parent of origin factor was not calculated algorithmically according to pedigrees and genotypical information but was determined previously to the estimation of the parameters and decided manually by expert-based review. They found that the penetrance was higher when the mutation was inherited from the mother than from the father (i.e. the survival curve was lower when the mutation was inherited from the mother than from the father). This difference was significant in the Portuguese families but not in the French sample. This results led them to the hypothesis of a genetically determined effect through an imprinting phenomenon.

Unfortunately, the gender of the transmitting parent is not known for most of carriers due to the fact that most of individuals are non-genotyped. And we should be skeptical as to the estimates made in this context. However, no simulation studies have been proposed to evaluate the quality of survival curve estimations in studying the parent-of-origin effect.

In this article, we propose to extend the semi-parametric method based on a Cox model described in [1] in order to determine the gender of the transmitting parent. The method for analysing the parent-of-origin effect in survival curve estimations of Mendelian disease when a part of genotypic information is unknown, is validated through simulations.

To do the link with the literature about the subject, we illustrate our results on Portuguese and French datasets of ATTRv already analyzed in [10] and compare our results with theirs.

### Methods

Survival curve are estimated with a semi-parametric method based on a Cox model and adapted for pedigree data. As mentioned before, the data available are families ascertained through at least one affected individual. For all individuals, we have information about his genotype, phenotype, age, gender, etc. Families data consist of individuals independents from each other conditionally to their genotype. The notion of family is taken into account exclusively in the estimation of the probability to carry the mutation. The method is described in detail in [1]. In this section, the method is briefly presented in the particular case of studying the parent-of-origin effect.

If we denote by ev the evidence which consists of all the available information: all time Ti and status δi, the individual genotype Xi and covariates Zi (possibly multidimensional: gender, comorbidities, ethnicity, etc.), and the genetic testing Gi (when available, NA is missing) we have:
P(X,ev)=∏iP(Ti,δi|Xi,Zi)P(Gi|Xi)1Gi≠NA︸ϕi(Xi)∏i∈FP(Xi)∏i∉FP(Xi|Xpati,Xmati)
where F is the set of founders and where Xpati and Xmati indicate respectively the genotype of the father and the mother on individual i. Thus, the genotype distribution among founders follows the Hardy-Weinberg equilibrium with disease allele frequency q, and the conditional distribution for non-founders follows the Mendelian transmission of alleles. If Gi = NA, then Xi ∈ {0, 1m, 1p, 2} corresponding respectively to a non-mutated individual, a mutated individual whose mutation comes from the father, a mutated individual whose mutation comes from the mother and finally, to a homozygous individual.

More precisely, we have:
P(Ti=t,δi=0|Xi=x,Zi)={exp(-Λ0(t)exp(Ziγ))ifx=1morx=2exp(-Λ0(t)exp(β+Ziγ))ifx=1p1ifXi=0
and
P(Ti=t,δi=1|Xi=x,Zi)={exp(-Λ0(t)exp(Ziγ))×(λ0(t)exp(Ziγ))ifx=1morx=2exp(-Λ0(t)exp(β+Ziγ))×(λ0(t)exp(β+Ziγ))ifx=1p0ifXi=0
which can be simplified as follows:
1λ0(t)P(Ti=t,δi=1|Xi=x,Zi)={exp(-Λ0(t)exp(Ziγ))exp(Ziγ)ifx=1morx=2exp(-Λ0(t)exp(β+Ziγ))exp(β+Ziγ)ifx=1p0ifXi=0
where λ0 is the baseline hazard and Λ0 the baseline cumulative hazard. Moreover, β is the Cox’s model parameter, to be estimated, for the parent-of-origin effect and γ the parameter vector corresponding to the other covariates. For the genetic testings, we have a simple model with false positive rate ε and false negative rate η:
P(Gi=1|Xi≠0)=1-εP(Gi=0|Xi=0)=1-η

Note that, in our case, these rates are probably very small (e.g. ε < 1/100 and η < 1/1000).

As Xi is either partially observed or not observed, we consider this variable as latent and use a classical Expectation-Maximization algorithm in order to maximize the log-likelihood model in parameter of interest.

Model parameters are: q, β, γ, ε, η, Λ0 and λ0. In order for the model to be identifiable, we assume that error rates ε and η are known as well as the disease allele frequency q. And since λ0 appears only has a proportional factor in the expression of ϕi(Xi) we can perform model inference without explicit value for this parameter. Our aim is therefore to estimate θ = (β, γ, Λ0) using a classical EM framework.

In this framework we alternate to steps until convergence:

E-Step give the weights wi using current parameter θold, computed for all i:
wipat=P(Xi=1p|ev;θold)andwimat=P(Xi=1morXi=2|ev;θold)

M-Step: create an artificial weighted dataset with the following 2n patients

and simply:

fit a (weighted) Cox’s proportional hazard model with Factor POO and the covariates to update β and γ

use non-parametric estimate (e.g. Kaplan-Meier or similar) of Λ0.

For computing the posterior weights wipat and wimat of the E-Step one has to integrate the model over all unobserved (even when partially observed) genotypes Xi. Since the total number of configurations for Xi is 4n it is clearly impossible to perform this summation with brute force. Fortunately, geneticist know for year that likelihood computations in genetic model can be performed efficiently using the Elston-Stewart algorithm [12]. Thanks to this algorithm, it is therefore possible to compute any posterior distribution as a simple likelihood ratio:
P(Xi=x|ev)=P(Xi=x,ev)P(ev)

In the probabilistic graphical model community, such computation can be dealt efficiently using the sum-product algorithm (also called: belief propagation, forward/backward, inward/outward, Kalman filter, etc.) in order to obtain in a single computational pass all
wipat and wimat. NB: in Totir [13], the authors suggest a generalization of the Elston-Stewart algorithm allowing to compute all posterior distribution in a single passe rather than repeating likelhood computations. Without surprise, Totir’s algorithm is in fact the exact reformulation of the classical sum-product algorithm.

Practically, initialization is performed by affecting random weights wi (ex: drawn from a uniform distribution on [0, 1]). EM iterations are stopped when we observe convergence on test survival estimates (ex: baseline survival at age 20, 40, 60, 80).

Research approved by ethic committee #00011558 (Institutional Review Board, Mondor Hospital, Creteil, France) To access the data, please contact the Institutional Review Board, Mondor Hospital, Creteil, France irb.mondor@aphp.fr.

### The model

If we denote by ev the evidence which consists of all the available information: all time Ti and status δi, the individual genotype Xi and covariates Zi (possibly multidimensional: gender, comorbidities, ethnicity, etc.), and the genetic testing Gi (when available, NA is missing) we have:
P(X,ev)=∏iP(Ti,δi|Xi,Zi)P(Gi|Xi)1Gi≠NA︸ϕi(Xi)∏i∈FP(Xi)∏i∉FP(Xi|Xpati,Xmati)
where F is the set of founders and where Xpati and Xmati indicate respectively the genotype of the father and the mother on individual i. Thus, the genotype distribution among founders follows the Hardy-Weinberg equilibrium with disease allele frequency q, and the conditional distribution for non-founders follows the Mendelian transmission of alleles. If Gi = NA, then Xi ∈ {0, 1m, 1p, 2} corresponding respectively to a non-mutated individual, a mutated individual whose mutation comes from the father, a mutated individual whose mutation comes from the mother and finally, to a homozygous individual.

More precisely, we have:
P(Ti=t,δi=0|Xi=x,Zi)={exp(-Λ0(t)exp(Ziγ))ifx=1morx=2exp(-Λ0(t)exp(β+Ziγ))ifx=1p1ifXi=0
and
P(Ti=t,δi=1|Xi=x,Zi)={exp(-Λ0(t)exp(Ziγ))×(λ0(t)exp(Ziγ))ifx=1morx=2exp(-Λ0(t)exp(β+Ziγ))×(λ0(t)exp(β+Ziγ))ifx=1p0ifXi=0
which can be simplified as follows:
1λ0(t)P(Ti=t,δi=1|Xi=x,Zi)={exp(-Λ0(t)exp(Ziγ))exp(Ziγ)ifx=1morx=2exp(-Λ0(t)exp(β+Ziγ))exp(β+Ziγ)ifx=1p0ifXi=0
where λ0 is the baseline hazard and Λ0 the baseline cumulative hazard. Moreover, β is the Cox’s model parameter, to be estimated, for the parent-of-origin effect and γ the parameter vector corresponding to the other covariates. For the genetic testings, we have a simple model with false positive rate ε and false negative rate η:
P(Gi=1|Xi≠0)=1-εP(Gi=0|Xi=0)=1-η

Note that, in our case, these rates are probably very small (e.g. ε < 1/100 and η < 1/1000).

### EM framework

As Xi is either partially observed or not observed, we consider this variable as latent and use a classical Expectation-Maximization algorithm in order to maximize the log-likelihood model in parameter of interest.

Model parameters are: q, β, γ, ε, η, Λ0 and λ0. In order for the model to be identifiable, we assume that error rates ε and η are known as well as the disease allele frequency q. And since λ0 appears only has a proportional factor in the expression of ϕi(Xi) we can perform model inference without explicit value for this parameter. Our aim is therefore to estimate θ = (β, γ, Λ0) using a classical EM framework.

In this framework we alternate to steps until convergence:

E-Step give the weights wi using current parameter θold, computed for all i:
wipat=P(Xi=1p|ev;θold)andwimat=P(Xi=1morXi=2|ev;θold)

M-Step: create an artificial weighted dataset with the following 2n patients

and simply:

fit a (weighted) Cox’s proportional hazard model with Factor POO and the covariates to update β and γ

use non-parametric estimate (e.g. Kaplan-Meier or similar) of Λ0.

### Posterior distributions

For computing the posterior weights wipat and wimat of the E-Step one has to integrate the model over all unobserved (even when partially observed) genotypes Xi. Since the total number of configurations for Xi is 4n it is clearly impossible to perform this summation with brute force. Fortunately, geneticist know for year that likelihood computations in genetic model can be performed efficiently using the Elston-Stewart algorithm [12]. Thanks to this algorithm, it is therefore possible to compute any posterior distribution as a simple likelihood ratio:
P(Xi=x|ev)=P(Xi=x,ev)P(ev)

In the probabilistic graphical model community, such computation can be dealt efficiently using the sum-product algorithm (also called: belief propagation, forward/backward, inward/outward, Kalman filter, etc.) in order to obtain in a single computational pass all
wipat and wimat. NB: in Totir [13], the authors suggest a generalization of the Elston-Stewart algorithm allowing to compute all posterior distribution in a single passe rather than repeating likelhood computations. Without surprise, Totir’s algorithm is in fact the exact reformulation of the classical sum-product algorithm.

Practically, initialization is performed by affecting random weights wi (ex: drawn from a uniform distribution on [0, 1]). EM iterations are stopped when we observe convergence on test survival estimates (ex: baseline survival at age 20, 40, 60, 80).

Research approved by ethic committee #00011558 (Institutional Review Board, Mondor Hospital, Creteil, France) To access the data, please contact the Institutional Review Board, Mondor Hospital, Creteil, France irb.mondor@aphp.fr.

### Results

We have simulated n families (n = 100 or 400) with 10 individuals as shown in Fig 1.

Genotypes were assigned respecting Mendelian transmission and heterozygous genotypes were ordered according to the gender (pat or mat) of the parent who transmitted the mutation. So, a new binary covariate (named POO for parent of origin) are introduced such as POO = mat if the mother transmitted the disease mutation and POO = pat if disease mutation was transmitted by the father. The disease allele frequency was set to q = 0.20 in our simulated dataset for the sake of speed (without simulating any ascertainment process). The age at event was simulated according to a piecewise constant hazard rate function. Moreover a censoring variable was added that follows a uniform distribution U[15,80].

To study the parent-of-origin effect, the POO was simulated in a proportional hazard model. The risk when the mutation is transmitted by the mother is given by λ0 as follows and the coefficient parameters in the hazard model is noted β (β = −0.60 or β = −1.2 in simulations). Thus, as β < 0, that’s mean that survival is upper when mutation provides on the father.
λ0(t)={0ift∈[0,20]0.02ift∈[20,40]0.10ift∈[40,60]0.05ift>60

To study the behavior of our method, we will consider different scenarios for simulations. In the first scenario (S0), all genotypes are set to unobserved, in the second more realistic scenario (S1), 80% of genotypes are observed in affected individuals and only 10% are observed in non affected individuals and in the third one (S2), all genotypes are observed. Moreover, we consider an other scenario where, in addition to the genotypes, we know, for each individual, the sex of the parent transmitting the mutation. We call this scenario “Oracle”. In order to compare the β estimate in the different scenarios, 200 replications are performed.

Fig 2 shows the violin plots according to the scenario for the percentage of observed genotypes, which allows to show the full distribution of the data. Thus, each color corresponds to a genotypic scenario for each of the cases A (n = 100 families was simulated with β = −0.6), B (n = 400 families was simulated with β = −0.6) and C (n = 100 families was simulated with β = −1.2). For each scenario, the mean and standard deviation are shown and the horizontal black line represents the true β value. In all cases the β estimate are unbiased whatever the scenario and the standard deviation increases with the percentage of unobserved genotypes. As expected, the oracle varies little. It is also interesting to look at the number of iterations needed before convergence of the EM algorithm, depending on the different scenarios.

A: n = 100 families; β = −0.6. B: n = 400 families; β = −0.6. C: n = 100 families; β = −1.2.

Fig 3 shows the violin plot including boxplots of the number of iterations required before the convergence of the EM algorithm according to scenarios (S0), (S1) and (S2) in the cases A and B. The violin plot shows that this number increases with the amount of unknown genotypes. Thus, when all genotypes are observed, the algorithm converges in 65 iterations on average while it converges in 89 iterations on average when no genotypes are observed. We also observe a higher variability of the number of iterations when the genotypes are not observed.

A: n = 100 families; β = −0.6. C: n = 100 families; β = −1.2.

Our method has been applied to two datasets of Portuguese families and French families already analyzed in [10]. All these analysed families had the ATTRV30M variant. TTRVal30Met (TTRV30M) was the only detected heterozygous mutation in the Portuguese kindreds. In contrast, 12 different TTR variants were detected in the French kindreds (42%), including TTRV30M in 58%. The clinical presentation of Portuguese patients was exclusively neurological and in French patients the phenotype was mainly mixed (neurologic and cardiac). For data analysis, the disease allele frequency was set to q = 0.04 as in [4]. The ascertainment bias was corrected by a classical method that consists in simply removing the phenotypic information of the proband, as done in [1]. The parent of origin gender parameter was not significantly different from zero in French families (β^=0.114 with p-value = 0.467). However, the β estimate was significantly different from zero in Portuguese families (β^=-0.999 and a p-value = 5.10−10) showing that the risk of being affected is higher when the mutation is transmitted by the mother than by the father. This difference in survival curves estimated is shown in Fig 4 and is totally consistent with the results of [10].

### Simulations study

We have simulated n families (n = 100 or 400) with 10 individuals as shown in Fig 1.

Genotypes were assigned respecting Mendelian transmission and heterozygous genotypes were ordered according to the gender (pat or mat) of the parent who transmitted the mutation. So, a new binary covariate (named POO for parent of origin) are introduced such as POO = mat if the mother transmitted the disease mutation and POO = pat if disease mutation was transmitted by the father. The disease allele frequency was set to q = 0.20 in our simulated dataset for the sake of speed (without simulating any ascertainment process). The age at event was simulated according to a piecewise constant hazard rate function. Moreover a censoring variable was added that follows a uniform distribution U[15,80].

To study the parent-of-origin effect, the POO was simulated in a proportional hazard model. The risk when the mutation is transmitted by the mother is given by λ0 as follows and the coefficient parameters in the hazard model is noted β (β = −0.60 or β = −1.2 in simulations). Thus, as β < 0, that’s mean that survival is upper when mutation provides on the father.
λ0(t)={0ift∈[0,20]0.02ift∈[20,40]0.10ift∈[40,60]0.05ift>60

To study the behavior of our method, we will consider different scenarios for simulations. In the first scenario (S0), all genotypes are set to unobserved, in the second more realistic scenario (S1), 80% of genotypes are observed in affected individuals and only 10% are observed in non affected individuals and in the third one (S2), all genotypes are observed. Moreover, we consider an other scenario where, in addition to the genotypes, we know, for each individual, the sex of the parent transmitting the mutation. We call this scenario “Oracle”. In order to compare the β estimate in the different scenarios, 200 replications are performed.

Fig 2 shows the violin plots according to the scenario for the percentage of observed genotypes, which allows to show the full distribution of the data. Thus, each color corresponds to a genotypic scenario for each of the cases A (n = 100 families was simulated with β = −0.6), B (n = 400 families was simulated with β = −0.6) and C (n = 100 families was simulated with β = −1.2). For each scenario, the mean and standard deviation are shown and the horizontal black line represents the true β value. In all cases the β estimate are unbiased whatever the scenario and the standard deviation increases with the percentage of unobserved genotypes. As expected, the oracle varies little. It is also interesting to look at the number of iterations needed before convergence of the EM algorithm, depending on the different scenarios.

A: n = 100 families; β = −0.6. B: n = 400 families; β = −0.6. C: n = 100 families; β = −1.2.

Fig 3 shows the violin plot including boxplots of the number of iterations required before the convergence of the EM algorithm according to scenarios (S0), (S1) and (S2) in the cases A and B. The violin plot shows that this number increases with the amount of unknown genotypes. Thus, when all genotypes are observed, the algorithm converges in 65 iterations on average while it converges in 89 iterations on average when no genotypes are observed. We also observe a higher variability of the number of iterations when the genotypes are not observed.

A: n = 100 families; β = −0.6. C: n = 100 families; β = −1.2.

### Application

Our method has been applied to two datasets of Portuguese families and French families already analyzed in [10]. All these analysed families had the ATTRV30M variant. TTRVal30Met (TTRV30M) was the only detected heterozygous mutation in the Portuguese kindreds. In contrast, 12 different TTR variants were detected in the French kindreds (42%), including TTRV30M in 58%. The clinical presentation of Portuguese patients was exclusively neurological and in French patients the phenotype was mainly mixed (neurologic and cardiac). For data analysis, the disease allele frequency was set to q = 0.04 as in [4]. The ascertainment bias was corrected by a classical method that consists in simply removing the phenotypic information of the proband, as done in [1]. The parent of origin gender parameter was not significantly different from zero in French families (β^=0.114 with p-value = 0.467). However, the β estimate was significantly different from zero in Portuguese families (β^=-0.999 and a p-value = 5.10−10) showing that the risk of being affected is higher when the mutation is transmitted by the mother than by the father. This difference in survival curves estimated is shown in Fig 4 and is totally consistent with the results of [10].

### Conclusion

In this article, we propose an extension of our method in order to take into account the gender of the transmitting parent in the estimation of the survival function from familial data in cases of age-dependent genetic diseases. The probability of the transmitting parent gender is calculated and its effect is estimated through a proportional hazard model. Our extension is assessed with simulations and applied on a French and a Portuguese dataset of ATTRv families. Obviously, the method provides confidence intervals for each estimation.

The simulation study shows that the method is unbiased and the increase in variance of the estimate with the number of unobserved genotypes.

We have applied this method to a French and Portuguese dataset and shows a significant difference of the penetrance function according to the gender of the transmitted parent. These results are consistent with the analyses already done in [10], even if the method used in [10] was less precise since the genotypic probabilities were not calculated. Results are also consistent with recent works in the Portuguese and Swedish ATTRV30M families showing that the disease’s risk and anticipation in the age of onset increase among offspring of affected mothers [14, 15]. The hypothesis of a role for mithochondrial DNA could explain this increased risk when the mutation is transmitted by the mother. The advantage here is that we have proposed an unified method able to take into account precisely these probabilities. Moreover, the method allows additional variables to be taken into account when estimating survival curves.

Finally, in this work, genotypic probabilities are calculated using a sum-product algorithm which is a very general method which can deal efficiently with very complex pedigree structure (ex: 2000 individuals with 50 loops). Unlike Elston-Stewart algorithm, the sum-product algorithm does not use loop breaking approaches to deal with loop pedigrees. Instead, the sum-product algorithm use an auxiliary tree called the junction tree (JT) which basically is a clique decomposition of the moral graph corresponding to the pedigree problem. JT and BP are well known is the graph theory (ex: JT can be used to solve a graph coloring problem) and in the mathematical field of probabilistic graphical models (Bayesian network, hidden Markov model, decision trees, Markov networks, etc.).



# SUPPLEMENTAL FILE 1: pone.0288958.pdf

# Preparing to download ...

[HHS Vulnerability Disclosure](https://www.hhs.gov/vulnerability-disclosure-policy/index.html)