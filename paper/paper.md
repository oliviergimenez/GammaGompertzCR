---
title: Fitting a Gamma-Gompertz survival model to capture-recapture data collected
  on free-ranging animal populations
authors:
- affiliation: 1
  name: Gilbert Marzolin
  orcid: 0000-0002-7819-5739
- affiliation: 1
  name: Olivier Gimenez
  orcid: 0000-0001-7001-5142
date: "12 January 2018"
output:
  html_document:
    fig_caption: yes
    keep_md: yes
  pdf_document: default
bibliography: paper.bib
tags:
- gamma-gompertz
- survival
- capture-recapture
- MCMC
- data cloning
- dipper
affiliations:
- index: 1
  name: Centre National de la Recherche Scientifique
---

# Introduction
  
`GammaGompertzCR` is an R package [@R] that allows estimating survival in free-ranging animal populations using a Gompertz capture-recapture model with a Gamma frailty to deal with individual heterogeneity [@MarzolinCharmantierGimenez2011; @Missov2013]. To use the package, users should be familiar with Bayesian MCMC techniques and in particular how to interpret convergence diagnostics. We refer to @robert2010introducing for an introduction to Bayesian MCMC techniques with R and to @king2009bayesian for an introduction for ecologists. In this paper, we introduce the theory underlying the model we implement in `GammaGompertzCR`, and illustrate the approach using a real example. 

# Theory

## The Gamma-Gompertz model

We consider a random variable $T \geq 0$ called time to event. When the event is death of some organism, $T$ is usually associated with a survivor function $S$ such as, for $t\geq 0$, $S(t) = P(T>t)$. Denoting $f$ the pdf of $T$, we have $S(t)= \int_t^{+ \infty }{f(t)dt}$ or $f(t)=-dS(t)/dt$. The hazard function or mortality rate $h$ yields, for $t \geq 0$, the rate of death $h(t)$ given the animal survived up to time $t$, that is $h(t) = f(t)/S(t)$ or $h(t) = -d\log{S(t)}/dt$. In case of discrete time steps - for instance age $j$ in years - annual survival $s$ by age is obtained as $s(j,j+1) = P(T > j+1|T>j) = S(j+1)/S(j)$. 

To deal with individual heterogeneity in survival, we can use individual random-effect models [@MarzolinCharmantierGimenez2011] or multiply the baseline mortality rate $h$ by a unit-specific random variable $u$ named frailty. When $h$ is Gompertz with $h(t) = a e^{b t}$ and $u$ is distributed as a $\Gamma(k,\lambda)$ (with mean $k/\lambda$), we have a multiplicative Gamma-Gompertz frailty model. Typically a $\Gamma$ with mean one is adopted, hence $k = \lambda$ and variance is $1/k$. In a Gamma-Gompertz model, the population survival function is obtained through marginalization. When two parameters are used, $S(t)= (1+(a/b\lambda)(e^{bt}-1))^{-k}$ [@Missov2013], hence the hazard equals $h(t)=a(k/\lambda)e^{bt}/(1+a(e^{bt}-1)/b\lambda)$. In our case, we use $\lambda = k$. For the detection probability, we used a yearly random effect normally distributed with mean $\psi$ and standard deviation $\sigma_{\eta}$ on the logit scale.

## Parameter estimation

To get maximum likelihood parameter estimates, we use data cloning in the Bayesian framework [@LeleDennisLutscher2007] implemented in the R package `dclone` [@Solymos2010]. Data cloning uses multiple copies of the data to produce prior-invariant inferences and converge towards a normal distribution centered at the maximum likelihood estimates. In addition, this method allows detecting non-identifiable parameters [@LeleNadeemSchmuland2010]. 

# Illustration

We now illustrate the estimation of Gamma-Gompertz model parameters using data cloning. We analysed capture-recapture data collected on 75 breeding females over 9 years in a Dipper (*Cinclus cinclus*) population [@MarzolinCharmantierGimenez2011]. These data are just a subset of the complete dataset (>1000 individuals, 35 years) that is provided with the package. We used 1, 10, 50 and 100 clones with 3 parallel MCMC chains for a total of 5000 iterates with a burn-in period of 1000. The model fitting took less than 5 minutes.

Parameter estimates are provided in Table 1. 

| Parameter  | Mean | SD | Credible interval |
| :---------: | :---: | :-: | :----------------: |
| $a$ | 0.12 | 0.02 | [0.09; 0.15] |
| $b$ | 0.48 | 0.03 | [0.43; 0.54] |
| $k$ | 4.90 | 0.03 | [4.84; 4.96] |
| $\psi$ | 1.40 | 0.03 | [1.34; 1.46] |
| $\sigma_{\eta}$ | 0.67 | 0.02 | [0.63; 0.71] |

Table 1: Parameter estimates of the Gamma-Gompertz model using the Dipper dataset. We refer to Theory section for more details about these parameters. We provide the posterior mean and standard deviation as well as the $95\%$ credible interval.

We also provide a plot of the relationship between estimated survival and age in Figure \ref{fig:survage}.

![The relationship between Dipper survival and age as estimated by the Gamma-Gompertz model.\label{fig:survage}](fig_agesurvival.png)

More details on how to use the package, including conducting convergence diagnostics, performing a sensitivity analysis and checking parameter identifiability, are provided at the package development repository (see next section). 

# Availability

The latest version of the package is available at 'https://github.com/oliviergimenez/GammaGompertzCR'.

# References
