---
title: 'Fitting a Gamma-Gompertz survival model to capture-recapture data collected on free-ranging animal populations'
bibliography: paper.bib
date: "24 February 2017"
tags:
- gamma-gompertz
- survival
- capture-recapture
- MCMC
- data cloning
- dipper
authors:
- affiliation: 1
  name: Gilbert Marzolin
  orcid: 0000-0002-7819-5739
- affiliation: 1
  name: Olivier Gimenez
  orcid: 0000-0001-7001-5142
affiliations:
- name: Centre National de la Recherche Scientifique
  index: 1

---

# Summary
  
`GammaGompertzCR` is an R package (@R) that allows estimating survival in free-ranging animal populations using a Gompertz capture-recapture model with a Gamma frailty to deal with individual heterogeneity (@MarzolinCharmantierGimenez2011, @Missov2013). 

In details, we consider a random variable $T \geq 0$ called time to event. When the event is death of some organism, $T$ is usually associated with a survivor function $S$ such as, for $t\geq 0$, $S(t) = P(T>t)$. Denoting $f$ the pdf of $T$, we have $S(t)= \int_t^{+ \infty }{f(t)dt}$ or $f(t)=-dS(t)/dt$. The hazard function or mortality rate $h$ yields, for $t \geq 0$, the rate of death $h(t)$ given the animal survived up to time $t$, that is $h(t) = f(t)/S(t)$ or $h(t) = -d\log{S(t)}/dt$. In case of discrete time steps - for instance age $j$ in years - annual survival $s$ by age is obtained as $s(j,j+1) = P(T > j+1|T>j) = S(j+1)/S(j)$. 

To deal with individual heterogeneity in survival, we use individual random-effect models (Marzolin et al. 2011) and multiply the baseline mortality rate $h$ by a unit-specific random variable $u$ named frailty. When $h$ is Gompertz with $h(t) = a e^{b t}$ and $u$ is distribution as a $\Gamma(k,\lambda)$ (with mean $k/\lambda$), we have a multiplicative Gamma-Gompertz frailty model. Typically a $\Gamma$ with mean one is adopted, hence $k = \lambda$ and variance is $1/k$. In a Gamma-Gompertz model, the population survival function is obtained through marginalization. When two parameters are used, $S(t)= (1+(a/b\lambda)(e^{bt}-1))^{-k}$ (Missov 2013), hence the hazard equals $h(t)=a(k/\lambda)e^{bt}/(1+a(e^{bt}-1)/b\lambda)$. In our case, we use $\lambda = k$.

To get maximum likelihood parameter estimates, we use data cloning in the Bayesian framework (@LeleDennisLutscher2007). Data cloning uses multiple copies of the data to produce prior-invariant inferences and a normal distribution centered at the maximum likelihood estimates. In addition, this method allows detecting non-identifiable parameters (@LeleNadeemSchmuland2010). 

The latest version of the package is available at 'https://github.com/oliviergimenez/GammaGompertzCR'.

# References
