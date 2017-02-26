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
  
`GammaGompertzCR` is an R package (@R) that allows estimating survival in free-ranging animal populations using a Gompertz capture-recapture model with a Gamma frailty to deal with individual heterogeneity (@MarzolinCharmantierGimenez2011, @Missov2013). It uses data cloning in the Bayesian framework to get maximum likelihood parameter estimates (@LeleDennisLutscher2007). Data cloning uses multiple copies of the data to produce prior-invariant inferences and a normal distribution centered at the maximum likelihood estimates. In addition, this method allows detecting non-identifiable parameters (@LeleNadeemSchmuland2010). The latest version of the package is available at 'https://github.com/oliviergimenez/GammaGompertzCR'.

# References
