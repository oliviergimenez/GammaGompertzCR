---
title: "Fitting a Gamma-Gompertz survival model to capture-recapture data collected on free-ranging animal populations"
author: "Gilbert Marzolin and Olivier Gimenez"
date: "25 février 2017"
output: html_document
---

# In brief
  
`GammaGompertzCR` is an R package that allows estimating survival in free-ranging animal populations using a Gompertz capture-recapture model with a Gamma frailty to deal with individual heterogeneity. 
It uses data cloning in the Bayesian framework to get maximum likelihood parameter estimates (Lele et al. 2007). Data cloning uses multiple copies of the data to produce prior-invariant inferences and a normal distribution centered at the maximum likelihood estimates. In addition, this method allows detecting non-identifiable parameters (Lele et al. 2010).

This repository hosts the development version of the package. Assuming JAGS is already on your system (otherwise: <http://mcmc-jags.sourceforge.net/>), you can install `GammaGompertzCR` as follows:
```{r}
if(!require(devtools)) install.packages("devtools")
library("devtools")
install_github('oliviergimenez/GammaGompertzCR')
```

# A session example

First, we need to load the package:
```{r}
library(GammaGompertzCR)
```

Then, we read in some data. Below we use a fragment of the dipper dataset used by Marzolin et al. (2011). Detections are the 1's, non-detections the 2's.
```{r}
mydat <- matrix(c(1,1,1,1,1,1,1,1,2,
1,1,1,1,1,1,1,2,2,
1,1,1,1,1,1,2,2,2,
1,1,1,1,1,2,2,2,2,
1,1,1,2,2,2,2,2,2,
1,1,1,2,2,2,2,2,2,
1,1,1,2,2,2,2,2,2,
1,1,2,1,2,2,2,2,2,
1,1,2,2,2,2,2,2,2,
1,1,2,2,2,2,2,2,2,
1,1,2,2,2,2,2,2,2,
1,2,1,1,1,1,1,2,2,
1,2,1,1,1,2,2,2,2,
1,2,1,1,2,1,1,2,2,
1,2,1,1,2,2,2,2,2,
1,2,1,2,2,2,2,2,2,
1,2,1,2,2,2,2,2,2,
1,2,1,2,2,2,2,2,2,
1,2,1,2,2,2,2,2,2,
1,2,1,2,2,2,2,2,2,
1,2,1,2,2,2,2,2,2,
1,2,1,2,2,2,2,2,2,
1,2,1,2,2,2,2,2,2,
1,2,1,2,2,2,2,2,2,
2,1,1,1,1,1,1,1,1,
2,1,1,1,1,1,1,2,2,
2,1,1,1,1,1,1,2,2,
2,1,1,1,1,2,2,2,2,
2,1,1,1,1,2,2,2,2,
2,1,1,1,1,2,2,2,2,
2,1,1,1,2,2,2,2,2,
2,1,1,1,2,2,2,2,2,
2,1,1,1,2,2,2,2,2,
2,1,1,1,2,2,2,2,2,
2,1,1,1,2,2,2,2,2,
2,1,1,1,2,2,2,2,2,
2,1,1,2,2,2,2,2,2,
2,1,1,2,2,2,2,2,2,
2,1,1,2,2,2,2,2,2,
2,1,1,2,2,2,2,2,2,
2,1,1,2,2,2,2,2,2,
2,1,1,2,2,2,2,2,2,
2,1,1,2,2,2,2,2,2,
2,1,1,2,2,2,2,2,2,
2,1,1,2,2,2,2,2,2,
2,1,1,2,2,2,2,2,2,
2,1,1,2,2,2,2,2,2,
2,1,1,2,2,2,2,2,2,
2,1,1,2,2,2,2,2,2,
2,1,1,2,2,2,2,2,2,
2,1,1,2,2,2,2,2,2,
2,1,2,1,1,1,1,1,1,
2,1,2,1,1,1,1,2,2,
2,1,2,1,2,2,2,2,2,
2,1,2,1,2,2,2,2,2,
2,1,2,1,2,2,2,2,2,
2,1,2,1,2,2,2,2,2,
2,1,2,1,2,2,2,2,2,
2,1,2,1,2,2,2,2,2,
2,1,2,1,2,2,2,2,2,
2,1,2,1,2,2,2,2,2,
2,1,2,1,2,2,2,2,2,
2,1,2,1,2,2,2,2,2,
2,1,2,1,2,2,2,2,2,
2,1,2,2,1,1,2,2,2,
2,2,1,1,1,1,1,1,2,
2,2,1,1,1,1,1,2,2,
2,2,1,1,1,1,1,2,2,
2,2,1,1,1,1,2,2,2,
2,2,1,1,1,1,2,2,2,
2,2,1,1,1,1,2,2,2,
2,2,1,1,1,2,2,2,2,
2,2,1,1,1,2,2,2,2,
2,2,1,1,2,2,2,2,2,
2,2,1,1,2,2,2,2,2),ncol=9,byrow=T)
head(mydat)
```

Fit Gamma-Gompertz model to the capture-recapture data:
```{r}
clo = c(2,5) # number of clones
nu = 1000 # number of updates
ni = 5000 # number of iterations
nt = 50 # thinning
nc = 2 # number of chains
# fitting procedure
post_inf = fit_ggcr(mydat,clo,nu,ni,nt,nc)
```

Let's display the parameter estimates:      
```{r}
summary(post_inf)
```

Get traceplots and posterior distributions for all parameters:
```{r}
plot(post_inf)
```

Let's represent survival as a function of age at the population level:
```{r}
# Get values from posterior distributions, chain 1
a = mean(post_inf[[1]][,'a'])
b = mean(post_inf[[1]][,'b'])
k = mean(post_inf[[1]][,'k'])
grid_age = seq(0,8,1)
S = exp(-k*log(1+(exp(b*grid_age)-1)*a/(b*k)))
plot(grid_age,S,xlab='age',ylab='estimated survival',lwd=2,type='l')
```

Data cloning allows testing identifiability rather easily (Lele 2010). The convergence of the data cloning algorithm is achieved when the largest eigenvalue of the post variance matrix tends to 0 at the approximate rate of 1 / `nclones`: hence we get a threshold indicating the minimal number of clones to agree with the desired results. When some parameters are not estimable, then the largest eigenvalue does not tend to 0. Two other statistics `ms.error` and `r.squared` are used to test the normality of the limit.
  
```{r}
dct <- dctable(post_inf)
plot(dct)
plot(dct, type= "log.var")
dcdiag(post_inf)
```

# References

Lele, SR, D Dennis, and F Lutscher. 2007. “Data Cloning: Easy Maximum Likelihood Estimation for Complex Ecological Models Using Bayesian Markov Chain Monte Carlo Methods.” Ecology Letters 10: 551–63.

Lele, SR, K Nadeem, and B Schmuland. 2010. “Estimability and Likelihood Inference for Generalized Linear Mixed Models Using Data Cloning.” Journal of the American Statistical Association 105: 1617–25.

Marzolin, G, A Charmantier, and O Gimenez. 2011. “Frailty in State-Space Models: Application to Actuarial Senescence in the Dipper.” Ecology 92: 562–67.

Missov, TI. 2013. “Gamma-Gompertz Life Expectancy at Birth.” Demographic Research 28: 259–70.

Solymos, P. 2010. “dclone: Data Cloning in R.” The R Journal 2: 29-37.
