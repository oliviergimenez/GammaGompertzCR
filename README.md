---
title: "Fitting a Gamma-Gompertz survival model to capture-recapture data collected on free-ranging animal populations"
author: "Gilbert Marzolin and Olivier Gimenez"
date: "25 février 2017"
output: html_document
---

# In a nutshell
  
`GammaGompertzCR` is an R package that allows estimating survival in free-ranging animal populations using a Gompertz capture-recapture model with a Gamma frailty to deal with individual heterogeneity. 
It uses data cloning in the Bayesian framework to get maximum likelihood parameter estimates (Lele et al. 2007). Data cloning uses multiple copies of the data to produce prior-invariant inferences and a normal distribution centered at the maximum likelihood estimates. In addition, this method allows detecting non-identifiable parameters (Lele et al. 2010). 

This repository hosts the development version of the package. Assuming JAGS is already on your system (otherwise: <http://mcmc-jags.sourceforge.net/>), you can install `GammaGompertzCR` as follows:
```{r}
if(!require(devtools)) install.packages("devtools")
library("devtools")
install_github('oliviergimenez/GammaGompertzCR')
```

# Model

We consider a random variable $T \geq 0$ called time to event. When the event is death of some organism, $T$ is usually associated with a survivor function $S$ such as, for $t\geq 0$, $S(t) = P(T>t)$. Denoting $f$ the pdf of $T$, we have $S(t)= \int_t^{+ \infty }{f(t)dt}$ or $f(t)=-dS(t)/dt$. The hazard function or mortality rate $h$ yields, for $t \geq 0$, the rate of death $h(t)$ given the animal survived up to time $t$, that is $h(t) = f(t)/S(t)$ or $h(t) = -d\log{S(t)}/dt$. In case of discrete time steps - for instance age $j$ in years - annual survival $s$ by age is obtained as $s(j,j+1) = P(T > j+1|T>j) = S(j+1)/S(j)$.

When studying individual heterogeneity in survival in a population, we can use individual random-effect models (Marzolin et al. 2011) and multiply the baseline mortality rate $h$ by a unit-specific random variable $u$ named frailty. When $h: t \rightarrow a e^{b t}$ is a Gompertz function and $u$ is a $\Gamma(k,\lambda)$ distribution (with mean $k/\lambda$), we have a multiplicative Gamma-Gompertz frailty model. Typically a $\Gamma$ with mean one is adopted, hence $k = \lambda$ and variance is $1/k$. In a Gamma-Gompertz model, the population survival function is obtained through marginalization. When two parameters are used, $S(t)= (1+(a/b\lambda)(e^{bt}-1))^{-k}$ (Missov 2013), hence the hazard equals $h(t)=a(k/\lambda)e^{bt}/(1+a(e^{bt}-1)/b\lambda)$. In our case, we use $\lambda = k$.

# A session example

First, we need to load the package:
```{r}
library(GammaGompertzCR)
```

Then, we read in some data. Below we use a fragment of the dipper dataset used by Marzolin et al. (2011). Detections are the 2's, non-detection the 1's.
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

Let's represent survival as a function of age at the population level: TO BE DONE
```{r}
post_inf[[1]]
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
