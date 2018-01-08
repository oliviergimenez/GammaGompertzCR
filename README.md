---
title: "Fitting a Gamma-Gompertz survival model to capture-recapture data collected on free-ranging animal populations"
author: "Gilbert Marzolin and Olivier Gimenez"
date: "12 January 2018"
output: 
  html_document: 
    keep_md: yes
---




# In brief
  
`GammaGompertzCR` is an R package that allows estimating survival in free-ranging animal populations using a Gompertz capture-recapture model with a Gamma frailty to deal with individual heterogeneity. It uses data cloning in the Bayesian framework to get maximum likelihood parameter estimates (Lele et al. 2007). Data cloning uses multiple copies of the data to produce prior-invariant inferences and a normal distribution centered at the maximum likelihood estimates. In addition, this method allows detecting non-identifiable parameters (Lele et al. 2010). To use the package, users should be familiar with Bayesian MCMC techniques and in particular how to interpret convergence diagnostics. We refer to Robert and Casella (2010) for an introduction to Bayesian MCMC techniques with R and to King et al. (2009) for an introduction for ecologists.

# Installation

This repository hosts the development version of the package. Assuming `JAGS` is already on your system (otherwise it'll take a minute to [get it](http://mcmc-jags.sourceforge.net/), you can install `GammaGompertzCR` as follows:

```r
if(!require(devtools)) install.packages("devtools")
library("devtools")
install_github('oliviergimenez/GammaGompertzCR')
```

# Getting help, reporting an issue or contribute

To report bugs/issues/feature requests, please file an [issue](https://github.com/oliviergimenez/GammaGompertzCR/issues). These are very welcome!

Feel free to contribute to the repository by forking and submitting a pull request. If you are new to GitHub, you might start with a [basic tutorial](https://help.github.com/articles/set-up-git) and check out [a more detailed guide to pull requests](https://help.github.com/articles/about-pull-requests/).

If you prefer to email, feel free to drop a message to the Gilbert and Olivier (see the DESCRIPTION file of this repository for their contact details).

# A session example

## Preliminary steps

First, we need to load the package:

```r
library(GammaGompertzCR)
```

```
## Thanks for using package GammaGompertzCR. 
##     Check out https://github.com/oliviergimenez/GammaGompertzCR for an example.
```

Then, we read in some data. Below we use a fragment of the dipper dataset `dipper_tiny` used by Marzolin et al. (2011) and provided with the package. The complete dataset `dipper_complete` is also provided. Detections are the 1's, non-detections the 2's.

```r
data(dipper_tiny)
head(dipper_tiny)
```

```
##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
## [1,]    1    1    1    1    1    1    1    1    2
## [2,]    1    1    1    1    1    1    1    2    2
## [3,]    1    1    1    1    1    1    2    2    2
## [4,]    1    1    1    1    1    2    2    2    2
## [5,]    1    1    1    2    2    2    2    2    2
## [6,]    1    1    1    2    2    2    2    2    2
```

## Model fitting

Now let's fit the Gamma-Gompertz model to these capture-recapture data:

```r
clo = c(1,10,50,100) # numbers of clones
nu = 1000 # number of updates
ni = 5000 # number of iterations
nt = 5 # thinning
nc = 3 # number of chains
initmeans = c(0.15,0.5,4.9,1.4,-0.4) # means of normally distributed priors for parameters
initprec = 1000 # common precision of the normal priors for the first number of clones

# fitting procedure (this may take a few minutes)
post_inf = fit_ggcr(dipper_tiny,clo,nu,ni,nt,nc,initmeans,initprec)
```

```
## 
## Parallel computation in progress
## 
## 
## Parallel computation in progress
## 
## 
## Parallel computation in progress
## 
## 
## Parallel computation in progress
```

Let's display the parameter estimates from posterior distributions :      

```r
summary(post_inf[[length(clo)]])
```

```
## 
## Iterations = 2005:7000
## Thinning interval = 5 
## Number of chains = 3 
## Sample size per chain = 1000 
## Number of clones = 100
## 
## 1. Empirical mean and standard deviation for each variable,
##    plus standard error of the mean:
## 
##          Mean      SD  DC SD  Naive SE Time-series SE  R hat
## a      0.1185 0.01723 0.1723 0.0003145      0.0003225 0.9995
## b      0.4868 0.02826 0.2826 0.0005160      0.0005081 1.0010
## k      4.8995 0.03179 0.3179 0.0005805      0.0005688 1.0013
## psi    1.4014 0.03241 0.3241 0.0005917      0.0005918 0.9996
## sigeta 0.6724 0.02091 0.2091 0.0003818      0.0003763 1.0022
## 
## 2. Quantiles for each variable:
## 
##           2.5%    25%    50%    75%  97.5%
## a      0.08706 0.1067 0.1174 0.1297 0.1542
## b      0.43182 0.4675 0.4863 0.5062 0.5424
## k      4.83570 4.8774 4.8999 4.9214 4.9612
## psi    1.33556 1.3805 1.4011 1.4230 1.4648
## sigeta 0.63374 0.6579 0.6718 0.6869 0.7139
```

Let's represent survival as a function of age at the population level. 


```r
# get values from posterior distributions
a<- coef(post_inf[[length(clo)]])[1]
b<- coef(post_inf[[length(clo)]])[2] 
k<- coef(post_inf[[length(clo)]])[3]
grid_age = seq(0,8,1)
S = exp(-k*log(1+(exp(b*grid_age)-1)*a/(b*k)))
plot(grid_age,S,xlab='age',ylab='estimated survival',lwd=2,type='l')
```

![](README_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

## Convergence diagnostics

First have a look at the traceplots for all parameter estimates that we obtained with the last (ascending order) number of clones:

```r
plot(post_inf[[length(clo)]], ask = TRUE)
```

![](README_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

![](README_files/figure-html/unnamed-chunk-8-2.png)<!-- -->

Now checking some correlations:

```r
autocorr.diag(post_inf[[length(clo)]])    
```

```
##                     a            b            k          psi      sigeta
## Lag 0    1.0000000000  1.000000000  1.000000000  1.000000000 1.000000000
## Lag 5    0.0344650219 -0.027349395 -0.015760295  0.001375407 0.001969149
## Lag 25   0.0097737538  0.004886749  0.021157276 -0.004250552 0.011828727
## Lag 50  -0.0004097896 -0.029179821 -0.013871670  0.010121597 0.004734142
## Lag 250  0.0022794338  0.007089348  0.001216246 -0.005261328 0.018063070
```

```r
crosscorr.plot(post_inf[[length(clo)]]) 
```

![](README_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
gelman.diag(post_inf[[length(clo)]])
```

```
## Potential scale reduction factors:
## 
##        Point est. Upper C.I.
## a               1       1.00
## b               1       1.00
## k               1       1.01
## psi             1       1.00
## sigeta          1       1.01
## 
## Multivariate psrf
## 
## 1
```

Parameters a and b seem to be slightly auto- and cross-correlated. 

## Sensitivity to priors

It is good advice to check if changing the priors leads to similar posterior results. In our case, the choice of `initmeans` is made through trial and error, because the models we're dealing with have likelihood with multiple spikes. For instance, we might begin with vague priors: `initmeans = c(1,1,5,2,-0.5), initprec=10, clo= c(1,50,100)`, which yields `c(0.11,0.51,4.39,1.44,0.41)` with quite large standard deviation. Modifying the prior means can affect the posterior estimates. Hence, additional tests are required. For each number of clones, the posterior variances of the parameters are divided by that obtained with `min(clo)` to get scaled variances. When increasing the number of clones, the means of the parameters converge to the maximum likelihood estimates while their scaled variances decrease to reach a lower bound. We get a visual threshold indicating the minimal number of clones to agree with
the desired results (Solymos 2010) by plotting the scaled variances or their logs against the number of clones.

First, create a list by hand or via a few lines of code when there are several lots of clones: `post_inf[[1]], ..., post_inf[[length(clo)]]`. Say we have 4 lots of clones:

```r
clo = c(1,10,50,100) 
lu = NULL   
for (i in 1:length(clo)) lu = c(lu,paste(c("post_inf[[",i,"]],"),sep="",collapse=""))
lu
```

```
## [1] "post_inf[[1]]," "post_inf[[2]]," "post_inf[[3]]," "post_inf[[4]],"
```

Second, use function `dct` to do the job:

```r
dct <- dctable(post_inf[[1]],post_inf[[2]],post_inf[[3]],post_inf[[4]]) 
plot(dct)        
```

![](README_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

```r
#log of scaled variances against number of clones.
plot(dct, type= "log.var")  
```

![](README_files/figure-html/unnamed-chunk-11-2.png)<!-- -->

## Checking parameter identifiability

Data cloning allows testing identifiability rather easily (Lele 2010). When some parameters are not estimable, the largest eigenvalue of the posterior covariance matrix `lambda.max` does not tend to 0. When it tends to 0, the parameters are estimable and the results are independent of prior choice. Two other statistics `ms.error` and `r.squared` are used to test the normality of the limit. Notice that the statistical accuracy is a function of the sample size but not of the number of cloned copies (Lele et al.2010). 

Get the necessary elements and checking:


```r
u<-v<-w<-x<-y<-numeric(0)
for (i in 1:length(clo)){
  u[i]<-clo[i]
  v[i]<-dcdiag(post_inf[[i]])$lambda.max
  w[i]<-dcdiag(post_inf[[i]])$ms.error
  x[i]<-dcdiag(post_inf[[i]])$r.squared
  y[i]<-dcdiag(post_inf[[i]])$r.hat
}
mat.dcdiag <- matrix(c(u,v,w,x,y),c(length(clo),5)) 
dimnames(mat.dcdiag)= list(NULL,c("clones","lambda.max","ms.error","r.squared","r.hat"))               
mat.dcdiag
```

```
##      clones   lambda.max    ms.error    r.squared    r.hat
## [1,]      1 0.0010388885 0.013733190 0.0009974044 1.001534
## [2,]     10 0.0009929232 0.008450080 0.0007712377 1.004060
## [3,]     50 0.0010165349 0.013395569 0.0010058643 1.001171
## [4,]    100 0.0010561432 0.007795492 0.0006619290 1.003630
```

# References

King R, Morgan BJT, Gimenez O and Brooks S (2009). Bayesian Analysis for Population Ecology. Chapman & Hall/CRC Interdisciplinary Statistics.

Lele, SR, D Dennis, and F Lutscher. 2007. Data Cloning: Easy Maximum Likelihood Estimation for Complex Ecological Models Using Bayesian Markov Chain Monte Carlo Methods. Ecology Letters 10: 551–63.

Lele, SR, K Nadeem, and B Schmuland. 2010. Estimability and Likelihood Inference for Generalized Linear Mixed Models Using Data Cloning. Journal of the American Statistical Association 105: 1617–25.

Marzolin, G, A Charmantier, and O Gimenez. 2011. Frailty in State-Space Models: Application to Actuarial Senescence in the Dipper. Ecology 92: 562–67.

Missov, TI. 2013. Gamma-Gompertz Life Expectancy at Birth. Demographic Research 28: 259–270.

R Core Team. 2016. R: A Language and Environment for Statistical Computing. Vienna, Austria: R Foundation for Statistical Computing. http://www.R-project.org/.

Robert CP and Casella G (2010). Introducing Monte Carlo Methods with R. Springer.

Solymos, P. 2010. dclone: Data Cloning in R. The R Journal 2: 29-37.


