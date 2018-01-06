---
title: "Fitting a Gamma-Gompertz survival model to capture-recapture data collected on free-ranging animal populations"
author: "Gilbert Marzolin and Olivier Gimenez"
date: "12 January 2018"
output: 
  html_document: 
    keep_md: yes
---

# In brief
  
`GammaGompertzCR` is an R package that allows estimating survival in free-ranging animal populations using a Gompertz capture-recapture model with a Gamma frailty to deal with individual heterogeneity. 
It uses data cloning in the Bayesian framework to get maximum likelihood parameter estimates (Lele et al. 2007). Data cloning uses multiple copies of the data to produce prior-invariant inferences and a normal distribution centered at the maximum likelihood estimates. In addition, this method allows detecting non-identifiable parameters (Lele et al. 2010).

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
## a      0.1180 0.01696 0.1696 0.0003096      0.0003163 0.9994
## b      0.4874 0.02824 0.2824 0.0005155      0.0005375 1.0011
## k      4.8995 0.03127 0.3127 0.0005709      0.0005718 1.0002
## psi    1.4014 0.03194 0.3194 0.0005831      0.0005657 1.0009
## sigeta 0.6722 0.02146 0.2146 0.0003919      0.0003691 1.0005
## 
## 2. Quantiles for each variable:
## 
##           2.5%    25%    50%    75%  97.5%
## a      0.08734 0.1059 0.1175 0.1293 0.1538
## b      0.43155 0.4686 0.4870 0.5065 0.5429
## k      4.83922 4.8784 4.9000 4.9203 4.9608
## psi    1.33744 1.3791 1.4019 1.4234 1.4607
## sigeta 0.63069 0.6572 0.6719 0.6868 0.7147
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

![](README2_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

## Convergence diagnostics

First have a look at the traceplots for all parameter estimates got with the last (ascending order) number of clones:

```r
plot(post_inf[[length(clo)]], ask = FALSE)
```

![](README2_files/figure-html/unnamed-chunk-7-1.png)<!-- -->![](README2_files/figure-html/unnamed-chunk-7-2.png)<!-- -->

Now checking some correlations:

```r
coda::autocorr.diag(post_inf[[length(clo)]])    
```

```
##                    a           b            k          psi        sigeta
## Lag 0    1.000000000  1.00000000  1.000000000  1.000000000  1.0000000000
## Lag 5    0.033694913  0.01362996  0.009879386 -0.032736556 -0.0031494722
## Lag 25   0.006880259  0.01077981 -0.016577855 -0.009879566  0.0001708971
## Lag 50  -0.003244604 -0.01037031 -0.031829974 -0.006790114  0.0302519217
## Lag 250 -0.011886721  0.01660218  0.007103807  0.021240631 -0.0032894622
```

```r
coda::crosscorr.plot(post_inf[[length(clo)]]) 
```

![](README2_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```r
coda::gelman.diag(post_inf[[length(clo)]])
```

```
## Potential scale reduction factors:
## 
##        Point est. Upper C.I.
## a           0.999          1
## b           1.001          1
## k           1.000          1
## psi         1.001          1
## sigeta      1.000          1
## 
## Multivariate psrf
## 
## 1
```

Parameters $a$ and $b$ seem to be slightly auto- and cross-correlated. 

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
dct <- dclone::dctable(post_inf[[1]],post_inf[[2]],post_inf[[3]],post_inf[[4]]) 
plot(dct)        
```

![](README2_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

```r
#log of scaled variances against number of clones.
plot(dct, type= "log.var")  
```

![](README2_files/figure-html/unnamed-chunk-10-2.png)<!-- -->

## Checking parameter identifiability

Data cloning allows testing identifiability rather easily (Lele 2010). When some parameters are not estimable, the largest eigenvalue of the posterior covariance matrix `lambda.max` does not tend to 0. When it tends to 0, the parameters are estimable and the results are independent of prior choice. Two other statistics `ms.error` and `r.squared` are used to test the normality of the limit. Notice that the statistical accuracy is a function of the sample size but not of the number of cloned copies (Lele et al.2010). 

Get the necessary elements and checking:


```r
u<-v<-w<-x<-y<-numeric(0)
for (i in 1:length(clo)){
  u[i]<-clo[i]
  v[i]<-dclone::dcdiag(post_inf[[i]])$lambda.max
  w[i]<-dclone::dcdiag(post_inf[[i]])$ms.error
  x[i]<-dclone::dcdiag(post_inf[[i]])$r.squared
  y[i]<-dclone::dcdiag(post_inf[[i]])$r.hat
}
mat.dcdiag <- matrix(c(u,v,w,x,y),c(length(clo),5)) 
dimnames(mat.dcdiag)= list(NULL,c("clones","lambda.max","ms.error","r.squared","r.hat"))               
mat.dcdiag
```

```
##      clones   lambda.max    ms.error    r.squared    r.hat
## [1,]      1 0.0009692336 0.009517570 0.0009526500 1.003052
## [2,]     10 0.0010381936 0.002541202 0.0002382847 1.002265
## [3,]     50 0.0010211120 0.006163615 0.0005751915 1.001143
## [4,]    100 0.0010262805 0.035804208 0.0026032228 1.001674
```

# References

Lele, SR, D Dennis, and F Lutscher. 2007. Data Cloning: Easy Maximum Likelihood Estimation for Complex Ecological Models Using Bayesian Markov Chain Monte Carlo Methods. Ecology Letters 10: 551–63.

Lele, SR, K Nadeem, and B Schmuland. 2010. Estimability and Likelihood Inference for Generalized Linear Mixed Models Using Data Cloning. Journal of the American Statistical Association 105: 1617–25.

Marzolin, G, A Charmantier, and O Gimenez. 2011. Frailty in State-Space Models: Application to Actuarial Senescence in the Dipper. Ecology 92: 562–67.

Missov, TI. 2013. Gamma-Gompertz Life Expectancy at Birth. Demographic Research 28: 259–270.

R Core Team. 2016. R: A Language and Environment for Statistical Computing. Vienna, Austria: R Foundation for Statistical Computing. http://www.R-project.org/.

Solymos, P. 2010. dclone: Data Cloning in R. The R Journal 2: 29-37.
