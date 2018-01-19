[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1154931.svg)](https://doi.org/10.5281/zenodo.1154931)

[![DOI](http://joss.theoj.org/papers/10.21105/joss.00216/status.svg)](https://doi.org/10.21105/joss.00216)

# Fitting a Gamma-Gompertz survival model to capture-recapture data collected on free-ranging animal populations

## In brief
  
`GammaGompertzCR` is an R package that allows estimating survival in free-ranging animal populations using a Gompertz capture-recapture model with a Gamma frailty to deal with individual heterogeneity. It uses data cloning in the Bayesian framework to get maximum likelihood parameter estimates (Lele et al. 2007). Data cloning uses multiple copies of the data to produce prior-invariant inferences and a normal distribution centered at the maximum likelihood estimates. In addition, this method allows detecting non-identifiable parameters (Lele et al. 2010). To use the package, users should be familiar with Bayesian MCMC techniques and in particular how to interpret convergence diagnostics. We refer to Robert and Casella (2010) for an introduction to Bayesian MCMC techniques with R and to King et al. (2009) for an introduction for ecologists.

## Installation

This repository hosts the development version of the package. Assuming `JAGS` is already on your system (otherwise it'll take a minute to [get it](http://mcmc-jags.sourceforge.net/), you can install `GammaGompertzCR` as follows:

```r
if(!require(devtools)) install.packages("devtools")
library("devtools")
install_github('oliviergimenez/GammaGompertzCR')
```

## Getting help, reporting an issue or contribute

To report bugs/issues/feature requests, please file an [issue](https://github.com/oliviergimenez/GammaGompertzCR/issues). These are very welcome!

Feel free to contribute to the repository by forking and submitting a pull request. If you are new to GitHub, you might start with a [basic tutorial](https://help.github.com/articles/set-up-git) and check out [a more detailed guide to pull requests](https://help.github.com/articles/about-pull-requests/).

If you prefer to email, feel free to drop a message to the Gilbert and Olivier (see the DESCRIPTION file of this repository for their contact details).

## A session example

### Preliminary steps

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

### Model fitting

Now let's fit the Gamma-Gompertz model to these capture-recapture data:

```r
clo <- c(1,10,50,100) # numbers of clones
nu <- 1000 # number of updates
ni <- 5000 # number of iterations
nt <- 5 # thinning
nc <- 3 # number of chains
initmeans <- c(0.15,0.5,4.9,1.4,-0.4) # means of normally distributed priors for parameters
initprec <- 1000 # common precision of the normal priors for the first number of clones

# fitting procedure (this may take a few minutes)
post_inf <- fit_ggcr(dipper_tiny,clo,nu,ni,nt,nc,initmeans,initprec)
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
## a      0.1187 0.01721 0.1721 0.0003142      0.0003184 0.9996
## b      0.4858 0.02840 0.2840 0.0005184      0.0005186 0.9998
## k      4.8999 0.03184 0.3184 0.0005814      0.0005729 1.0010
## psi    1.4031 0.03177 0.3177 0.0005801      0.0005802 1.0020
## sigeta 0.6723 0.02120 0.2120 0.0003871      0.0003871 0.9999
## 
## 2. Quantiles for each variable:
## 
##           2.5%    25%    50%    75%  97.5%
## a      0.08708 0.1068 0.1181 0.1298 0.1546
## b      0.43072 0.4671 0.4856 0.5046 0.5423
## k      4.83788 4.8780 4.9004 4.9218 4.9615
## psi    1.34099 1.3821 1.4027 1.4248 1.4677
## sigeta 0.63147 0.6576 0.6718 0.6864 0.7150
```

Let's represent survival as a function of age at the population level. 


```r
# get values from posterior distributions
a <- coef(post_inf[[length(clo)]])[1]
b <- coef(post_inf[[length(clo)]])[2] 
k <- coef(post_inf[[length(clo)]])[3]
grid_age <- seq(0,8,1)
S <- exp(-k*log(1+(exp(b*grid_age)-1)*a/(b*k)))
plot(grid_age,S,xlab='age',ylab='estimated survival',lwd=2,type='l')
```

![](README_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

### How to choose the precision, number of iterations and clones?

For the Gamma-Gompertz model, all parameters are positive. The user begins with a quite moderate precision, say 10, for the normal priors. This precision is only valid during the run with the first element of vector clone. For the second element of clone and the following ones, mean and precision of a parameter prior are chosen as being the posterior mean and precision of this same parameter obtained at the previous run (function `update`). Several chains run in parallel may converge to different values according to their starting points, while the standard deviation (sd) of the parameter estimates may be quite large. Then, one increases the new prior initial precision `initprec` from 100 or 1000, associated with new `initmeans` picked in the neighborhood of the previous mean estimates, and start the iterations again. To choose among the different values, we look at the corresponding tests `lambda.max`, `ms.error` and `r.squared`. In the above example, with `prec = 10`, `lamba.max` varies around 0.12 according to the number of clones, while with `prec = 100` or `prec = 1000`, it decreases to 0.002, which indicates a more relevant choice. Again, the statistical accuracy is a function of the sample size but not of the number of cloned copies (Lele et al. 2010). Regarding the number of iterations, `ni`, the plots of the different chains usually help. For each parameter, a display by chain of the iterations can be obtained with `densityplot(post_inf[[1]][,"a"])`, the `'[[1]]'` being for the first element of vector clone. 

### Convergence diagnostics

First have a look at the traceplots for all parameter estimates that we obtained with the last (ascending order) number of clones:

```r
plot(post_inf[[length(clo)]], ask = TRUE)
```

![](README_files/figure-html/unnamed-chunk-7-1.png)<!-- -->![](README_files/figure-html/unnamed-chunk-7-2.png)<!-- -->

Now checking some correlations:

```r
library(coda)
autocorr.diag(post_inf[[length(clo)]])    
```

```
##                   a           b            k           psi       sigeta
## Lag 0   1.000000000 1.000000000  1.000000000  1.0000000000  1.000000000
## Lag 5   0.004885454 0.007062062 -0.007529550  0.0003981668  0.020619596
## Lag 25  0.028170386 0.039102982  0.007367834 -0.0093045888  0.018273569
## Lag 50  0.003122581 0.003448430 -0.009122560  0.0285703985 -0.004180051
## Lag 250 0.013335939 0.012153872 -0.006019243 -0.0068175773 -0.025717130
```

```r
crosscorr.plot(post_inf[[length(clo)]]) 
```

![](README_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```r
gelman.diag(post_inf[[length(clo)]])
```

```
## Potential scale reduction factors:
## 
##        Point est. Upper C.I.
## a               1       1.00
## b               1       1.00
## k               1       1.00
## psi             1       1.01
## sigeta          1       1.00
## 
## Multivariate psrf
## 
## 1
```

Parameters a and b seem to be slightly auto- and cross-correlated. 

### Sensitivity to priors

It is good advice to check if changing the priors leads to similar posterior results. In our case, the choice of `initmeans` is made through trial and error, because the models we're dealing with have likelihood with multiple spikes. For instance, we might begin with vague priors: `initmeans = c(1,1,5,2,-0.5), initprec=10, clo= c(1,50,100)`, which yields `c(0.11,0.51,4.39,1.44,0.41)` with quite large standard deviation. Modifying the prior means can affect the posterior estimates. Hence, additional tests are required. For each number of clones, the posterior variances of the parameters are divided by that obtained with `min(clo)` to get scaled variances. When increasing the number of clones, the means of the parameters converge to the maximum likelihood estimates while their scaled variances decrease to reach a lower bound. We get a visual threshold indicating the minimal number of clones to agree with
the desired results (Solymos 2010) by plotting the scaled variances or their logs against the number of clones.

Say we have 4 lots of clones:

```r
clo <- c(1,10,50,100) 
```

We use function `dctable` to do the job:

```r
library(dclone)
```

```
## Loading required package: parallel
```

```
## Loading required package: Matrix
```

```
## dclone 2.1-2 	 2016-03-11
```

```r
dct <- do.call(dctable,post_inf) 
plot(dct)        
```

![](README_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

```r
#log of scaled variances against number of clones.
plot(dct, type= "log.var")  
```

![](README_files/figure-html/unnamed-chunk-10-2.png)<!-- -->

### Checking parameter identifiability

Data cloning allows testing identifiability rather easily (Lele 2010). When some parameters are not estimable, the largest eigenvalue of the posterior covariance matrix `lambda.max` does not tend to 0. When it tends to 0, the parameters are estimable and the results are independent of prior choice. Two other statistics `ms.error` and `r.squared` are used to test the normality of the limit. Notice that the statistical accuracy is a function of the sample size but not of the number of cloned copies (Lele et al.2010). 

Get the necessary elements and checking:


```r
u<-v<-w<-x<-y<-numeric(0)
for (i in 1:length(clo)){
  u[i] <- clo[i]
  v[i] <- dcdiag(post_inf[[i]])$lambda.max
  w[i] <- dcdiag(post_inf[[i]])$ms.error
  x[i] <- dcdiag(post_inf[[i]])$r.squared
  y[i] <- dcdiag(post_inf[[i]])$r.hat
}
mat.dcdiag <- matrix(c(u,v,w,x,y),c(length(clo),5)) 
dimnames(mat.dcdiag) <- list(NULL,c("clones","lambda.max","ms.error","r.squared","r.hat"))               
mat.dcdiag
```

```
##      clones   lambda.max    ms.error    r.squared    r.hat
## [1,]      1 0.0010108009 0.007055683 0.0006908974 1.002480
## [2,]     10 0.0009812864 0.018414260 0.0016752113 1.003646
## [3,]     50 0.0010256574 0.011075997 0.0010300549 1.000944
## [4,]    100 0.0010446431 0.008833334 0.0008841003 1.001552
```

## References

King R, Morgan BJT, Gimenez O and Brooks S (2009). Bayesian Analysis for Population Ecology. Chapman & Hall/CRC Interdisciplinary Statistics.

Lele, SR, D Dennis, and F Lutscher. 2007. Data Cloning: Easy Maximum Likelihood Estimation for Complex Ecological Models Using Bayesian Markov Chain Monte Carlo Methods. Ecology Letters 10: 551–63.

Lele, SR, K Nadeem, and B Schmuland. 2010. Estimability and Likelihood Inference for Generalized Linear Mixed Models Using Data Cloning. Journal of the American Statistical Association 105: 1617–25.

Marzolin, G, A Charmantier, and O Gimenez. 2011. Frailty in State-Space Models: Application to Actuarial Senescence in the Dipper. Ecology 92: 562–67.

Missov, TI. 2013. Gamma-Gompertz Life Expectancy at Birth. Demographic Research 28: 259–270.

R Core Team. 2016. R: A Language and Environment for Statistical Computing. Vienna, Austria: R Foundation for Statistical Computing. http://www.R-project.org/.

Robert CP and Casella G (2010). Introducing Monte Carlo Methods with R. Springer.

Solymos, P. 2010. dclone: Data Cloning in R. The R Journal 2: 29-37.


