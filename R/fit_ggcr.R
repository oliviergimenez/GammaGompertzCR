#' fit_ggcr
#'
#' This function fits a Gamma-Gompertz model to capture-recapture data using data cloning. To improve parameter estimates clones of the data are used. To speed up convergence at each new number of clones the means and precisions of the priors supposed normally distributed (the last parameter in alphabetical order is logged) are updated. Only for the first number of clones the means must be put in vector initmeans, and a common precision is set in initprec. The specifics of the Markov chains along with the different series of clones are defined
#' @param mydata a matrix of 1's and 2's for detections and non-detections (the first 1 being the birth year)
#' @param clo vector of numbers of clones in each lot: increasing order is advisable and start with 1 clone is not compulsory. Default is c(1,10,50,100)
#' @param nu number of updates in the Markov chains. Default is 1000.
#' @param ni number of iterations in the Markov chains. Default is 5000.
#' @param nt number of thinning (see Link and Eaton 2012). Default is 5.
#' @param nc number of Markov chains. Default is 3.
#' @param initmeans for the first lot of clones, vector of the normally distributed means of the parameter priors, the last parameter being logged. Default is c(0.15,0.5,4.9,1.4,-0.4).
#' @param initprec a common precision of the prior normal distributions is set for all the parameters. Default is 1000.
#' @return This function returns values from the posterior distribution of the parameters of a Gamma-Gompertz model.
#' \itemize{
#'   \item a is the coefficient in the baseline mortality rate
#'   \item b is the coefficient in the exponent of the mortality rate
#'   \item k is the inverse of the Gamma variance
#'   \item psi is the mean of the detection rate on the logit scale
#'   \item sigeta is the standard deviation of the detection rate on the logit scale
#' }
#' @references Link, W. A. and Eaton, M. J. (2012), On thinning of chains in MCMC. Methods in Ecology and Evolution, 3: 112–115. doi:10.1111/j.2041-210X.2011.00131.x
#' @author Gilbert Marzolin, Olivier Gimenez
#' @keywords package
#' @importFrom stats sd
#' @import dcmle
#' @import dclone
#' @import parallel

#' @export

fit_ggcr <- function(mydata, clo = c(1,10,50,100), nu = 1000, ni = 5000, nt = 5, nc = 3, initmeans = c(0.15,0.5,4.9,1.4,-0.4), initprec = 1000){

# parameters to be estimated (alphabetical order)
params <- c("a","b","k","psi","sigeta")

# get descriptors of dataset
dat_features <- prep_data(mydata)

# updating the priors from one clone to the next to speed up convergence
# length(par)=5:the last param is sigeta but updating works on pri
# the last of which being log.sigeta
# dcsd is datacloning sd ie product of standard deviation by square root of nclones
upfun <- function(x) {
	if (missing(x)) {
		np <- 5
            return(cbind(initmeans, rep(initprec,np)))
 } 	else {
 		ncl <- nclones(x)
 		if (is.null(ncl))
 		ncl <- 1
 		par <- coef(x)
 		prec <- exp(-2*log(dcsd(x)))
 		log.sigeta <- mcmcapply(x[,"sigeta"],log)
 		par[length(par)] <- mean(log.sigeta)
 		prec[length(prec)] <- exp(-2*log(sd(log.sigeta)* sqrt(ncl)))
 		return(cbind(par, prec))
 	}
}

known_state <- function(ms,notseen){
state <- ms
state[state == notseen] <- NA
for (i in 1:dim(ms)[1]){
	n1 <- min(which(state[i,]==1))
	n2 <- max(which(state[i,]==1))
	state[i,n1:n2]<-1
	state[i,n1]<- NA
}
return(state)
}

# list of data
dat <- list(mydata = mydata, nind = dat_features$nind, nyears = dat_features$nyears, first = dat_features$first, age = dat_features$age, pri = upfun(), clo = 1, st = known_state(mydata,2))

cl <- makePSOCKcluster(3)

# Run the MCMC analysis
    m <- lapply(clo, function(z)
         jags.parfit( cl = cl, data = dclone(dat, n.clones = z,attrib = TRUE, unchanged = c("nyears", "nind")),
         params = params, model = ggcr_model,
         n.update = nu, n.iter = ni, thin = nt,
         n.chains = nc, multiply = NULL, update = "pri", updatefun = upfun,
         initsfun = NULL, flavour = c("jags"), partype = c("parchains")))

stopCluster(cl)

return(m)

}
