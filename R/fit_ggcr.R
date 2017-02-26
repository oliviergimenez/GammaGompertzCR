#' fit_ggcr
#'
#' This function fits a Gamma-Gompertz model to capture-recapture data using data cloning. 
#' @param mydata a matrix of 1's and 2's for non-detections and detection
#' @param clo number of clones
#' @param nu number of updates
#' @param ni number of iterations
#' @param nt number of thinning iterations
#' @param nc number of chains
#' @return This function returns values from the posterior distribution of the parameters of a Gamma-Gompertz model.
#' @author Gilbert Marzolin, Olivier Gimenez
#' @keywords package
#' @export

fit_ggcr <- function(mydat,clo,nu,ni,nt,nc){

# parameters to be estimated (alphabetical order)
params <- c("a","b","k","psi","sigeta")

# get descriptors of dataset
dat_features <- prep_data(mydat)

# updating the priors from one clone to the next to speed up convergence 
# length(par)=5:the last param is sigeta but updating works on pri
# the last of which being log.sigeta : 
upfun <- function(x) {
	if (missing(x)) {
		np <- 5
		return(cbind(c(0.3,0.49,5.3,1.5,0.5),rep(100,np)))  
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
dat <- list(mydat = mydat, nind = dat_features$nind, nyears = dat_features$nyears, first = dat_features$first, age = dat_features$age, pri = upfun(), clo = 1, st = known_state(mydat,2))

cl <- makePSOCKcluster(3)

# Run the MCMC analysis 
m <- dc.parfit(cl, data = dat, params = params, model = ggcr_model, n.clones = clo, n.update = nu, n.iter = ni, thin = nt, n.chains = nc, multiply = "clo", unchanged = c("nyears","nind","pri"), update = "pri", updatefun = upfun, partype = "parchains", flavour = "jags")  

stopCluster(cl)

return(m)

}   
