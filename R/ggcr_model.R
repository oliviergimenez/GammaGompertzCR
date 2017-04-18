#' ggcr_model: Gamma-Gompertz capture-recapture model
#'
#' This function implements a gamma-gompertz capture-recapture model in the BUGS syntax
#' @author Gilbert Marzolin, Olivier Gimenez
#' @keywords package
#' @export

ggcr_model <- function(){

# mean demographic parameters prepared for updating priors
a ~ dnorm(pri[1,1],pri[1,2])              
b ~ dnorm(pri[2,1],pri[2,2])               
k ~ dnorm(pri[3,1],pri[3,2])               
psi ~ dnorm(pri[4,1],pri[4,2])            
log.sigeta ~ dnorm(pri[5,1],pri[5,2])

for (i in 1:nind){  # for each unit
	st[i,first[i]] <- mydat[i,first[i]]

	for (j in first[i]:(nyears-1)){
		
        # phi[i,j] is survival for unit i from year j to j+1
		phi[i,j] <- exp(k*(log(abs(b*k-a+a*exp(b*age[i,j])))-log(abs(b*k-a+a*exp(b*(age[i,j]+1))))))

        ## STATE EQUATIONS ##
		# draw st[i,j+1] given st[i,j]: st[i,j] is state 1 if alive and 2 if dead in year j
		st[i,j+1]~ dcat((c(phi[i,j],1-phi[i,j])*equals(st[i,j],1))+c(0,1)*equals(st[i,j],2))
        
        ## OBSERVATION EQUATIONS ##
		# draw Obs[i,j+1] given st[i,j+1]
		mydat[i,j+1]~ dcat(c(p[i,j],1-p[i,j])*equals(st[i,j+1],1)+ c(0,1)*equals(st[i,j+1],2))
		
		# p[i,j] is detection prob at j+1
		logit(p[i,j]) <- psi + eta[j]  #same detection for all units with time random effect.

     }	#j

} #i

for (j in 1:(nyears-1)) { eta[j] ~ dnorm(0,taueta) }
taueta <- pow(sigeta,-2)
sigeta <- exp(log.sigeta)

} #model
