#' ggcr_model: Gamma-Gompertz capture-recapture model 
#'
#' This function implements a gamma-gompertz capture-recapture model in the BUGS syntax. The survival of unit 'i' at age 'j' is expressed through the parameters of a Gompertz baseline mortality and a Gamma frailty. Then it is used along with the detection probability to tie the animal real states (alive or dead when unobserved) to the observed capture-recapture matrix, hence to deduce the survival parameters by maximum likelihood through data-cloning. The model is parameterized with a (coefficient of the baseline mortality rate), b (coefficient in the exponent of the baseline mortality rate) and k (inverse of variance of the Gamma function) for the Gamma-Gompertz survival along with the nuisance detection parameters psi (mean of logit detection) and sigeta (standard deviation of logit detection). Through marginalization (Missov 2013) the population survival function enables the plot of survival by ages. 
#' @author Gilbert Marzolin, Olivier Gimenez
#' @references Missov, TI. 2013. Gamma-Gompertz Life Expectancy at Birth. Demographic Research 28:259-270. doi:10.4054/DemRes.2013.28.9
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
	st[i,first[i]] <- mydata[i,first[i]]

	for (j in first[i]:(nyears-1)){
		
        # phi[i,j] is survival for unit i from year j to j+1
		phi[i,j] <- exp(k*(log(abs(b*k-a+a*exp(b*age[i,j])))-log(abs(b*k-a+a*exp(b*(age[i,j]+1))))))

        ## STATE EQUATIONS ##
		# draw st[i,j+1] given st[i,j]: st[i,j] is state 1 if alive and 2 if dead in year j
		st[i,j+1]~ dcat((c(phi[i,j],1-phi[i,j])*equals(st[i,j],1))+c(0,1)*equals(st[i,j],2))
        
        ## OBSERVATION EQUATIONS ##
		# draw Obs[i,j+1] given st[i,j+1]
		mydata[i,j+1]~ dcat(c(p[i,j],1-p[i,j])*equals(st[i,j+1],1)+ c(0,1)*equals(st[i,j+1],2))
		
		# p[i,j] is detection prob at j+1
		logit(p[i,j]) <- psi + eta[j]  #same detection for all units with time random effect.

     }	#j

} #i

for (j in 1:(nyears-1)) { eta[j] ~ dnorm(0,taueta) }
taueta <- pow(sigeta,-2)
sigeta <- exp(log.sigeta)

} #model
