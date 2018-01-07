#' prep_data
#'
#' This function gets number of individuals, number of capture occasions, occasions of first capture and age from a capture-recapture dataset. Given a matrix of capture-recapture data (1 for detection, 2 for non-detection) at regular time periods, for instance years, when the first capture of an animal occurs during its birth year, at each year the age of that animal is known. Quantities nind and nyears are respectively the numbers of animals and of time units in this matrix.
#' @param mydata a matrix of 1's and 2's for detections and non-detections
#' @return This function returns a list with components the number of individuals, the number of capture occasions, the occasion of first capture for each individual and the age of each individual over time.
#' \itemize{
#'   \item nind is the number of rows in mydata
#'   \item nyears is the number of columns in mydata
#'   \item first is a vector gathering for each row the rank of the first '1' denoting first detection
#'   \item age is defined as follows: as the first '1' in row 'i' denotes the time of birth of the 'i'th animal, 'age[i,j]' for j in (first[i],nyears-1) is the age at time 'j'
#' }
#' @author Gilbert Marzolin, Olivier Gimenez
#' @references King R, Morgan BJT, Gimenez O and Brooks S (2009). Bayesian Analysis for Population Ecology. Chapman & Hall/CRC Interdisciplinary Statistics.
#' @references Robert CP and Casella G (2010). Introducing Monte Carlo Methods with R. Springer.
#' @keywords package
#' @export

prep_data <- function(mydata){

# number of individuals and capture occasions
nind <- dim(mydata)[[1]]
nyears <- dim(mydata)[[2]]

# vector with year of first captures
get_first <- function(x) min(which(x==1))
first <- apply(mydata,1,get_first)

# define age (age is not used for j = nyears)
age <- array(NA,c(nind,nyears-1))
for (i in 1:nind){
	for (j in first[i]:(nyears-1)){
		age[i,j]<-j-first[i]
	}
}

res = list(nind=nind,nyears=nyears,first=first,age=age)
return(res)

}   
