#' prep_data
#'
#' This function gets number of individuals, number of capture occasions, occasions of first capture and age from a capture-recapture dataset
#' @param mydata a matrix of 1's and 2's for non-detections and detection
#' @return This function returns a list with components the number of individuals, the number of capture occasions, the occasions of first capture and age.
#' @author Gilbert Marzolin, Olivier Gimenez
#' @keywords package
#' @export

prep_data <- function(mydata){

# number of individuals and capture occasions
nind <- dim(mydat)[[1]]
nyears <- dim(mydat)[[2]]

# vector with year of first captures
get_first <- function(x) min(which(x==1))
first <- apply(mydat,1,get_first)

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
