context("fit_ggcr")

test_that("fit_ggcr output", {
	data(dipper_tiny)
	clo <- c(10) # numbers of clones
	nu <- 1000 # number of updates
	ni <- 5000 # number of iterations
	nt <- 5 # thinning
	nc <- 2 # number of chains
	initmeans <- c(0.15,0.5,4.9,1.4,-0.4) # means of normally distributed priors for parameters
	initprec <- 1000 # common precision of the normal priors for the first number of clones
	post_inf <- fit_ggcr(dipper_tiny,clo,nu,ni,nt,nc,initmeans,initprec)
	a <- coef(post_inf[[length(clo)]])[1]
	b <- coef(post_inf[[length(clo)]])[2] 
	k <- coef(post_inf[[length(clo)]])[3]
	rounded_res <- round(c(a,b,k),2)
	expect_equivalent(rounded_res, c(0.12,0.49,4.90))
})

