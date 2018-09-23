estfun_observed <- function(gen, r) {
	data.frame(
		mean=mean(gen),
		RR=1/mean(exp(-r*gen))
	)
}

estfun_nonparametric <- function(gen, r) {
	data.frame(
		mean=weighted.mean(gen, exp(r*gen)),
		RR=mean(exp(r*gen))
	)
}

bootfun <- function(gen, 
					logr, logr.sd,
					estfun,
					nsim=1000) {
	bootlist <- vector('list', nsim)
	
	for (i in 1:nsim) {
		r <- exp(rnorm(1, logr, logr.sd))
		
		bootgen <- sample(gen, length(gen), replace=TRUE)
		
		bootlist[[i]] <- estfun(bootgen, r/7)
	}
	
	bootdata <- do.call("rbind", bootlist) 
	
	est <- rbind(estfun(gen, exp(logr)/7), 
				 apply(bootdata, 2, quantile, c(0.025, 0.975), na.rm=TRUE))

	rownames(est) <- c("estimate", "lwr", "upr")
	
	est
}
