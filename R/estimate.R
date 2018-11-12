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

minuslogl.gamma <- function(log.mean, log.shape, gen) {
	mean <- exp(log.mean)
	shape <- exp(log.shape)
	rate <- shape/mean
	
	-sum(dgamma(gen, shape=shape, rate=rate, log=TRUE))
}

populationfun <- function(mean, shape, gen, r) {
	m <- mle2(minuslogl.gamma, 
			  start=list(log.mean=log(mean), log.shape=log(shape)),
			  method="Nelder-Mead",
			  optimizer="optim",
			  data=list(gen=gen),
			  control=list(maxit=10000)
	)
	
	cc <- exp(coef(m))
	
	mean <- cc[1]
	shape <- cc[2]
	rate <- shape/mean
	adj.rate <- rate - r
	adj.mean <- shape/adj.rate
	
	data.frame(
		estimate=unname(c((1+r*adj.mean/shape)^(shape), adj.mean, shape)),
		type=c("RR", "mean", "shape")
	) 
}

minuslogl.full <- function(log.R, log.mean, log.shape, data, tmax) {
	R <- exp(log.R)
	mean <- exp(log.mean)
	shape <- exp(log.shape)
	scale <- mean/shape
	
	t_censor <- tmax-data$t_infected
	
	nll <- -(sum(
		log.R +
			dgamma(data$generation, shape=shape, scale=scale, log=TRUE),
		na.rm=TRUE
	) + sum(
		-R * pgamma(t_censor, shape=shape, scale=scale),
		na.rm=TRUE
	))
	
	if (nll==0) return(Inf) else return(nll)
}

parametricfun <- function(R, mean, shape,
						  data,
						  tmax) {
	m <- mle2(minuslogl.full, 
		 start=list(log.R=log(R), log.mean=log(mean), log.shape=log(shape)),
		 method="Nelder-Mead",
		 optimizer="optim",
		 data=list(data=data, tmax=tmax),
		 control=list(maxit=10000)
	)
	
	pp <- profile(m, 1:2)
	
	ci <- exp(confint(pp))
	
	est <- rbind(
		exp(coef(m)[1:2]),
		t(ci)
	)
	
	rownames(est) <- c("estimate", "lwr", "upr")
	
	cbind(est[,2], est[,1])
}
