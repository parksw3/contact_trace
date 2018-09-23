source("../sim/ebola_param_seminr.R")
source("../R/generation.R")
source("../R/estimate.R")

load("../sim/seminr_sim.rda")

load("seminr_growth.rda")

estlist <- vector('list', length(reslist))

for (i in 1:length(reslist)) {
	print(i)
	sim <- reslist[[i]]
	fit <- fitlist[[i]]
	
	logr <- fit@mle2@coef[1]
	logr.sd <- sqrt(diag(fit@mle2@vcov)[1])
	
	if (is.na(logr.sd)) {
		logr.sd <- max((log(fit@growthRate[3]) - log(fit@growthRate[1]))/1.96,
					   (log(fit@growthRate[1]) - log(fit@growthRate[2]))/1.96)	
	} 
	
	gen <- network.generation(sim, interval=c(0, 7*15), plot=FALSE)
	
	est <- as.data.frame(t(bootfun(gen, logr, logr.sd, estfun_nonparametric)))
	
	est$type <- c("mean", "RR")
	
	est[,4] <- c(est[1,2] < true.mean && true.mean < est[1,3], 
				 est[2,2] < true.R && true.R < est[2,3])
	
	colnames(est)[4] <- "coverage"
	
	estlist[[i]] <- est
}

save("estlist", file="seminr_nonparametric.rda")
