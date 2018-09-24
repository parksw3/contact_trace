library(bbmle)
source("../sim/ebola_param_seminr.R")
source("../R/generation.R")
source("../R/estimate.R")

load("../sim/seminr_sim.rda")

weeks <- 3:5*3

estlist <- vector('list', length(weeks))

for (j in 1:length(weeks)) {
	sublist <- vector('list', length(reslist))
	
	for (i in 1:length(reslist)) {
		print(i)
		sim <- reslist[[i]]
		
		gendata <- generation.data(sim, tmax=7*weeks[j])
		
		est <- as.data.frame(t(parametricfun(true.R, true.mean, true.shape, gendata, 7*weeks[j])))
		
		est[,4] <- c(est[1,2] < true.mean && true.mean < est[1,3], 
					 est[2,2] < true.R && true.R < est[2,3])
		
		colnames(est)[4] <- "coverage"
		
		est$type <- c("mean", "RR")
		
		est$sim <- i
		
		sublist[[i]] <- est
	}
	
	subdf <- do.call("rbind", sublist)
	
	rownames(subdf) <- NULL
	
	subdf$weeks <- weeks[j]
	
	estlist[[j]] <- subdf
}

save("estlist", file="seminr_parametric.rda")
