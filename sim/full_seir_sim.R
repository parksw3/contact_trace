library(igraph)
source("ebola_param_seminr.R")
source("../R/seminr.R")

nsim <- 100

set.seed(101)
reslist <- vector('list', nsim)
i <- 1
while (i <= nsim) {
	print(i)
	rr <- seminr.full(N, beta/N, sigma, 1, gamma, 1, I0 = 10, imax=1010, tmax=200,
					  keep.intrinsic=TRUE)
	if(nrow(rr$data) > 1000) {
		reslist[[i]] <- rr
		i <- i+1
	}
	
	save("reslist", file="full_seir_sim.rda")
}

save("reslist", file="full_seir_sim.rda")
