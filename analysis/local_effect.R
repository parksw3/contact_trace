library(igraph)
source("ebola_param_seminr.R")
source("../R/seminr.R")
source("../R/generation.R")

beta <- 1

N <- 5

nsim <- 1000

set.seed(101)
reslist_full <- vector('list', nsim)
i <- 1
while (i <= nsim) {
	rr <- seminr.full(N, beta/(N-1), sigma, m.sigma, gamma, n.gamma, I0 = 1, imax=10,
					  keep.intrinsic = TRUE)
	reslist_full[[i]] <- rr
	i <- i+1
}

g <- graph.tree(5, children=4)

set.seed(101)
reslist_local <- vector('list', nsim)
i <- 1
while (i <= nsim) {
	rr <- seminr(g, beta/(N-1), sigma, m.sigma, gamma, n.gamma, initial_infected=1, imax=10,
				 keep.intrinsic = TRUE)
	
	if (nrow(rr$data) > 1) {
		reslist_local[[i]] <- rr
		i <- i+1
	}
}

gendata_full <- lapply(reslist_full, network.generation, type="backward", plot=FALSE)
gendata_full <- unlist(gendata_full)
gendata_local <- lapply(reslist_local, network.generation, type="backward", plot=FALSE)
gendata_local <- unlist(gendata_local)
mean(gendata_full)
mean(gendata_local)

t <- seq(0, 110, 1)
genden <- hist(gen, breaks=t, plot=FALSE)$density


hist(gendata_local, freq=FALSE, breaks=30)
lines(head(t, -1), genden)

true.mean <- 1/sigma + 1/gamma




