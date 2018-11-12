library(igraph)
source("ebola_param_seminr.R")
source("../R/seminr.R")

cmpGraph <- read.table("../data/CA-CondMat.txt")
cmpGraph <- graph.data.frame(cmpGraph)
cmpGraph <- as.undirected(cmpGraph)
cmpGraph <- delete_vertex_attr(cmpGraph, "name")

nsim <- 50

beta <- 0.4/5

set.seed(101)
reslist <- vector('list', nsim)
i <- 1
while (i <= nsim) {
	print(i)
	rr <- seminr(cmpGraph, beta, sigma, m.sigma, gamma, n.gamma, I0 = 10, tmax=2000)
	if(nrow(rr$data) > 1000) {
		reslist[[i]] <- rr
		i <- i+1
	}
	
	save("reslist", file="cmp_seminr_sim.rda")
}

save("reslist", file="cmp_seminr_sim.rda")
