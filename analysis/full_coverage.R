library(igraph)
library(dplyr)
library(bbmle)
source("../sim/ebola_param_seminr.R")
source("../R/generation.R")
source("../R/estimate.R")
source("../R/empirical.R")

load("../sim/full_seir_sim.rda")

intrinsic_fun <- function(tau) {
	sigma*gamma/(gamma-sigma) * (exp(-sigma*tau)-exp(-gamma*tau))
}

true.var <- integrate(function(x) (x-true.mean)^2 * intrinsic_fun(x), 0, Inf)[[1]]

true.shape <- true.mean^2/true.var

individual.est <- lapply(reslist, function(x){
	tmax <- x$data[1001,1]
	gendata <- generation.data(x, tmax=tmax)
	
	est <- as.data.frame(t(parametricfun(true.R, true.mean, true.shape, gendata, tmax)))
	
	est$type <- c("RR", "mean", "shape")
	
	est$coverage <- c(
		est$lwr[1] < true.R && true.R < est$upr[1],
		est$lwr[2] < true.mean && true.mean < est$upr[2],
		est$lwr[3] < true.shape && true.shape < est$upr[3]
	)
	
	print(est)
	
	est
}) %>%
	bind_rows(.id="sim")

save("individual.est", file="full_coverage.rda")
