library(igraph)
library(dplyr)
library(bbmle)
library(epigrowthfit)
source("../sim/ebola_param_seminr.R")
source("../R/generation.R")
source("../R/estimate.R")
source("../R/empirical.R")

load("../sim/cmp_seminr_sim.rda")

cmpGraph <- read.table("../data/CA-CondMat.txt")
cmpGraph <- graph.data.frame(cmpGraph)
cmpGraph <- as.undirected(cmpGraph)
cmpGraph <- delete_vertex_attr(cmpGraph, "name")

d <- degree(cmpGraph)

kappa <- var(d)/mean(d) + mean(d) -1

beta <- gamma/7

## This is SEIR
## R0fun.network <- function(r) (gamma+r)/(gamma*sigma/(sigma+r)+r/kappa)

## This is SE2IR
R0fun.network <- function(r) {
	lambda <- (2 * sigma+r)^2 * (r + gamma)/((kappa-1) * 4 * sigma^2 - 4 * sigma * r - r^2)
	
	kappa * lambda/(lambda + gamma)
}



censor.gi <- lapply(
	reslist
	, network.generation
	, plot=FALSE
	, interval=c(10, 400)
	, interval.type="cases"
)

growth <- lapply(reslist, function(x){
	day <- floor(x$data$time) + 1
	
	fitdata <- data.frame(
		day=unique(day),
		cases=as.vector(table(day))
	)
	
	fitdata <- fitdata[cumsum(fitdata$cases) <= 1000,]
	
	ee <- epigrowthfit(fitdata$day, fitdata$cases,
				 model="logistic")
	
	dd <- as.data.frame(t(ee@growthRate))
	
	dd
}) %>%
	bind_rows(.id="sim")

individual.est <- lapply(reslist, function(x){
	tmax <- x$data[1000-9,1]
	gendata <- generation.data(x, tmax=tmax)
	
	est <- as.data.frame(t(parametricfun(true.R, true.mean, true.shape, gendata, tmax)))
	
	est$type <- c("mean", "RR")
	
	est
}) %>%
	bind_rows(.id="sim")

population.est <- lapply(1:50, function(x){
	cg <- censor.gi[[x]]
	gg <- growth[x,2]
	
	data.frame(
		estimate=c(weighted.mean(cg, exp(gg*cg)), mean(exp(gg*cg))),
		type=c("mean", "RR")
	)
}) %>%
	bind_rows(.id="sim")

population.est2 <- lapply(1:50, function(x){
	cg <- censor.gi[[x]]
	gg <- growth[x,2]
	
	populationfun(true.mean, true.shape, cg, gg)
}) %>%
	bind_rows(.id="sim")

observed.est <- lapply(1:50, function(x){
	cg <- censor.gi[[x]]
	gg <- growth[x,2]
	
	data.frame(
		estimate=c(mean(cg), 1/mean(exp(-gg*cg))),
		type=c("mean", "RR")
	)
}) %>%
	bind_rows(.id="sim")

intrinsic.est <- lapply(growth$value, function(x) 1/mean(exp(-x*gen)))

network.est <- lapply(growth$value, R0fun.network)


#plot(density(unlist(intrinsic.est)), ylim=c(0, 0.8))
#lines(density(unlist(network.est)), col=4)
#lines(density(sapply(individual.est, function(x) x[2,1])), col=2)
#lines(density((observed.est %>% filter(type=="RR"))$estimate), col=3)




