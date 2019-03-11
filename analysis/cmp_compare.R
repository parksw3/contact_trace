library(igraph)
library(dplyr)
library(bbmle)
library(epigrowthfit)
source("../sim/ebola_param_seminr.R")
source("../R/generation.R")
source("../R/estimate.R")
source("../R/empirical.R")

load("../sim/cmp_seir_sim.rda")

intrinsic_fun <- function(tau) {
	sigma*gamma/(gamma-sigma) * (exp(-sigma*tau)-exp(-gamma*tau))
}

conditional_fun <- function(tau) {
	sigma*(gamma+lambda)/(gamma+lambda-sigma) * (exp(-sigma*tau)-exp(-(gamma+lambda)*tau))
}

true.var <- integrate(function(x) (x-true.mean)^2 * intrinsic_fun(x), 0, Inf)[[1]]

true.shape <- true.mean^2/true.var

cmpGraph <- read.table("../data/CA-CondMat.txt")
cmpGraph <- graph.data.frame(cmpGraph)
cmpGraph <- as.undirected(cmpGraph)
cmpGraph <- delete_vertex_attr(cmpGraph, "name")

d <- degree(cmpGraph)

kappa <- var(d)/mean(d) + mean(d) -1

beta <- 0.13/5

## This is SEIR
R0fun.network <- function(r) (gamma+r)/(gamma*sigma/(sigma+r)+r/kappa)

## This is SE2IR
## R0fun.network <- function(r) {
##	lambda <- (2 * sigma+r)^2 * (r + gamma)/((kappa-1) * 4 * sigma^2 - 4 * sigma * r - r^2)
##	
##	kappa * lambda/(lambda + gamma)
## }

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
	
	ee <- MASS::glm.nb(cases~day, data=fitdata)
	
	dd <- data.frame(r=coef(ee)[2])
	
	dd
}) %>%
	bind_rows(.id="sim")

individual.est <- lapply(reslist, function(x){
	tmax <- x$data[1000-9,1]
	gendata <- generation.data(x, tmax=tmax)
	
	est <- as.data.frame(t(parametricfun(5, true.mean, true.shape, gendata, tmax)))
	
	est$type <- c("RR", "mean", "shape")
	
	print(est)
	
	est
}) %>%
	bind_rows(.id="sim")

individual.RR <- lapply(1:100, function(x){
	gg <- growth[x,2]
	
	mean <- filter(individual.est, type=="mean")[x,2]
	shape <- filter(individual.est, type=="shape")[x,2]
	
	(1 + gg * 1/shape * mean)^(shape)
})

population.est <- lapply(1:100, function(x){
	cg <- censor.gi[[x]]
	gg <- growth[x,2]
	
	data.frame(
		estimate=c(weighted.mean(cg, exp(gg*cg)), mean(exp(gg*cg))),
		type=c("mean", "RR")
	)
}) %>%
	bind_rows(.id="sim")

population.est2 <- lapply(1:100, function(x){
	cg <- censor.gi[[x]]
	gg <- growth[x,2]
	
	populationfun(true.mean, true.shape, cg, gg)
}) %>%
	bind_rows(.id="sim")

observed.est <- lapply(1:100, function(x){
	cg <- censor.gi[[x]]
	gg <- growth[x,2]
	
	data.frame(
		estimate=c(mean(cg), 1/mean(exp(-gg*cg))),
		type=c("mean", "RR")
	)
}) %>%
	bind_rows(.id="sim")

intrinsic.est <- lapply(growth$r, function(x) 1/mean(exp(-x*gen)))

network.est <- lapply(growth$r, R0fun.network)

trapman.est <- lapply(reslist, trapman.R0, 100)

empirical.est <- lapply(reslist, empirical.R0, 100)

RRdata <- data.frame(
	observed=observed.est$estimate[observed.est$type=="RR"],
	individual=unlist(individual.RR),
	population=population.est$estimate[population.est$type=="RR"],
	intrinsic=unlist(intrinsic.est),
	network=unlist(network.est),
	trapman=unlist(trapman.est),
	empirical=unlist(empirical.est)
)

gendata <- data.frame(
	observed=observed.est$estimate[observed.est$type=="mean"],
	individual=individual.est$estimate[individual.est$type=="mean"],
	population=population.est$estimate[population.est$type=="mean"]
)

save("RRdata", "gendata", file="cmp_compare.rda")
