library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
source("../sim/ebola_param.R")
source("../R/seir.R")
source("../R/generation.R")

beta <- 1
N <- 5
lambda <- beta/(N-1)

nsim <- 5000

set.seed(101)
reslist_full <- vector('list', nsim)
i <- 1
while (i <= nsim) {
	rr <- seir.full(N, beta/(N-1), sigma, gamma, I0 = 1, imax=10)
	reslist_full[[i]] <- rr
	i <- i+1
}

set.seed(101)
reslist_local <- vector('list', nsim)
i <- 1
while (i <= nsim) {
	rr <- seir.local(4, beta/(N-1), sigma, gamma)
	
	reslist_local[[i]] <- rr
	i <- i+1
}

gendata_full <- reslist_full %>%
	lapply(network.generation, type="backward", plot=FALSE) %>%
	unlist

gendata_cond <- reslist_local %>%
	lapply("[[", "conditional") %>%
	unlist

gendata_int <- reslist_local %>%
	lapply("[[", "intrinsic") %>%
	unlist
	
gendata_mean <- data.frame(
	type=c("intrinsic", "conditional", "homogeneous"),
	value=c(mean(gendata_int), mean(gendata_cond), mean(gendata_full))
)

true_mean <- data.frame(
	key=c("intrinsic", "conditional"),
	value=c(1/sigma+1/gamma, 1/sigma + 1/(gamma+lambda))
)
	
gendata_total <- list(
	intrinsic=data.frame(
		value=gendata_int
	),
	conditional=data.frame(
		value=gendata_cond
	),
	homogeneous=data.frame(
		value=gendata_full
	) 
) %>%
	bind_rows(.id="type") %>%
	mutate(type=factor(type, levels=c("intrinsic", "conditional", "homogeneous")))

tvec <- seq(0, 80, by=0.1)

dendata <- list(
	data.frame(
		key="intrinsic",
		x=tvec,
		density=intrinsic_fun(tvec)
	),
	data.frame(
		key="conditional",
		x=tvec,
		density=conditional_fun(tvec)
	) 
) %>%
	bind_rows

gglocal <- ggplot(gendata_total) +
	geom_histogram(aes(value, y=..density..), bins=30, col="black", fill=NA) +
	geom_line(data=dendata, aes(x, density, col=key), lwd=1.2) +
	geom_vline(data=true_mean, aes(xintercept=value, col=key), lwd=1.2) +
	geom_vline(data=gendata_mean, aes(xintercept=value), lty=2, lwd=1.2) +
	geom_text(data=gendata_mean, aes(55, 0.075, label=paste0("observed mean: ", round(value, 1)))) +
	scale_x_continuous("generation interval (days)",limits=c(0, 80)) +
	scale_y_continuous(expand=c(0, 0), limits=c(0, 0.08)) +
	facet_wrap(~type) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank(),
		axis.title.y = element_blank(),
		axis.text.y = element_blank(),
		axis.ticks.y = element_blank(),
		panel.border = element_blank(),
		axis.line.x = element_line(),
		legend.title = element_blank(),
		legend.position="top"
	)

ggsave("local_effect.pdf", gglocal, width=8, height=3.5)

