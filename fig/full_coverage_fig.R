library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
source("../sim/ebola_param_seminr.R")

load("../analysis/full_coverage.rda")

intrinsic_fun <- function(tau) {
	sigma*gamma/(gamma-sigma) * (exp(-sigma*tau)-exp(-gamma*tau))
}

true.var <- integrate(function(x) (x-true.mean)^2 * intrinsic_fun(x), 0, Inf)[[1]]

true.shape <- true.mean^2/true.var

coverdata <- individual.est %>%
	group_by(type) %>%
	summarize(
		coverage=mean(coverage)
	) %>%
	filter(
		type!="RR"
	) %>%
	mutate(
		type=factor(type, labels=c("Mean", "CV"))
	)

g1 <- ggplot(filter(individual.est, type=="mean")) +
	geom_boxplot(aes(y=estimate)) +
	geom_hline(yintercept=true.mean, lty=2) +
	scale_y_continuous("Mean generation time (days)") +
	ggtitle("A") +
	theme(
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank()
	)

g2 <- ggplot(filter(individual.est, type=="shape")) +
	geom_boxplot(aes(y=sqrt(1/estimate))) +
	geom_hline(yintercept=sqrt(1/true.shape), lty=2)  +
	scale_y_continuous("Coefficient of variation (CV)")  +
	ggtitle("B") +
	theme(
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank()
	)

g3 <- ggplot(coverdata) +
	geom_point(aes(type, coverage), shape=1, size=3) +
	geom_hline(yintercept=0.95, lty=2) +
	scale_y_continuous("Coverage probability", limits=c(0, 1)) +
	scale_x_discrete("Parameters")  +
	ggtitle("C")

xvec <- seq(0, 51, by=0.1)

realdensity <- data.frame(
	x=xvec,
	density=intrinsic_fun(xvec)
)

estdensity <- lapply(1:100, function(x) {
	ii <- filter(individual.est, sim==x)
	
	est.mean <- ii$estimate[2]
	est.shape <- ii$estimate[3]
	
	data.frame(
		x=xvec,
		density=dgamma(xvec, est.shape, scale=est.mean/est.shape)
	)
}) %>%
	bind_rows(.id="sim")

g4 <- ggplot(estdensity) +
	geom_line(aes(x, density, group=sim, col="estimated"), alpha=0.2) +
	geom_line(data=realdensity, aes(x, density, col="true"), lwd=1.2, lty=2) +
	scale_x_continuous("Generation time (days)", expand=c(0, 0)) +
	scale_y_continuous("Density", expand=c(0, 0)) +
	scale_color_manual(values=c("black", "red")) +
	ggtitle("D") +
	theme(
		panel.grid = element_blank(),
		legend.position=c(0.75, 0.69),
		legend.title=element_blank()
	)
	
rvec <- seq(0, 0.11, by=0.001)

realdata <- data.frame(
	r=rvec,
	R0=(1 + rvec/sigma) * (1 + rvec/gamma)
)

estdata <- lapply(1:100, function(x) {
	ii <- filter(individual.est, sim==x)
	
	est.mean <- ii$estimate[2]
	est.shape <- ii$estimate[3]
	
	data.frame(
		r=rvec,
		R0=(1 + est.mean * 1/est.shape * rvec)^est.shape
	)
}) %>%
	bind_rows(.id="sim") %>%
	group_by(r) %>%
	summarize(
		median=median(R0),
		lwr=quantile(R0, 0.025),
		upr=quantile(R0, 0.975)
	)

g5 <- ggplot(estdata) +
	geom_ribbon(aes(r, ymin=lwr, ymax=upr), alpha=0.2, col="black") +
	geom_line(aes(r, median), lwd=1.2) +
	geom_line(data=realdata, aes(r, R0), lwd=1.2, col="red", lty=2) +
	scale_x_continuous(expression(Growth~rate~(day^{-1})), expand=c(0, 0)) +
	scale_y_continuous("Reproductive number") +
	ggtitle("E")

gtop <- arrangeGrob(g1, g2, g3, nrow=1)
gbottom <- arrangeGrob(g4, g5, nrow=1)

gall <- arrangeGrob(gtop, gbottom, nrow=2, heights=c(5.3, 4.5))

ggsave("full_coverage_fig.pdf", gall, width=6, height=4)
