library(epigrowthfit)
library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
source("../sim/ebola_param_seminr.R")
source("../R/seminr.R")
source("../R/generation.R")
source("../R/estimate.R")

weekmax <- 28
weekmin <- 9

set.seed(101)
rr <- seminr.full(N, beta/(N-1), sigma, m.sigma, gamma, n.gamma, I0 = 10, tmax=7*weekmax)

weeks <- weekmin:weekmax

reslist <- vector('list', length(weeks))

for (i in 1:length(weeks)) {
	print(i)
	week <- weeks[i]
	
	obsdata <- rr$data %>%
		filter(time < 7*week)
 	
	fitdata <- data.frame(
		week=1:week,
		cases=as.vector(table(floor(obsdata$time/7)))
	)
	
	epifit <- epigrowthfit(fitdata$week, fitdata$cases,
						   model="logistic")
	
	logr <- epifit@mle2@coef[1]
	logr.sd <- sqrt(diag(epifit@mle2@vcov)[1])
	
	if (is.na(logr.sd) || logr.sd <= 0) {
		logr.sd <- max((log(epifit@growthRate[3]) - log(epifit@growthRate[1]))/1.96,
					   (log(epifit@growthRate[1]) - log(epifit@growthRate[2]))/1.96)	
	} 
	
	gen <- network.generation(rr, plot=FALSE, interval=c(0, 7*week))
	
	gendata <- generation.data(rr, tmax=7*week)
	
	nonparaest <- as.data.frame(t(bootfun(gen, logr, logr.sd, estfun_nonparametric)))
	
	nonparaest$type <- c("mean", "RR")
	nonparaest$method <- "nonparametric"
	
	obsest <- as.data.frame(t(bootfun(gen, logr, logr.sd, estfun_observed)))
	
	obsest$type <- c("mean", "RR")
	obsest$method <- "observed"
	
	paraest <- as.data.frame(t(parametricfun(true.R, true.mean, true.shape, gendata, 7*week)))
	
	paraest$type <- c("mean", "RR")
	paraest$method <- "parametric"
	
	growthest <- data.frame(
		estimate=epifit@growthRate[1],
		lwr=epifit@growthRate[2],
		upr=epifit@growthRate[3],
		type="rr",
		method="epigrowthfit"
	)
	
 	allest <- rbind(obsest, nonparaest, paraest, growthest)
	
 	rownames(allest) <- NULL
 	
 	allest$time <- week
 	
 	reslist[[i]] <- allest
}

resdata <- reslist %>%
	bind_rows

rdata <- resdata %>%
	filter(type=="rr")

meandata <- resdata %>%
	filter(type=="mean") %>%
	mutate(method=factor(method, levels=c("observed", "nonparametric", "parametric")))

RRdata <- resdata %>%
	filter(type=="RR") %>%
	mutate(method=factor(method, levels=c("observed", "nonparametric", "parametric")))

obsdata <- rr$data

fitdata <- data.frame(
	week=1:max(floor(obsdata$time/7)),
	cases=head(as.vector(table(floor(obsdata$time/7))), -1)
)

gincidence <- ggplot(fitdata) +
	geom_line(aes(week, cases), lwd=1.1) +
	scale_x_continuous(limits=c(1, 28), expand=c(0,0)) +
	scale_y_log10("weekly incidence") +
	theme(
		axis.title.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.text.x = element_blank()
	)

ggrowth <- ggplot(rdata) +
	geom_line(aes(time, estimate), lwd=1.1) +
	geom_ribbon(aes(time, ymin=lwr, ymax=upr), alpha=0.1) +
	scale_x_continuous(limits=c(1, 28), expand=c(0,0)) +
	scale_y_continuous("growth rate") +
	theme(
		axis.title.y = element_text(margin = margin(t = 0, r = 10.2, b = 0, l = 0)),
		axis.title.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.text.x = element_blank()
	)

gmean <- ggplot(meandata) +
	geom_line(aes(time, estimate, col=method), lwd=1.1) +
	geom_ribbon(aes(time, ymin=lwr, ymax=upr, fill=method), alpha=0.1) +
	geom_hline(yintercept = true.mean, lty=2) +
	scale_x_continuous(limits=c(1, 28), expand=c(0,0)) +
	scale_y_continuous("mean generation") +
	theme(
		axis.title.y = element_text(margin = margin(t = 0, r = 17.8, b = 0, l = 0)),
		axis.title.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.text.x = element_blank(),
		legend.position = "none"
	)

grepro <- ggplot(RRdata) +
	geom_line(aes(time, estimate, col=method), lwd=1.1) +
	geom_ribbon(aes(time, ymin=lwr, ymax=upr, fill=method), alpha=0.1) +
	geom_hline(yintercept = true.R, lty=2) +
	scale_x_continuous("week", limits=c(1, 28), expand=c(0,0)) +
	scale_y_continuous("reproductive number") +
	theme(
		axis.title.y = element_text(margin = margin(t = 0, r = 15.5, b = 0, l = 0)),
		legend.position = "bottom",
		legend.title = element_blank()
	)

gfinal <- arrangeGrob(gincidence, ggrowth, gmean, grepro, ncol=1,
			 heights=c(2, 2, 2, 3))

ggsave("example.pdf", gfinal, width=6, height=8)
