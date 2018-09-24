library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

source("../sim/ebola_param_seminr.R")

dir <- c("../analysis/")
files <- c("seminr_observed.rda", "seminr_nonparametric.rda", "seminr_parametric.rda")

comblist <- vector('list', 3)

for (i in 1:length(files)) {
	load(paste0(dir, files[i]))
	
	comblist[[i]] <- estlist
}

names(comblist) <- c("observed", "nonparametric", "parametric")

combdata <- comblist %>%
	lapply(bind_rows) %>%
	bind_rows(.id="method") %>%
	mutate(method=factor(method, levels=c("observed", "nonparametric", "parametric")),
		   weeks=factor(weeks, levels=3:5*3),
		   type=factor(type, labels=c("mean generation interval", "reproductive number")))

coverdata <- combdata %>%
	group_by(method, type, weeks) %>%
	summarize(coverage=mean(coverage, na.rm=TRUE))

estdata <- combdata %>%
	group_by(method, type, weeks) %>%
	summarize(
		median=median(estimate),
		lwr=quantile(estimate, 0.025),
		upr=quantile(estimate, 0.975)
	) %>%
	filter(!(method=="parametric" & weeks==3))

truedata <- data.frame(
	type=c("mean generation interval", "reproductive number"),
	value=c(true.mean, true.R)
)

g1 <- ggplot(estdata) +
	geom_point(aes(weeks, median, col=method, shape=method), position=position_dodge(width=0.2)) +
	geom_errorbar(aes(weeks, ymin=lwr, ymax=upr, col=method), width=0.2, position=position_dodge(width=0.2)) + 
	geom_hline(data=truedata, aes(yintercept=value), lty=2) +
	scale_y_continuous("estimate") +
	facet_wrap(~type, scale="free_y") +
	theme(
		axis.title.x = element_blank(),
		strip.background = element_blank(),
		legend.position = "top",
		legend.title = element_blank()
	)

g2 <- ggplot(coverdata) +
	geom_point(aes(weeks, coverage, col=method, shape=method)) +
	geom_path(aes(as.numeric(weeks), coverage, col=method)) +
	geom_hline(yintercept=0.95, lty=2) +
	# scale_x_continuous("weeks") +
	scale_y_continuous("coverage probability",limits=c(0, 1)) +
	facet_wrap(~type, scale="free_y") +
	theme(
		strip.background = element_blank(),
		legend.position = "none",
		strip.text = element_blank()
	)

gcomp <- arrangeGrob(g1, g2, heights=c(6, 5))

ggsave("compare_methods.pdf", gcomp, width=8, height=6)
