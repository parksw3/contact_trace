library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

source("../sim/ebola_param_seminr.R")

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

load("../analysis/cmp_compare.rda")

RRgroup <- data.frame(
	key=c("contact\ntracing", "population\ncorrection", "individual\ncorrection", "initial\nforward", "empirical", "egocentric", "intrinsic"),
	group=factor(c("tracing based", "tracing based", "tracing based", "empirical", "empirical", "individual based", "individual based"),
				 levels=c("tracing based", "empirical", "individual based"))
)

RRdata2 <- RRdata %>%
	gather(key, value) %>%
	mutate(key=factor(key, 
					  levels=c("observed", "population", "individual", "forward", "empirical", "network", "intrinsic"),
					  labels=c("contact\ntracing", "population\ncorrection", "individual\ncorrection", 
					  		 "initial\nforward", "empirical", "egocentric", "intrinsic"))) %>%
	merge(RRgroup)

ggR <- ggplot(RRdata2) +
	geom_boxplot(aes(key, value, fill=key), alpha=0.7) +
	scale_y_log10("Reproductive number", breaks=c(1, 2, 4, 8, 16)) +
	facet_grid(~group, scale="free", space="free_x") +
	theme(
		panel.spacing.x = unit(0, "lines"),
		strip.background = element_blank(),
		legend.position = "none",
		axis.title.x=element_blank()
	)

g2 <- ggplot(RRdata) +
  geom_point(aes(forward, individual), shape=1) +
  geom_abline(lty=2) +
  scale_x_continuous("initial forward") +
  scale_y_continuous("individual correction") +
  theme(
    panel.grid = element_blank()
  )

g3 <- ggplot(RRdata) +
  geom_point(aes(forward, population), shape=1) +
  geom_abline(lty=2) +
  scale_x_continuous("initial forward") +
  scale_y_continuous("population correction") +
  theme(
    panel.grid = element_blank()
  )

g4 <- ggplot(RRdata) +
  geom_point(aes(forward, empirical), shape=1) +
  geom_abline(lty=2) +
  scale_x_continuous("initial forward") +
  scale_y_continuous("empirical") +
  theme(
    panel.grid = element_blank()
  )

gtot <- arrangeGrob(ggR, g2, g3, g4, layout_matrix = matrix(c(1, 1, 1, 2, 3, 4), nrow=3), widths=c(2, 1))

ggsave("cmp_reproductive.pdf", ggR, width=6, height=3)
ggsave("cmp_reproductive2.pdf", gtot, width=10, height=6)
