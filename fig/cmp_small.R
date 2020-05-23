library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

source("../sim/ebola_param_seminr.R")

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

load("../analysis/cmp_compare_small.rda")

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

ggsave("cmp_reproductive_small.pdf", ggR, width=6, height=3)
