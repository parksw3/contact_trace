library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

source("../sim/ebola_param_seminr.R")

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

genlist <- list()
RRlist <- list()

load("../analysis/cmp_compare.rda")

genlist$big <- gendata
RRlist$big <- RRdata

load("../analysis/cmp_compare_2.rda")

genlist$small <- gendata
RRlist$small <- RRdata

RRgroup <- data.frame(
	key=c("contact\ntracing", "population\ncorrection", "individual\ncorrection", "tree\nbased", "empirical", "local\ncorrection", "intrinsic"),
	group=factor(c("tracing based", "tracing based", "tracing based", "empirical", "empirical", "individual based", "individual based"),
				 levels=c("tracing based", "empirical", "individual based"))
)

RRdata2 <- RRlist %>%
	bind_rows(.id="type") %>%
	gather(key, value, -type) %>%
	mutate(key=factor(key, 
					  levels=c("observed", "population", "individual", "trapman", "empirical", "network", "intrinsic"),
					  labels=c("contact\ntracing", "population\ncorrection", "individual\ncorrection", 
					  		 "tree\nbased", "empirical", "local\ncorrection", "intrinsic"))) %>%
	merge(RRgroup) %>%
	filter(key != "tree\nbased") %>%
	mutate(type=factor(type, levels=c("small", "big"), labels=c("Smaller epidemic", "Larger epidemic")))

gendata2 <- genlist %>%
	bind_rows(.id="type") %>%
	gather(key, value, -type) %>%
	mutate(key=factor(key, 
					  levels=c("observed", "population", "individual", "trapman", "empirical", "network", "intrinsic"),
					  labels=c("contact\ntracing", "population\ncorrection", "individual\ncorrection", 
					  		 "tree\nbased", "empirical", "local\ncorrection", "intrinsic"))) %>%
	merge(RRgroup) %>%
	filter(key != "tree\nbased")

ggR <- ggplot(RRdata2 %>% filter(type=="Larger epidemic")) +
	geom_boxplot(aes(key, value, fill=key), alpha=0.7) +
	ylab("Reproductive number") +
	facet_grid(~group, scale="free", space="free_x") +
	theme(
		panel.spacing.x = unit(0, "lines"),
		strip.background = element_blank(),
		legend.position = "none",
		axis.title.x=element_blank()
	)

ggsave("cmp_reproductive.pdf", ggR, width=6, height=3)
