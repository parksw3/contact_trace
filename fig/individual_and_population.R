library(ggplot2); theme_set(theme_bw())
library(gridExtra)
source("../sim/ebola_param.R")

lat <- 1/sigma
inf <- 1/gamma

t <- seq(0, 35, by=0.001)

dd_individual <- data.frame(
	time=t,
	density=ifelse(t>lat & t < lat+inf, 2/inf, 0),
	type="Individual-level kernel"
)

dd_population <- data.frame(
	time=t,
	density=2 * intrinsic_fun(t),
	type="Population-level kernel"
)

dd_population_poly <- data.frame(
	time=c(t, 35),
	density=c(2 * intrinsic_fun(t), 0)
)

g1 <- ggplot(dd_individual) +
	geom_line(aes(time, density), lwd=1.2) +
	scale_x_continuous("Time (days)", expand=c(0, 0)) +
	scale_y_continuous(expression(K[a](tau)), expand=c(0,0), limit=c(0, 0.45)) +
	facet_wrap(~type) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank()
	)

g2 <- g1 %+% dd_population

grid.arrange(
	g1 + geom_polygon(aes(time, density), alpha=0.3),
	g2 + geom_polygon(data=dd_population_poly, aes(time, density), alpha=0.3),
	nrow=1
)



