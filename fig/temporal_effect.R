library(deSolve)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

source("../sim/ebola_param.R")

## ode

seir.ode <- function(t, y, parms) {
    with(as.list(c(y, parms)), {
        dS <- -beta/N*S*I
        dE <- beta/N*S*I - sigma*E
        dI <- sigma*E - gamma*I
        list(c(dS, dE, dI))
    })
}

sumfun <- function(out,
				   tmax,
				   t,
				   tauby=5) {
	incidence <- beta*out$S*out$I/N
	c.incidence <- cumsum(incidence)
	
	out$incidence <- incidence
	
	genden <- intrinsic_fun(t)
	
	censor <- function(tau, tmax) {
		m <- which(t==tmax)
		total <- c.incidence[m]
		
		i <- which(t==tau)
		
		beta/gamma*sum(incidence[1:(m-i+1)] * out$S[i:m]/N)/total*genden[i]
	}
	
	backward <- function(tau, tmax) {
		m <- which(zapsmall((tmax-tau)-t)==0)
		i <- which(t==tau)
		
		incidence[m] * genden[i]
	}
	
	mean.censor <- function(tmax) {
		if (tmax<=0){
			return(0)
		} else{
			t <- t[t <= tmax]
			
			dc <- sapply(t, censor, tmax=tmax)
			return(sum(dc*t)/sum(dc))
		}
	}
	
	mean.backward <- function(tmax) {
		if (tmax<=0){
			return(0)
		} else{
			t <- t[t <= tmax]
			
			dc <- sapply(t, backward, tmax=tmax)
			return(sum(dc*t)/sum(dc))
		}
	}
	
	tau <- seq(0, tmax, tauby)
	cm <- sapply(tau, mean.censor)
	bm <- sapply(tau, mean.backward)
	
	censor.df <- data.frame(
		t=tau
		, aggregated=cm
		, backward=bm
	)
}

parms <- c(beta=beta, sigma=sigma, gamma=gamma, N=N)
y <- c(S = N-10, E=0, I=10)

dt <- 0.1
tmax <- 305
t <- seq(0, tmax, dt)
out1 <- as.data.frame(ode(y, t, seir.ode, parms))
names(out1)[1] <- "t"

gensum <- sumfun(out1, tmax=tmax, t=t, tauby=5)

outdata <- out1 %>%
	mutate(
		incidence=beta*S*I/N
	)

g1 <- ggplot(outdata, aes(t, incidence)) +
    geom_line(lwd=1.2) +
    scale_x_continuous(expand=c(0,0), limits=c(0, 300), breaks=seq(0, 250, by=50)) +
    scale_y_continuous("daily incidence", expand=c(0,0), limits=c(-5, 700)) +
	# facet_wrap(~type) +
    theme(
        panel.grid=element_blank(),
        axis.title.x=element_blank(),
        plot.margin=unit(c(5.5, 10, 5.5, 5.5), "points"),
        strip.background = element_blank(),
        panel.spacing = unit(0, "cm")
    )

gendata <- gensum %>%
	gather(key, value, -t)

g2 <- ggplot(gendata, aes(t, value, lty=key)) +
    geom_line(lwd=1.2) +
    geom_hline(yintercept = 1/sigma + 1/gamma, lty=3) +
    scale_x_continuous("time (days)", expand=c(0,0), limits=c(0, 300), breaks=seq(0, 250, by=50)) +
    scale_y_continuous("mean GI (days)", expand=c(0,0), limits=c(0, 30)) +
	scale_linetype_manual(values=c(1, 2)) +
	# facet_wrap(~type) +
	theme(
        panel.grid=element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 8.1, b = 0, l = 0)),
        plot.margin=unit(c(5.5, 10, 5.5, 5.5), "points"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(0, "cm"),
        legend.title = element_blank(),
        legend.position = c(0.16, 0.8),
        legend.key.width = grid::unit(2, "cm")
    )

gg_temporal <- arrangeGrob(g1, g2, nrow=2, heights=c(0.4, 0.6))

ggsave("temporal_effect.pdf", gg_temporal, width=6, height=4)
