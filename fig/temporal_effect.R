library(deSolve)
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

## we want to simulate equilibrium state!
seirs.ode <- function(t, y, parms) {
	with(as.list(c(y, parms)), {
		dS <- gamma * I -beta/N*S*I
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
		
		beta*sum(incidence[1:(m-i+1)] * out$S[i:m]/N)/total*genden[i]
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
	
	tau <- seq(0, tmax, tauby)
	cm <- sapply(tau, mean.censor)
	
	censor.df <- data.frame(
		t=tau
		, mean=cm
	)
}

parms <- c(beta=beta, sigma=sigma, gamma=gamma, N=N)
y <- c(S = N-10, E=0, I=10)

dt <- 0.1
tmax <- 305
t <- seq(0, tmax, dt)
out1 <- as.data.frame(ode(y, t, seir.ode, parms))
out2 <- as.data.frame(ode(y, t, seirs.ode, parms))
out3 <- as.data.frame(ode(unlist(tail(out2, 1))[-1] , t, seirs.ode, parms))
out3 <- as.data.frame(ode(unlist(tail(out3, 1))[-1] , t, seirs.ode, parms))
names(out1)[1] <- names(out2)[1] <- names(out3)[1] <- "t"

outlist <- list(out1, out3)

gensum <- lapply(outlist, sumfun, tmax=305, t=t, tauby=5)

# gensum_long <- sumfun(out4, tmax=10000, t=0:10000, tauby=1000)

outdata <- bind_rows(out1, out3, .id="type") %>%
	mutate(
		incidence=beta*S*I/N,
		type=factor(type, label=c("epidemic", "equilibrium"))
	)

## We probably don't need to show the endemic case.. seems trivial...

g1 <- ggplot(outdata %>% filter(type=="epidemic"), aes(t, incidence)) +
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
	bind_rows(.id="type")

g2 <- ggplot(gendata %>% filter(type==1), aes(t, mean)) +
    geom_line(lwd=1.2) +
    geom_hline(yintercept = 1/sigma + 1/gamma, lty=2) +
    scale_x_continuous("time (days)", expand=c(0,0), limits=c(0, 300), breaks=seq(0, 250, by=50)) +
    scale_y_continuous("mean observed GI", expand=c(0,0), limits=c(0, 1/sigma+1/gamma+4)) +
	# facet_wrap(~type) +
	theme(
        panel.grid=element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 8.1, b = 0, l = 0)),
        plot.margin=unit(c(5.5, 10, 5.5, 5.5), "points"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(0, "cm")
    )

gg_temporal <- arrangeGrob(g1, g2, nrow=2)

ggsave("temporal_effect.pdf", gg_temporal, width=6, height=4)
