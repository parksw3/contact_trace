library(deSolve)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

intrinsic_fun <- function(tau, parm) {
	with(as.list(parm), {
		sigma*gamma/(gamma-sigma) * (exp(-sigma*tau)-exp(-gamma*tau))
	})
}

## ode

measles.ode <- function(t, y, parms) {
	with(as.list(c(y, parms)), {
		
		beta <- b0 * (1 + b1 * cos(2 * pi * t/period))
		
		incidence <- beta * S * I/N
		
		dS <- mu * (N-S) - incidence
		dE <- incidence - (sigma+mu)*E
		dI <- sigma*E - (gamma+mu)*I
		dR <- (gamma+mu)*I - mu*R
		list(c(dS, dE, dI, dR))
	})
}

sumfun <- function(out,
				   tmax,
				   t,
				   parm,
				   tauby=5) {
	incidence <- out$incidence
	c.incidence <- cumsum(incidence)
	
	genden <- intrinsic_fun(t, parm)
	
	censor <- function(tau, tmax) {
		m <- which(t==tmax)
		total <- c.incidence[m]
		
		i <- which(t==tau)
		
		sum(incidence[1:(m-i+1)] * out$S[i:m]/parm[["N"]])*genden[i]
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
	
	censor.df
}

parms <- c(
	b0=1250/365,
	b1=0.15,
	sigma=1/8,
	gamma=1/5,
	mu=1/50/365,
	period=365,
	N=5e7
)

y <- c(S=0.05, E=0.0001, I=0.0001, R=1-0.05-0.0001-0.0001) * 5e7

dt <- 1
dt2 <- 0.2
tmax <- 43 * 365
tmax2 <- 4 * 365
t <- seq(0, tmax, dt)
t2 <- seq(0, tmax2, dt2)
out1 <- as.data.frame(rk(y, t, measles.ode, parms))

y2 <- unlist(tail(out1, 1)[-1])
out2 <- as.data.frame(rk(y2, t2, measles.ode, parms))

out2$incidence <- with(as.list(parms), {
	b0 * (1 + b1 * cos(2 * pi * t2/period)) * out2$S * out2$I/N
})

gensum <- sumfun(out2, tmax2, t2, parms, tauby=10)

g1 <- ggplot(out2, aes(time, S)) +
	geom_line(lwd=1) +
	scale_y_continuous("susceptible") +
	scale_x_continuous(expand=c(0,0)) +
	ggtitle("number of susceptible") +
	theme(
		plot.title = element_text(hjust = 0.5),
		panel.grid=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		plot.margin=unit(c(5.5, 10, 5.5, 5.5), "points"),
		strip.background = element_blank(),
		panel.spacing = unit(0, "cm")
	)

g2 <- ggplot(out2, aes(time, incidence)) +
	geom_line(lwd=1) +
	scale_y_continuous("daily incidence") +
	scale_x_continuous(expand=c(0,0)) +
	ggtitle("daily incidence") +
	theme(
		plot.title = element_text(hjust = 0.5),
		panel.grid=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		plot.margin=unit(c(5.5, 10, 5.5, 15.5), "points"),
		strip.background = element_blank(),
		panel.spacing = unit(0, "cm")
	)

g3 <- ggplot(gensum, aes(t, mean)) +
	geom_line(lwd=1) +
	geom_hline(yintercept=13, lty=2) +
	scale_y_continuous("mean observed GI", limit=c(0, 16), expand=c(0,0)) +
	scale_x_continuous("time (days)", expand=c(0,0)) +
	ggtitle("mean observed GI") +
	theme(
		plot.title = element_text(hjust = 0.5),
		panel.grid=element_blank(),
		axis.title.y=element_blank(),
		plot.margin=unit(c(5.5, 10, 5.5, 30.5), "points"),
		strip.background = element_blank(),
		panel.spacing = unit(0, "cm")
	)
	

gg_temporal <- arrangeGrob(g1, g2, g3, nrow=3)

ggsave("measles.pdf", gg_temporal, width=6, height=4)
