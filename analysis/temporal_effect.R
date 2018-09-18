source("ebola_param_seminr.R")
library(deSolve)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

intrinsic_fun <- function(tau) {
    sigma*gamma/(gamma-sigma) * (exp(-sigma*tau)-exp(-gamma*tau))
}

## ode

seir.ode <- function(t, y, parms) {
    with(as.list(c(y, parms)), {
        dS <- -beta/N*S*I
        dE1 <- beta/N*S*I - m.sigma*sigma*E1
        dE2 <- m.sigma*sigma*E1 - m.sigma*sigma*E2
        dI <- m.sigma*sigma*E2 - gamma*I
        list(c(dS, dE1, dE2, dI))
    })
}

parms <- c(beta=beta, sigma=sigma, m.sigma=m.sigma, gamma=gamma, N=N)
y <- c(S = N-1, E1=0, E2=0, I=10)

dt <- 1
tmax <- 300
t <- seq(0, tmax, dt)
out <- as.data.frame(ode(y, t, seir.ode, parms))
names(out)[1] <- "t"

incidence <- beta*out$S*out$I/N
c.incidence <- cumsum(incidence)

out$incidence <- incidence

genden <- hist(gen, breaks=t, plot=FALSE)$density

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

tau <- seq(0, tmax, 5)
cm <- sapply(tau, mean.censor)

censor.df <- data.frame(
    t=tau
    , mean=cm
)

exp.df <- data.frame(
    t=t[which(t==100)],
    incidence=incidence[which(t==100)],
    mean=cm[which(tau==100)]
)

out.shade <- rbind(out[out$t <= 100, c("t", "incidence")], c(100, 0), c(0,0))
colnames(out.shade)[1] <- "t"

g1 <- ggplot(out, aes(t, incidence)) +
    geom_line(lwd=1.2) +
    geom_polygon(data=out.shade, alpha=0.2) +
    geom_segment(data=data.frame(t=100, incidence=0, tend=100, iend=out$incidence[which(t==100)]), aes(xend=tend, yend=iend),
                 lty=1) +
    geom_point(data=exp.df, size=4) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous("daily incidence", expand=c(0,0), limits=c(0, 650)) +
    theme(
        panel.grid=element_blank(),
        axis.title.x=element_blank(),
        plot.margin=unit(c(5.5, 10, 5.5, 5.5), "points")
    )

g2 <- ggplot(censor.df, aes(t, mean)) +
    geom_line(lwd=1.2) +
    geom_point(data=exp.df, size=4) +
    geom_segment(data=data.frame(t=100, mean=0, tend=100, iend=censor.df$mean[which(censor.df$t==100)]), aes(xend=tend, yend=iend),
                 lty=1) +
    geom_hline(yintercept = 1/sigma + 1/gamma, lty=2) +
    scale_x_continuous("time (days)", expand=c(0,0)) +
    scale_y_continuous("mean observed GI", expand=c(0,0), limits=c(0, 1/sigma+1/gamma+1)) +
    theme(
        panel.grid=element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.margin=unit(c(5.5, 10, 5.5, 5.5), "points")
    )

gg_temporal <- arrangeGrob(g1, g2)

ggsave("temporal_effect.pdf", gg_temporal, width=6, height=4)
