set.seed(101)
beta <- 2/5
sigma <- 1/11.4
m.sigma <- round(1.75)
gamma <- 1/5
n.gamma <- round((5/4.7)^2)
N <- 100000

true.mean <- 1/sigma + 1/gamma

nsamp <- 1e5

lat <- rgamma(nsamp, m.sigma, sigma*m.sigma)
inf <- rgamma(nsamp, n.gamma, gamma*n.gamma)

gen <- lat + runif(nsamp, 0, sample(inf, prob=inf, replace=TRUE)) 

true.var <- var(gen)

true.shape <- true.mean^2/true.var

true.R <- beta/gamma