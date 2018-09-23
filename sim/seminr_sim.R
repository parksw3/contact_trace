source("../sim/ebola_param_seminr.R")
source("../R/seminr.R")
source("../R/generation.R")

nsim <- 100

set.seed(101)
reslist <- vector('list', nsim)
i <- 1
while (i <= nsim) {
    print(i)
    rr <- seminr.full(N, beta/(N-1), sigma, m.sigma, gamma, n.gamma, I0 = 10, tmax=7*15)
    if(nrow(rr$data) > 100) {
        reslist[[i]] <- rr
        i <- i+1
    }
}

save("reslist", file="seminr_sim.rda")
