##' @param x a vector of nodes
##' @param size number of nodes to pick at random
sample2 <- function(x, size) {
    if(length(x)==1) {
        rep(x, size)
    } else{
        sample(x, size, replace=TRUE)
    }
}

##' Simulate an epidemic until done or imax reached
##' @param g igraph object
##' @param beta contact rate
##' @param sigma 1/(latent period)
##' @param gamma 1/(infectious period)
##' @param imax maximum number of infected individuals
seir <- function(g,
                 beta,
                 sigma,
                 gamma,
                 I0,
                 initial_infected,
                 seed=NULL,
                 imax,
                 keep.intrinsic=FALSE){
    if(class(g) != "igraph"){
        stop("g must be an igraph object")
    }
    
    if (!is.null(seed)) set.seed(seed)
    
    V <- as.vector(V(g))
    
    if (missing(initial_infected)) {
        if (missing(I0)) stop("specify the initial conditions")
        
        initial_infected <- sample(V, I0)    
    }

    I0 <- length(initial_infected)
    
    if (missing(imax)) imax <- length(V)
    
    t <- 0
    
    queue_v <- queue_t <- queue_infector <- rep(NA, I0)
    
    queue_v[1:I0] <- initial_infected 
    queue_t[1:I0] <- 0
    
    t_infected <- t_infectious <- t_recovered <- rep(NA, length(V))
    t_infected[initial_infected] <- 0
    
    t_gillespie <- NULL
    c_infected <- 0
    
    if (keep.intrinsic) {
        intrinsic_generation <- vector('list', length(V))
    } else {
        intrinsic_generation <- NULL
    }
    
    done <- rep(FALSE, length(V))
    infected_by <- rep(NA, length(V))
    
    stop <- FALSE
    
    while (!stop) {
        j.index <- which.min(queue_t)
        j <- queue_v[j.index]
        
        infected_by[j] <- queue_infector[j.index]
        t_infected[j] <- queue_t[j.index]
        
        t <- queue_t[j.index]; t_gillespie <- c(t_gillespie, t)
        
        latent <- rexp(1, sigma)
        t_infectious[j] <- t + latent
        
        c_infected <- c_infected +1
        
        n <- as.vector(neighbors(g, j))
        if (!keep.intrinsic) n <- n[!done[n]]
        
        rate <- length(n) * beta + gamma
        
        prob <- gamma/rate
        
        ncontact <- rnbinom(1, size=1, prob=prob)
        
        time_between <- rexp(ncontact+1, rate=rate)
        c_time <- latent + cumsum(time_between)
        if (keep.intrinsic) intrinsic_generation[[j]] <- c_time[1:ncontact]
        
        if (ncontact > 0) {
            contact <- sample2(n, ncontact)
            filter <- which(!duplicated(contact))
            
            queue_v <- c(queue_v, contact[filter])
            queue_infector <- c(queue_infector, rep(j, length(contact[filter])))
            queue_t <- c(queue_t, t + c_time[filter])
        }
        
        t_recovered[j] <- t+tail(c_time,1)
        done[j] <- TRUE
        
        filter2 <- !done[queue_v]
        queue_v <- queue_v[filter2]
        queue_infector <- queue_infector[filter2]
        queue_t <- queue_t[filter2]
        
        stop <- (c_infected == length(V) || all(done[queue_v]) || c_infected == imax)
    }
    
    return(
        list(
            data=data.frame(
                time=t_gillespie[(I0):c_infected],
                infected=(I0):c_infected
            ),
            intrinsic_generation=intrinsic_generation,
            t_infected=t_infected,
            t_infectious=t_infectious,
            t_recovered=t_recovered,
            infected_by=infected_by
        )
    )
}

## creating a full graph takes up a lot of memory
## so here's an optimized version that doesn't rely on igraph but does the same thing
seir.full <- function(size,
                      beta,
                      sigma,
                      gamma,
                      I0,
                      seed = NULL,
                      imax,
                      keep.intrinsic=FALSE){
    
    if (!is.null(seed)) set.seed(seed)
    
    V <- 1:size
    
    if (missing(I0)) {
        if (missing(I0)) stop("specify the initial conditions")
        
    }
    initial_infected <- 1:I0
    
    if (missing(imax)) imax <- size
    
    queue_v <- queue_t <- queue_infector <- rep(NA, I0)
    
    queue_v[1:I0] <- initial_infected 
    queue_t[1:I0] <- 0
    
    t_infected <- t_infectious <- t_recovered <- rep(NA, size)
    t_infected[initial_infected] <- 0
    
    t_gillespie <- NULL
    c_infected <- 0
    
    if (keep.intrinsic) {
        intrinsic_generation <- vector('list', length(V))
    } else {
        intrinsic_generation <- NULL
    }
    
    done <- rep(FALSE, size)
    infected_by <- rep(NA, size)
    
    stop <- FALSE
    
    while (!stop) {
        j.index <- which.min(queue_t)
        j <- queue_v[j.index]
        
        infected_by[j] <- queue_infector[j.index]
        t_infected[j] <- queue_t[j.index]
        
        t <- queue_t[j.index]; t_gillespie <- c(t_gillespie, t)
        
        latent <- rexp(1, sigma)
        t_infectious[j] <- t + latent
        
        c_infected <- c_infected +1
        
        n <- V[V != j]
        if (!keep.intrinsic) n <- n[!done[n]]
        
        rate <- length(n) * beta + gamma
        
        prob <- gamma/rate
        
        ncontact <- rnbinom(1, size=1, prob=prob)
        if (ncontact > 0) {
            queue_v <- c(queue_v, sample2(n, ncontact))
            queue_infector <- c(queue_infector, rep(j, ncontact))
        }
        
        time_between <- rexp(ncontact+1, rate=rate)
        c_time <- latent + cumsum(time_between)
        
        if (keep.intrinsic) intrinsic_generation[[j]] <- c_time[1:ncontact]
        
        if (ncontact > 0) {
            queue_t <- c(queue_t, t + c_time[1:ncontact])
        }
        
        t_recovered[j] <- t+tail(c_time,1)
        done[j] <- TRUE
        
        filter2 <- !done[queue_v]
        queue_v <- queue_v[filter2]
        queue_infector <- queue_infector[filter2]
        queue_t <- queue_t[filter2]
        
        stop <- (c_infected == length(V) || all(done[queue_v]) || c_infected == imax)
    }
    
    return(
        list(
            data=data.frame(
                time=t_gillespie[(I0):c_infected],
                infected=(I0):c_infected
            ),
            intrinsic_generation=intrinsic_generation,
            t_infected=t_infected,
            t_infectious=t_infectious,
            t_recovered=t_recovered,
            infected_by=infected_by
        )
    )
}

##' SEIR on a local tree (one node connected to n branches)
##' @param n number of branches
seir.local <- function(n,
                       beta,
                       sigma,
                       gamma,
                       seed = NULL){
    latent <- rexp(1, sigma)
    
    rate <- n * beta + gamma
    
    prob <- gamma/rate
    
    ncontact <- rnbinom(1, size=1, prob=prob)
    time_between <- rexp(ncontact+1, rate=rate)
    c_time <- latent + cumsum(time_between)
    
    v <- sample2(1:n, ncontact)
    list(
    	conditional=c(head(c_time, -1)[which(!duplicated(v))]),
    	intrinsic=head(c_time,-1)
    )
}
