intrinsic.generation <- function(x, plot=TRUE, ...) {
    ig <- unlist(x$intrinsic_generation)
    if(plot) hist(ig, freq=FALSE, ...)
    invisible(ig)
}

network.generation <- function(x, plot=TRUE,
                               type=c("backward", "forward"),
                               interval,
                               interval.type=c("time", "cases"),
                               ...) {
    type <- match.arg(type)
    interval.type <- match.arg(interval.type)
    
    if (missing(interval)) {
        interval <- c(0, max(x$data$time))
    } else if (interval.type=="cases") {
        interval <- x$data$time[match(interval, x$data$infected)]
    }
    
    cohort <- which(x$t_infected >= interval[1] & x$t_infected <= interval[2])
    
    if (type == "backward") {
        infector <- x$infected_by[cohort]
        ng <- x$t_infected[cohort] -  x$t_infected[infector]
        ng <- ng[!is.na(ng)]
    } else {
        ng <- unlist(lapply(cohort, function(z) {
            infectee <- which(z==x$infected_by)
            x$t_infected[infectee] - x$t_infected[z]
        }))
    }
    
    if (plot) hist(ng, freq=FALSE, ...)
    
    invisible(ng)
}

generation.data <- function(x,
                            tmax,
                            method=c("backward", "forward")
                            ) {
    method <- match.arg(method)
    
    if (method=="backward") {
        index <- which(x$t_infected < tmax)
        
        index <- index[order(x$t_infected[index])]
        
        data.frame(
            index=index,
            t_infected=x$t_infected[index],
            infected_by=x$infected_by[index],
            generation=x$t_infected[index]-x$t_infected[x$infected_by[index]]
        )
    } else {
        index <- which(x$t_infected < tmax)
        
        index <- index[order(x$t_infected[index])]
        
        data.frame(
            index=x$infected_by[index],
            t_infected=x$t_infected[x$infected_by[index]],
            infected=index,
            generation=x$t_infected[index]-x$t_infected[x$infected_by[index]]
        )
    }
}
