beta <- 2/5
sigma <- 1/11.4
gamma <- 1/5
N <- 60000
lambda <- beta/(N-1)

true.mean <- 1/sigma + 1/gamma

intrinsic_fun <- function(tau) {
	sigma*gamma/(gamma-sigma) * (exp(-sigma*tau)-exp(-gamma*tau))
}

conditional_fun <- function(tau) {
	sigma*(gamma+lambda)/(gamma+lambda-sigma) * (exp(-sigma*tau)-exp(-(gamma+lambda)*tau))
}

true.var <- integrate(function(x) (x-true.mean)^2 * intrinsic_fun(x), 0, Inf)[[1]]

true.cv <- sqrt(true.var)/true.mean

true.R <- beta/gamma