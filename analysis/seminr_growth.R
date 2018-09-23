library(epigrowthfit)

load("../sim/seminr_sim.rda")

fitlist <- vector('list', length(reslist))

for (i in 1:length(reslist)) {
	print(i)
	sim <- reslist[[i]]

	week <- floor(sim$data$time/7) + 1
	
	fitdata <- data.frame(
		week=1:15,
		cases=as.vector(table(head(week, -1)))
	)
	
	fitlist[[i]] <- epigrowthfit(fitdata$week, fitdata$cases,
						model="logistic")
}

save(fitlist, file="seminr_growth.rda")
