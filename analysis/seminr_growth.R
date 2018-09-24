library(epigrowthfit)

load("../sim/seminr_sim.rda")

weeks <- 3:5*3

fitlist <- vector('list', length(weeks))

for (j in 1:length(weeks)) {
	
	sublist <- vector('list', length(reslist))
	
	for (i in 1:length(reslist)) {
		print(i)
		sim <- reslist[[i]]
		
		week <- floor(sim$data$time/7) + 1
		
		fitdata <- data.frame(
			week=1:15,
			cases=as.vector(table(head(week, -1)))
		)
		
		fitdata <- fitdata[1:weeks[j],]
		
		sublist[[i]] <- epigrowthfit(fitdata$week, fitdata$cases,
									 model="logistic")
	}
	
	fitlist[[j]] <- sublist
}

save(fitlist, file="seminr_growth.rda")
