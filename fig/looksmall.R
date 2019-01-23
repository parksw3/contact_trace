library(dplyr)
library(readr)

## line list
ll <- read_csv(input_files[[1]])

## merge 
ml <- (ll
	%>% left_join(
		ll
		, by=c("ID"="IDSource")
		, suffix=c("one", "the_other")
	) %>% mutate(gen = Exposure_dateone-Exposure_datethe_other)
)

print(summary(ml))
print((ml$gen))

quit()

## data
dd <- data.frame(
	time=sp2$exposure,
	gen=sp2$exposure - sp2$exposure[match(sp2$IDSource, sp2$ID)],
	serial=sp2$symptom - sp2$symptom[match(sp2$IDSource, sp2$ID)],
	source=sp2$IDSource,
	recipient=sp2$ID,
	type=sp2$Type
) %>%
	filter(!is.na(gen), gen > 0) %>% ## there is one case where gen < 0 so I get rid of it
	arrange(time)

## All contacts (as opposed to case)
dd_contact <- dd %>%
	group_by(time) %>%
	summarize(
		sum=sum(gen),
		weight=length(gen)
	) %>%
	mutate(
		mean=cumsum(sum)/cumsum(weight)
	)

dd_case <- dd %>%
	filter(type=="Case") %>%
	group_by(time) %>%
	summarize(
		sum=sum(gen),
		weight=length(gen)
	) %>%
	mutate(
		mean=cumsum(sum)/cumsum(weight)
	)

dd_serial <- dd %>%
	filter(type=="Case", !is.na(serial)) %>%
	group_by(time) %>%
	summarize(
		sum=sum(serial),
		weight=length(serial)
	) %>%
	mutate(
		mean=cumsum(sum)/cumsum(weight)
	)

pdf("smallpox.pdf", width=10, height=3)
par(mfrow=c(1, 3))
plot(dd_contact$time, dd_contact$mean, type="l",
	 xlab="Time (days since the first case)",
	 ylab="Mean generation interval",
	 ylim=c(13, 18.5))
title("All contacts")
abline(h=17.7, lty=2)

plot(dd_case$time, dd_case$mean, type="l",
	 xlab="Time (days since the first case)",
	 ylab="Mean generation interval",
	 ylim=c(13, 18.5))
title("Symptomatic cases")
abline(h=17.7, lty=2)

plot(dd_serial$time, dd_serial$mean, type="l",
	 xlab="Time (days since the first case)",
	 ylab="Mean serial interval")
title("Symptomatic cases")
abline(h=17.7, lty=2)
dev.off()
