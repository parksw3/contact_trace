library(dplyr)
library(outbreaks)
library(ggplot2); theme_set(theme_bw())

ll <- mers_korea_2015$linelist
cc <- mers_korea_2015$contacts

toDate <- as.numeric(as.Date(ll$dt_onset[match(cc$to, ll$id)]))
toDate <- toDate - 16569
fromDate <- as.numeric(as.Date(ll$dt_onset[match(cc$from, ll$id)]))
fromDate <- fromDate - 16569

serial_data <- data.frame(
	time=toDate,
	gen=toDate - fromDate,
	within=ll$loc_hosp[match(cc$to, ll$id)] == ll$loc_hosp[match(cc$from, ll$id)],
	hospital=ll$loc_hosp[match(cc$to, ll$id)]
)

boxplot(gen~within, data=serial_data)
t.test(gen~within, data=serial_data)

ggplot(serial_data) +
	geom_boxplot(aes(within, gen))

dd <- data.frame(
	time=ll$dt_onset[match(cc$to, ll$id)],
	gen=toDate - fromDate
) %>%
	arrange(time) %>%
	group_by(time) %>%
	summarize(
		sum=sum(gen),
		weight=length(gen)
	) %>%
	mutate(mean=cumsum(sum)/cumsum(weight))

plot(dd$time, dd$mean)

plot(table(ll$dt_onset))


