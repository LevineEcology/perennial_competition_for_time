
weather <- read.csv("weather.csv")
flux <- read.csv("flux.csv")

for (i in 1:nrow(flux)) {

  flux[i, "time"] <- strsplit(flux[i, "datetime"], "T")[[1]][2]
  flux[i, "date"] <- strsplit(flux[i, "datetime"], "T")[[1]][1]

}


flux <- flux[flux$time == "13:00",]
flux$date <- as.Date(flux$date)
flux.2014 <- flux[flux$year == 2015,]

weather$date <- as.Date(weather$date)
weather.2015 <- weather[weather$date < as.Date("2016-01-01") & weather$date > as.Date("2014-12-31"),]
plot(flux.2014$date, flux.2014$net.rad)
points(weather.2015$date, weather.2015$prec, pch = ".", color = "red")
