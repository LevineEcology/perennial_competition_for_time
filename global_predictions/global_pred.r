
library(ncdf4)
library(RColorBrewer)
library(raster)
library(ggplot2)

present_data <- nc_open("~/Downloads/pr_3hr_GFDL-CM3_rcp45_r1i1p1_2016010100-2020123123.nc")
present_data

lat <- ncvar_get(present_data, "lat")
lon <- ncvar_get(present_data, "lon")
time <- ncvar_get(present_data, "time")
tunits <- ncatt_get(present_data, "time", "units")
fillvalue <- ncatt_get(present_data, "pr", "_FillValue")$value

precip <- ncvar_get(present_data, "pr")
precip[precip == fillvalue] <- NA
precip <- precip*10800
precip[precip < 0.3] <- 0

output <- precip[,,1]
output[] <- NA
freq <- output
sd_freq <- output
size <- output
sd_size <- output

for (i in 1:nrow(output)) {

  for (j in 1:ncol(output)) {

    data <- precip[i,j,]
    stormsize <- numeric(0)
    gaps <- numeric(0)
    gap <- 0
    ss <- 0

    for (k in 1:length(data)) {

      if (k == 1) {
        if (data[k] == 0) gaps <- gaps + 0.125
        else ss <- ss + data[k]
      }
      else {

        if (data[k] > 0 & data[k-1] == 0) {
          ss <- ss + data[k]
          gaps <- c(gaps, gap)
          gap <- 0
        }
        else if (data[k] == 0 & data[k-1] > 0) {
          gap <- gap + 0.125
          stormsize <- c(stormsize, ss)
          ss <- 0
        }
        else if (data[k] == 0) gap <- gap + 0.125
        else if (data[k] > 0) ss <- ss + data[k]

      }

    }

    freq[i,j] <- mean(gaps)
    sd_freq[i,j] <- sd(gaps)
    size[i,j] <- mean(data[data > 0])
    sd_size[i,j] <- sd(data[data > 0])

  }

}

future_data <- nc_open("~/Downloads/pr_3hr_GFDL-CM3_rcp45_r1i1p1_2096010100-2100123123.nc")

lat <- ncvar_get(future_data, "lat")
lon <- ncvar_get(future_data, "lon")
time <- ncvar_get(future_data, "time")
tunits <- ncatt_get(future_data, "time", "units")
fillvalue <- ncatt_get(future_data, "pr", "_FillValue")$value

precip <- ncvar_get(future_data, "pr")
precip[precip == fillvalue] <- NA
precip <- precip*10800
precip[precip < 0.05] <- 0

output <- precip[,,1]
output[] <- NA
future_freq <- output
future_sd_freq <- output
future_size <- output
future_sd_size <- output

for (i in 1:nrow(output)) {

  for (j in 1:ncol(output)) {

    data <- precip[i,j,]
    stormsize <- numeric(0)
    gaps <- numeric(0)
    gap <- 0
    ss <- 0

    for (k in 1:length(data)) {

      if (k == 1) {
        if (data[k] == 0) gaps <- gaps + 0.125
        else ss <- ss + data[k]
      }
      else {

        if (data[k] > 0 & data[k-1] == 0) {
          ss <- ss + data[k]
          gaps <- c(gaps, gap)
          gap <- 0
        }
        else if (data[k] == 0 & data[k-1] > 0) {
          gap <- gap + 0.125
          stormsize <- c(stormsize, ss)
          ss <- 0
        }
        else if (data[k] == 0) gap <- gap + 0.125
        else if (data[k] > 0) ss <- ss + data[k]

      }

    }

    future_freq[i,j] <- mean(gaps)
    future_sd_freq[i,j] <- sd(gaps)
    future_size[i,j] <- mean(data[data > 0])
    future_sd_size[i,j] <- sd(data[data > 0])

  }

}

x <- raster(t(log(future_freq)), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat),
            crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
x <- flip(x, direction = 'y')

test_spdf <- as(x, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")

ggplot() +
  geom_tile(data = test_df, aes(x=x, y=y, fill=value)) +
  scale_fill_gradient2()



future_freq <- rbind(future_freq[73:144,],future_freq[1:72,])
freq <- rbind(freq[73:144,],freq[1:72,])

x <- raster(t(log(future_freq/freq)), xmn= -180, xmx= 180, ymn=min(lat), ymx=max(lat),
            crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
x <- flip(x, direction = 'y')

test_spdf <- as(x, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")

countries <- map_data("world")
ggplot() +
  geom_tile(data = test_df, aes(x=x, y=y, fill=value)) +
  geom_polygon(data = countries, aes(x= long, y = lat, group = group), fill = NA, color = "black") +
  scale_fill_gradient2()


future_size <- rbind(future_size[73:144,],future_size[1:72,])
size <- rbind(size[73:144,],size[1:72,])

x <- raster(t(log(future_size/size)), xmn=-180, xmx=180, ymn=min(lat), ymx=max(lat),
            crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
x <- flip(x, direction = 'y')

test_spdf <- as(x, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")

countries <- map_data("world")
ggplot() +
  geom_tile(data = test_df, aes(x=x, y=y, fill=value)) +
  geom_polygon(data = countries, aes(x= long, y = lat, group = group), fill = NA, color = "black") +
  scale_fill_gradient2()
