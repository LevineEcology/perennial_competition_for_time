library(raster)

aspect <- readRDS("aspect.rds")
aspect[aspect > 180] <- abs(aspect[aspect > 180] - 360)
writeRaster(aspect, "aspect.tif", format = "GeoTiff", overwrite = TRUE)

dem <- readRDS("dem.rds")
writeRaster(dem, "dem.tif", format = "GeoTiff", overwrite = TRUE)
