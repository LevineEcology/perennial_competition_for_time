library(ggplot2)
library(lme4)
library(brms)
library(reticulate)
setwd("~/Documents/Science/Princeton/water_comp_perennials/meta_analysis/")

data <- read.csv("chase_data_div_cor.csv")

cdsapi <- import("cdsapi")
source_python("reanalysis_pull.py")
pull_cds_data(c('2006', '2007'), '1', '1', 'test.nc')

for (study in unique(data$study)) {

  for (site in unique(data[data$study == study, "site"])) {

  start_date <- as.Date(data[data$study == study, "treatment_start_date"][1], format = "%m/%d/%y")
  end_date <- as.Date(data[data$study == study, "treatment_end_date"][1], format = "%m/%d/%y")
  start_year <- format(start_date,"%Y")
  end_year <- format(end_date,"%Y")
  years <- as.character(seq(as.numeric(start_year), as.numeric(end_year)))

  pull_cds_data(years,
                data[data$study == study & data$site == site, "latitude"][1],
                data[data$study == study & data$site == site, "longitude"][1],
                paste0(study, "_", site, ".nc"))

  }

}


study <- unique(data$study)[1]

start_date <- as.Date(data[data$study == study, "treatment_start_date"][1], format = "%m/%d/%y")
end_date <- as.Date(data[data$study == study, "treatment_start_date"][2], format = "%m/%d/%y")

start_year <- format(start_date,"%Y")
end_year <- format(start_date,"%Y")

years <- as.character(seq(as.numeric(start_year), as.numeric(end_year)))


as.Date(start_date, format = "%m/%d/%y")
