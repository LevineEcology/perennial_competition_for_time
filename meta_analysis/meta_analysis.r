
library(ggplot2)
library(tidyverse)
library(brms)
library(sf)

source("../SedgwickWaterCompetition/code/utility/plot_utility.R")

data <- read.csv("meta_analysis/meta_analysis_timing_trim.csv")
data <- data[1:91,]

for (i in which(data$first_author == "zhang")) {

  data[i, "uncertainty"] <- (((as.numeric(strsplit("1.05;2.05", ";")[[1]])[2] - as.numeric(strsplit("1.05;2.05", ";")[[1]])[1]) / 50 ) * 42.5) / 1.96

}
data$uncertainty <- as.numeric(data$uncertainty)


for (i in unique(data$unique_comp_id)) {

  data[data$unique_comp_id == i, "diversity_change"] <- (log(data[data$unique_comp_id == i, "diversity"]) -
    log(data[data$unique_comp_id == i & data$timing_trt_pct == 0.0, "diversity"]))

}

for (i in unique(data$unique_comp_id)) {

  data[data$unique_comp_id == i, "day_dif"] <- (data[data$unique_comp_id == i, "timing_trt_days"] -
    data[data$unique_comp_id == i & data$timing_trt_days == min(data[data$unique_comp_id == i, "timing_trt_days"]), "timing_trt_days"])

}

data[data$uncertainty == 0, "uncertainty"] <- 0.01


data$first_author <- factor(data$first_author)
contrasts(data$first_author) <- "contr.sum"

data$treatment <- "control"
data[data$timing_trt_pct != 0, "treatment"] <- "increase_time"

meta_trt <- brm(diversity | se(uncertainty) ~ treatment + (1|first_author),
                data = data,
                chains = 2,
                cores = 2,
                iter = 4000,
                warmup = 2000)
summary(meta_trt)
conditional_effects(meta_trt, points = TRUE, re_formula = NA)


meta <- brm(diversity_change ~ timing_trt_pct + (1|first_author),
            data = data[data$treatment == "increase_time",],
            chains = 2,
            cores = 2,
            iter = 4000,
            warmup = 2000)
summary(meta)
conditional_effects(meta)

plot(data[data$treatment == "increase_time", "timing_trt_days"], data[data$treatment == "increase_time", "diversity_change"])

meta <- brm(diversity | se(uncertainty) ~ day_dif + (1|first_author),
                data = data,
                chains = 2,
                cores = 2,
                iter = 4000,
                warmup = 2000)
summary(meta)

newdata <- data.frame(day_dif = seq(min(data$day_dif), max(data$day_dif), length.out = 40),
                      year = rep(4, times = 40), uncertainty = rep(0.1, times = 40))
newdata$pred <- predict(meta, newdata, re_formula = NA)[,1]
newdata$lower <- predict(meta, newdata, re_formula = NA)[,3]
newdata$upper <- predict(meta, newdata, re_formula = NA)[,4]

data$upper <- data$diversity + data$uncertainty
data$lower <- data$diversity - data$uncertainty

for (i in 1:nrow(data)) {

  data[i, "study_name"] <- tools::toTitleCase(paste0(data[i, "first_author"], " ", data[i, "publ_year"]))

}

second_axis(ggplot(data = data, aes(x = day_dif)) +
  geom_point(aes(color = study_name, y = diversity),
             size = 4, position = position_dodge(width = 1)) +
  #geom_errorbar(aes(color = first_author, ymin = lower, ymax = upper),
  #              position = position_dodge(width = 1)) +
  geom_line(data = newdata, aes(x = day_dif, y = pred), linewidth = 2) +
  geom_ribbon(data = newdata, aes(x = day_dif, ymin = lower, ymax = upper),
              fill = "gray", alpha = 0.3) +
  scale_color_brewer(type = "qual", palette = 3) +
  scale_x_continuous(expand = c(0.01,0), limits = c(NA, NA)) +
  scale_y_continuous(expand = c(0,0), limits = c(NA, 15)) +
  scale_color_discrete(name = "Study") +
  ylab("Species Richness") +
  xlab("Change in Inter-Rain Interval Length") +
  theme_jabo())



meta <- brm(diversity | se(uncertainty) ~ timing_trt_pct + (1|first_author),
                data = data,
                chains = 2,
                cores = 2,
                iter = 4000,
                warmup = 2000)
summary(meta)


newdata <- data.frame(timing_trt_pct = seq(min(data$timing_trt_pct), max(data$timing_trt_pct), length.out = 40),
                      year = rep(4, times = 40), uncertainty = rep(0.1, times = 40))
newdata$pred <- predict(meta, newdata, re_formula = NA)[,1]
newdata$lower <- predict(meta, newdata, re_formula = NA)[,3]
newdata$upper <- predict(meta, newdata, re_formula = NA)[,4]

data$upper <- data$diversity + data$uncertainty
data$lower <- data$diversity - data$uncertainty

second_axis(ggplot(data = data, aes(x = timing_trt_pct)) +
  geom_point(aes(color = study_name, y = diversity),
             size = 4, position = position_dodge(width = 1)) +
  #geom_errorbar(aes(color = first_author, ymin = lower, ymax = upper),
  #              position = position_dodge(width = 1)) +
  geom_line(data = newdata, aes(x = timing_trt_pct, y = pred), linewidth = 2) +
  geom_ribbon(data = newdata, aes(x = timing_trt_pct, ymin = lower, ymax = upper),
              fill = "gray", alpha = 0.3) +
  scale_color_brewer(type = "qual", palette = 3) +
  scale_x_continuous(expand = c(0.0,0), limits = c(-0.8, 4.9), labels = scales::percent) +
  scale_y_continuous(expand = c(0,0), limits = c(NA, 15)) +
  scale_color_discrete(name = "Study") +
  ylab("Species Richness") +
  xlab("Change in Inter-Rain Interval Length") +
  theme_jabo())





ggplot(data = data, aes(x = treatment, y = diversity)) +
  geom_jitter(aes(color = first_author)) +
  geom_boxplot()


ggplot(data = data, aes(x = year, y = diversity)) +
  geom_jitter(aes(color = first_author))


locations <- st_as_sf(data, coords = c("lon", "lat"))

ggplot(locations) +
  geom_map(data = map_data("world"), map = map_data("world"), aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.1) +
  geom_sf(aes(color = study_name), size = 3) +
  scale_color_brewer(type = "qual", palette = 3)
