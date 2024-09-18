pacman::p_load(tidyverse)

secveg <- read_csv("~/research/data/degradation/secondary_veg_coverage_basins.csv")


secveg %>%
  dplyr::filter(bandName == "age_2021") %>%
  pluck("meanValue") %>%
  #hist(breaks = 1000)
  summary()


secveg2 <- read_csv("~/research/data/degradation/secondary_veg_age_basins.csv")


secveg2 %>%
  dplyr::filter(bandName == "age_2021") %>%
  pluck("meanValue") %>%
  #hist(breaks = 1000)
  summary()
