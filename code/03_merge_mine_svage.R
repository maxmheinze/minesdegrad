
# Header ------------------------------------------------------------------

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("00_prelims.R")

pacman::p_load(
  sf,
  dplyr,
  readr
)

mibas <- read_csv(pd("degradation/downstream_upstream_distance.csv"))
svage <- read_csv(pd("degradation/age.csv"))
svveg <- read_csv(pd("degradation/veg_cleaned.csv"))

# Merge -------------------------------------------------------------------

svage <- svage %>%
  dplyr::filter(bandName != "age_2011") %>%
  mutate(year = as.integer(substr(bandName, 26, 29)))

mbage <- left_join(mibas, svage)

mbveg <- left_join(mbage, dplyr::filter(svveg, defclass == 3))

mbsvg <- mbveg %>%
  dplyr::select(-bandName, -defclass, -share) %>%
  dplyr::rename(mean_svage = meanValue,
                median_svage = medianValue,
                min_svage = minValue,
                max_svage = maxValue,
                svshare = rel_area)


# Export CSV --------------------------------------------------------------

write_csv(mbsvg, pd("degradation/mine_basins_secveg.csv"))  
