# Header ------------------------------------------------------------------

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("00_prelims.R")

pacman::p_load(
  dplyr,
  sf,
  readxl,
  readr
)


# Read in Data ------------------------------------------------------------

# outcome
mbsvg <- read_csv(pd("degradation/mine_basins_secveg.csv"))

# treatment
dup <- read_csv(pd("degradation/downstream_upstream_distance_ordered.csv"))

# geographical controls
geo_controls <- readRDS(pd("degradation/geo_data_agg.RDS"))

# meteorological controls
met_controls <- readRDS(pd("degradation/met_data_agg.RDS"))

# population 
pop <- readRDS(pd("degradation/population.RDS"))


# Prepare Data for Regression ---------------------------------------------

df_reg <- full_join(dup, mbsvg) |> relocate(year, .after = iso3c) |> 
  left_join(geo_controls, by = "HYBAS_ID") |>
  left_join(met_controls, by = c("HYBAS_ID", "year")) |>
  left_join(pop, by = c("HYBAS_ID", "year")) |> 
  mutate(t.trend = year - 2000, 
         distance_bin = cut(distance, breaks = c(-Inf, 10, 20, 30, 40, 50, Inf)),
         distance_centroid_bin = cut(distance_centroid, breaks = c(-Inf, 10, 20, 30, 40, 50, Inf)))

# assign a downstream basin to be unique for a mine
df_reg <- df_reg %>%
  filter(!is.na(eco_id), year < 2024) |>
  mutate(downstream = ifelse(distance == 0, 1, downstream)) %>%
  group_by(HYBAS_ID, year) %>%
  arrange(distance) %>%
  slice_head(n = 1) %>%
  ungroup()

saveRDS(df_reg, file = pd("degradation/df_reg.RDS"))
