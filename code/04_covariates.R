# Load Packages -----------------------------------------------------------

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("00_prelims.R")

pacman::p_load(
  sf,
  dplyr,
  tmap,
  tidyverse,
  countrycode,
  geosphere,
  terra,
  ncdf4,
  magrittr,
  parallel
)

mbsvg <- read_csv(pd("degradation/mine_basins_secveg.csv"))  
relevant_basins <- st_read(pd("degradation/relevant_basins.gpkg"))


# Population --------------------------------------------------------------

rstlist <- paste0(pd("jde/worldpop/"), list.files(path = pd("jde/worldpop/"))) 

worldpop <- mclapply(
  rstlist,
  function(x) {
    current_rstr <- terra::rast(x)
    
    extracted <- terra::extract(terra::subset(current_rstr, 1), relevant_basins, fun = sum, exact = T, na.rm = T)
    
    return(extracted)
  },
  mc.cores = 12
)

worldpop1 <- worldpop

for (i in seq_along(worldpop)) {
  worldpop[[i]] <- dplyr::mutate(worldpop[[i]], year = 2010 + i,
                                 HYBAS_ID = relevant_basins$HYBAS_ID)
  colnames(worldpop[[i]]) <- c("id", "pop", "year", "HYBAS_ID")
}

worldpop %<>%
  bind_rows() %>%
  select(HYBAS_ID, year, pop)

worldpop_2015 <- worldpop |> filter(year == 2015) |> 
  transmute(HYBAS_ID, pop_2015 = pop)

worldpop_pred <- worldpop |> rename(pop_pred = pop)
for(i in unique(worldpop_pred$HYBAS_ID)) {
  if(all(is.na(worldpop_pred |> filter(HYBAS_ID == i) |> pull(pop_pred)))) {
    worldpop_pred <- worldpop_pred |> rows_append(
      tibble(HYBAS_ID = i, year = 2021:2023, pop_pred = NA)
    )
  } else {
    lm_temp <- lm(pop_pred ~ year, data = worldpop_pred |> filter(HYBAS_ID == i))
    worldpop_pred <- worldpop_pred |> rows_append(
      tibble(HYBAS_ID = i, year = 2021:2023, 
             pop_pred = predict(lm_temp, 
                                newdata = tibble(year = 2021:2023)))
    )
  }
}

worldpop <- full_join(worldpop_pred, worldpop) |> 
  left_join(worldpop_2015, by = "HYBAS_ID") |> 
  group_by(HYBAS_ID) |> 
  mutate(pop_const = pop, 
         pop_const = replace(pop_const, year %in% 2021:2023, 
                             pop_const[year == 2020]))


saveRDS(worldpop, file = pd("degradation/population.RDS"))



# Elevation ---------------------------------------------------------------

elevation <- terra::rast(pd("redd/grid_pnas/elevation.tif"))

elevation_df <- terra::extract(x = elevation, y = relevant_basins,
                               weights = TRUE, exact = TRUE) |>
  filter(!is.na(elevation)) |>
  group_by(ID) |> transmute(weighted = elevation * weight / sum(weight)) |>
  summarise(elevation = sum(weighted, na.rm = T)) |>
  left_join(tibble(ID = seq_len(nrow(relevant_basins)), HYBAS_ID = relevant_basins$HYBAS_ID), by = "ID") |>
  relocate(HYBAS_ID, .before = ID) |> dplyr::select(-ID)

saveRDS(elevation_df, file = pd("degradation/elevation.RDS"))

rm(elevation); gc()



# Slope -------------------------------------------------------------------

slope <- terra::rast(pd("redd/grid_pnas/slope.tif"))

slope_df <- terra::extract(x = slope, y = relevant_basins,
                           weights = TRUE, exact = TRUE) |>
  filter(!is.na(slope)) |>
  group_by(ID) |> transmute(weighted = slope * weight / sum(weight)) |>
  summarise(slope = sum(weighted, na.rm = T)) |>
  left_join(tibble(ID = seq_len(nrow(relevant_basins)), HYBAS_ID = relevant_basins$HYBAS_ID), by = "ID") |>
  relocate(HYBAS_ID, .before = ID) |> dplyr::select(-ID)

saveRDS(slope_df, file = pd("degradation/slope.RDS"))

rm(slope); gc()



# SoilGrid ----------------------------------------------------------------

# taking the value with the highest weight

soilgrid <- terra::rast(pd("redd/grid_pnas/soilgrid.tif"))

soil_prim_df <- terra::extract(x = soilgrid, y = relevant_basins,
                               weights = TRUE, exact = TRUE) |>
  filter(!is.na(soilgrid)) |>
  rename(soil_prim = soilgrid) |>
  group_by(ID) |> mutate(weight = weight / sum(weight)) |>
  left_join(tibble(ID = seq_len(nrow(relevant_basins)), HYBAS_ID = relevant_basins$HYBAS_ID), by = "ID") |>
  relocate(HYBAS_ID, .before = ID) |> dplyr::select(-ID) |>
  group_by(HYBAS_ID, soil_prim) |>
  summarise(weight = sum(weight, na.rm = T)) |>
  filter(!is.na(soil_prim)) |>
  group_by(HYBAS_ID) |>
  slice_max(weight, n = 1) |>
  dplyr::select(-weight)

saveRDS(soil_prim_df, file = pd("degradation/soil_primary.RDS"))

groups <- read.csv(pm("soilgrid_grouped.csv"))

soil_prim_df$soilgrid_grouped <- droplevels(as.factor(
  groups[match(soil_prim_df$soil_prim, groups$Number), "WRB_group"]))

saveRDS(soil_prim_df, file = pd("degradation/soil_primary-grouped.RDS"))

rm(soilgrid); gc()


# Biome -------------------------------------------------------------------

# assigning eco-region with highest weight per basin

ecoregions <- terra::rast(pd("redd/grid_pnas/ecoregions_2017.tif"))
ecoregions_concordance <- read_csv(pd("redd/grid_pnas/ecoregions_2017_concordance_tbl.csv"))
names(ecoregions_concordance) <- tolower(names(ecoregions_concordance))

ecoregions_df <- terra::extract(x = ecoregions, y = relevant_basins,
                                weights = TRUE, exact = TRUE) |>
  filter(!is.na(ecoregions_2017)) |>
  rename(eco_id = ecoregions_2017) |>
  group_by(ID) |> mutate(weight = weight / sum(weight)) |>
  left_join(tibble(ID = seq_len(nrow(relevant_basins)), HYBAS_ID = relevant_basins$HYBAS_ID), by = "ID") |>
  dplyr::select(-ID) |> relocate(HYBAS_ID, .before = eco_id) |>
  group_by(HYBAS_ID, eco_id) |>
  summarise(ecoregion_weight = sum(weight, na.rm = T)) |>
  filter(!is.na(eco_id)) |>
  group_by(HYBAS_ID) |>
  slice_max(ecoregion_weight, n = 1) |>
  dplyr::select(-ecoregion_weight) |>
  left_join(ecoregions_concordance, by = "eco_id")

ecoregions_df <- ecoregions_df |> 
  mutate(biome_group = "Forests", 
         biome_group = replace(biome_group, biome_num %in% c(7, 9, 10), "Grasslands"), 
         biome_group = replace(biome_group, biome_num %in% 13, "Deserts"))

saveRDS(ecoregions_df, file = pd("degradation/ecoregion.RDS"))

rm(ecoregions, ecoregions_concordance); gc()



# Accessibility to Cities 2015 --------------------------------------------

accessibility_to_cities_2015 <- terra::rast(pd("redd/grid_pnas/accessibility_to_cities_2015.tif"))

accessibility_to_cities_2015_df <- terra::extract(x = accessibility_to_cities_2015,
                                                  y = relevant_basins,
                                                  weights = TRUE, exact = TRUE) |>
  filter(accessibility_to_cities_2015 > 0, !is.na(accessibility_to_cities_2015)) |> 
  group_by(ID) |> transmute(weighted = accessibility_to_cities_2015 * weight / sum(weight)) |>
  summarise(accessibility_to_cities_2015 = sum(weighted, na.rm = T)) |>
  left_join(tibble(ID = seq_len(nrow(relevant_basins)), HYBAS_ID = relevant_basins$HYBAS_ID), by = "ID") |>
  relocate(HYBAS_ID, .before = ID) |> dplyr::select(-ID)

saveRDS(accessibility_to_cities_2015_df, file = pd("degradation/accessibility_to_cities_2015.RDS"))

rm(accessibility_to_cities_2015); gc()


# Merging time-invariant geo variables in one dataframe -------------------

elevation_df <- readRDS(pd("degradation/elevation.RDS"))
slope_df <- readRDS(pd("degradation/slope.RDS"))
soil_prim_df <- readRDS(pd("degradation/soil_primary-grouped.RDS"))
ecoregions_df <- readRDS(pd("degradation/ecoregion.RDS"))
accessibility_to_cities_2015_df <- readRDS(pd("degradation/accessibility_to_cities_2015.RDS"))

geo_data_df <- left_join(ecoregions_df, elevation_df, by = "HYBAS_ID") |>
  left_join(slope_df, by = "HYBAS_ID") |>
  left_join(soil_prim_df, by = "HYBAS_ID") |>
  left_join(accessibility_to_cities_2015_df, by = "HYBAS_ID")

saveRDS(geo_data_df, file = pd("degradation/geo_data_agg.RDS"))



# Time-Varying Climate Variables ------------------------------------------

# Weather variables from Climatic Research Unit (CRU)
# https://www.nature.com/articles/s41597-020-0453-3

dates <- paste0(rep(2001:2023, each = 12), "-", stringr::str_pad(1:12, 2, pad = "0"))


# Mean Temperature --------------------------------------------------------

tmp_data <- terra::rast(pd("AP_FP/CRU_unzipped/cru_ts4.08.1901.2023.tmp.dat.nc"), 
                        subds = "tmp") 
pos <- which(substr(time(tmp_data), 1, 7) %in% dates)

start <- Sys.time()
tmp_df <- terra::extract(tmp_data[[pos]], relevant_basins, weights = T, exact = T) |>
  filter(if_all(contains("tmp"), ~ !is.na(.x))) |> 
  group_by(ID) |>
  mutate(weight = weight / sum(weight),
         across(contains("tmp"), function(x) x * weight)) |>
  dplyr::select(-weight) |>
  summarise(across(contains("tmp"), function(x) sum(x, na.rm = T))) |>
  left_join(tibble(ID = seq_len(nrow(relevant_basins)), HYBAS_ID = relevant_basins$HYBAS_ID), by = "ID") |>
  relocate(HYBAS_ID, .before = ID) |> dplyr::select(-ID) |>
  tidyr::pivot_longer(contains("tmp")) |>
  left_join(tibble(name = names(tmp_data),
                   time = substr(time(tmp_data), 1, 7)),
            by = "name") |>
  transmute(HYBAS_ID,
            year = as.numeric(substr(time, 1, 4)),
            month = as.numeric(substr(time, 6, 7)),
            tmp_mean = value) |>
  arrange(year, month, HYBAS_ID)
cat("Calculation finished after", format(Sys.time() - start), "\n")

saveRDS(tmp_df, file = pd("degradation/cru_tmp.RDS"))



# Max Temperature ---------------------------------------------------------

tmx_data <- terra::rast(pd("AP_FP/CRU_unzipped/cru_ts4.08.1901.2023.tmx.dat.nc"),
                        subds = "tmx")
pos <- which(substr(time(tmx_data), 1, 7) %in% dates)

start <- Sys.time()
tmx_df <- terra::extract(tmx_data[[pos]], relevant_basins, weights = T, exact = T) |>
  filter(if_all(contains("tmx"), ~ !is.na(.x))) |> 
  group_by(ID) |>
  mutate(weight = weight / sum(weight),
         across(contains("tmx"), function(x) x * weight)) |>
  dplyr::select(-weight) |>
  summarise(across(contains("tmx"), function(x) sum(x, na.rm = T))) |>
  left_join(tibble(ID = seq_len(nrow(relevant_basins)), HYBAS_ID = relevant_basins$HYBAS_ID), by = "ID") |>
  relocate(HYBAS_ID, .before = ID) |> dplyr::select(-ID) |>
  tidyr::pivot_longer(contains("tmx")) |>
  left_join(tibble(name = names(tmx_data),
                   time = substr(time(tmx_data), 1, 7)),
            by = "name") |>
  transmute(HYBAS_ID,
            year = as.numeric(substr(time, 1, 4)),
            month = as.numeric(substr(time, 6, 7)),
            tmp_max = value) |>
  arrange(year, month, HYBAS_ID)
cat("Calculation finished after", format(Sys.time() - start), "\n")

saveRDS(tmx_df, file = pd("degradation/cru_tmx.RDS"))



# Min Temperature ---------------------------------------------------------

tmn_data <- terra::rast(pd("AP_FP/CRU_unzipped/cru_ts4.08.1901.2023.tmn.dat.nc"),
                        subds = "tmn")
pos <- which(substr(time(tmn_data), 1, 7) %in% dates)

start <- Sys.time()
tmn_df <- terra::extract(tmn_data[[pos]], relevant_basins, weights = T, exact = T) |>
  filter(if_all(contains("tmn"), ~ !is.na(.x))) |> 
  group_by(ID) |>
  mutate(weight = weight / sum(weight),
         across(contains("tmn"), function(x) x * weight)) |>
  dplyr::select(-weight) |>
  summarise(across(contains("tmn"), function(x) sum(x, na.rm = T))) |>
  left_join(tibble(ID = seq_len(nrow(relevant_basins)), HYBAS_ID = relevant_basins$HYBAS_ID), by = "ID") |>
  relocate(HYBAS_ID, .before = ID) |> dplyr::select(-ID) |>
  tidyr::pivot_longer(contains("tmn")) |>
  left_join(tibble(name = names(tmn_data),
                   time = substr(time(tmn_data), 1, 7)),
            by = "name") |>
  transmute(HYBAS_ID,
            year = as.numeric(substr(time, 1, 4)),
            month = as.numeric(substr(time, 6, 7)),
            tmp_min = value) |>
  arrange(year, month, HYBAS_ID)
cat("Calculation finished after", format(Sys.time() - start), "\n")

saveRDS(tmn_df, file = pd("degradation/cru_tmn.RDS"))


# Precipitation -----------------------------------------------------------

pre_data <- terra::rast(pd("AP_FP/CRU_unzipped/cru_ts4.08.1901.2023.pre.dat.nc"),
                        subds = "pre")
pos <- which(substr(time(pre_data), 1, 7) %in% dates)

start <- Sys.time()
pre_df <- terra::extract(pre_data[[pos]], relevant_basins, weights = T, exact = T) |>
  filter(if_all(contains("pre"), ~ !is.na(.x))) |> 
  group_by(ID) |>
  mutate(weight = weight / sum(weight),
         across(contains("pre"), function(x) x * weight)) |>
  dplyr::select(-weight) |>
  summarise(across(contains("pre"), function(x) sum(x, na.rm = T))) |>
  left_join(tibble(ID = seq_len(nrow(relevant_basins)), HYBAS_ID = relevant_basins$HYBAS_ID), by = "ID") |>
  relocate(HYBAS_ID, .before = ID) |> dplyr::select(-ID) |>
  tidyr::pivot_longer(contains("pre")) |>
  left_join(tibble(name = names(pre_data),
                   time = substr(time(pre_data), 1, 7)),
            by = "name") |>
  transmute(HYBAS_ID,
            year = as.numeric(substr(time, 1, 4)),
            month = as.numeric(substr(time, 6, 7)),
            precipitation = value) |>
  arrange(year, month, HYBAS_ID)
cat("Calculation finished after", format(Sys.time() - start), "\n")

saveRDS(pre_df, file = pd("degradation/cru_pre.RDS"))



# Cloud Cover -------------------------------------------------------------

cld_data <- terra::rast(pd("AP_FP/CRU_unzipped/cru_ts4.08.1901.2023.cld.dat.nc"),
                        subds = "cld")
pos <- which(substr(time(cld_data), 1, 7) %in% dates)

start <- Sys.time()
cld_df <- terra::extract(cld_data[[pos]], relevant_basins, weights = T, exact = T) |>
  filter(if_all(contains("cld"), ~ !is.na(.x))) |> 
  group_by(ID) |>
  mutate(weight = weight / sum(weight),
         across(contains("cld"), function(x) x * weight)) |>
  dplyr::select(-weight) |>
  summarise(across(contains("cld"), function(x) sum(x, na.rm = T))) |>
  left_join(tibble(ID = seq_len(nrow(relevant_basins)), HYBAS_ID = relevant_basins$HYBAS_ID), by = "ID") |>
  relocate(HYBAS_ID, .before = ID) |> dplyr::select(-ID) |>
  tidyr::pivot_longer(contains("cld")) |>
  left_join(tibble(name = names(cld_data),
                   time = substr(time(cld_data), 1, 7)),
            by = "name") |>
  transmute(HYBAS_ID,
            year = as.numeric(substr(time, 1, 4)),
            month = as.numeric(substr(time, 6, 7)),
            cloud_cover = value) |>
  arrange(year, month, HYBAS_ID)
cat("Calculation finished after", format(Sys.time() - start), "\n")

saveRDS(cld_df, file = pd("degradation/cru_cld.RDS"))



# Rainy Days --------------------------------------------------------------

wet_data <- terra::rast(pd("AP_FP/CRU_unzipped/cru_ts4.08.1901.2023.wet.dat.nc"),
                        subds = "wet")
pos <- which(substr(time(wet_data), 1, 7) %in% dates)

start <- Sys.time()
wet_df <- terra::extract(wet_data[[pos]], relevant_basins, weights = T, exact = T) |>
  filter(if_all(contains("wet"), ~ !is.na(.x))) |> 
  group_by(ID) |>
  mutate(weight = weight / sum(weight),
         across(contains("wet"), function(x) x * weight)) |>
  dplyr::select(-weight) |>
  summarise(across(contains("wet"), function(x) sum(x, na.rm = T))) |>
  left_join(tibble(ID = seq_len(nrow(relevant_basins)), HYBAS_ID = relevant_basins$HYBAS_ID), by = "ID") |>
  relocate(HYBAS_ID, .before = ID) |> dplyr::select(-ID) |>
  tidyr::pivot_longer(contains("wet")) |>
  left_join(tibble(name = names(wet_data),
                   time = substr(time(wet_data), 1, 7)),
            by = "name") |>
  transmute(HYBAS_ID,
            year = as.numeric(substr(time, 1, 4)),
            month = as.numeric(substr(time, 6, 7)),
            wet_days = value) |>
  arrange(year, month, HYBAS_ID)
cat("Calculation finished after", format(Sys.time() - start), "\n")

saveRDS(wet_df, file = pd("degradation/cru_wet.RDS"))



# Frost Days --------------------------------------------------------------

frs_data <- terra::rast(pd("AP_FP/CRU_unzipped/cru_ts4.08.1901.2023.frs.dat.nc"),
                        subds = "frs")
pos <- which(substr(time(frs_data), 1, 7) %in% dates)

start <- Sys.time()
frs_df <- terra::extract(frs_data[[pos]], relevant_basins, weights = T, exact = T) |>
  filter(if_all(contains("frs"), ~ !is.na(.x))) |> 
  group_by(ID) |>
  mutate(weight = weight / sum(weight),
         across(contains("frs"), function(x) x * weight)) |>
  dplyr::select(-weight) |>
  summarise(across(contains("frs"), function(x) sum(x, na.rm = T))) |>
  left_join(tibble(ID = seq_len(nrow(relevant_basins)), HYBAS_ID = relevant_basins$HYBAS_ID), by = "ID") |>
  relocate(HYBAS_ID, .before = ID) |> dplyr::select(-ID) |>
  tidyr::pivot_longer(contains("frs")) |>
  left_join(tibble(name = names(frs_data),
                   time = substr(time(frs_data), 1, 7)),
            by = "name") |>
  transmute(HYBAS_ID,
            year = as.numeric(substr(time, 1, 4)),
            month = as.numeric(substr(time, 6, 7)),
            frs_days = value) |>
  arrange(year, month, HYBAS_ID)
cat("Calculation finished after", format(Sys.time() - start), "\n")

saveRDS(frs_df, file = pd("degradation/cru_frs.RDS"))



# Aggregating monthly frequency of time-varying climate variables  --------

tmp_df <- readRDS(pd("degradation/cru_tmp.RDS"))
tmx_df <- readRDS(pd("degradation/cru_tmx.RDS"))
tmn_df <- readRDS(pd("degradation/cru_tmn.RDS"))
pre_df <- readRDS(pd("degradation/cru_pre.RDS"))
cld_df <- readRDS(pd("degradation/cru_cld.RDS"))
wet_df <- readRDS(pd("degradation/cru_wet.RDS"))
frs_df <- readRDS(pd("degradation/cru_frs.RDS"))

met_data_df <- left_join(tmp_df, tmx_df, by = c("HYBAS_ID", "year", "month")) |>
  left_join(tmn_df, by = c("HYBAS_ID", "year", "month")) |>
  left_join(pre_df, by = c("HYBAS_ID", "year", "month")) |>
  left_join(cld_df, by = c("HYBAS_ID", "year", "month")) |>
  left_join(wet_df, by = c("HYBAS_ID", "year", "month")) |>
  left_join(frs_df, by = c("HYBAS_ID", "year", "month")) |>
  group_by(HYBAS_ID, year) |>
  summarise(tmp_mean = mean(tmp_mean, na.rm = T),
            tmp_max = max(tmp_max, na.rm = T),
            tmp_min = min(tmp_min, na.rm = T),
            precipitation = sum(precipitation, na.rm = T),
            cloud_cover = mean(cloud_cover, na.rm = T),
            wet_days = sum(wet_days, na.rm = T),
            frs_days = sum(frs_days, na.rm = T))

saveRDS(met_data_df, file = pd("degradation/met_data_agg.RDS"))
