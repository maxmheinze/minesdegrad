

# Setup -------------------------------------------------------------------

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("00_prelims.R")

pacman::p_load(tidyverse,
               sf)


# Read Data ---------------------------------------------------------------

veg <- read_csv(pd("degradation/veg.csv"))



# Clean veg ---------------------------------------------------------------

# Isolate year and land use values
veg <- veg %>%
  mutate(year = as.numeric(substr(bandName, 16, 20))) %>%
  dplyr::select(-bandName) %>%
  filter(histogram != "{}") %>% 
  mutate(histogram = str_replace_all(histogram, "[{}]", "")) %>%
  separate_rows(histogram, sep = ",") %>% 
  separate(histogram, into = c("landuse", "share"), sep = "=") %>%  
  mutate(landuse = as.numeric(landuse), 
         share = as.numeric(share))

# Calculate deforestation class shares
veg <- veg %>%
  mutate(defclass = floor(landuse/100)) %>%
  group_by(HYBAS_ID, year, defclass) %>%
  summarize(NEXT_DOWN = first(NEXT_DOWN),
            NEXT_SINK = first(NEXT_SINK),
            MAIN_BAS=first(MAIN_BAS),
            DIST_SINK = first(DIST_SINK),
            DIST_MAIN = first(DIST_MAIN),
            SUB_AREA = first(SUB_AREA),
            UP_AREA = first(UP_AREA),
            PFAF_ID = first(PFAF_ID),
            SIDE = first(SIDE),
            LAKE = first(LAKE),
            ENDO = first(ENDO),
            COAST = first(COAST),
            ORDER = first(ORDER),
            SORT = first(SORT),
            share = sum(share)) %>%
  ungroup() %>%
  group_by(HYBAS_ID, year) %>%
  mutate(rel_area = share / sum(share)) %>%
  ungroup()

write_csv(veg, pd("degradation/veg_cleaned.csv"))

