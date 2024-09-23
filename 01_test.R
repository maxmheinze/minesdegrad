

# Setup -------------------------------------------------------------------

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("00_prelims.R")

pacman::p_load(tidyverse,
               sf)


# Read Data ---------------------------------------------------------------

svcov <- read_csv(pd("degradation/secondary_veg_coverage_basins.csv"))

svage <- read_csv(pd("/degradation/secondary_veg_age_basins.csv"))


# Test --------------------------------------------------------------------


