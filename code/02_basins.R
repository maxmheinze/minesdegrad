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
  geosphere
)


# Load and Prepare Mine Data ----------------------------------------------

mines <- st_read(pd("jde/mines/global_mining_polygons_v2.gpkg"))

africa_codes <- codelist |>
  filter(continent == "Americas") |>
  pluck("iso3c")

m <- mines |> filter(ISO3_CODE %in% africa_codes) |> st_centroid()



# Load and Prepare HydroBASINS Data ---------------------------------------

# https://data.hydrosheds.org/file/hydrobasins/standard/hybas_af_lev01-12_v1c.zip
lvl <- "12"
s <- st_read(pd("hydrobasins_sa_lakes/hybas_lake_sa_lev", lvl, "_v1c.shp"))
s <- st_make_valid(s)
s <- s |> mutate(NEXT_DOWN = ifelse(DIST_SINK == 0, 0, NEXT_DOWN))
lake_id <- s |> st_drop_geometry() |> filter(LAKE != 0) |> pull(HYBAS_ID)
s <- dplyr::filter(s, LAKE == 0)
d <- s |> select(HYBAS_ID, NEXT_DOWN) |> st_drop_geometry()

int <- st_intersects(s, m)

mine_ids <- int[lengths(int) > 0]
names(mine_ids) <- s[["HYBAS_ID"]][lengths(int) > 0]
basin_mines <- data.frame(
  "mine_basin" = rep(names(mine_ids), lengths(mine_ids)),
  "mine_id" = mine_ids |> purrr::map_dfr(as_tibble) |> pull(value))
basin_mines <- basin_mines |>
  left_join(m |> st_drop_geometry() |> transmute(
    mine_id = seq_len(nrow(m)), iso3c = ISO3_CODE, mine_area_km2 = AREA),
    by = "mine_id")
# assign country by larger area of mine
basin_mines_iso <- basin_mines |>
  group_by(mine_basin, iso3c) |>
  summarise(mine_area_km2 = sum(mine_area_km2)) |>
  group_by(mine_basin) |>
  slice_max(mine_area_km2, n = 1) |>
  select(-mine_area_km2)
basin_mines <- basin_mines |>
  group_by(mine_basin) |>
  summarise(mine_area_km2 = sum(mine_area_km2), 
            mine_number = n()) |>
  mutate(mine_avg_area_km2 = mine_area_km2 / mine_number) |> 
  left_join(basin_mines_iso) |>
  left_join(s |> st_drop_geometry() |> 
              transmute(mine_basin = as.character(HYBAS_ID), mine_basin_pfaf_id = PFAF_ID)) |> 
  mutate(mine_basin = as.double(mine_basin), 
         iso3c = replace(iso3c, iso3c == "ESH", "MAR")) |>
  relocate(iso3c:mine_basin_pfaf_id, .after = mine_basin)

# Get the IDs of the basins that contain mines
treated_id <- s[["HYBAS_ID"]][lengths(int) > 0]


# Get Downstream/Upstream Basins ------------------------------------------

# Function that determines up- and downstream basins from NEXT_DOWN field

stream <- \(id, n = 1L, max = 11L, down = TRUE) {
  if(n >= max) return()
  if(down) {
    id_next <- d[["NEXT_DOWN"]][d[["HYBAS_ID"]] == id]
  } else { # Upstream
    id_next <- d[["HYBAS_ID"]][d[["NEXT_DOWN"]] %in% id]
  }
  if(all(id_next == 0) || length(id_next) == 0) return()
  return(c(id_next, Recall(id_next, n = n + 1L, max = max, down = down)))
}

stream_ordered <- function(id, n = 1L, max = 11L, down = TRUE) {
  if (n >= max) return(NULL)
  if (down) {
    id_next <- d[["NEXT_DOWN"]][d[["HYBAS_ID"]] == id]
  } else {
    id_next <- d[["HYBAS_ID"]][d[["NEXT_DOWN"]] %in% id]
  }
  if (all(id_next == 0) || length(id_next) == 0) return()
  results <- cbind(id_next, n)
  id_next <- unique(id_next)
  recursive_results <- Recall(id_next, n = n + 1L, max = max, down = down)
  if (!is.null(recursive_results)) {
    results <- rbind(results, recursive_results)
  }
  return(results)
}

# Get basins downstream and upstream of mines
downstream_ids <- lapply(treated_id, stream, down = TRUE)
upstream_ids <- lapply(treated_id, stream, down = FALSE)

# Tracking orders of basins
downstream_ids_ordered <- lapply(treated_id, stream_ordered, down = TRUE)
downstream_ids_ordered <- lapply(downstream_ids_ordered,
                                 function(x) x[which(!x[,1] %in% lake_id), ])
upstream_ids_ordered <- lapply(treated_id, stream_ordered, down = FALSE)



# Relevant Basins Shapefile -----------------------------------------------

# Extracting the treated and untreated basins to prepare in GEE
downstream_ids_vector <- unlist(downstream_ids)
downstream_ids_vector <- downstream_ids_vector[which(!downstream_ids_vector %in% lake_id)]
upstream_ids_vector <- unlist(upstream_ids)

su <- s |> mutate(
  status = ifelse(HYBAS_ID %in% treated_id, "mine",
                  ifelse(HYBAS_ID %in% downstream_ids_vector, "downstream",
                         ifelse(HYBAS_ID %in% upstream_ids_vector, "upstream",
                                NA_character_)))
)
su |> pull(status) |> table()

# Plot areas
# t <- st_crop(s, st_bbox(s |> filter(!is.na(status)))) |>
#   tm_shape() +
#   tm_fill("status", colorNA = "transparent") +
#   tm_borders() +
#   tm_shape(m) +
#   tm_dots(col = "black", border.col = "white", size = .3, shape = 23)
# tmap_save(t, paste0("basins_mines-l", lvl, ".png"))


# extracting treated and untreated polygons
relevant_basins <- su |>
  filter(!is.na(status))
relevant_basins$basin_area_km2 <- units::drop_units(st_area(relevant_basins) / 10^6)
basin_area <- data.frame(HYBAS_ID = relevant_basins$HYBAS_ID,
                         basin_area_km2 = units::drop_units(st_area(relevant_basins) / 10^6))

write_sf(relevant_basins, pd("degradation/relevant_basins.gpkg"))

# Necessary for Earth Engine
write_sf(relevant_basins, pd("degradation/relevant_basins.shp"))


# Calculating Elaborated Distances ----------------------------------------

basin_centroids <- relevant_basins %>%
  st_centroid %>%
  dplyr::select(HYBAS_ID, NEXT_DOWN, DIST_SINK) %>%
  mutate(NEXT_DOWN = ifelse(DIST_SINK == 0, 0, NEXT_DOWN))

# Centroid Distances
basin_distances_centroids <- basin_centroids %>%
  st_drop_geometry() %>%
  left_join(basin_centroids, by = join_by("NEXT_DOWN" == "HYBAS_ID")) %>%
  dplyr::select(-NEXT_DOWN.y) %>%
  mutate(geometry.x = basin_centroids$geometry) %>%
  rename(geometry.y = geometry) %>%
  relocate(geometry.y, .after = geometry.x) %>%
  mutate(distance = ifelse(NEXT_DOWN == 0, 0,
                           as.numeric(st_distance(st_sfc(geometry.x), # distance in km
                                                  st_sfc(geometry.y), 
                                                  by_element = TRUE)) / 10^3)) %>%
  dplyr::select(HYBAS_ID, NEXT_DOWN, distance)

# River Flow Distances
basin_distances_river <- relevant_basins %>%
  st_drop_geometry() %>%
  dplyr::select(HYBAS_ID, NEXT_DOWN, NEXT_SINK, DIST_SINK, DIST_MAIN) %>%
  mutate(NEXT_DOWN = ifelse(DIST_SINK == 0, 0, NEXT_DOWN)) %>%
  left_join(., ., by = join_by("NEXT_DOWN" == "HYBAS_ID")) %>%
  mutate(distance = ifelse(NEXT_DOWN ==0, DIST_SINK.x, DIST_SINK.x - DIST_SINK.y)) %>%
  dplyr::select(HYBAS_ID, NEXT_DOWN, distance)


# Elaborated Downstream Distances -----------------------------------------

for (i in seq_along(downstream_ids_ordered)) {
  downstream_ids_ordered[[i]] <- rbind(c(treated_id[i], 0), downstream_ids_ordered[[i]])
}

get_distance <- function(id1, id2, distances_df) {
  if (is.na(id1) | is.na(id2) | id1 == id2) return(0)
  row <- distances_df[distances_df$HYBAS_ID == id2 & distances_df$NEXT_DOWN == id1, ]
  if (nrow(row) == 0) return(NA)
  return(row$distance)
}

for (i in seq_along(downstream_ids_ordered)) {
  current_list <- downstream_ids_ordered[[i]]
  distances <- c(0)
  distances_centroid <- c(0)
  if (nrow(current_list) > 1) {  
    for (j in 2:nrow(current_list)) {
      id_prev <- current_list[j-1, 1]
      id_curr <- current_list[j, 1]
      distances <- c(distances, get_distance(id_curr, id_prev, basin_distances_river))
      distances_centroid <- c(distances_centroid, get_distance(id_curr, id_prev, basin_distances_centroids))
    }
  }
  
  distances <- cumsum(distances)
  distances_centroid <- cumsum(distances_centroid)
  downstream_ids_ordered[[i]] <- cbind(current_list, distances, distances_centroid)
}

for(i in seq_along(downstream_ids_ordered)) {
  downstream_ids_ordered[[i]] <- data.frame(downstream_ids_ordered[[i]])
  downstream_ids_ordered[[i]] <- cbind(downstream_ids_ordered[[i]], treated_id[i])
  colnames(downstream_ids_ordered[[i]]) <- c("basin_id", "dist_n", "dist_km", "dist_km_centroid", "mine_basin")
}

downstream_distances_list <- downstream_ids_ordered

downstream_distances_df <- bind_rows(downstream_distances_list)


# Elaborated Upstream Distances -------------------------------------------

for (i in seq_along(upstream_ids_ordered)) {
  upstream_ids_ordered[[i]] <- rbind(c(treated_id[i], 0), upstream_ids_ordered[[i]])
  colnames(upstream_ids_ordered[[i]]) <- c("id_next", "n")
  upstream_ids_ordered[[i]] <- dplyr::left_join(data.frame(upstream_ids_ordered[[i]]), data.frame(basin_distances_river[,1:2]), by = join_by("id_next" == "HYBAS_ID"))
  upstream_ids_ordered[[i]][1,3] <- NA
}

get_distance_upstream <- function(id1, id2, distances_df) {
  if (is.na(id1) | is.na(id2)) return(0)
  row <- distances_df[distances_df$HYBAS_ID == id1 & distances_df$NEXT_DOWN == id2, ]
  if (nrow(row) == 0) return(NA)
  return(row$distance)
}

for (i in seq_along(upstream_ids_ordered)) {
  current_list <- upstream_ids_ordered[[i]]
  
  distances <- numeric(nrow(current_list))
  distances[1] <- 0
  
  distances_centroid <- numeric(nrow(current_list))
  distances_centroid[1] <- 0
  
  if (nrow(current_list) > 1) {
    for (j in 2:nrow(current_list)) {
      id1 <- current_list[j, 1] 
      id2 <- current_list[j, 3]  
      distances[j] <- get_distance_upstream(id1, id2, basin_distances_river)
      distances_centroid[j] <- get_distance_upstream(id1, id2, basin_distances_centroids)
    }
  }
  
  accumulated_distances <- numeric(nrow(current_list))
  accumulated_distances_centroid <- numeric(nrow(current_list))
  for (j in 1:nrow(current_list)) {
    current_distance <- distances[j]
    current_distance_centroid <- distances_centroid[j]
    id2 <- current_list[j, 3]
    while (!is.na(id2) && id2 != current_list[1, 1]) {  # Traverse upstream until the first row
      idx <- which(current_list[, 1] == id2)
      if (length(idx) == 0) break
      current_distance <- current_distance + distances[idx]
      current_distance_centroid <- current_distance_centroid + distances_centroid[idx]
      id2 <- current_list[idx, 3]
    }
    accumulated_distances[j] <- current_distance
    accumulated_distances_centroid[j] <- current_distance_centroid
  }
  
  upstream_ids_ordered[[i]] <- cbind(current_list, distances, distances_centroid, accumulated_distances, accumulated_distances_centroid)
}


for(i in seq_along(upstream_ids_ordered)) {
  upstream_ids_ordered[[i]] <- data.frame(upstream_ids_ordered[[i]][,c(1,2,6,7)])
  upstream_ids_ordered[[i]] <- cbind(upstream_ids_ordered[[i]], treated_id[i])
  colnames(upstream_ids_ordered[[i]]) <- c("basin_id", "dist_n", "dist_km", "dist_km_centroid", "mine_basin")
}

upstream_distances_list <- upstream_ids_ordered

upstream_distances_df <- bind_rows(upstream_distances_list)


# Calculating Distances (Legacy)  -----------------------------------------

# relevant_basins <- read_sf("~/minesfood/data/relevant_basins.gpkg")

# downstream_ids_with_name <- downstream_ids
# names(downstream_ids_with_name) <- treated_id
# downstream_ids_with_name <- lapply(downstream_ids_with_name, function(x)
#   x[which(!x %in% lake_id)])
# 
# upstream_ids_with_name <- upstream_ids
# names(upstream_ids_with_name) <- treated_id
# 
# relevant_basins_centroid <- st_centroid(relevant_basins)
# 
# basins_df <- relevant_basins_centroid
# 
# # Function to calculate the distance between two points
# calculate_distance <- function(id1, id2, basins_df) {
#   basin1 <- basins_df[basins_df$HYBAS_ID == id1,]
#   basin2 <- basins_df[basins_df$HYBAS_ID == id2,]
# 
#   if (nrow(basin1) == 0 || nrow(basin2) == 0) {
#     return(NA)  # If either basin is not found, return NA
#   }
# 
#   # Extract coordinates
#   coords1 <- st_coordinates(st_centroid(basin1))
#   coords2 <- st_coordinates(st_centroid(basin2))
# 
#   # Calculate the distance
#   dist <- distHaversine(coords1, coords2)
#   return(dist)
# }


# Downstream Distances (Legacy) -------------------------------------------

# # Initialize an empty list to store the results
# downstream_distances_list <- list()
# 
# # Loop through each basin and calculate distances to its downstream basins
# for (basin_id in names(downstream_ids_with_name)) {
#   downstream_ids <- downstream_ids_with_name[[basin_id]]
# 
#   if (!is.null(downstream_ids)) {
#     for (downstream_id in downstream_ids) {
#       distance <- calculate_distance(as.numeric(basin_id), as.numeric(downstream_id), basins_df)
#       downstream_distances_list[[length(downstream_distances_list) + 1]] <- data.frame(
#         basin_id = as.double(basin_id),
#         downstream_id = downstream_id,
#         distance = distance
#       )
#     }
#   }
# }
# 
# downstream_distances_df <- bind_rows(downstream_distances_list)


# Upstream Distances (Legacy) ---------------------------------------------

# # Initialize an empty list to store the results
# upstream_distances_list <- list()
# 
# # Loop through each basin and calculate distances to its upstream basins
# for (basin_id in names(upstream_ids_with_name)) {
#   upstream_ids <- upstream_ids_with_name[[basin_id]]
# 
#   if (!is.null(upstream_ids)) {
#     for (upstream_id in upstream_ids) {
#       distance <- calculate_distance(as.numeric(basin_id), as.numeric(upstream_id), basins_df)
#       upstream_distances_list[[length(upstream_distances_list) + 1]] <- data.frame(
#         basin_id = as.double(basin_id),
#         upstream_id = upstream_id,
#         distance = distance
#       )
#     }
#   }
# }
# 
# # Combine the list into a data frame
# upstream_distances_df <- bind_rows(upstream_distances_list)


# Merge Upstream Distances/Downstream Distances DFs -----------------------

downstream_distances_df <- downstream_distances_df |>
  mutate(downstream = 1) |>
  rename(HYBAS_ID = basin_id)

upstream_distances_df <- upstream_distances_df |>
  mutate(downstream = 0) |>
  rename(HYBAS_ID = basin_id)

downstream_upstream_distance <- rbind(downstream_distances_df, upstream_distances_df)

downstream_upstream_distance <- downstream_upstream_distance %>%
  dplyr::filter(!(HYBAS_ID == mine_basin & downstream == 0)) %>%
  rename(order = dist_n, distance = dist_km, distance_centroid = dist_km_centroid)

downstream_upstream_distance <- left_join(downstream_upstream_distance, basin_area) |>
  left_join(basin_mines) |>
  left_join(s |> st_drop_geometry() |> transmute(HYBAS_ID = HYBAS_ID, basin_pfaf_id = PFAF_ID)) |> 
  mutate(mine_area_km2 = replace(mine_area_km2, mine_basin != HYBAS_ID, 0)) |> 
  transmute(HYBAS_ID, HYBAS_PFAF_ID = basin_pfaf_id, mine_basin, mine_basin_pfaf_id, 
            iso3c, downstream, order, distance, distance_centroid, basin_area_km2, mine_area_km2,
            mine_number, mine_avg_area_km2)

write.csv(downstream_upstream_distance,
          pd("degradation/downstream_upstream_distance.csv"), row.names = FALSE)


# Lookup DF with Order ----------------------------------------------------

downstream_df_ordered <- downstream_distances_df %>%
  mutate(status = "downstream") %>%
  dplyr::select(-dist_km, -dist_km_centroid, -downstream) %>%
  `colnames<-`(c("basin", "order", "mine_basin", "status"))

upstream_df_ordered <- upstream_distances_df %>%
  mutate(status = "upstream") %>%
  dplyr::select(-dist_km, -dist_km_centroid, -downstream) %>%
  `colnames<-`(c("basin", "order", "mine_basin", "status"))

mine_df_ordered <- data.frame(treated_id, treated_id) |>
  mutate(status = "mine",
         order = 0) |>
  `colnames<-`(c("basin", "mine_basin", "status", "order")) %>%
  dplyr::relocate(basin, order, mine_basin, status)

basins_ordered <- rbind(downstream_df_ordered, upstream_df_ordered, mine_df_ordered)


# Relevant Basins Shapefile with Order ------------------------------------

# If a basin is down- or upstream of multiple mines, it is considered to
# be only down- or upstream of the closest mine.

basins_ordered_unique <- basins_ordered |>
  group_by(basin) |>
  arrange(order, .by_group = T) |>
  slice_head(n = 1)

so <- s |>
  left_join(basins_ordered_unique, by = join_by("HYBAS_ID" == "basin")) |>
  dplyr::filter(!is.na(status)) |>
  dplyr::select(HYBAS_ID, mine_basin, status, order) |>
  left_join(basin_area) |>
  left_join(basin_mines) |>
  left_join(s |> st_drop_geometry() |> transmute(HYBAS_ID, basin_pfaf_id = PFAF_ID)) |> 
  relocate(iso3c, .after = mine_basin) |>
  relocate(geometry, .after = mine_area_km2)

write_sf(so, pd("degradation/relevant_basins_ordered.gpkg"))
write_sf(so, pd("degradation/relevant_basins_ordered.shp"))


# Add Order to Downstream/Upstream Distance DF ----------------------------

# This is legacy code that is redundant (I think) but I'm not sure whether
# we use this some time later so I just ensured consistency and left it
# in --Max

downstream_upstream_distance_ordered <- downstream_upstream_distance |>
  transmute(HYBAS_ID = as.double(HYBAS_ID), mine_basin = as.double(mine_basin), 
            downstream, distance, distance_centroid) |>
  full_join(basins_ordered_unique, by = join_by("HYBAS_ID" == "basin", "mine_basin")) |>
  drop_na(status) |>
  mutate_at(vars(downstream, distance, distance_centroid), ~replace_na(., 0)) |>
  left_join(basin_area) |>
  left_join(basin_mines) |>
  left_join(s |> st_drop_geometry() |> transmute(HYBAS_ID, basin_pfaf_id = PFAF_ID)) |> 
  mutate(mine_area_km2 = replace(mine_area_km2, mine_basin != HYBAS_ID, 0)) |> 
  transmute(HYBAS_ID, HYBAS_PFAF_ID = basin_pfaf_id, mine_basin, mine_basin_pfaf_id, 
            iso3c, downstream, status, order, distance, distance_centroid, basin_area_km2,
            mine_area_km2, mine_number, mine_avg_area_km2)

write.csv(downstream_upstream_distance_ordered,
          pd("degradation/downstream_upstream_distance_ordered.csv"), row.names = FALSE)
