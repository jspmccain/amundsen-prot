
## This script makes the sampling map
library(ggOceanMaps)
library(stringr)

## load the lat and long
station_lat_longs <- read.csv("data/station_lat_longs_manual_edit.csv")

## format the lat and long
all_lats <- str_split(station_lat_longs$lat,pattern = "  ")
all_lats_unlist <- -1*(unlist(lapply(all_lats, `[[`, 1)) %>% as.numeric())
all_longs <- str_split(station_lat_longs$long, pattern = "  ")
all_longs_unlist <- unlist(lapply(all_longs, `[[`, 1)) %>% as.numeric()

## make them into a df
stations_formatted <- data.frame(station = station_lat_longs$man_sel.station,
                                 long = all_longs_unlist,
                                 lat = all_lats_unlist)

## remove the stations we didn't do proteomics on
stations_formatted_all_msz <- stations_formatted %>% dplyr::filter(station != 4,
                                                                   station != 1,
                                                                   station != 10,
                                                                   station != 11)

## format the dataframe to work with ggOceanMaps
dt <- stations_formatted_all_msz %>% 
  dplyr::rename(lon = long) %>% 
  # dplyr::select(-station) %>% 
  dplyr::mutate(long = lon - 2*180) %>% 
  dplyr::select(-lon)

## make the map
map_of_sampling <- basemap(data = dt, bathymetry = TRUE, rotate = TRUE, bathy.style = "rcb", 
        limits = c(-128, -112, -75,  -70)) + 
  ggspatial::geom_spatial_point(data = dt %>% inner_join(mean_taxa_by_station, by = 'station'), 
                                aes(x = long, y = lat, colour = mean_diatom_by_station), 
                                size = 7, 
                                alpha = 1) +
  scale_colour_gradient(high = "darkgreen", low = "white", limits = c(0.1, 50)) +
  ggspatial::geom_spatial_text(data = dt, size = 3,
                                aes(x = long, y = lat,
                                    label = as.character(station))) + #;map_of_sampling #+
  theme(legend.position = 'side');map_of_sampling


map_of_sampling_w_depth <- basemap(data = dt, bathymetry = TRUE, rotate = TRUE, bathy.style = "rcb", 
        limits = c(-128, -112, -75,  -70)) + 
  ggspatial::geom_spatial_point(data = dt, aes(x = long, y = lat), 
                                color = "black", size = 7, 
                                pch = 21, 
                                fill= 'bisque', 
                                alpha = 1) +
  ggspatial::geom_spatial_text(data = dt, size = 3,
                               aes(x = long, y = lat,
                                   label = as.character(station))) 

ggsave(map_of_sampling_w_depth, filename = 'figures/map_of_sampling_w_depth.pdf', height = 10, width = 10)
