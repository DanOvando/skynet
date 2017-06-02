# run_skynet -------
# Author: Dan Ovando
# Project: skynet
# Summary: Master script to run project skynet


# load libraries ------------------------------
# Summary: Load things you need
rm(list = ls())
set.seed(123)
library(tidyverse)
library(forcats)
library(modelr)
library(stringr)
library(rstan)
library(rstanarm)
library(car)
library(AER)
library(broom)
library(bigrquery)
library(FishData)
library(tmap)
library(rgdal)
library(leaflet)
library(tmaptools)
library(raster)
library(viridis)

demons::load_functions()

# set options ------------------------------
# Summary: set options for model run

run_name <- 0.5

run_dir <- paste('results', run_name, sep = '/')

run_description <-
  'Just a bit of scratch paper'

if (dir.exists(run_dir) == F) {
  dir.create(run_dir, recursive = T)
}

write(run_description,
      file = paste(run_dir, 'RUN_DESCRIPTION.txt', sep = '/'))

# Run Options ------------------------------
# Summary: Set model options

fishing_threshold <-  0.75 # Set a prob. cutoff for what is considered "fishing"


# Load in data ------------------------------
# Summary: Load data

project <- "ucsb-gfw"

fishing_connection <-
  src_bigquery(project, "skynet") # This function initiliazes a connection with a BQ dataset


ebs <- fishing_connection %>%
  tbl("ebs_w_vessels") %>%
  collect(n = Inf)


  mutate(y_m_d = date(timestamp)) %>%
  mutate(
    obs_year = year(y_m_d),
    obs_month = month(y_m_d),
    round_lat = round(lat, as.integer(2)),
    round_lon = round(lon, as.integer(2))
  ) %>%
  filter(obs_month %in% summer) %>%
  group_by(round_lat, round_lon, obs_year,mmsi) %>%
  summarise(
    fishing_hours = sum(hours * measure_new_score * as.numeric(measure_new_score >= fishing_threshold)),
    mean_dist_from_shore = mean(distance_from_shore),
    mean_dist_from_port = mean(distance_from_port)
  ) %>%
  collect() %>%
  rename(year = obs_year) %>%
  mutate(fishing_hours = ifelse(fishing_hours == 0, NA, fishing_hours)) %>%
  mutate(log_fishing_hours = log(fishing_hours))

ebs_fleet <- vessel_connection %>%
  tbl('complete_fishing_fleet_characteristics_2015') %>%
  collect() %>%
  filter(mmsi %in% unique(ebs$mmsi))

ebs <- ebs %>%
  group_by(round_lat, round_lon, year) %>%
  summarise(
    fishing_hours = sum(fishing_hours),
    mean_dist_from_shore = mean(mean_dist_from_shore),
    mean_dist_from_port = mean(mean_dist_from_port)
  ) %>%
  mutate(fishing_hours = ifelse(fishing_hours == 0, NA, fishing_hours)) %>%
  mutate(log_fishing_hours = log(fishing_hours))


ebs_trawl <-
  FishData::download_catch_rates(survey = "EBSBTS", species_set = 50) %>%
  as_data_frame() %>%
  set_names(tolower(colnames(.)))

ebs_trawl_summary <- ebs_trawl %>%
  set_names(tolower(colnames(.))) %>%
  mutate(round_lat = round(lat, 2),
         round_lon = round(long, 2)) %>%
  group_by(year, round_lat, round_lon,sci) %>%
  summarise(sci_mean_cpue = mean(wt, na.rm = T)) %>%
  group_by(year,round_lat, round_lon) %>%
  summarise(total_wt = sum(sci_mean_cpue, na.rm = T),
            mean_cpue = mean(sci_mean_cpue, na.rm = T)) %>%
  mutate(log_mean_cpue = log(mean_cpue),
         log_total_wt = log(total_wt))

ram_stocks <- readxl::read_excel("~/Box Sync/Databases/RAM/RAM v3.80/DB Files With Assessment Data/RLSADB v3.80 (assessment data only).xlsx",
                   sheet = 'stock') %>%
  mutate(scientificname = tolower(scientificname))

ram_stocks$scientificname[ram_stocks$scientificname == 'theragra chalcogramma'] <- 'gadus chalcogrammus'

alaska_ram_stocks <- ram_stocks %>%
  filter(str_detect(region,'Alaska')) %>%
  dplyr::select(scientificname) %>%
  unique() %>%
  unlist() %>%
  as.character() %>%
  tolower()

ram_ebs_trawl <- ebs_trawl %>%
  set_names(tolower(colnames(.))) %>%
  mutate(round_lat = round(lat, 2),
         round_lon = round(long, 2),
         scientificname = str_replace(sci,'_',' ') %>% tolower()) %>%
  filter(scientificname %in% alaska_ram_stocks)

ram_ebs_trawl_summary <- ram_ebs_trawl %>%
  set_names(tolower(colnames(.))) %>%
  mutate(round_lat = round(lat, 2),
         round_lon = round(long, 2)) %>%
  group_by(year, round_lat, round_lon,sci) %>%
  summarise(sci_mean_cpue = mean(wt, na.rm = T)) %>%
  group_by(year,round_lat, round_lon) %>%
  summarise(total_wt = sum(sci_mean_cpue, na.rm = T),
            mean_cpue = mean(sci_mean_cpue, na.rm = T)) %>%
  mutate(log_mean_cpue = log(mean_cpue),
         log_total_wt = log(total_wt))
# scratch ------------------------------
# Summary: just getting started comparing trawls and gfw

joint_ebs <- ebs %>%
  right_join(ram_ebs_trawl_summary, by = c('year', 'round_lat', 'round_lon')) %>%
  filter(is.na(total_wt) == F)

gfw_vs_trawl_plot <- joint_ebs %>%
  ggplot(aes(
    log(fishing_hours),
    log(mean_cpue),
    fill = log(mean_dist_from_port)
  )) +
  geom_point(shape = 21, size = 5, alpha = 0.75) +
  geom_smooth(method = 'lm', se = F, color = 'red',
              linetype = 2) +
  xlab('GFW Fishing Hours') +
  ylab('Survey Biomass') +
  scale_fill_viridis(name = 'Distance from Port',
                     direction = -1) +
  theme_ipsum(base_size = 18,
              axis_title_size = 24) +
  theme(legend.position = 'bottom',
        legend.title = element_text(size = 16),
        legend.key.width = unit(3, 'cm'))

ggsave(filename = paste0(run_dir,'/correlation.png'),gfw_vs_trawl_plot,
       width = 9, height = 5)

# Test some regressions ------------------------------
# Summary:

joint_ebs <- joint_ebs %>%
  ungroup() %>%
  mutate(slog_fishing_hours = (log_fishing_hours - mean(log_fishing_hours, na.rm = T)) / (2*sd(log_fishing_hours, na.rm = T)),
         smean_dist_from_port = (mean_dist_from_port - mean(mean_dist_from_port, na.rm = T)) / (2 * sd(mean_dist_from_port, na.rm = T)))

model <- lm(log_mean_cpue ~ slog_fishing_hours + slog_fishing_hours^2 + smean_dist_from_port, data = joint_ebs)

summary(model)

joint_ebs <- joint_ebs %>%
  add_predictions(model)

joint_ebs %>%
  ggplot(aes(log_mean_cpue, pred)) +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_smooth(method = 'lm')

# Make some Maps ------------------------------
# Summary: Make some Maps

coords <- SpatialPoints(ebs %>% dplyr::select(round_lon, round_lat) %>% ungroup())

gfw_ebs <- SpatialPointsDataFrame(coords, ebs %>% ungroup())

proj4string(gfw_ebs) <- CRS("+proj=longlat +datum=WGS84")

alaska <- read_shape('~/Box Sync/Databases/alaska_coastline/alaska_coastline.shp')

alaska <-  crop_shape(alaska, gfw_ebs)

gfw_ebs <- spTransform(gfw_ebs,crs(alaska))

gfw_ebs_map <- tm_shape(alaska, title = 'test')  +
  tm_lines(col = 'black') +
  tm_shape(gfw_ebs, title = 'test') +
  tm_bubbles(size = 'fishing_hours', col = 'log_fishing_hours',
             palette = plasma(20),
             legend.size.show = F,
             alpha = 0.5) +
  tm_legend(legend.outside.position = c("right"),
            legend.outside = TRUE
            ) #format legend

save_tmap(tm = gfw_ebs_map, filename = paste(run_dir,'/gfw_effort_map.png', sep = ''),
          width = 8, height = 5, units = 'in')

gfw_port_ebs_map <- tm_shape(alaska)  +
  tm_lines(col = 'black') +
  tm_shape(gfw_ebs) +
  tm_bubbles(size = 'mean_dist_from_port', col = 'mean_dist_from_port',
             palette = viridis(10),
             legend.size.show = F) +
  tm_legend(legend.outside.position = c("right"),
            legend.outside = TRUE) #format legend



coords <- SpatialPoints(ram_ebs_trawl_summary %>% ungroup() %>%  dplyr::select(round_lon, round_lat))

fviz_ebs <- SpatialPointsDataFrame(coords, ram_ebs_trawl_summary %>% ungroup())

proj4string(fviz_ebs) <- CRS("+proj=longlat +datum=WGS84")



fviz_ebs_map <- tm_shape(alaska, fill = 'grey')  +
  tm_lines(col = 'black') +
  tm_shape(fviz_ebs) +
  tm_bubbles(size = .05, col = 'total_wt',
             palette = plasma(20),
             legend.size.show = F,
             alpha = 0.5) +
  tm_legend(legend.outside.position = c("right"),
            legend.outside = TRUE) #format legend

ebs_trawl_comp <- ebs_trawl %>%
  group_by(year, long, lat) %>%
  summarise(total_wt = sum(wt, na.rm = T)) %>%
  mutate(log_total_wt = log(total_wt)) %>%
  filter(year == 2016)

coords <- SpatialPoints(ebs_trawl_comp %>% ungroup() %>%  dplyr::select(long, lat))

fviz_ebs <- SpatialPointsDataFrame(coords, ebs_trawl_comp %>% ungroup())

proj4string(fviz_ebs) <- CRS("+proj=longlat +datum=WGS84")

fviz_ebs_map <- tm_shape(alaska, fill = 'grey')  +
  tm_lines(col = 'black') +
  tm_shape(fviz_ebs) +
  tm_bubbles(size = .05, col = 'log_total_wt',
             palette = plasma(20),
             legend.size.show = F,
             alpha = 0.5) +
  tm_legend(legend.outside.position = c("right"),
            legend.outside = TRUE) #format legend

save_tmap(tm = fviz_ebs_map, filename = paste(run_dir,'/fviz_trawl_map.png', sep = ''),
          width = 8, height = 5, units = 'in')


# Save Results ------------------------------
# Summary:

plots <- ls()[str_detect(ls(), '_plot')]

plot_list <- purrr::map(plots, get) %>%
  set_names(plots)

maps <- ls()[str_detect(ls(), '_map')]

map_list <- purrr::map(maps, get) %>%
  set_names(maps)


save(file = paste0(run_dir,'/ebs.Rdata'), ebs)

save(file = paste0(run_dir,'/ebs_trawl.Rdata'), ebs_trawl)

save(file = paste0(run_dir,'/ebs_ram_trawl.Rdata'), ram_ebs_trawl)

save(file = paste0(run_dir,'/plots.Rdata'), plot_list)

save(file = paste0(run_dir,'/maps.Rdata'), map_list)

demons::save_plots(plot_dir = run_dir)

demons::save_tmaps(plot_dir = run_dir)

tmap_mode('plot')
map_list$gfw_ebs_map
