# run-skynet -------
# Author: Dan Ovando
# Project: skynet
# Summary: A wrapper script to fit and analyze models predicting abundance
# using GFW data

# prep run ----------------------------------------------------------

set.seed(42)
library(bigrquery)
library(lubridate)
library(leaflet)
library(rgdal)
library(FishData)
library(ggmap)
library(stringr)
library(viridis)
library(purrr)
library(purrrlyr)
library(modelr)
library(VAST)
library(TMB)
library(trelliscopejs)
library(ggjoy)
library(modelr)
library(caret)
library(taxize)
library(tidyverse)

demons::load_functions('functions')

run_name <- 'testing'

run_description <- 'development of full paper'

run_dir <- file.path('results', run_name, '')

if (dir.exists(run_dir) == F) {
  dir.create(run_dir)
}

write(run_description, file = paste0(run_dir, 'description.txt'))


# set run options ---------------------------------------------------------

query_fishdata <-  F

query_gfw <- F

query_environmentals <-  F

query_mpas <-  F

vasterize <-  T

fished_only <-  T

hack_zeros <-  T

gfw_project <- "ucsb-gfw"

gfw_dataset <- 'skynet'

lat_lon_res <-
  0.25 # round data to intervels of 0.25 degrees lat lon, as in 56.25

num_knots <-  100

species_list <- c('Gadus_chalcogrammus') #, 'Hippoglossoides_elassodon','Gadus_macrocephalus')

survey_list <- 'ebsbts'

# get fishdata ------------------------------------------------------------

if (query_fishdata == T) {
  fishdata_names <- tribble(
    ~ survey,
    ~ survey_name,
    #--/--
    'EBSBTS',
    'eastern bering sea bottom trawl survey',
    # 'WCGBTS','west coast groundfish bottom trawl survey',
    # 'WCGHL', 'West_coast_groundfish_hook_and_line',
    'GOABTS',
    'gulf of alaska bottom trawl survey',
    'AIBTS',
    'aluetian islands bottom trawl survey'
  )
  fish_data <- fishdata_names %>%
    mutate(survey_data = map(
      survey,
      ~ download_catch_rates(.x, species_set = 50,
                             error_tol = 1e-2)
    ))

  fish_data <- unnest(fish_data) %>%
    mutate(survey = tolower(survey))

  save(file = 'data/fish_data.Rdata', fish_data)

} else {
  load("data/fish_data.Rdata")

}


survey_bbox <- fish_data %>%
  filter(Long < 0) %>%
  group_by(survey) %>%
  summarise(
    min_year = min(Year, na.rm = T),
    max_year = max(Year, na.rm = T),
    min_lat = min(Lat, na.rm = T),
    max_lat = max(Lat, na.rm = T),
    min_lon = min(Long, na.rm = T),
    max_lon = max(Long, na.rm = T)
  ) %>%
  mutate(survey = tolower(survey))

survey_bbox

survey_names <- paste0(fish_data$survey %>% tolower(), '_gfw') %>% unique()

# query_gfw ---------------------------------------------------------------


if (query_gfw == T) {
  gfw_data <- data_frame(survey = survey_names)

  fishing_connection <- DBI::dbConnect(bigrquery::dbi_driver(),
                                       project = gfw_project,
                                       dataset = gfw_dataset)

  gfw_data <- gfw_data %>%
    mutate(
      gfw_data = map(
        survey,
        ~ fishing_connection %>%
          tbl(.x) %>%
          collect(n = Inf) %>%
          select(-b_mmsi,-b_year) %>%
          set_names(str_replace_all(colnames(.), '(a_)|(b_)', '')) %>%
          filter(is.na(on_fishing_list_nn) == F) %>%
          arrange(year, mmsi)

      )
    )

  save(file = 'data/gfw_data.Rdata', gfw_data)

} else {
  load(file = 'data/gfw_data.Rdata')
}

gfw_bbox <- gfw_data %>%
  unnest() %>%
  group_by(survey) %>%
  summarise(
    min_year = min(year, na.rm = T),
    max_year = max(year, na.rm = T),
    min_lat = min(rounded_lat, na.rm = T),
    max_lat = max(rounded_lat, na.rm = T),
    min_lon = min(rounded_lon, na.rm = T),
    max_lon = max(rounded_lon, na.rm = T)
  )


# query environmental data ------------------------------------------------


if (query_environmentals == T) {
  chl_data <- gfw_bbox %>%
    mutate(
      chl_a = pmap(
        list(
          min_year = pmax(2012, min_year),
          max_year = max_year,
          min_lat = min_lat,
          max_lat = max_lat,
          min_lon = min_lon,
          max_lon = max_lon
        ),
        query_erddap,
        desired_data = 'chl-a',
        date_interval = 1,
        space_interval = 1,
        runit = lat_lon_res
      )
    )

  chl_plot_foo <- function(dat) {
    qmplot(
      x = rlon,
      y = rlat,
      data = dat,
      color = log(mean_chlorophyll),
      maptype = 'toner-lite',
      source = 'google'
    ) +
      scale_color_viridis(guide = F) +
      facet_wrap( ~ year) +
      theme_classic() +
      labs(caption = 'SST(c)')


  }

  chl_data <-  chl_data %>%
    mutate(plots = map_plot(chl_a, chl_plot_foo))

  save(file = 'data/chl_env_data.Rdata', chl_data)


  sst_data <- gfw_bbox %>%
    mutate(
      sst = pmap(
        list(
          min_year = pmax(2012, min_year),
          max_year = max_year,
          min_lat = min_lat,
          max_lat = max_lat,
          min_lon = min_lon,
          max_lon = max_lon
        ),
        query_erddap,
        desired_data = 'sst',
        date_interval = 14,
        space_interval = 25,
        runit = lat_lon_res
      )
    )

  sst_plot_foo <- function(dat) {
    qmplot(
      x = rlon,
      y = rlat,
      data = dat,
      color = log(mean_analysed_sst),
      maptype = 'toner-lite',
      source = 'google'
    ) +
      scale_color_viridis(guide = F) +
      facet_wrap( ~ year) +
      theme_classic() +
      labs(caption = 'SST(c)')


  }

  sst_data <-  sst_data %>%
    mutate(plots = map_plot(sst, sst_plot_foo))

  save(file = 'data/sst_env_data.Rdata', sst_data)

  wave_data <- gfw_bbox %>%
    mutate(
      wave_height = pmap(
        list(
          min_year = pmax(2012, min_year),
          max_year = max_year,
          min_lat = min_lat,
          max_lat = max_lat,
          min_lon = min_lon,
          max_lon = max_lon
        ),
        query_erddap,
        desired_data = 'waves',
        date_interval = 24 * 14,
        #1 hour interval
        space_interval = 0.5,
        runit = lat_lon_res
      )
    )


  wave_data <- wave_data %>%
    mutate(plots = map_plot(
      wave_height,
      ~ quick_map(
        .x,
        lat_var = quo(rlat),
        lon_var = quo(rlon),
        plot_var = quo(mean_Thgt),
        facet_var = quo(year)
      )
    ))

  save(file = 'data/wave_env_data.Rdata', wave_data)


  wind_data <- gfw_bbox %>%
    mutate(
      wind_speed = pmap(
        list(
          min_year = pmax(2012, min_year),
          max_year = max_year,
          min_lat = min_lat,
          max_lat = max_lat,
          min_lon = min_lon,
          max_lon = max_lon
        ),
        query_erddap,
        desired_data = 'wind',
        date_interval = 1,
        #1 hour interval
        space_interval = 1,
        runit = lat_lon_res
      )
    )


  wind_data <- wind_data %>%
    mutate(wind_vectors = map(
      wind_speed,
      ~ .x %>%
        spread(wind_direction, mean_wind_speed) %>%
        mutate(
          wind_vector = map2(x_wind, y_wind, ~ vectorize_wind(x_wind = .x, y_wind = .y)),
          wind_speed = map_dbl(wind_vector, "wind_speed"),
          wind_angle = map_dbl(wind_vector, "wind_angle")
        ) %>%
        select(-wind_vector, -x_wind, -y_wind)
    )) %>%
    select(-wind_speed)



  wind_data <- wind_data %>%
    mutate(plots = map_plot(
      wind_vectors,
      ~ quick_map(
        .x,
        lat_var = quo(rlat),
        lon_var = quo(rlon),
        plot_var = quo(wind_speed),
        facet_var = quo(year)
      )
    ))

  save(file = 'data/wind_env_data.Rdata', wind_data)

  topo_data <- gfw_bbox %>%
    mutate(
      altitude = pmap(
        list(
          min_year = pmax(2012, min_year),
          max_year = max_year,
          min_lat = min_lat,
          max_lat = max_lat,
          min_lon = min_lon,
          max_lon = max_lon
        ),
        query_erddap,
        desired_data = 'topography',
        date_interval = 1,
        #1 hour interval
        space_interval = 1,
        runit = lat_lon_res
      )
    )

  topo_data  <- topo_data %>%
    mutate(altitude = map(
      altitude,
      ~ .x %>%
        mutate(m_below_sea_level = -1 * mean_altitude)
    ))



  topo_data <- topo_data %>%
    mutate(plots = map_plot(
      altitude,
      ~ quick_map(
        .x,
        lat_var = quo(rlat),
        lon_var = quo(rlon),
        plot_var = quo(m_below_sea_level),
        facet_var = quo(mean_altitude_units)
      )
    ))

  save(file = 'data/topo_env_data.Rdata', topo_data)

} else {
  data_files <- list.files("data/")

  data_files <-
    paste0('data/', data_files[str_detect(data_files, 'env_data.Rdata')])

  for (i in seq_along(data_files)) {
    load(data_files[i])

  }

}


# query mpas --------------------------------------------------------------

if (query_mpas == T) {
  mpa_ids <-
    data_frame(mpa_type = c('no_take_mpas_id', 'restricted_use_mpas_id'))

  fishing_connection <- DBI::dbConnect(bigrquery::dbi_driver(),
                                       project = project,
                                       dataset = 'regions_ids')

  mpa_ids <- mpa_ids %>%
    mutate(data = map(mpa_type, ~ fishing_connection %>%
                        tbl(.x) %>%
                        collect(n = Inf))) %>%
    unnest()

  mpa_ids <- mpa_ids %>%
    mutate(
      mpa_type = str_split(region_id, ':', simplify = T)[, 1],
      mpa_id = str_split(region_id, ':', simplify = T)[, 2]
    )

  save(file = 'data/mpa_ids.Rdata', mpa_ids)
} else {
  load(file = 'data/mpa_ids.Rdata')

}



# perform filters and data processing ---------------------------------------------------------

# decide species that go in
# for now, query FAO for species fished... this is hacky though since means that depends on static data...

fao <- data.table::fread("~/Box Sync/Databases/FAO/FAO_Capture_1950to2012.csv", header = T) %>%
  as_data_frame() %>%
  filter(isscaap_group < 60) # get fao data

fao_species <- unique(fao$sciname) # unique fao species

species_seen <- fish_data$Sci %>% unique() %>%
  str_replace('_',' ') # species in fishdata

safe_syns <- safely(taxize::synonyms)

species_data <- data_frame(sci_name = species_seen) %>%
  filter(is.na(sci_name) == F) %>%
  mutate(worms = map(sci_name, taxize::get_wormsid, ask = F)) %>%
  mutate(worms_id = map(worms, `[[`,1) %>% as.character()) %>% # extract the worms id
  filter(is.na(worms_id) == F) %>%
  mutate(known_synonyms = map(worms_id %>% as.numeric(), safe_syns, db = 'worms'))

species_data <- species_data %>%
  mutate(syn_results = map(known_synonyms,'result')) %>%
  mutate(syns = map(syn_results,`[[`,1)) %>%
  mutate(sci_syns = map(syns, "scientificname")) # extract synonyms


check_fao <- function(sci, syns, faos){

  in_fao <- sci %in% faos | any(syns %in% faos)

}

fished_species <- species_data %>%
  mutate(in_fao = map2_lgl(sci_name, sci_syns, check_fao, fao_species)) %>%
  select(sci_name, sci_syns, in_fao) %>%
  filter(in_fao == T) # find species whose primary name or synonyms are in fao database

# if (fished_only == T){ # include only fished species
#
#   fish_data <- fish_data %>%
#     filter(str_replace(Sci, '_',' ') %in% fished_species$sci_name)
#
# }
# sci2comm(fished_species$sci_name)


# MAJOR HACK TO DEAL WITH PERFECTLY OBSERVED SPECIES FIX FIX FIX FIX ONCE VAST IS FIXED
# can't run species that are perfectly observed at this point, so for now going to do a
# MAJOR hack and for species without any zeros, convert values below some quantile threshold to zeros

if (hack_zeros == T) {

hack_foo <- function(x){

  min_wt <- min(x$Wt)

  if (min_wt > 0){

    min_quantile <- quantile(x$Wt,0.025) %>% as.numeric()

    x$Wt[x$Wt <=min_quantile ] <-  0

  }

  return(x)
}

fish_data <- fish_data %>%
  nest(-survey,-Sci,- Year) %>%
  mutate(hackdat = map(data,hack_foo )) %>%
  select(-data) %>%
  unnest()
} # close hack zeros


fish_data <- fish_data %>%
  filter(is.na(Sci) == F)

# vasterize ---------------------------------------------------------------

species_list <- unique(fish_data$Sci)

# species_list <- c('Gadus_macrocephalus')

subset_fish_data <- fish_data %>%
  filter(survey %in% survey_list , Sci %in% species_list) %>%  dplyr::rename(Lon = Long,
                                                                              Catch_KG = Wt,
                                                                              spp  = Sci) %>%
  mutate(AreaSwept_km2 = 1 , Vessel = 'missing') %>%
  select(Year, Lat, Lon, Vessel, AreaSwept_km2, Catch_KG, spp) %>%
  set_names(tolower(colnames(.))) %>%
  filter(year >= 2012) %>%
  group_by(year, spp) %>%
  mutate(seen_types = length(unique(catch_kg > 0))) %>%
  filter(seen_types == 2,
         spp != 'Paguridae') %>%
  select(-seen_types) %>%
  ungroup() %>%
  mutate(vessel = as.factor(vessel),
         spp = as.factor(spp)) %>%
  as.data.frame() #note that this breaks if it's a tibble, should make a bit more flexible


# arg <- subset_fish_data %>%
#   group_by(year,spp) %>%
#   summarise(nobs = mean(catch_kg > 0)) %>%
#   ggplot(aes(spp, nobs)) +
#   geom_col() +
#   coord_flip() +
#   facet_wrap(~year)

set.seed(42)

if (vasterize == T) {
  vast_fish <-
    vasterize_index(
      raw_data = subset_fish_data,
      run_dir = run_dir,
      nstart = 100,
      region = 'Eastern_Bering_Sea',
      n_x = num_knots,
      obs_model = c(2,0)
    )

  # vast_fish$spatial_list$loc_x_lat_long
  #
  # vast_fish$spatial_densities
#
#   ggmap::qmplot(
#     x = approx_long,
#     y = approx_lat,
#     data = vast_fish$spatial_densities,
#     color = density
#   ) +
#     facet_wrap(~species) +
#     scale_color_viridis()

  save(file = paste0(run_dir, 'vast_fish.Rdata'),
       vast_fish)
} else {
  load(file = paste0(run_dir, 'vast_fish.Rdata'))

}
# build database ----------------------------------------------------------

# pre-process gfw data

ebs_trawlers <- gfw_data %>%
  filter(survey == "ebsbts_gfw") %>%
  unnest() %>%
  filter(inferred_label_allyears == 'trawlers')

skynet_data <- ebs_trawlers %>%
  group_by(year, rounded_lat, rounded_lon) %>%
  summarise(
    num_vessels = length(unique(mmsi)),
    total_hours = sum(total_hours),
    total_engine_hours = sum(total_hours * inferred_engine_power),
    dist_from_port = mean(mean_distance_from_port),
    dist_from_shore = mean(mean_distance_from_shore),
    mean_vessel_length = mean(inferred_length),
    no_take_mpa = any(is.na(mpant) == F),
    restricted_use_mpa =  any(is.na(mparu) == F)
  ) %>%
  ungroup() %>%
  mutate(vessel_hours = total_engine_hours * num_vessels,
         any_fishing = total_engine_hours > 0) %>%
  group_by(rounded_lat, rounded_lon) %>%
  arrange(year) %>%
  mutate(
    total_engine_hours_lag1 = lag(total_engine_hours, 1),
    total_engine_hours_lag2 = lag(total_engine_hours, 2),
    total_engine_hours_lag3 = lag(total_engine_hours, 3),
    port_engine_hours = dist_from_port * total_engine_hours,
    port_hours = dist_from_port * total_hours,
    port_numbers = dist_from_port * num_vessels,
    shore_numbers = dist_from_shore * num_vessels,
    shore_hours = dist_from_shore * total_hours,
    shore_engine_hours = dist_from_shore * total_engine_hours,
    large_vessels = num_vessels * mean_vessel_length
  )


# pre process fishdata

if (fished_only == T){ # include only fished species

  total_fish_data <- vast_fish$spatial_densities  %>%
    filter(str_replace(species, '_',' ') %in% fished_species$sci_name)

} else {
  total_fish_data <- vast_fish$spatial_densities
}


total_fish_data <- total_fish_data %>%
  group_by(knot, year) %>%
  summarise(density = sum(density))


# aggregate data

skynet_data <-
  create_skynet_data(
    gfw_data = skynet_data,
    fish_data = total_fish_data,
    fish_knots = vast_fish$spatial_list$loc_x_lat_long,
    surveys = 'ebsbts_gfw',
    topo_data = topo_data,
    wind_data = wind_data,
    chl_data = chl_data,
    wave_data = wave_data,
    sst_data = sst_data
  )


# skynet_data %>%
#   group_by(rounded_lat, rounded_lon) %>%
#   summarise(percent_missing = mean(is.na(mean_altitude))) %>%
#   ggplot(aes(rounded_lon, rounded_lat, fill = factor(percent_missing))) +
#   geom_tile()

# qmplot(rounded_lon, rounded_lat, color = log(density), data = skynet_data) +
#   scale_color_viridis()

# prepare models ----------------------------------------------------------

skynet_names <- colnames(skynet_data)

do_not_scale <-
  c(
    'year',
    'rounded_lat',
    'rounded_lon',
    'no_take_mpa',
    'restricted_use_mpa',
    'knot',
    'distance',
    'any_fishing',
    'log_density',
    'density'
  ) # variables that should not be centered and scaled

do_scale <- skynet_names[!skynet_names %in% do_not_scale]

skynet_data <- skynet_data %>%
  ungroup() %>%
  purrrlyr::dmap_at(do_scale, ~ (.x - mean(.x, na.rm = T)) / (2 * sd(.x, na.rm = T))) #center and scale variables that should be

dep_var <- c('log_density')

ind_vars <-
  c(
    paste(skynet_names[!skynet_names %in% c(dep_var,
                                            'year',
                                            'rounded_lat',
                                            'rounded_lon',
                                            'knot',
                                            'distance')], collapse = '+'),
    paste(
      c(
        'total_hours' ,
        'large_vessels' ,
        'port_hours' ,
        'port_engine_hours' ,
        'port_numbers' ,
        'shore_hours' ,
        'shore_engine_hours' ,
        'shore_numbers' ,
        'dist_from_port' ,
        'dist_from_shore' ,
        'num_vessels' ,
        'mean_vessel_length'
      ),
      collapse = '+'
    ),
    paste(
      c(
        'large_vessels' ,
        'port_numbers' ,
        'dist_from_port' ,
        'dist_from_shore' ,
        'mean_vessel_length'
      ),
      collapse = '+'
    ),
    paste(
      c(
        'large_vessels' ,
        'port_numbers' ,
        'dist_from_port' ,
        'dist_from_shore' ,
        'mean_vessel_length',
        'm_below_sea_level',
        'mean_analysed_sst',
        'mean_chlorophyll',
        'wind_speed',
        'wind_angle'
      ),
      collapse = '+'
    ),
    paste(
      c(
        'large_vessels' ,
        'port_numbers' ,
        'total_engine_hours',
        'restricted_use_mpa',
        'any_fishing',
        'port_engine_hours',
        'dist_from_port' ,
        'dist_from_shore' ,
        'mean_vessel_length',
        'm_below_sea_level',
        'mean_analysed_sst',
        'mean_chlorophyll',
        'wind_speed',
        'wind_angle',
        'random_var'
      ),
      collapse = '+'
    ),
    paste(
      c(
        'port_numbers' ,
        'dist_from_port' ,
        'dist_from_shore' ,
        'mean_vessel_length',
        'm_below_sea_level',
        'mean_analysed_sst',
        'mean_chlorophyll',
        'wind_speed',
        'wind_angle',
        'random_var'
      ),
      collapse = '+'
    )
  )

models <- 'rf'

test_train_data <- modelr::crossv_kfold(skynet_data, k = 2)

skynet_models <-  purrr::cross_df(list(
  dep_var = dep_var,
  ind_vars = list(ind_vars[[6]]),
  temp = list(test_train_data),
  model = models
)) %>%
  unnest(temp, .drop = F) %>%
  mutate(ind_vars = map(ind_vars, unlist))

# run models --------------------------------------------------------------

sfm <- safely(fit_skynet)

skynet_models <- skynet_models %>%
  slice(1) %>%
  mutate(fitted_model = pmap(
    list(
      dep_var = dep_var,
      ind_vars = ind_vars,
      model = model,
      training = train,
      testing = test
    ),
    sfm,
    fitcontrol_number = 5,
    fitcontrol_repeats = 2
  ))

# diagnose models ---------------------------------------------------------

wee <- map(skynet_models$fitted_model,'result')

a <- wee[[1]]

a$model$finalModel

varImpPlot(a$model$finalModel)

arg <- skynet_models$test[[1]] %>%
  as_data_frame() %>%
  na.omit() %>%
  ungroup() %>%
  mutate(ld_hat = predict(a$model$finalModel, newdata = .))

arg %>%
  ggplot(aes(log_density, ld_hat)) +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1))


huh <- a$model$trainingData %>%
  mutate(ld_hat = a$model$finalModel$predicted)

huh %>%
  ggplot(aes(.outcome, ld_hat)) +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1))
