# run-skynet -------
# Author: Dan Ovando
# Project: skynet
# Summary: A wrapper script to fit and analyze models predicting abundance
# using GFW data

# prep run ----------------------------------------------------------

set.seed(42)
rm(list = ls())
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

run_name <- 'full-test'

run_description <- 'development of full paper'

run_dir <- file.path('results', run_name, '')

if (dir.exists(run_dir) == F) {
  dir.create(run_dir)
}

write(run_description, file = paste0(run_dir, 'description.txt'))


# set section options (what to run) ---------------------------------------------------------

run_models <-  F # fit statistical models to data

vasterize <-  F # run vast or load saved object

plot_data <-  F # plot covariates and maps

query_fishdata <-  F # get trawl survey data or load saved

query_gfw <- F # get gfw data or load saved

query_environmentals <-  F #get environmental data or load saved

query_mpas <-  F # get mpa data or load saved

query_synonyms <-  F # look up synonyms for each species

query_coastal_map <-  F # get coastal map or load


# set run options (how to run) ---------------------------------------------------------

fished_only <-  T # only include fished species in predictive model fitting

hack_zeros <-  T # hack zeros into perfectly observed species

survey_years_only <-  T #only include years that were actually surveyed

include_all_species <-  T #include all species in VAST

bottom_gears_only <- F # only include bottom gear

clip_gfw  <-  T # clip GFW data to extent of trawl surveys

gfw_project <- "ucsb-gfw"

gfw_dataset <- 'skynet'

max_percent_missing <-  0.25

lat_lon_res <-
  0.25 # round data to intervels of 0.25 degrees lat lon, as in 56.25

min_year <-  2012

species_list <-
  c('Gadus_chalcogrammus',
    'Hippoglossoides_elassodon',
    'Gadus_macrocephalus')

survey_list <- c('ebsbts', 'goabts', 'aibts', 'wcgbts', 'wcghl')

knots_per_survey <- tribble(~ survey,
                            ~ knots,
                            #--/--
                            'ebsbts',
                            200,
                            'goabts',
                            200,
                            'aibts',
                            200,
                            'wcgbts',
                            1000,
                            'wcghl',
                            100)

# get fishdata ------------------------------------------------------------

if (query_fishdata == T) {
  fishdata_names <- tribble(
    ~ survey,
    ~ survey_name,
    ~ survey_region,
    #--/--
    'ebsbts',
    'eastern bering sea bottom trawl survey',
    'eastern_bering_sea',
    'wcgbts',
    'west coast groundfish bottom trawl survey',
    'california_current',
    'wcghl',
    'West_coast_groundfish_hook_and_line',
    'wcghl_domain',
    'goabts',
    'gulf of alaska bottom trawl survey',
    'gulf_of_alaska',
    'aibts',
    'aluetian islands bottom trawl survey',
    'aleutian_islands'
  )

  sdcr  <- safely(download_catch_rates)

  fish_data <- fishdata_names %>%
    mutate(survey_data = map(
      survey,
      ~ download_catch_rates(.x, species_set = 50,
                             error_tol = 1e-2)
    ))

  fish_data <- unnest(fish_data) %>%
    mutate(survey = tolower(survey))

  save(file = 'data/fish_data.Rdata', fish_data)
  #
  #   a = fish_data %>%
  #     mutate(rlat = round(Lat,2),
  #            rlong = round(Long,2)) %>%
  #     group_by(rlat,rlong) %>%
  #     summarise(tw = sum(Wt, na.rm = T)) %>%
  #     filter(rlong > -160, rlong < 0)
  #
  #   qmplot(rlong, rlat, color = log(tw), data = a) +
  #     scale_color_viridis() +
  #     theme_minimal()
  #
} else {
  load("data/fish_data.Rdata")

}

fish_data <- fish_data %>%
  left_join(knots_per_survey, by = 'survey')

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

survey_names <- fish_data$survey %>% tolower() %>% unique()



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
          select(-b_mmsi, -b_year) %>%
          set_names(str_replace_all(colnames(.), '(a_)|(b_)', '')) %>%
          filter(is.na(on_fishing_list_nn) == F) %>%
          arrange(year, mmsi)

      )
    )
  save(file = 'data/gfw_data.Rdata', gfw_data)

} else {
  load(file = 'data/gfw_data.Rdata')
}

gfw_bbox <- gfw_data %>% # calculate bounding box for GFW data
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

# process gfw data

merge_known_and_inferred <- function(known, inferred, cap = NA) {
  if (is.na(cap)) {
    cap <-  1.1 * max(inferred, na.rm = T)
  }

  out <- ifelse(is.na(known), pmin(inferred, cap), known)

}

gfw_data <- gfw_data %>% # merge known and inferred vessel characteristics
  unnest() %>%
  mutate(
    engine_power = map2_dbl(
      known_engine_power,
      inferred_engine_power,
      merge_known_and_inferred,
      cap = 1.1 * max(inferred_engine_power, na.rm = T)
    ),
    vessel_length = map2_dbl(known_length,
                             inferred_length,
                             merge_known_and_inferred,
                             cap = 150),
    tonnage = map2_dbl(
      known_tonnage,
      inferred_tonnage,
      merge_known_and_inferred,
      cap = 1.1 * max(inferred_tonnage, na.rm = T)
    )
  ) %>%
  nest(-survey)


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
      color = log(mean_chlorophyll)
    ) +
      scale_color_viridis(guide = F) +
      facet_wrap(~ year) +
      theme_classic() +
      labs(caption = 'SST(c)')


  }

  chl_data <-  chl_data %>%
    mutate(plots = map_plot(chl_a, chl_plot_foo))

  a <- chl_data$chl_a[[4]]

  qmplot(rlon, rlat, color = mean_chlorophyll, data = a)

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
      color = log(mean_analysed_sst)
    ) +
      scale_color_viridis(guide = F) +
      facet_wrap(~ year) +
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


  # wave_data <- wave_data %>%
  #   mutate(plots = map_plot(
  #     wave_height,
  #     ~ quick_map(
  #       .x,
  #       lat_var = quo(rlat),
  #       lon_var = quo(rlon),
  #       plot_var = quo(mean_Thgt),
  #       facet_var = quo(year)
  #     )
  #   ))

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
        select(-wind_vector,-x_wind,-y_wind)
    )) %>%
    select(-wind_speed)



  # wind_data <- wind_data %>%
  #   mutate(plots = map_plot(
  #     wind_vectors,
  #     ~ quick_map(
  #       .x,
  #       lat_var = quo(rlat),
  #       lon_var = quo(rlon),
  #       plot_var = quo(wind_speed),
  #       facet_var = quo(year)
  #     )
  #   ))

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



  # topo_data <- topo_data %>%
  #   mutate(plots = map_plot(
  #     altitude,
  #     ~ quick_map(
  #       .x,
  #       lat_var = quo(rlat),
  #       lon_var = quo(rlon),
  #       plot_var = quo(m_below_sea_level),
  #       facet_var = quo(mean_altitude_units)
  #     )
  #   ))

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
                                       project = gfw_project,
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

if (query_synonyms == T) {
  fao <-
    data.table::fread("~/Box Sync/Databases/FAO/FAO_Capture_1950to2012.csv",
                      header = T) %>%
    as_data_frame() %>%
    filter(isscaap_group < 60) # get fao data

  fao_species <- unique(fao$sciname) # unique fao species

  species_seen <- fish_data$Sci %>% unique() %>%
    str_replace('_', ' ') # species in fishdata

  safe_syns <- safely(taxize::synonyms)

  species_data <- data_frame(sci_name = species_seen) %>%
    filter(is.na(sci_name) == F) %>%
    mutate(worms = map(sci_name, taxize::get_wormsid, ask = F)) %>%
    mutate(worms_id = map(worms, `[[`, 1) %>% as.character()) %>% # extract the worms id
    filter(is.na(worms_id) == F) %>%
    mutate(known_synonyms = map(worms_id %>% as.numeric(), safe_syns, db = 'worms'))

  species_data <- species_data %>%
    mutate(syn_results = map(known_synonyms, 'result')) %>%
    mutate(syns = map(syn_results, `[[`, 1)) %>%
    mutate(sci_syns = map(syns, "scientificname")) # extract synonyms


  check_fao <- function(sci, syns, faos) {
    in_fao <- sci %in% faos | any(syns %in% faos)

  }

  fished_species <- species_data %>%
    mutate(in_fao = map2_lgl(sci_name, sci_syns, check_fao, fao_species)) %>%
    select(sci_name, sci_syns, in_fao) %>%
    filter(in_fao == T) # find species whose primary name or synonyms are in fao database

  save(file = 'data/fished_species.Rdata', fished_species)

} else {
  load(file = 'data/fished_species.Rdata')

}
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
  hack_foo <- function(x) {
    min_wt <- min(x$Wt)

    if (min_wt > 0) {
      min_quantile <- quantile(x$Wt, 0.025) %>% as.numeric()

      x$Wt[x$Wt <= min_quantile] <-  0

    }

    return(x)
  }

  fish_data <- fish_data %>%
    nest(-survey, -Sci, -Year) %>%
    mutate(hackdat = map(data, hack_foo)) %>%
    select(-data) %>%
    unnest()
} # close hack zeros


fish_data <- fish_data %>%
  filter(is.na(Sci) == F)

if (bottom_gears_only == T) {
  gfw_data <- gfw_data %>%
    unnest() %>%
    filter(
      inferred_label_allyears == 'trawlers' |
        inferred_sublabel_allyears == 'pots_and_traps'
    ) %>%
    nest(-survey)
}

gfw_data <- gfw_data %>%
  filter(survey %in%  survey_list)

if (include_all_species == T) {
  species_list <- unique(fish_data$Sci)
}


subset_fish_data <- fish_data %>%
  filter(survey %in% survey_list , Sci %in% species_list) %>%  dplyr::rename(Lon = Long,
                                                                             Catch_KG = Wt,
                                                                             spp  = Sci) %>%
  mutate(AreaSwept_km2 = 1 , Vessel = 'missing') %>%
  select(survey,
         survey_region,
         knots,
         Year,
         Lat,
         Lon,
         Vessel,
         AreaSwept_km2,
         Catch_KG,
         spp) %>%
  set_names(tolower(colnames(.))) %>%
  filter(year >= min_year) %>%
  group_by(survey, year, spp) %>%
  mutate(seen_types = length(unique(catch_kg > 0))) %>%
  group_by(survey, spp) %>%
  mutate(seen_every_year = any(seen_types < 2)) %>%
  filter(seen_every_year == F) %>%
  select(-seen_types, -seen_every_year) %>%
  ungroup() %>%
  # mutate(vessel = as.factor(vessel),
  #        spp = as.factor(spp)) %>%
  as.data.frame() %>%
  nest(-survey, -survey_region, -knots) #note that this breaks if it's a tibble, should make a bit more flexible


# vasterize ---------------------------------------------------------------

set.seed(42)

if (vasterize == T) {
  # a = subset_fish_data$data[[1]] %>%
  #   group_by(year, spp) %>%
  #   summarise(ps = mean(catch_kg > 0)) %>%
  #   ggplot(aes(spp, ps)) +
  #   geom_col() +
  #   facet_wrap(~year) +
  #   coord_flip()

  vast_fish <- subset_fish_data %>%
    mutate(
      vasterized_data = purrr::pmap(
        list(
          region = survey_region,
          raw_data = data,
          n_x = knots
        ),
        vasterize_index,
        run_dir = run_dir,
        nstart = 100,
        obs_model = c(2, 0)
      )
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

skynet_data <- gfw_data %>%
  unnest() %>%
  group_by(survey, year, rounded_lat, rounded_lon) %>%
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
    port_engine_hours_interaction = dist_from_port * total_engine_hours,
    port_hours_interaction = dist_from_port * total_hours,
    port_numbers_interaction = dist_from_port * num_vessels,
    shore_numbers_interaction = dist_from_shore * num_vessels,
    shore_hours_interaction = dist_from_shore * total_hours,
    shore_engine_hours_interaction = dist_from_shore * total_engine_hours,
    large_vessels_interaction = num_vessels * mean_vessel_length
  ) %>%
  ungroup()


# pre process fishdata

knots  <- vast_fish %>%
  select(survey, vasterized_data) %>%
  mutate(knots = map(vasterized_data, c('spatial_list', 'loc_x_lat_long'))) %>%
  select(-vasterized_data)

# mutate(knot_plot =  map(knots, ~ quick_map(
#   .x,
#   lat_var = quo(approx_lat),
#   lon_var = quo(approx_long),
#   plot_var = quo(knot))))


# quick_map(knots$knots[[4]],lat_var = quo(approx_lat),
#           lon_var = quo(approx_long),
#           plot_var = quo(knot), min_lon = -160)

qmplot(approx_long,
       approx_lat,
       color = knot,
       data = knots$knots[[1]])


if (fished_only == T) {
  # include only fished species

  total_fish_data <- vast_fish %>%
    mutate(spatial_densities = map(vasterized_data, 'spatial_densities')) %>%
    select(survey, spatial_densities) %>%
    unnest() %>%
    filter(str_replace(species, '_', ' ') %in% fished_species$sci_name)

} else {
  total_fish_data <- vast_fish %>%
    mutate(spatial_densities = map(vasterized_data, 'spatial_densities')) %>%
    select(survey, spatial_densities) %>%
    unnest()
}


total_fish_data <- total_fish_data %>%
  group_by(survey, knot, year) %>%
  summarise(density = sum(density)) %>%
  ungroup() %>%
  nest(-survey, .key = fish_data)

if (survey_years_only == T) {

  survey_years <- fish_data %>%
    group_by(survey, Year) %>%
    summarise(year_surveyed = T)  %>%
    set_names(tolower)

total_fish_data <- total_fish_data %>%
  unnest() %>%
  left_join(survey_years, by = c('survey', 'year')) %>%
  filter(year_surveyed == T) %>%
  nest(-survey, .key = fish_data)
}

# aggregate data


skynet_data <- skynet_data %>%
  nest(-survey, .key = gfw_data) %>%
  left_join(total_fish_data, by = 'survey') %>%
  left_join(knots, by = 'survey')


if (clip_gfw == T) {
  skynet_data <- skynet_data %>%
    mutate(gfw_data = map2(knots, gfw_data, clip_gfw_to_fishdata))
}

skynet_data <- skynet_data %>%
  mutate(
    combo_data = pmap(
      list(
        gfw_data = gfw_data,
        fish_data = fish_data,
        fish_knots = knots,
        surveys = survey
      ),
      create_skynet_data,
      topo_data = topo_data,
      wind_data = wind_data,
      chl_data = chl_data,
      wave_data = wave_data,
      sst_data = sst_data
    )
  ) %>%
  select(survey, combo_data) %>%
  unnest()

missing_too_many <- skynet_data %>%
  map_df( ~ mean(is.na(.x))) %>%
  gather(variable, percent_missing) %>%
  arrange(desc(percent_missing)) %>%
  filter(percent_missing > max_percent_missing,
         !str_detect(variable, 'lag'), !str_detect(variable, 'density')) %>%
         {
           .$variable
         }

skynet_data <- skynet_data %>%
  select_(paste0('-', missing_too_many))

# skynet_data %>%
#   group_by(rounded_lat, rounded_lon) %>%
#   summarise(percent_missing = mean(is.na(mean_altitude))) %>%
#   ggplot(aes(rounded_lon, rounded_lat, fill = factor(percent_missing))) +
#   geom_tile()

# qmplot(
#   rounded_lon,
#   rounded_lat,
#   color = knot,
#   data = skynet_data %>% filter(rounded_lon > -175, is.na(density) == F, survey == 'wcghl_gfw')
# ) +
#   scale_color_viridis() +
#   facet_grid(survey~year, scales = 'free')



# plot plot plot plot -----------------------------------------------------

# plot distributions of data by year and survey

if (plot_data == T) {
  plot_vars <- skynet_data %>%
    select(-survey,-year,-rounded_lat,-rounded_lon) %>%
    colnames()

  plot_covariates <- function(dat, variable, run_dir) {
    p <- dat %>%
      ggplot() +
      geom_joy(aes_(
        x = as.name(variable),
        y = ~ year,
        group = ~ year,
        fill = ~ survey
      )) +
      facet_wrap( ~ survey, scales = 'free_x') +
      scale_fill_viridis_d(guide = F) +
      labs(title = variable)

    ggsave(filename = paste0(run_dir, variable, '_plot.pdf'), p, height = 8, width = 8)

  }

  walk(plot_vars,
       ~ plot_covariates(
         dat = skynet_data,
         variable = .x,
         run_dir = run_dir
       ))



  temp_maps <- skynet_data %>%
    nest(-survey)

  map_plot_foo <- function(survey, dat, run_dir, query_coastal_map = F) {
    a <-  dat %>%
      group_by(rounded_lat, rounded_lon, year) %>%
      summarise(
        density = unique(density),
        gfw_hours = unique(total_hours),
        knot = unique(knot)
      ) %>%
      gather('variable', 'value', density, gfw_hours) %>%
      group_by(variable, year) %>%
      mutate(value = value / max(value))


    a <-  a %>%
      dplyr::mutate(geometry = purrr::map2(rounded_lon, rounded_lat, ~ sf::st_point(x = c(.x, .y), dim = 'XY'))) %>%
      ungroup() %>%
      mutate(geometry = sf::st_sfc(geometry, crs =
                                     "+proj=longlat +datum=WGS84 +no_defs")) %>%
      sf::st_sf() %>%
      select(year, variable, value, geometry)


    if (query_coastal_map == T){

    # coast_map <-
    #   rnaturalearth::ne_download(scale = 'medium',
    #                              category = 'physical',
    #                              type = 'coastline') %>% sf::st_as_sf()


    global_map <-
      rnaturalearth::ne_download(scale = 'medium',
                                 category = 'physical',
                                 type = 'land') %>% sf::st_as_sf()


    save(file = "data/global_map.Rdata", global_map)
    } else {

     load(file = "data/global_map.Rdata")
    }

    bbox <- sf::st_bbox(a)

    var_map <- a %>%
      ggplot() +
      geom_sf(data = global_map, fill = 'grey60') +
      geom_sf(aes(color = value), size = 0.01, alpha = 0.5) +
      facet_grid(variable ~ year) +
      coord_sf(xlim = c(bbox['xmin'], bbox['xmax']),
               ylim = c(bbox['ymin'], bbox['ymax'])) +
      scale_color_viridis()

    ggsave(
      filename = paste0(run_dir, survey, '_map.pdf'),
      var_map,
      height = 10,
      width = 10
    )

  }

  walk2(
    temp_maps$survey,
    temp_maps$data,
    ~ map_plot_foo(
      survey = .x,
      dat = .y,
      run_dir = run_dir
    ),
    query_coastal_map = query_coastal_map
  )

  # map_plot_foo(temp_maps$survey[[1]], dat = temp_maps$data[[1]],
  #              run_dir = run_dir)


  # a %>%
  #   filter(survey == 'ebsbts_gfw') %>%
  #   ggplot() +
  #   geom_raster(aes(rounded_lon, rounded_lat, fill = density)) +
  #   geom_contour(aes(rounded_lon, rounded_lat, z = gfw_hours, colour = ..level..)) +
  #   facet_grid(. ~ year, scales = 'free') +
  #   scale_fill_gradient(low = 'white', high = 'red') +
  #   scale_color_viridis(option = 'D') +
  #   theme_classic()


}
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
    'density',
    'survey'
  ) # variables that should not be centered and scaled

do_scale <- skynet_names[!skynet_names %in% do_not_scale]

skynet_data <- skynet_data %>%
  ungroup() %>%
  purrrlyr::dmap_at(do_scale, ~ (.x - mean(.x, na.rm = T)) / (2 * sd(.x, na.rm = T))) #center and scale variables that should be

dep_var <- c('log_density')

never_ind_vars <-
  c(
    'log_density',
    'density',
    'survey',
    'year',
    'rounded_lat',
    'rounded_lon',
    'knot',
    'distance',
    'in_poly'
  )



tree_candidate_vars <-  skynet_names[!skynet_names %in% c(dep_var,
                                                                never_ind_vars) & str_detect(skynet_names, '_interaction')
 == F]


lm_candidate_vars <-  skynet_names[!skynet_names %in% c(dep_var,
                                                          never_ind_vars)]
models <- 'rf'

data_sources <- tibble(
lag_0_skynet_data = list(skynet_data %>%
  select(-contains('lag')) %>%
  na.omit()),
lag_1_skynet_data = list(skynet_data %>%
  select(-total_engine_hours_lag2, -total_engine_hours_lag3) %>%
  na.omit())) %>%
  gather(data_subset, data)

# lag_2_skynet_data <- skynet_data %>%
#   select(-total_engine_hours_lag3) %>%
#   na.omit()
#
# lag_3_skynet_data <- skynet_data %>%
#   select() %>%
#   na.omit()


generate_test_training <- function(dat, test_set, kfolds = 2, cut_year = 2014) {

  if (test_set == 'random') {
    out <- modelr::crossv_kfold(dat, k = kfolds, id = 'test_set') %>%
      mutate(test_set = 'random') %>%
      slice(1)

  }
  if (test_set == 'west_coast') {
    train <- dat %>%
      filter(!survey %in% c('wcgbts', 'wcghl'))

    test <- dat %>%
      filter(survey %in% c('wcgbts', 'wcghl'))

    out <-
      data_frame(
        train = list(resample(train, 1:nrow(train))),
        test = list(resample(test, 1:nrow(test))),
        test_set = 'west_coast'
      )

  }

  if (test_set == 'alaska') {
    train <- dat %>%
      filter(survey %in% c('wcgbts', 'wcghl'))

    test <- dat %>%
      filter(!survey %in% c('wcgbts', 'wcghl'))

    out <-
      data_frame(
        train = list(resample(train, 1:nrow(train))),
        test = list(resample(test, 1:nrow(test))),
        test_set = 'alaska'
      )

  }
  if (test_set == 'historic') {
    train <- dat %>%
      filter(year <= cut_year)

    test <- dat %>%
      filter(year > cut_year)

    out <-
      data_frame(
        train = list(resample(train, 1:nrow(train))),
        test = list(resample(test, 1:nrow(test))),
        test_set = paste0('year_gr_than_', cut_year)
      )

  }

  return(out)
}


test_train_data <- purrr::cross_df(list(
  test_sets = c('random','alaska','west_coast','historic'),
  data_subset = c('lag_0_skynet_data', 'lag_1_skynet_data')
)) %>%
  left_join(data_sources, by = 'data_subset')


test_train_data <- test_train_data %>%
  mutate(resampled = map2(data, test_sets, ~generate_test_training(dat = .x, test_set = .y))) %>%
  select(-data) %>%
  unnest()


# test_train_data <- modelr::crossv_kfold(lag_0_skynet_data, k = 1)

skynet_models <-  purrr::cross_df(list(
  dep_var = dep_var,
  # ind_vars = list(ind_vars[[6]]),
  # never_ind_vars = list(never_ind_vars),
  temp = list(test_train_data),
  model = models
)) %>%
  unnest(temp, .drop = F) #%>%
  # mutate(ind_vars = map(ind_vars, unlist))

# run models --------------------------------------------------------------


if (run_models == T) {

  sfm <- safely(fit_skynet)

  skynet_models <- skynet_models %>%
    mutate(fitted_model = pmap(
      list(
        dep_var = dep_var,
        model = model,
        training = train,
        testing = test
      ),
      sfm,
      fitcontrol_number = 10,
      fitcontrol_repeats = 10,
      never_ind_vars = never_ind_vars,
      tree_candidate_vars = tree_candidate_vars,
      lm_candidate_vars = lm_candidate_vars
    ))

  save(file = paste0(run_dir, 'skynet_models.Rdata'),
       skynet_models)

} else {
  load(file = paste0(run_dir, 'skynet_models.Rdata'))
}

# diagnose models ---------------------------------------------------------

mse_foo <- function(dat, dep_var, pred_var) {
  obs <- dat[,dep_var]

  pred <- dat[,pred_var]

  mse <- mean((obs - pred)^2)
}

diagnostic_plot_foo <- function(data,
                                r2,
                                dep_var,
                                test_region,
                                data_set){

  data %>%
    ggplot(aes_(as.name(dep_var), ~ pred)) +
    geom_abline(aes(intercept = 0, slope = 1), color = 'red', linetype = 2) +
    geom_point(alpha = 0.75) +
    labs(title = paste0('Test Region is ', test_region,'; R2 = ', r2),
         subtitle = paste0('Data Subset ', data_set))


}


skynet_models <- skynet_models %>%
  mutate(test_data = map(fitted_model, c('result', 'test_predictions'))) %>%
  mutate(
    mse = map2_dbl(test_data, dep_var, mse_foo, pred_var = 'pred'),
    var_y = map2_dbl(test_data, dep_var, ~ var(.x[, .y])),
    psuedo_r2 = round(1 - mse / var_y, 2)
  ) %>%
  mutate(test_plot = pmap(list(data = test_data,
                               r2 = psuedo_r2,
                               test_region = test_sets,
                               data_set = data_subset,
                               dep_var = dep_var),diagnostic_plot_foo))

save_foo <- function(test_plot,
                     test_region,
                     data_set,
                     run_dir) {

  ggsave(filename = paste0(run_dir,test_region,'-',data_set,'.pdf'), test_plot,
         height = 8, width = 8)


}

pwalk(list(test_plot = skynet_models$test_plot,
  test_region = skynet_models$test_sets,
                          data_set = skynet_models$data_subset),save_foo, run_dir = run_dir)

# wee <- map(skynet_models$fitted_model, c('result'))
#
# a <- wee[[3]]
#
# a$model$finalModel
#
# varImpPlot(a$model$finalModel)
#
# arg <- skynet_models$test[[3]] %>%
#   as_data_frame() %>%
#   na.omit() %>%
#   ungroup() %>%
#   mutate(ld_hat = predict(a$model$finalModel, newdata = .))
#
# arg %>%
#   ggplot(aes(log_density, ld_hat)) +
#   geom_point() +
#   geom_abline(aes(intercept = 0, slope = 1)) +
#   labs(x = 'Observed Multi-Species Density',
#        y = 'Predicted Multi-Species Density', title = 'Testing data omitted from model fitting') +
#   hrbrthemes::theme_ipsum()
#
#
# huh <- a$model$trainingData %>%
#   mutate(ld_hat = a$model$finalModel$predicted)
#
# huh %>%
#   ggplot(aes(.outcome, ld_hat)) +
#   geom_point(aes(color = mean_analysed_sst)) +
#   geom_abline(aes(intercept = 0, slope = 1))
