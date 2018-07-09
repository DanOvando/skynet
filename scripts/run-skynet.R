# run-skynet -------
# Author: Dan Ovando
# Project: skynet
# Summary: A wrapper script to fit and analyze models predicting abundance
# using GFW data

# prep run ----------------------------------------------------------

set.seed(42)
library(bigrquery)
library(lubridate)
library(FishData)
library(ggmap)
library(stringr)
library(viridis)
library(purrr)
library(purrrlyr)
library(modelr)
library(VAST)
library(TMB)
library(ggridges)
library(caret)
library(taxize)
library(gbm)
library(recipes)
library(sf)
library(rstan)
library(furrr)
library(tidyverse)


rstan::rstan_options(auto_write = TRUE)

filter <- dplyr::filter

functions <- list.files(here::here("functions"))

walk(functions, ~ here::here("functions", .x) %>% source()) # load local functions

run_name <- "d1.0"

run_description <-
  "dissertation default results"

run_dir <- file.path("results", run_name, "")

if (dir.exists(run_dir) == F) {
  dir.create(run_dir)
}

write(run_description, file = paste0(run_dir, "description.txt"))


# set section options (what to run) ---------------------------------------------------------

num_cores <- 3

run_models <- TRUE # fit gfw models to fishdata

tune_pars <- TRUE # pre-tune machine learning models

models <-
  c("ranger",
    "bagged_mars",
    "mars",
    "gbm",
    "structural",
    "engine_power",
    "hours")

# models <- c("mars")

vasterize <- F # run vast or load saved object

impute_missing <- F

plot_data <- F # plot covariates and maps

query_fishdata <- F # get trawl survey data or load saved

query_gfw <- F # get gfw data or load saved

query_environmentals <- F # get environmental data or load saved

query_mpas <- F # get mpa data or load saved

query_synonyms <- F # look up synonyms for each species

query_coastal_map <- F # get coastal map or load

query_prices <- F # get exvessel price data for observed species

round_environmentals <- T


# set run options (how to run) ---------------------------------------------------------

raw_fish <- F

dep_vars <- c("density","lag_density", "lag_economic_density","lag_economic_density")

gfw_vars_only <- T

fished_only <-
  T # only include fished species in predictive model fitting

unfished_only <- F # trial option fitting to unfished only

no_wcghl <- TRUE # get rid of WCGHL in model fitting

hack_zeros <- T # hack zeros into perfectly observed species

survey_years_only <-
  T # only include years that were actually surveyed

include_all_species <- T # include all species in VAST

bottom_gears_only <- T # only include bottom gear

clip_gfw <- T # clip GFW data to extent of trawl surveys

vars_to_drop <- c("wind_speed", "wind_angle")

gfw_project <- "ucsb-gfw"

gfw_dataset <- "skynet"

max_percent_missing <- 0.2

lat_lon_res <-
  0.25 # round data to intervels of 0.25 degrees lat lon, as in 56.25

res <- 1 / lat_lon_res

min_times_seen <- 9

min_dist_from_shore <- 500

min_speed <- 0.01

max_speed <- 20

min_year <- 2010

species_list <-
  c("Gadus_chalcogrammus",
    "Hippoglossoides_elassodon",
    "Gadus_macrocephalus")

# survey_list <- c('ebsbts', 'goabts', 'aibts', 'wcgbts', 'wcghl')
survey_list <- c("ebsbts", "goabts", "aibts", "wcgbts", "wcghl")

knots_per_survey <- tribble( ~ survey,
                             ~ knots,
                             # --/--
                             "ebsbts",
                             300,
                             "goabts",
                             400,
                             "aibts",
                             200,
                             "wcgbts",
                             1000,
                             "wcghl",
                             100)

# get fishdata ------------------------------------------------------------

if (query_fishdata == T) {
  fishdata_names <- tribble(
    ~ survey,
    ~ survey_name,
    ~ survey_region,
    # --/--
    "ebsbts",
    "eastern bering sea bottom trawl survey",
    "eastern_bering_sea",
    "wcgbts",
    "west coast groundfish bottom trawl survey",
    "california_current",
    "wcghl",
    "West_coast_groundfish_hook_and_line",
    "wcghl_domain",
    "goabts",
    "gulf of alaska bottom trawl survey",
    "gulf_of_alaska",
    "aibts",
    "aluetian islands bottom trawl survey",
    "aleutian_islands"
  )

  sdcr <- safely(download_catch_rates)

  fishdata_names$survey <- fishdata_names$survey %>% toupper()

  fish_data <- fishdata_names %>%
    mutate(survey_data = map(
      survey,
      ~ download_catch_rates(.x, species_set = 50,
                             error_tol = 1e-2)
    ))

  fish_data$survey_data <-
    fish_data$survey_data %>% map(`[`, c("Sci", "Year", "TowID", "Lat", "Long", "Wt", "Num"))

  fish_data <- unnest(fish_data) %>%
    mutate(survey = tolower(survey))

  save(file = "data/fish_data.Rdata", fish_data)
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
  left_join(knots_per_survey, by = "survey")

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

survey_names <- fish_data$survey %>% tolower() %>% unique()

survey_years <- fish_data %>%
  group_by(survey, Year) %>%
  summarise() %>%
  mutate(surveyed_year = TRUE) %>%
  set_names(tolower) %>%
  ungroup()


# query_gfw ---------------------------------------------------------------


if (query_gfw == T) {
  fishing_connection <- DBI::dbConnect(
    bigrquery::dbi_driver(),
    project = gfw_project,
    dataset = gfw_dataset,
    allowLargeResults = T
  )



  temp_gfw_data <- list()

  for (i in 1:nrow(survey_bbox)) {
    bbox <- survey_bbox %>%
      slice(i) %>%
      select(contains("lat"), contains("lon")) %>%
      gather(var, val) %>%
      {
        .$val
      }

    gfw_query <- glue::glue_sql(
      "SELECT
      *
      FROM (
      SELECT
      YEAR(timestamp) year,
      MONTH(timestamp) month,
      mmsi,
      FLOOR(lat*{res})/{res} rounded_lat,
      FLOOR(lon*{res})/{res} rounded_lon,
      SUM(IF( nnet_score == 1, hours, 0)) AS total_hours,
      AVG(distance_from_shore ) AS mean_distance_from_shore,
      AVG(distance_from_port) AS mean_distance_from_port,
      INTEGER(REGEXP_REPLACE( IF(REGEXP_EXTRACT(regions,'\"(mparu:.*?)\"') CONTAINS \".\", LEFT(REGEXP_EXTRACT(regions,'\"(mparu:.*?)\"'),INSTR(REGEXP_EXTRACT(regions,'\"(mparu:.*?)\"'),\".\")-1),REGEXP_EXTRACT(regions,'\"(mparu:.*?)\"')), '[^0-9 ]','')) mparu,
      INTEGER(REGEXP_REPLACE( IF(REGEXP_EXTRACT(regions,'\"(mpant:.*?)\"') CONTAINS \".\", LEFT(REGEXP_EXTRACT(regions,'\"(mpant:.*?)\"'),INSTR(REGEXP_EXTRACT(regions,'\"(mpant:.*?)\"'),\".\")-1),REGEXP_EXTRACT(regions,'\"(mpant:.*?)\"')), '[^0-9 ]','')) mpant
      FROM
      [world-fishing-827:gfw_research.nn]
      WHERE
      (distance_from_shore > ({min_dist_from_shore})
      AND (implied_speed > ({min_speed}) AND implied_speed < ({max_speed})))
      AND lat >= ({bbox[1]})
      AND lat <= ({bbox[2]})
      AND lon > ({bbox[3]})
      AND lon <= ({bbox[4]})
      AND _PARTITIONTIME BETWEEN TIMESTAMP('2010-01-01')
      AND TIMESTAMP('2017-12-31')
      AND seg_id IN (
      SELECT
      seg_id
      FROM
      [world-fishing-827:gfw_research.good_segments])
      GROUP BY
      mmsi,
      year,
      month,
      rounded_lat,
      rounded_lon,
      mparu,
      mpant) a
      LEFT JOIN (
      SELECT
      year,
      mmsi,
      shipname,
      inferred_engine_power,
      inferred_length,
      inferred_tonnage,
      known_engine_power,
      known_length,
      known_tonnage,
      best_label,
      known_label,
      inferred_label_allyears,
      inferred_sublabel_allyears,
      good_shiptype_msg,
      on_fishing_list_nn
      FROM
      [world-fishing-827:gfw_research.vessel_info]
      WHERE
      (on_fishing_list_nn IS TRUE OR on_fishing_list IS TRUE)
      AND spoofing_factor < 1.01
      AND offsetting IS NULL) b
      ON
      a.mmsi = b.mmsi
      AND a.year = b.year",
      bbox = bbox,
      min_dist_from_shore = min_dist_from_shore,
      min_speed = min_speed,
      max_speed = max_speed,
      .con = fishing_connection
      )

    temp_gfw_data[[i]] <-
      bigrquery::query_exec(gfw_query, project = "ucsb-gfw", max_pages = Inf)
  }

  gfw_data <- data_frame(survey = survey_bbox$survey) %>%
    mutate(gfw_data = temp_gfw_data) %>%
    unnest() %>%
    select(-b_mmsi, -b_year) %>%
    set_names(str_replace_all(colnames(.), "(a_)|(b_)", "")) %>%
    nest(-survey)


  save(file = "data/gfw_data.Rdata", gfw_data)
} else {
  load(file = "data/gfw_data.Rdata")
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
    cap <- 1.1 * max(inferred, na.rm = T)
  }

  out <- ifelse(is.na(known), pmin(inferred, cap), known)
}

gfw_data <-
  gfw_data %>% # merge known and inferred vessel characteristics
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
        desired_data = "chl-a",
        date_interval = 1,
        space_interval = 1,
        runit = 0.25
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
      labs(caption = "SST(c)")
  }

  # chl_data <- chl_data %>%
  #   mutate(plots = trelliscopejs::map_plot(chl_a, chl_plot_foo))

  # a <- chl_data$chl_a[[4]]

  # qmplot(rlon, rlat, color = mean_chlorophyll, data = a)

  save(file = "data/chl_env_data.Rdata", chl_data)


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
        desired_data = "sst",
        date_interval = 14,
        space_interval = 25,
        runit = 0.25
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
      labs(caption = "SST(c)")
  }

  # sst_data <- sst_data %>%
  #   mutate(plots = trelliscopejs::map_plot(sst, sst_plot_foo))

  save(file = "data/sst_env_data.Rdata", sst_data)

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
        desired_data = "waves",
        date_interval = 24 * 14,
        # 1 hour interval
        space_interval = 0.5,
        runit = 0.25
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

  save(file = "data/wave_env_data.Rdata", wave_data)


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
        desired_data = "wind",
        date_interval = 1,
        # 1 hour interval
        space_interval = 1,
        runit = 0.25
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

  save(file = "data/wind_env_data.Rdata", wind_data)

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
        desired_data = "topography",
        date_interval = 1,
        # 1 hour interval
        space_interval = 1,
        runit = 0.25
      )
    )

  topo_data <- topo_data %>%
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

  save(file = "data/topo_env_data.Rdata", topo_data)
} else {
  data_files <- list.files("data/")

  data_files <-
    paste0("data/", data_files[str_detect(data_files, "env_data.Rdata")])

  for (i in seq_along(data_files)) {
    load(data_files[i])
  }
}

if (round_environmentals == T) {
  data_files <- list.files("data/")

  data_files <-
    paste0("data/", data_files[str_detect(data_files, "env_data.Rdata")])

  enviro_types <- str_extract(data_files, "(?<=/).*(?=_)")

  enviro_types <- str_replace(enviro_types, "env", "data")

  for (i in enviro_types) {
    desired_data <- str_extract(i, ".*(?=_)")

    temp <- round_enviro(get(i),
                         desired_data = desired_data, runit = lat_lon_res)

    assign(i, temp)
  }

  topo_data <- topo_data %>%
    mutate(altitude = map(
      altitude,
      ~ .x %>%
        mutate(m_below_sea_level = -1 * mean_altitude)
    ))
}


# query mpas --------------------------------------------------------------

if (query_mpas == T) {
  mpa_ids <-
    data_frame(mpa_type = c("no_take_mpas_id", "restricted_use_mpas_id"))

  fishing_connection <- DBI::dbConnect(bigrquery::dbi_driver(),
                                       project = gfw_project,
                                       dataset = "regions_ids")

  mpa_ids <- mpa_ids %>%
    mutate(data = map(mpa_type, ~ fishing_connection %>%
                        tbl(.x) %>%
                        collect(n = Inf))) %>%
    unnest()

  mpa_ids <- mpa_ids %>%
    mutate(
      mpa_type = str_split(region_id, ":", simplify = T)[, 1],
      mpa_id = str_split(region_id, ":", simplify = T)[, 2]
    )

  save(file = "data/mpa_ids.Rdata", mpa_ids)
} else {
  load(file = "data/mpa_ids.Rdata")
}



# data processing ---------------------------------------------------------

# decide species that go in
# for now, query FAO for species fished... this is hacky though since means that depends on static data...

if (query_synonyms == T) {
  fao <-
    data.table::fread("~/Box Sync/Databases/FAO/FAO_Capture_1950to2012.csv",
                      header = T) %>%
    as_data_frame() %>%
    filter(isscaap_group < 60) # get fao data

  fao_species <- unique(fao$sciname) # unique fao species

  species_seen <- fish_data$Sci %>%
    unique() %>%
    str_replace("_", " ") # species in fishdata

  safe_syns <- safely(taxize::synonyms)

  species_data <- data_frame(sci_name = species_seen) %>%
    filter(is.na(sci_name) == F) %>%
    mutate(worms = map(sci_name, taxize::get_wormsid, ask = F)) %>%
    mutate(worms_id = map(worms, `[[`, 1) %>% as.character()) %>% # extract the worms id
    filter(is.na(worms_id) == F) %>%
    mutate(known_synonyms = map(worms_id %>% as.numeric(), safe_syns, db = "worms"))

  species_data <- species_data %>%
    mutate(syn_results = map(known_synonyms, "result")) %>%
    mutate(syns = map(syn_results, `[[`, 1)) %>%
    mutate(sci_syns = map(syns, "scientificname")) # extract synonyms


  check_fao <- function(sci, syns, faos) {
    in_fao <- sci %in% faos | any(syns %in% faos)
  }

  fished_species <- species_data %>%
    mutate(in_fao = map2_lgl(sci_name, sci_syns, check_fao, fao_species)) %>%
    select(sci_name, sci_syns, in_fao) %>%
    filter(in_fao == T) # find species whose primary name or synonyms are in fao database

  save(file = "data/fished_species.Rdata", fished_species)
} else {
  load(file = "data/fished_species.Rdata")
}


if (query_prices == T) {
  # process exvessel price data

  fish_prices <-
    readr::read_csv("~/Box Sync/Databases/Exvessel Prices/Exvessel Price Database.csv") %>%
    filter(Year > 2010) %>%
    group_by(scientific_name) %>%
    summarise(mean_exvessel_price = mean(exvessel, na.rm = T))

  fish_prices <- fish_prices %>%
    rename(species = scientific_name) %>%
    mutate(species = str_replace(species, " ", "_"))

  species_seen <- fish_data$Sci %>%
    unique() %>%
    str_replace("_", " ") # species in fishdata

  safe_syns <- safely(taxize::synonyms)

  species_synonyms <- data_frame(sci_name = species_seen) %>%
    filter(is.na(sci_name) == F) %>%
    mutate(worms = map(sci_name, taxize::get_wormsid, ask = F)) %>%
    mutate(worms_id = map(worms, `[[`, 1) %>% as.character()) %>% # extract the worms id
    filter(is.na(worms_id) == F) %>%
    mutate(known_synonyms = map(worms_id %>% as.numeric(), safe_syns, db = "worms"))

  species_synonyms <- species_synonyms %>%
    mutate(syn_results = map(known_synonyms, "result")) %>%
    mutate(syns = map(syn_results, `[[`, 1)) %>%
    mutate(sci_syns = map(syns, "scientificname")) %>% # extract synonyms
    select(sci_name, sci_syns) %>%
    mutate(species = str_replace(sci_name, " ", "_")) %>%
    select(-sci_name)


  price_foo <- function(sci_name, sci_syns, fish_prices) {
    joint_names <- data_frame(species = c(sci_name, sci_syns)) %>%
      mutate(species = stringr::str_replace_all(species, " ", "_"))

    syn_prices <- joint_names %>%
      left_join(fish_prices, by = "species")

    mean_price <- mean(syn_prices$mean_exvessel_price, na.rm = T)
  }

  species_prices <- species_synonyms %>%
    mutate(mean_exvessel_price = map2_dbl(species, sci_syns, price_foo, fish_prices = fish_prices)) %>%
    select(species, mean_exvessel_price)

  save(file = "data/fish_prices.Rdata", species_prices)
} else {
  load(file = "data/fish_prices.Rdata")
}

# MAJOR HACK TO DEAL WITH PERFECTLY OBSERVED SPECIES FIX FIX FIX FIX ONCE VAST IS FIXED
# can't run species that are perfectly observed at this point, so for now going to do a
# MAJOR hack and for species without any zeros, convert values below some quantile threshold to zeros

if (hack_zeros == T) {
  hack_foo <- function(x) {
    min_wt <- min(x$Wt)

    if (min_wt > 0) {
      min_quantile <- quantile(x$Wt, 0.025) %>% as.numeric()

      x$Wt[x$Wt <= min_quantile] <- 0
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
      (
        inferred_label_allyears == "trawlers" |
          inferred_sublabel_allyears == "pots_and_traps" |
          str_detect(inferred_sublabel_allyears, "set_") |
          inferred_sublabel_allyears == "trawlers"
      ) &
        !inferred_label_allyears %in% c("cargo_or_tanker", "passenger")
    ) %>%
    nest(-survey)
}

gfw_data <- gfw_data %>%
  filter(survey %in% survey_list)

if (include_all_species == T) {
  species_list <- unique(fish_data$Sci)
}


subset_fish_data <- fish_data %>%
  filter(survey %in% survey_list, Sci %in% species_list) %>%
  dplyr::rename(Lon = Long,
                Catch_KG = Wt,
                spp = Sci) %>%
  mutate(AreaSwept_km2 = 1, Vessel = "missing") %>%
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
  mutate(times_seen = sum(catch_kg > 0)) %>%
  group_by(survey, spp) %>%
  mutate(seen_every_year = all(times_seen > min_times_seen)) %>%
  filter(seen_every_year == T) %>%
  select(-times_seen, -seen_every_year) %>%
  ungroup() %>%
  as.data.frame() %>%
  nest(-survey, -survey_region, -knots,-spp) # note that this breaks if it's a tibble, should make a bit more flexible



# vasterize ---------------------------------------------------------------

set.seed(42)

if (vasterize == T) {

  vast_fish <- subset_fish_data %>%
    mutate(
      vasterized_data = purrr::pmap(
        list(
          region = survey_region,
          raw_data = data,
          n_x = knots,
          spp = spp
        ),
        safely(vasterize_index),
        version = "VAST_v4_1_0",
        run_dir = run_dir,
        nstart = 200,
        obs_model = c(2, 0),
        epsilon1 = 1
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

  error <- map(vast_fish$vasterized_data, "error")

  vast_fish$vasterized_data <-
    map(vast_fish$vasterized_data, "result")

  save(file = here::here("data", "vast_fish.Rdata"),
       vast_fish)
} else {
  load(file = here::here("data", "vast_fish.Rdata"))
}

# build database ----------------------------------------------------------

# pre process fishdata

vast_worked <- !map_lgl(vast_fish$vasterized_data, is.null)

vast_fish <- vast_fish[vast_worked, ]

knots <- vast_fish %>%
  select(survey, vasterized_data) %>%
  mutate(knots = map(vasterized_data, c("spatial_list", "loc_x_lat_long"))) %>%
  select(-vasterized_data) %>%
  unnest() %>%
  unique() %>%
  nest(-survey, .key = "knots")

# prepare candidate data streams

candidate_data <-
  tribble(
    ~ fished_only,
    ~ unfished_only,
    ~ survey_months_only,
    ~ raw_fish,
    T,
    F,
    F,
    F,
    T,
    F,
    T,
    F,
    F,
    T,
    T,
    F
  ) %>%
  mutate(
    data = pmap(
      list(
        fished_only = fished_only,
        unfished_only = unfished_only,
        survey_months_only = survey_months_only,
        raw_fish = raw_fish
      ),
      prepare_data,
      vast_fish = vast_fish,
      species_prices = species_prices,
      vars_to_drop = vars_to_drop,
      survey_years = survey_years
    )
  ) %>%
  mutate(
    total_fish_data = map(data, "total_fish_data"),
    skynet_data = map(data, "skynet_data")
  ) %>%
  select(-data)


save(file = here::here("results", run_name, "skynet_data.Rdata"),
     candidate_data)

# plot data -----------------------------------------------------

# plot distributions of data by year and survey

if (plot_data == T) {
  plot_vars <- candidate_data$skynet_data[[1]] %>%
    select(-survey, -year, -rounded_lat, -rounded_lon) %>%
    colnames()

  plot_covariates <- function(dat, variable, run_dir) {
    p <- dat %>%
      ggplot() +
      geom_density_ridges(aes_(
        x = as.name(variable),
        y = ~ year,
        group = ~ year,
        fill = ~ survey
      )) +
      facet_wrap(~ survey, scales = "free_x") +
      scale_fill_viridis_d(guide = F) +
      labs(title = variable)

    ggsave(
      filename = paste0(run_dir, variable, "_plot.pdf"),
      p,
      height = 8,
      width = 8
    )
  }

  walk(
    plot_vars,
    ~ plot_covariates(
      dat = candidate_data$skynet_data[[1]],
      variable = .x,
      run_dir = run_dir
    )
  )


  temp_maps <- candidate_data$skynet_data[[1]] %>%
    nest(-survey)

  map_plot_foo <-
    function(survey,
             dat,
             run_dir,
             query_coastal_map = F) {
      a <- dat %>%
        group_by(rounded_lat, rounded_lon, year) %>%
        summarise(
          density = unique(density),
          gfw_hours = unique(total_hours),
          knot = unique(knot)
        ) %>%
        gather("variable", "value", density, gfw_hours) %>%
        group_by(variable, year) %>%
        mutate(value = value / max(value))


      a <- a %>%
        dplyr::mutate(geometry = purrr::map2(rounded_lon, rounded_lat, ~ sf::st_point(x = c(.x, .y), dim = "XY"))) %>%
        ungroup() %>%
        mutate(geometry = sf::st_sfc(geometry, crs =
                                       "+proj=longlat +datum=WGS84 +no_defs")) %>%
        sf::st_sf() %>%
        select(year, variable, value, geometry)


      if (query_coastal_map == T) {
        # coast_map <-
        #   rnaturalearth::ne_download(scale = 'medium',
        #                              category = 'physical',
        #                              type = 'coastline') %>% sf::st_as_sf()

        global_map <-
          rnaturalearth::ne_download(scale = "medium",
                                     category = "physical",
                                     type = "land") %>%
          sf::st_as_sf()


        save(file = "data/global_map.Rdata", global_map)
      } else {
        load(file = "data/global_map.Rdata")
      }

      bbox <- sf::st_bbox(a)

      var_map <- a %>%
        ggplot() +
        geom_sf(data = global_map, fill = "grey60") +
        geom_sf(aes(color = value), size = 0.01, alpha = 0.5) +
        facet_grid(variable ~ year) +
        coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
                 ylim = c(bbox["ymin"], bbox["ymax"])) +
        scale_color_viridis()

      ggsave(
        filename = paste0(run_dir, survey, "_map.pdf"),
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
} # close plot chunk

# prepare models ----------------------------------------------------------


never_ind_vars <-
  c(
    "biomass",
    "log_biomass",
    "mean_knot_area",
    "log_density",
    "density",
    "lag_density",
    "log_lag_density",
    "lag_economic_density",
    "log_lag_economic_density",
    "cs_log_density",
    "cs_density",
    "survey",
    "year",
    "rounded_lat",
    "rounded_lon",
    "knot",
    "distance",
    "in_poly",
    "mean_altitude",
    "aggregate_price",
    "vessel_hours",
    "total_engine_hours",
    "aggregate_price.x",
    "aggregate_price.y",
    "index",
    "surveyed_year",
    "relative_density",
    "log_economic_density",
    "economic_density"
  )


skynet_names <- colnames(candidate_data$skynet_data[[1]])

gfw_vars <-
  skynet_names[str_detect(skynet_names,
                          "(vessel)|(length)|(hours)|(engine)|(fishing)|(numbers)")]

tree_candidate_vars <- skynet_names[!skynet_names %in% c(dep_vars,
                                                         never_ind_vars) &
                                      str_detect(skynet_names, "_interaction")
                                    == F &
                                      str_detect(skynet_names, "lag") == F &
                                      !str_detect(skynet_names, "any_fishing")]

gfw_only_tree_candidate_vars <-
  c(tree_candidate_vars[tree_candidate_vars %in% gfw_vars], "random_var")

enviro_only_tree_candidate_vars <-
  c(tree_candidate_vars[!tree_candidate_vars %in% gfw_only_tree_candidate_vars &
                          !tree_candidate_vars %in% (never_ind_vars)], "random_var")

lm_candidate_vars <- skynet_names[!skynet_names %in% c(dep_vars,
                                                       never_ind_vars)]

basic_skynet_data <- candidate_data %>%
  filter(fished_only == T,
         unfished_only == FALSE,
         survey_months_only == T) %>% {
           .$skynet_data[[1]]
         }

basic_all_months_skynet_data <- candidate_data %>%
  filter(fished_only == T,
         unfished_only == FALSE,
         survey_months_only == F) %>% {
           .$skynet_data[[1]]
         }

unfished_skynet_data <- candidate_data %>%
  filter(fished_only == FALSE,
         unfished_only == TRUE,
         survey_months_only == TRUE) %>% {
           .$skynet_data[[1]]
         }

delta_skynet <- basic_skynet_data %>%
  select(-contains("lag")) %>%
  gather(
    variable,
    value,-survey,
    -year,
    -rounded_lat,
    -rounded_lon,
    -dist_from_port,
    -dist_from_shore,
    -no_take_mpa,
    -restricted_use_mpa,
    -any_fishing,
    -m_below_sea_level,
    -random_var,
    -knot,
    -distance,
    -aggregate_price
  )

delta_vars <- unique(delta_skynet$variable)

delta_candidate_vars <-
  c(tree_candidate_vars[tree_candidate_vars %in% delta_vars],
    "random_var")

delta_skynet <- delta_skynet %>%
  group_by(survey, variable, rounded_lat, rounded_lon) %>%
  arrange(year) %>%
  mutate(lag_value = lag(value, 1)) %>%
  mutate(delta_value = value - lag_value) %>%
  select(-value, -lag_value) %>%
  spread(variable, delta_value) %>%
  na.omit() %>%
  ungroup()


# merge candidate data sources
data_sources <- tibble(
  skynet = list(basic_skynet_data %>%
                  na.omit()),
  unfished_skynet = list(unfished_skynet_data %>%
                           na.omit()),
  all_months_skynet_data  = list(
    basic_all_months_skynet_data %>%
    na.omit()
  ),
  delta_skynet = list(delta_skynet),
  skynet_25km = list(
    rescale_data(basic_skynet_data, resolution = 25) %>%
      na.omit()
  ),
  skynet_100km = list(
    rescale_data(basic_skynet_data, resolution = 100) %>%
      na.omit()
  )
) %>%
  gather(data_subset, data)


test_train_data <- purrr::cross_df(list(
  test_sets = c(
    "random",
    "alaska",
    "west_coast",
    "historic",
    "random_alaska",
    "random_west_coast",
    "spatial_alaska",
    "spatial_west_coast",
    "historic_alaska",
    "historic_west_coast",
    "california",
    "wa-or"
  ),
  data_subset = data_sources$data_subset
)) %>%
  left_join(data_sources, by = "data_subset")


test_train_data <- test_train_data %>%
  mutate(resampled = map2(data, test_sets, ~ generate_test_training(dat = .x, test_set = .y))) %>%
  select(-data) %>%
  unnest()


# test_train_data <- modelr::crossv_kfold(lag_0_skynet_data, k = 1)


skynet_models <- purrr::cross_df(list(
  dep_var = (dep_vars),
  temp = list(test_train_data),
  model = models,
  weight_surveys = c(F),
  variables = c("gfw_only", "enviro_only", "gfw_and_enviro")
)) %>%
  unnest(temp, .drop = F) %>%
  filter(!(data_subset == 'delta_skynet' &
             model == 'structural')) %>%
  # filter(!(model == 'structural' & data_subset != 'uncentered_skynet')) %>%
  filter(!(
    model %in% c('ranger', 'gbm') &
      data_subset == 'uncentered_skynet'
  )) %>%
  filter(!(weight_surveys == T &
             !(test_sets %in% c(
               "random", "historic"
             ))))

# run models --------------------------------------------------------------

if (run_models == T) {
  if (tune_pars == T) {

    prepped_train <- skynet_models %>%
      filter(
        data_subset == "skynet",
        test_sets %in% c("california"),
        model %in% c("ranger", "gbm", "mars", "bagged_mars"),
        dep_var == "density"
      )  %>%
      group_by(model, variables,test_sets) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(candidate_vars = ifelse(
        str_detect(.$data_subset, "delta"),
        list(delta_candidate_vars),
        ifelse(
          variables == "gfw_only",
          list(gfw_only_tree_candidate_vars),
          ifelse(
            variables == "enviro_only",
            list(enviro_only_tree_candidate_vars),
            list(tree_candidate_vars)
          )
        )
      )) %>%
      mutate(
        fitted_model = pmap(
          list(
            data_subset = data_subset,
            dep_var = dep_var,
            model_name = model,
            training = train,
            testing = test,
            tree_candidate_vars = candidate_vars
          ),
          prep_train,
          fitcontrol_number = 10,
          fitcontrol_repeats = 2,
          never_ind_vars = never_ind_vars,
          tune_model = T,
          cores = num_cores,
          data_sources = data_sources
        )
      )

    tuned_pars <- prepped_train %>%
      ungroup() %>%
      select(model, test_sets, variables,fitted_model) %>%
      mutate(tuned_pars = map(fitted_model, "tuned_pars"))


    kfold_preds <- prepped_train %>%
      ungroup() %>%
      select(model, variables,fitted_model) %>%
      mutate(kfold_preds = map(fitted_model, "kfold_preds"))

    # kfold_preds <- map(prepped_train$fitted_model, "kfold_preds") %>%
    #   set_names(prepped_train$model)

    saveRDS(file = paste0(run_dir, "tuned_pars.RDS"),
            tuned_pars)

    saveRDS(file = paste0(run_dir, "kfold_preds.RDS"),
            kfold_preds)

  } else {
    tuned_pars <- readRDS(file = paste0(run_dir, "tuned_pars.RDS"))
  }

  sfm <- safely(fit_skynet)

  # future::plan(future::multiprocess, workers = num_cores)


  skynet_models <- skynet_models %>%
    # filter(data_subset == "skynet",
    #        test_sets == "random") %>%
    # group_by(model, dep_var) %>%
    # slice(1) %>%
    ungroup() %>%
    mutate(candidate_vars = ifelse(
      str_detect(.$data_subset, "delta"),
      list(delta_candidate_vars),
      ifelse(
        variables == "gfw_only",
        list(gfw_only_tree_candidate_vars),
        ifelse(
          variables == "enviro_only",
          list(enviro_only_tree_candidate_vars),
          list(tree_candidate_vars)
        )
      )
    )) %>%
    mutate(
      fitted_model = pmap(
        list(
          data_subset = data_subset,
          dep_var = dep_var,
          model_name = model,
          training = train,
          testing = test,
          tree_candidate_vars = candidate_vars,
          variable = variables
        ),
        sfm,
        fitcontrol_number = 1,
        fitcontrol_repeats = 1,
        never_ind_vars = never_ind_vars,
        tune_model = T,
        cores = num_cores,
        data_sources = data_sources,
        tuned_pars = tuned_pars
      )
    )

  print('ran models')

  save(file = paste0(run_dir, "skynet_models.Rdata"),
       skynet_models)

  print('saved models')

} else {
  load(file = paste0(run_dir, "skynet_models.Rdata"))
}

# diagnose models ---------------------------------------------------------



skynet_models <- skynet_models %>%
  mutate(index = 1:nrow(.)) %>%
  mutate(error = map(fitted_model, "error")) %>%
  mutate(no_error = map_lgl(error, is.null)) %>%
  filter(no_error)



diagnostic_plot_foo <- function(data,
                                r2,
                                dep_var,
                                test_region,
                                train_region,
                                data_set) {
  data %>%
    ggplot(aes_(as.name(dep_var), ~ pred)) +
    geom_abline(aes(intercept = 0, slope = 1),
                color = "red",
                linetype = 2) +
    geom_point(alpha = 0.75) +
    labs(
      title = paste0(
        "tested on ",
        test_region,
        "- trained on ",
        train_region,
        "; R2 = ",
        r2
      ),
      subtitle = paste0("data subset is: ", data_set)
    )
} # close diagnostic functions


skynet_models <- skynet_models %>%
  mutate(
    test_data = map(fitted_model, c("result", "test_predictions")),
    training_data = map(fitted_model, c("result", "training_predictions"))
  ) %>%
  select(-fitted_model) %>%
  mutate(
    var_y = map2_dbl(test_data, dep_var, ~ var(.x[, .y])),
    r2 = map2_dbl(test_data, dep_var, ~ yardstick::rsq(.x, .y, "pred")),
    rmse = map2_dbl(test_data, dep_var, ~ yardstick::rmse(.x, .y, "pred")),
    ccc = map2_dbl(test_data, dep_var, ~ yardstick::ccc(.x, .y, "pred")),
    r2_training = map2_dbl(training_data, dep_var, ~ yardstick::rsq(.x, .y, "pred")),
    rmse_training = map2_dbl(training_data, dep_var, ~ yardstick::rmse(.x, .y, "pred")),
    ccc_training = map2_dbl(training_data, dep_var, ~ yardstick::ccc(.x, .y, "pred"))
  ) %>%
  mutate(
    test_plot = pmap(
      list(
        data = test_data,
        r2 = r2,
        test_region = test_sets,
        train_region = train_set,
        data_set = data_subset,
        dep_var = dep_var
      ),
      diagnostic_plot_foo
    ),
    training_plot = pmap(
      list(
        data = training_data,
        r2 = r2_training,
        test_region = paste0(test_sets),
        train_region = paste0(train_set, "- training plot"),
        data_set = data_subset,
        dep_var = dep_var
      ),
      diagnostic_plot_foo
    )
  )


skynet_models %>%
  ggplot(aes(model, r2_training, color = dep_var)) +
  geom_point() +
  coord_flip()

print('processed models')


