# run-skynet -------
# Author: Dan Ovando
# Project: skynet
# Summary: A wrapper script to fit and analyze models predicting abundance
# using GFW data

# prep run ----------------------------------------------------------

set.seed(42)
library(bigrquery)
library(lubridate)
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
library(ggridges)
library(modelr)
library(caret)
library(taxize)
library(gbm)
library(recipes)
library(rlang)
library(tidyverse)

demons::load_functions("functions")

run_name <- "v1.0"

run_description <- "development of full paper"

run_dir <- file.path("results", run_name, "")

if (dir.exists(run_dir) == F) {
  dir.create(run_dir)
}

write(run_description, file = paste0(run_dir, "description.txt"))


# set section options (what to run) ---------------------------------------------------------

run_models <- F # fit statistical models to data

vasterize <- F # run vast or load saved object

raw_fish <- F

plot_data <- F # plot covariates and maps

query_fishdata <- F # get trawl survey data or load saved

query_gfw <- F # get gfw data or load saved

min_times_seen <- 9

min_dist_from_shore <- 500

min_speed <- 0.01

max_speed <- 20

query_environmentals <- F # get environmental data or load saved

query_mpas <- F # get mpa data or load saved

query_synonyms <- F # look up synonyms for each species

query_coastal_map <- F # get coastal map or load

query_prices <- F # get exvessel price data for observed species

# set run options (how to run) ---------------------------------------------------------

center_and_scale <- F # center and scale some data?

gfw_vars_only <- F

fished_only <-
  T # only include fished species in predictive model fitting

hack_zeros <- T # hack zeros into perfectly observed species

survey_years_only <-
  T # only include years that were actually surveyed

include_all_species <- T # include all species in VAST

bottom_gears_only <- T # only include bottom gear

clip_gfw <- T # clip GFW data to extent of trawl surveys

vars_to_drop <- c("wind_speed", "wind_angle")

gfw_project <- "ucsb-gfw"

gfw_dataset <- "skynet"

max_percent_missing <- 0.25

lat_lon_res <-
  0.25 # round data to intervels of 0.25 degrees lat lon, as in 56.25

min_year <- 2012

species_list <-
  c(
    "Gadus_chalcogrammus",
    "Hippoglossoides_elassodon",
    "Gadus_macrocephalus"
  )

# survey_list <- c('ebsbts', 'goabts', 'aibts', 'wcgbts', 'wcghl')
survey_list <- c("ebsbts", "goabts", "aibts", "wcgbts", "wcghl")

knots_per_survey <- tribble(
  ~ survey,
  ~ knots,
  # --/--
  "ebsbts",
  200,
  "goabts",
  200,
  "aibts",
  200,
  "wcgbts",
  1000,
  "wcghl",
  100
)

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

  fish_data <- fishdata_names %>%
    mutate(survey_data = map(
      survey,
      ~ download_catch_rates(
        .x, species_set = 50,
        error_tol = 1e-2
      )
    ))

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
    mmsi,
    ROUND(lat*4)/4 rounded_lat,
    ROUND(lon*4)/4 rounded_lon,
    SUM(IF( nnet_score == 1, hours, 0)) AS total_hours,
    AVG(distance_from_shore ) AS mean_distance_from_shore,
    AVG(distance_from_port) AS mean_distance_from_port,
    INTEGER(REGEXP_REPLACE( IF(REGEXP_EXTRACT(regions,'\"(mparu:.*?)\"') CONTAINS \".\", LEFT(REGEXP_EXTRACT(regions,'\"(mparu:.*?)\"'),INSTR(REGEXP_EXTRACT(regions,'\"(mparu:.*?)\"'),\".\")-1),REGEXP_EXTRACT(regions,'\"(mparu:.*?)\"')), '[^0-9 ]','')) mparu,
    INTEGER(REGEXP_REPLACE( IF(REGEXP_EXTRACT(regions,'\"(mpant:.*?)\"') CONTAINS \".\", LEFT(REGEXP_EXTRACT(regions,'\"(mpant:.*?)\"'),INSTR(REGEXP_EXTRACT(regions,'\"(mpant:.*?)\"'),\".\")-1),REGEXP_EXTRACT(regions,'\"(mpant:.*?)\"')), '[^0-9 ]','')) mpant
    FROM
    [world-fishing-827:gfw_research.nn]
    WHERE
    (distance_from_shore > ({min_dist_from_shore})
    OR (implied_speed > ({min_speed}) AND implied_speed < ({max_speed})))
    AND lat >= ({bbox[1]})
    AND lat <= ({bbox[2]})
    AND lon > ({bbox[3]})
    AND lon <= ({bbox[4]})
    AND _PARTITIONTIME BETWEEN TIMESTAMP('2010-01-01')
    AND TIMESTAMP('2016-12-31')
    AND seg_id IN (
    SELECT
    seg_id
    FROM
    [world-fishing-827:gfw_research.good_segments])
    GROUP BY
    mmsi,
    year,
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

    temp_gfw_data[[i]] <- bigrquery::query_exec(gfw_query, project = "ucsb-gfw", max_pages = Inf)
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
    vessel_length = map2_dbl(
      known_length,
      inferred_length,
      merge_known_and_inferred,
      cap = 150
    ),
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
      labs(caption = "SST(c)")
  }

  chl_data <- chl_data %>%
    mutate(plots = map_plot(chl_a, chl_plot_foo))

  a <- chl_data$chl_a[[4]]

  qmplot(rlon, rlat, color = mean_chlorophyll, data = a)

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
      labs(caption = "SST(c)")
  }

  sst_data <- sst_data %>%
    mutate(plots = map_plot(sst, sst_plot_foo))

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
        runit = lat_lon_res
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


# query mpas --------------------------------------------------------------

if (query_mpas == T) {
  mpa_ids <-
    data_frame(mpa_type = c("no_take_mpas_id", "restricted_use_mpas_id"))

  fishing_connection <- DBI::dbConnect(
    bigrquery::dbi_driver(),
    project = gfw_project,
    dataset = "regions_ids"
  )

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



# perform filters and data processing ---------------------------------------------------------

# decide species that go in
# for now, query FAO for species fished... this is hacky though since means that depends on static data...

if (query_synonyms == T) {
  fao <-
    data.table::fread(
      "~/Box Sync/Databases/FAO/FAO_Capture_1950to2012.csv",
      header = T
    ) %>%
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
      (inferred_label_allyears == "trawlers" |
        inferred_sublabel_allyears == "pots_and_traps" |
        str_detect(inferred_sublabel_allyears, "set_") |
        inferred_sublabel_allyears == "trawlers") &
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
  dplyr::rename(
    Lon = Long,
    Catch_KG = Wt,
    spp = Sci
  ) %>%
  mutate(AreaSwept_km2 = 1, Vessel = "missing") %>%
  select(
    survey,
    survey_region,
    knots,
    Year,
    Lat,
    Lon,
    Vessel,
    AreaSwept_km2,
    Catch_KG,
    spp
  ) %>%
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
  nest(-survey, -survey_region, -knots) # note that this breaks if it's a tibble, should make a bit more flexible


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
        safely(vasterize_index),
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

  vast_fish$vasterized_data <- map(vast_fish$vasterized_data, "result")

  save(
    file = paste0(run_dir, "vast_fish.Rdata"),
    vast_fish
  )
} else {
  load(file = paste0(run_dir, "vast_fish.Rdata"))
}

save(file = here::here("results", run_name, "gfw_data.Rdata"), gfw_data)

# build database ----------------------------------------------------------

# pre-process gfw data

skynet_data <- gfw_data %>%
  unnest() %>%
  group_by(survey, year, rounded_lat, rounded_lon) %>%
  summarise(
    num_vessels = length(unique(mmsi)),
    total_hours = sum(total_hours),
    total_engine_power = sum(inferred_engine_power),
    total_engine_hours = sum(total_hours * inferred_engine_power),
    dist_from_port = mean(mean_distance_from_port),
    dist_from_shore = mean(mean_distance_from_shore),
    mean_vessel_length = mean(inferred_length),
    no_take_mpa = any(is.na(mpant) == F),
    restricted_use_mpa = any(is.na(mparu) == F)
  ) %>%
  ungroup() %>%
  mutate(
    vessel_hours = total_engine_hours * num_vessels,
    any_fishing = total_engine_hours > 0
  ) %>%
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

knots <- vast_fish %>%
  select(survey, vasterized_data) %>%
  mutate(knots = map(vasterized_data, c("spatial_list", "loc_x_lat_long"))) %>%
  select(-vasterized_data)


if (fished_only == T) {
  # include only fished species

  total_fish_data <- vast_fish %>%
    mutate(spatial_densities = map(vasterized_data, "spatial_densities")) %>%
    select(survey, spatial_densities) %>%
    unnest() %>%
    filter(str_replace(species, "_", " ") %in% fished_species$sci_name)
} else {
  total_fish_data <- vast_fish %>%
    mutate(spatial_densities = map(vasterized_data, "spatial_densities")) %>%
    select(survey, spatial_densities) %>%
    unnest()
}


if (raw_fish == T) {
  # option to sub in raw survey densities instead of VAST outputs

  knot_foo <- function(dat) {
    dat <- dat %>%
      mutate(knot = paste(approx_long, approx_lat, sep = "-") %>% as.factor() %>% as.numeric())
  }


  total_fish_data <- fish_data %>%
    {
      if (fished_only == T) {
        filter(
          .,
          str_replace(.$Sci, "_", " ") %in% fished_species$sci_name
        )
      } else {
        .
      }
    } %>%
    select(survey, Sci, Year, Wt, Long, Lat) %>%
    rename(
      species = Sci,
      density = Wt,
      approx_long = Long,
      approx_lat = Lat
    ) %>%
    set_names(tolower) %>%
    mutate(
      approx_long = round(approx_long * (1 / lat_lon_res)) / (1 / lat_lon_res),
      approx_lat = round(approx_lat * (1 / lat_lon_res)) / (1 / lat_lon_res)
    ) %>%
    nest(-survey) %>%
    mutate(data = map(data, knot_foo)) %>%
    unnest() %>%
    select(survey, knot, species, year, density, approx_long, approx_lat) %>%
    group_by(year, survey, knot, species, approx_long, approx_lat) %>%
    summarise(density = sum(density)) %>%
    ungroup()


  knots <- total_fish_data %>%
    select(survey, approx_long, approx_lat, knot) %>%
    unique() %>%
    nest(-survey, .key = "knots")
}


mean_survey_prices <- total_fish_data %>%
  left_join(species_prices, by = "species") %>%
  group_by(survey) %>%
  summarise(aggregate_price = sum(unique(mean_exvessel_price), na.rm = T))

total_fish_data <- total_fish_data %>%
  group_by(survey, knot, year) %>%
  summarise(density = sum(density)) %>%
  ungroup() %>%
  nest(-survey, .key = fish_data)

if (survey_years_only == T) {
  survey_years <- fish_data %>%
    group_by(survey, Year) %>%
    summarise(year_surveyed = T) %>%
    set_names(tolower)

  total_fish_data <- total_fish_data %>%
    unnest() %>%
    left_join(survey_years, by = c("survey", "year")) %>%
    filter(year_surveyed == T) %>%
    nest(-survey, .key = fish_data)
}

# aggregate data -----------------------------------------------



skynet_data <- skynet_data %>%
  nest(-survey, .key = gfw_data) %>%
  left_join(total_fish_data, by = "survey") %>%
  left_join(knots, by = "survey")


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
  unnest() %>%
  filter(is.na(density) == F)

missing_some <- map_lgl(skynet_data, ~any(is.na(.x)))

impute_locations <- which(str_detect(colnames(skynet_data),'lag') == F &  map_lgl(skynet_data, is.numeric)
 == T & !colnames(skynet_data) %in% c("year", "rounded_lat", "rounded_lon") & missing_some)

skynet_data2 <- skynet_data %>%
  dmap_at(impute_locations, fill_enviros, data = skynet_data)

missing_too_many <- skynet_data %>%
  map_df(~ mean(is.na(.x))) %>%
  gather(variable, percent_missing) %>%
  arrange(desc(percent_missing)) %>%
  filter(
    percent_missing > max_percent_missing,
    !str_detect(variable, "lag"), !str_detect(variable, "density")
  ) %>%
  {
    .$variable
  }

if (length(missing_too_many) > 0){
skynet_data <- skynet_data %>%
  select_(paste0("-", missing_too_many)) %>%
  left_join(mean_survey_prices, by = "survey")
}

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
      dat = skynet_data,
      variable = .x,
      run_dir = run_dir
    )
  )



  temp_maps <- skynet_data %>%
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
        mutate(geometry = sf::st_sfc(
          geometry, crs =
            "+proj=longlat +datum=WGS84 +no_defs"
        )) %>%
        sf::st_sf() %>%
        select(year, variable, value, geometry)


      if (query_coastal_map == T) {
        # coast_map <-
        #   rnaturalearth::ne_download(scale = 'medium',
        #                              category = 'physical',
        #                              type = 'coastline') %>% sf::st_as_sf()


        global_map <-
          rnaturalearth::ne_download(
            scale = "medium",
            category = "physical",
            type = "land"
          ) %>%
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
        coord_sf(
          xlim = c(bbox["xmin"], bbox["xmax"]),
          ylim = c(bbox["ymin"], bbox["ymax"])
        ) +
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

skynet_data <- skynet_data %>%
  select(-one_of(vars_to_drop))

skynet_names <- colnames(skynet_data)

do_not_scale <-
  c(
    "year",
    "rounded_lat",
    "rounded_lon",
    "no_take_mpa",
    "restricted_use_mpa",
    "knot",
    "distance",
    "any_fishing",
    'log_density',
    'cs_log_density',
    'cs_density',
    'density',
    "survey",
    "aggregate_price"
  ) # variables that should not be centered and scaled

skynet_data <- skynet_data %>%
  ungroup() %>%
  mutate(cs_log_density = (log_density - mean(log_density)) / sd(log_density),
         cs_density = (density - mean(density)) / sd(density))

uc_skynet_data <- skynet_data

if (center_and_scale == T) {
  do_scale <- skynet_names[!skynet_names %in% do_not_scale]

  skynet_data <- skynet_data %>%
    group_by(survey) %>%
    purrrlyr::dmap_at(do_scale, ~ (.x - mean(.x, na.rm = T)) / (sd(.x, na.rm = T))) %>% # center and scale variables that should be
    ungroup()
}


save(
  file = here::here("results", run_name, "skynet_data.Rdata"),
  skynet_data, uc_skynet_data
)


dep_var <- c("cs_log_density")

never_ind_vars <-
  c(
    "log_density",
    "density",
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
    "total_engine_hours"
  )


gfw_vars <-
  skynet_names[str_detect(
    skynet_names,
    "(vessel)|(distance)|(length)|(hours)|(engine)|(port)|(shore)|(mpa)"
  )]

tree_candidate_vars <- skynet_names[!skynet_names %in% c(
  dep_var,
  never_ind_vars
) &
  str_detect(skynet_names, "_interaction")
  == F]

if (gfw_vars_only == T) {
  tree_candidate_vars <-
    c(tree_candidate_vars[tree_candidate_vars %in% gfw_vars], "random_var")
}



lm_candidate_vars <- skynet_names[!skynet_names %in% c(
  dep_var,
  never_ind_vars
)]
# models <- c('structural')
models <- c("ranger", "gbm", "structural")

delta_skynet <- skynet_data %>%
  select(-contains("lag")) %>%
  gather(
    variable,
    value, -survey,
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

delta_candidate_vars <- c(
  tree_candidate_vars[tree_candidate_vars %in% delta_vars],
  "random_var"
)

delta_skynet <- delta_skynet %>%
  group_by(survey, variable, rounded_lat, rounded_lon) %>%
  arrange(year) %>%
  mutate(lag_value = lag(value, 1)) %>%
  mutate(delta_value = value - lag_value) %>%
  select(-value, -lag_value) %>%
  spread(variable, delta_value) %>%
  na.omit() %>%
  ungroup()

data_sources <- tibble(
  skynet = list(skynet_data %>%
    select(-contains("lag")) %>%
    na.omit()),

  lag_1_skynet_data = list(
    skynet_data %>%
      select(-total_engine_hours_lag2, -total_engine_hours_lag3) %>%
      na.omit()
  ),
  delta_skynet = list(delta_skynet),
  uncentered_skynet = list(uc_skynet_data %>%
                             select(-contains("lag")) %>%
na.omit())
) %>%
  gather(data_subset, data)



test_train_data <- purrr::cross_df(list(
  test_sets = c(
    "random",
    "alaska",
    "west_coast",
    "historic",
    "goa-ai",
    "ebs-ai"
  ),
  data_subset = c("skynet", "delta_skynet","uncentered_skynet")
)) %>%
  left_join(data_sources, by = "data_subset")

test_train_data <- test_train_data %>%
  mutate(resampled = map2(data, test_sets, ~ generate_test_training(dat = .x, test_set = .y))) %>%
  select(-data) %>%
  unnest()


# test_train_data <- modelr::crossv_kfold(lag_0_skynet_data, k = 1)

dep_vars <- c('log_density','cs_density','density')

skynet_models <- purrr::cross_df(list(
  dep_var = (dep_vars),
  temp = list(test_train_data),
  model = models
)) %>%
  unnest(temp, .drop = F) %>%
  filter(!(data_subset == 'delta_skynet' & model == 'structural')) %>%
  filter(!(model == 'structural' & data_subset != 'uncentered_skynet')) %>%
  filter( !(model %in% c('ranger','gbm') & data_subset == 'uncentered_skynet')) %>%
  filter(!(model == 'structural' & dep_var != dep_vars[1]))

#
# check <- skynet_models %>%
#   filter(model == 'structural')

# run models --------------------------------------------------------------


if (run_models == T) {
  sfm <- safely(fit_skynet)

  skynet_models <- skynet_models %>%
    # filter(model == "ranger") %>%
    # slice(1) %>%
    # filter(train_set == 'random') %>%
    # slice(1:4) %>%
    # group_by(model, data_subset, dep_var) %>%
    # mutate(i = 1:length(train)) %>%
    # filter(i <=1) %>%
    # ungroup() %>%
    mutate(candidate_vars = ifelse(
      str_detect(.$data_subset, "delta"),
      list(delta_candidate_vars),
      list(tree_candidate_vars)
    )) %>%
    mutate(
      fitted_model = pmap(
        list(
          dep_var = dep_var,
          model_name = model,
          training = train,
          testing = test,
          tree_candidate_vars = candidate_vars
        ),
        sfm,
        fitcontrol_number = 2,
        fitcontrol_repeats = 1,
        never_ind_vars = never_ind_vars,
        tune_model = T
      )
    )

  print('ran models')
  save(
    file = paste0(run_dir, "skynet_models.Rdata"),
    skynet_models
  )
  print('saved models')

} else {
  load(file = paste0(run_dir, "skynet_models.Rdata"))
}

# diagnose models ---------------------------------------------------------

skynet_models <- skynet_models %>%
  mutate(error = map(fitted_model, "error")) %>%
  mutate(no_error = map_lgl(error, is.null)) %>%
  filter(no_error)

mse_foo <- function(dat, dep_var, pred_var) {
  obs <- dat[, dep_var]

  pred <- dat[, pred_var]

  mse <- mean((obs - pred) ^ 2)
}

diagnostic_plot_foo <- function(data,
                                r2,
                                dep_var,
                                test_region,
                                train_region,
                                data_set) {
  data %>%
    ggplot(aes_(as.name(dep_var), ~ pred)) +
    geom_abline(
      aes(intercept = 0, slope = 1),
      color = "red",
      linetype = 2
    ) +
    geom_point(alpha = 0.75) +
    labs(
      title = paste0("tested on ", test_region, "- trained on ", train_region, "; R2 = ", r2),
      subtitle = paste0("data subset is: ", data_set)
    )
}


skynet_models <- skynet_models %>%
  mutate(
    test_data = map(fitted_model, c("result", "test_predictions")),
    training_data = map(fitted_model, c("result", "training_predictions")),
    results = map(fitted_model, c("result", "model"))
  ) %>%
  mutate(
    mse = map2_dbl(test_data, dep_var, mse_foo, pred_var = "pred"),
    var_y = map2_dbl(test_data, dep_var, ~ var(.x[, .y])),
    psuedo_r2 = round(1 - mse / var_y, 2),
    mse_training = map2_dbl(training_data, dep_var, mse_foo, pred_var = "pred"),
    var_y_training = map2_dbl(training_data, dep_var, ~ var(.x[, .y])),
    psuedo_r2_training = round(1 - mse_training / var_y_training, 2)
  ) %>%
  mutate(
    test_plot = pmap(
      list(
        data = test_data,
        r2 = psuedo_r2,
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
        r2 = psuedo_r2_training,
        test_region = paste0(test_sets),
        train_region = paste0(train_set, "- training plot"),
        data_set = data_subset,
        dep_var = dep_var
      ),
      diagnostic_plot_foo
    )
  )

skynet_models <- skynet_models %>%
  mutate(resolution_test = map2(test_data, dep_var, resolution_effect)) %>%
  mutate(test_resolution_plot = map(resolution_test, plot_resolution_effect))


print('processed models')

save_foo <- function(test_plot,
                     model,
                     train_region,
                     test_region,
                     data_set,
                     dep_var,
                     run_dir) {
  ggsave(
    filename = paste0(run_dir, model, "-train_", train_region, "-test_", test_region, "-data_", data_set,'-depvar_',dep_var, ".pdf"),
    test_plot,
    height = 8,
    width = 8
  )
}

pwalk(
  list(
    model = skynet_models$model,
    train_region = skynet_models$train_set,
    test_plot = skynet_models$test_plot,
    test_region = skynet_models$test_sets,
    data_set = skynet_models$data_subset,
    dep_var = skynet_models$dep_var
  ),
  save_foo,
  run_dir = run_dir
)

pwalk(
  list(
    model = skynet_models$model,
    test_plot = skynet_models$training_plot,
    train_region = skynet_models$train_set,
    test_region = paste0(skynet_models$test_sets, "-training plot"),
    data_set = skynet_models$data_subset,
    dep_var = skynet_models$dep_var
  ),
  save_foo,
  run_dir = run_dir
)

pwalk(
  list(
    model = skynet_models$model,
    train_region = skynet_models$train_set,
    test_plot = skynet_models$test_resolution_plot,
    test_region = paste0(skynet_models$test_sets, "-resolution_plot"),
    data_set = skynet_models$data_subset,
    dep_var = skynet_models$dep_var
  ),
  save_foo,
  run_dir = run_dir
)


print('printed models')

save(file = here::here("results", run_name, "processed_skynet_models.Rdata"), skynet_models)

print('saved processed models')
