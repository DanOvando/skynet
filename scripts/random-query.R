library(bigrquery)
library(lubridate)
library(ggmap)
library(stringr)
library(viridis)
library(purrr)
library(purrrlyr)
library(modelr)
library(sf)
library(gganimate)
library(hrbrthemes)
library(tidyverse)

gfw_project <- "ucsb-gfw"

gfw_dataset <- "skynet"

max_percent_missing <- 0.2

lat_lon_res <-
  0.5 # round data to intervels of 0.25 degrees lat lon, as in 56.25

res <- 1 / lat_lon_res


min_year <- 2012

min_times_seen <- 9

min_dist_from_shore <- 1000

min_speed <- 0.01

max_speed <- 20


survey_bbox <-
  data_frame(
    min_year = 2012,
    max_year = 2016,
    min_lat = -37,
    max_lat = 36,
    min_lon = -20,
    max_lon = 16
  )


fishing_connection <- DBI::dbConnect(
  bigrquery::dbi_driver(),
  project = gfw_project,
  dataset = gfw_dataset,
  allowLargeResults = T
)

temp_gfw_data <- list()

  bbox <- survey_bbox %>%
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

  temp_gfw_data <-bigrquery::query_exec(gfw_query, project = "ucsb-gfw", max_pages = Inf)

  africa_data <- temp_gfw_data %>%
    select(-b_mmsi,-b_year) %>%
    set_names(str_replace_all(colnames(.), "(a_)|(b_)", "")) %>%
    group_by(year, rounded_lat, rounded_lon) %>%
    summarise(total_hours = sum(total_hours))


  global_map <-
    rnaturalearth::ne_download(scale = "medium",
                               category = "physical",
                               type = "land") %>%
    sf::st_as_sf()

  africa_fishing <-  africa_data %>%
    dplyr::mutate(geometry = purrr::map2(rounded_lon, rounded_lat, ~ sf::st_point(x = c(.x, .y), dim = 'XY'))) %>%
    ungroup() %>%
    mutate(geometry = sf::st_sfc(geometry, crs =
                                   "+proj=longlat +datum=WGS84 +no_defs")) %>%
    sf::st_sf()

  years <- expand.grid(id = 1:nrow(global_map), year = unique(africa_fishing$year))

  yearly_map <- global_map %>%
    mutate(id = 1:nrow(.)) %>%
    right_join(years, by = "id")

  bbox <- sf::st_bbox(africa_fishing)

  translate_variable <- c(
    total_abundance = "Fish",
    total_engine_hours = "Fishing"
  )

  africa_map <- africa_fishing %>%
    ggplot(aes(frame = year)) +
    geom_raster(aes(x = rounded_lon, y = rounded_lat, fill = log(total_hours)), show.legend = F) +
    geom_sf(data = yearly_map, fill = 'grey60') +
    scale_fill_viridis() +
    coord_sf(xlim = c(bbox['xmin'], bbox['xmax']),
             ylim = c(bbox['ymin'], bbox['ymax'])) +
    labs(x = "Longitude", y = "Latitude") +
    theme_ipsum(base_size = 18, axis_title_size = 20,subtitle_size = 18, axis_text_size = 18, plot_title_size = 20, strip_text_size = 20) +
    theme(plot.margin = unit(c(.25,.25,.25,.25), "in")) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.text = element_text(size = 40))


  gganimate(africa_map, filename = "africa_map.gif", ani.width = 500, ani.height = 600, title_frame = FALSE)


  naive_map <- ebs_fishing %>%
    gather(variable, value, total_abundance, total_engine_hours) %>%
    mutate(variable = fct_relevel(variable,c("total_engine_hours","total_abundance"))) %>%
    ggplot(aes(frame = year)) +
    geom_raster(aes(x = recenter_lon, y = rounded_lat, fill = value), show.legend = F) +
    geom_sf(data = yearly_map, fill = 'grey60') +
    scale_fill_viridis() +
    coord_sf(xlim = c(bbox['xmin'], bbox['xmax']),
             ylim = c(bbox['ymin'], bbox['ymax'])) +
    facet_grid(.~variable,labeller = labeller(variable = translate_variable), as.table = FALSE) +
    theme_ipsum(base_size = 18, axis_title_size = 20,subtitle_size = 18, axis_text_size = 18, plot_title_size = 20, strip_text_size = 20) +
    theme(plot.margin = unit(c(.25,.25,.25,.25), "in")) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.text = element_text(size = 40))

  gganimate(naive_map, filename = "alaska_map.gif", ani.width = 800, ani.height = 600)

