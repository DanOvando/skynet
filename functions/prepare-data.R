prepare_data <-
  function(vast_fish,
           fished_only,
           unfished_only,
           survey_months_only,
           raw_fish,
           no_wcghl = T,
           species_prices,
           vars_to_drop = vars_to_drop) {

    if (survey_months_only == T){
    survey_months <- tribble(~survey,~survey_months,
                             "wcgbts", c(5:10),
                             "ebsbts",c(7,8),
                             "wcgbts",c(7,8),
                             "goabts", c(7,8),
                             "wcghl", c(7,8),
                             "aibts", c(7,8)
    )


    gfw_data <- gfw_data %>%
      left_join(survey_months, by = "survey") %>%
      mutate(data = map2(data,survey_months, ~filter(.x, month %in% .y))) %>%
      select(-survey_months)
    }

    skynet_data <- gfw_data %>%
      unnest() %>%
      group_by(survey, year, rounded_lat, rounded_lon) %>%
      rename(te = total_hours) %>%
      summarise(
        num_vessels = length(unique(mmsi)),
        total_hours = sum(te),
        total_engine_power = sum(inferred_engine_power),
        total_engine_hours = sum(te * inferred_engine_power),
        positive_hours = sum(te > 0),
        median_hours = median(te),
        median_engine_power = median(inferred_engine_power),
        median_engine_hours = median(te * inferred_engine_power),
        dist_from_port = mean(mean_distance_from_port),
        dist_from_shore = mean(mean_distance_from_shore),
        mean_vessel_length = mean(inferred_length),
        no_take_mpa = any(is.na(mpant) == F),
        restricted_use_mpa = any(is.na(mparu) == F)
      ) %>%
      ungroup() %>%
      mutate(vessel_hours = total_engine_hours * num_vessels,
             any_fishing = total_engine_hours > 0) %>%
      group_by(rounded_lat, rounded_lon) %>%
      arrange(year) %>%
      mutate(
        cumulative_hours = cumsum(total_engine_hours),
        cumulative_numbers = cumsum(total_hours),
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


    if (fished_only == T) {
      # include only fished species

      total_fish_data <- vast_fish %>%
        mutate(spatial_densities = map(vasterized_data, "spatial_densities")) %>%
        select(survey, spatial_densities) %>%
        unnest() %>%
        filter(
          str_replace(species, "_", " ") %in% fished_species$sci_name |
            str_detect(species, "Sebastes")
        )
    } else if (unfished_only == T) {
      total_fish_data <- vast_fish %>%
        mutate(spatial_densities = map(vasterized_data, "spatial_densities")) %>%
        select(survey, spatial_densities) %>%
        unnest() %>%
        filter(!(
          str_replace(species, "_", " ") %in% fished_species$sci_name  |
            str_detect(species, "Sebastes")
        ))
    } else{
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
          filter(.,
                 str_replace(.$Sci, "_", " ") %in% fished_species$sci_name)
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
        select(survey,
               knot,
               species,
               year,
               density,
               approx_long,
               approx_lat) %>%
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
      left_join(species_prices, by = "species")

    if (all(is.na(total_fish_data$mean_exvessel_price)) |
        is.null(total_fish_data$mean_exvessel_price)) {
      total_fish_data$mean_exvessel_price <- 1

    }

    total_fish_data$mean_exvessel_price[is.na(total_fish_data$mean_exvessel_price)] <-
      mean(total_fish_data$mean_exvessel_price, na.rm = T)

    total_fish_data <- total_fish_data %>%
      group_by(survey, knot, year) %>%
      summarise(
        density = sum(density),
        mean_knot_area = mean(area),
        biomass = sum(biomass),
        economic_density = sum(density * (mean_exvessel_price * .001)),
        economic_biomass = sum(biomass * (mean_exvessel_price * .001))
      ) %>%
      ungroup() %>%
      nest(-survey, .key = fish_data)


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

    missing_some <- map_lgl(skynet_data, ~ any(is.na(.x)))

    if (impute_missing == T) {
      impute_locations <-
        which(
          str_detect(colnames(skynet_data), 'lag') == F &
            map_lgl(skynet_data, is.numeric)
          == T &
            !colnames(skynet_data) %in% c("year", "rounded_lat", "rounded_lon") &
            missing_some
        )

      skynet_data <- skynet_data %>%
        dmap_at(impute_locations, fill_enviros, data = skynet_data)
    }
    missing_too_many <- skynet_data %>%
      map_df(~ mean(is.na(.x))) %>%
      gather(variable, percent_missing) %>%
      arrange(desc(percent_missing)) %>%
      filter(
        percent_missing > max_percent_missing,
        !str_detect(variable, "lag"),!str_detect(variable, "density")
      ) %>%
      {
        .$variable
      }

    if (length(missing_too_many) > 0) {
      skynet_data <- skynet_data %>%
        select_(paste0("-", missing_too_many)) %>%
        left_join(mean_survey_prices, by = "survey")
    } else{
      skynet_data <- skynet_data %>%
        left_join(mean_survey_prices, by = "survey")
    }

    skynet_data <- skynet_data %>%
      select(-one_of(vars_to_drop))


    if (no_wcghl == T){
      skynet_data <-  skynet_data %>%
        filter(survey != "wcghl")

      total_fish_data <-  total_fish_data %>%
        filter(survey != "wcghl")

    }

    skynet_data <- skynet_data %>%
      ungroup() %>%
      group_by(survey) %>%
      mutate(
        cs_log_density = (log_density - mean(log_density)) / sd(log_density),
        cs_density = (density - mean(density)) / sd(density),
        relative_density = density / max(density)
      ) %>%
      ungroup()

    total_fish_data$index <- 1:nrow(total_fish_data)

    skynet_data$index <- 1:nrow(skynet_data)

    return(list(total_fish_data = total_fish_data,
                skynet_data = skynet_data))

  }