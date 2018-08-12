snap_to_grid <-
  function(old_grid,
           new_grid,
           old_lon_name,
           old_lat_name,
           new_lon_name,
           new_lat_name,
           dep_var = "log_density") {

    old_lon_name <- enquo(old_lon_name)

    old_lat_name <- enquo(old_lat_name)

    new_lon_name <- enquo(new_lon_name)

    new_lat_name <- enquo(new_lat_name)

    # old_grid <- test_resolution$data[[1]]
    #
    # new_grid <- test_resolution$new_grid[[1]]

    old_coords <- old_grid %>%
      select(!!old_lon_name, !!old_lat_name)

    new_coords <- new_grid %>%
      select(!!new_lon_name, !!new_lat_name)

    nearest_knot <-
      RANN::nn2(new_coords, old_coords, k = 1)

    if (!"id" %in% colnames(new_grid)){

      new_grid$id <- 1:nrow(new_grid)

    }

    old_grid$new_knot <- new_grid$id[nearest_knot$nn.idx]

    old_grid <- old_grid %>%
      left_join(new_grid, by = c("new_knot" = "id"))

    if (!"mean_knot_area" %in% colnames(old_grid)){

      old_grid$mean_knot_area <- 1

    }

    if (!"pred" %in% colnames(old_grid)){

      old_grid$pred <- NA

    }

    new_grid <- old_grid %>%
      mutate(survey_knot = paste(survey, knot, sep = "-")) %>%
      group_by(year,survey_knot) %>%
      mutate(spread_over = n_distinct(new_knot)) %>%
      mutate(dilute_total = 1/spread_over) %>%
      group_by(year, survey, knot, new_knot, lat, lon) %>%
      summarise(b = unique(biomass),
                mp = unique(aggregate_price),
                dt = unique(dilute_total),
                ka = unique(mean_knot_area),
                d = unique(density),
                th = sum(total_hours),
                teh = sum(total_engine_hours),
                tn = sum(num_vessels),
                mpred = mean(pred)
                ) %>%
      group_by(year, survey, new_knot, lat, lon) %>%
      summarise(total_biomass = sum(b * dt),
                total_revenue = sum(b * mp * dt),
                agg_mean_density = weighted.mean(d,ka),
                agg_mean_economic_density = weighted.mean(d * mp,ka),
                approx_survey = unique(survey)[1],
                old_obs = length(b),
                agg_total_hours = sum(th),
                agg_total_engine_hours = sum(teh),
                agg_total_num_vessels = sum(tn),
                agg_pred = ifelse(str_detect(dep_var, "density"), weighted.mean(mpred, w = ka, na.rm = T), sum(mpred * dt, na.rm = T)))


      # old_grid %>%
      #   group_by(survey,knot, rounded_lat, rounded_lon) %>%
      #   summarise(tb = sum(biomass)) %>%
      #   ggplot(aes(rounded_lon, rounded_lat, fill = tb)) +
      #   geom_point(shape = 21) +
      #   scale_fill_viridis()

      # new_fish %>%
      #   group_by(survey,new_knot, lat, lon) %>%
      #   summarise(tb = sum(total_biomass)) %>%
      #   ggplot(aes(lon, lat, fill = tb)) +
      #   geom_point(shape = 21) +
      #   scale_fill_viridis()

  }