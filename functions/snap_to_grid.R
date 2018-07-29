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
      group_by(year, new_knot) %>%
      summarise(
        old_obs = length(density),
        lat = unique(!!new_lat_name),
        lon = unique(!!new_lon_name),
        approx_survey = unique(survey)[1],
        agg_mean_density = sum(unique(density * mean_knot_area) / sum(unique(mean_knot_area))) ,
        total_biomass = sum(unique(density * mean_knot_area)),
        total_revenue = sum(unique(economic_density * mean_knot_area)),
        agg_mean_economic_density = weighted.mean(economic_density, w = mean_knot_area),
        agg_total_engine_hours = sum(total_engine_hours),
        agg_total_hours = sum(total_hours),
        agg_total_num_vessels = sum(num_vessels),
        agg_pred = ifelse(str_detect(dep_var, "density"), weighted.mean(pred, w = mean_knot_area, na.rm = T), sum(pred, na.rm = T))
      )

  }