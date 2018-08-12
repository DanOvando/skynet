#' downscale_data takes skynet data rescales it to a new spatial resolution
#'
#' @param data data to be rescaled
#' @param resolution km2 resolution of new data
#'
#' @return a rescaled dataframe
#' @export
#'
#' @examples
rescale_data <- function(data, resolution = 25, lon_name = rounded_lon, lat_name = rounded_lat) {

  orig_names <- colnames(data)

  lon_name <- enquo(lon_name)

  lat_name <- enquo(lat_name)

  old_coords <- data %>%
    select(!!lon_name, !!lat_name)

    new_grid <- create_grid(data = data, resolution = resolution, lon_name = rounded_lon, lat_name = rounded_lat)

    new_coords <- new_grid %>%
      select(lon, lat)

    nearest_knot <-
      RANN::nn2(new_coords, old_coords, k = 1)

    if (!"id" %in% colnames(new_grid)){

      new_grid$id <- 1:nrow(new_grid)

    }

    data$new_knot <- new_grid$id[nearest_knot$nn.idx]

    data <- data %>%
      left_join(new_grid, by = c("new_knot" = "id"))

    if (!"mean_knot_area" %in% colnames(data)){

      data$mean_knot_area <- 1

    }

    # test <- data %>%
    #   group_by(new_knot, year) %>%
    #   mutate(obs = length(density)) %>%
    #   ungroup() %>%
    #   filter(year == 2012) %>%
    #   filter(obs > 2)
    #
    # test %>%
    #   ggplot(aes(rounded_lon, rounded_lat, fill =factor(knot))) +
    #   geom_point(shape = 21)

    # assign new lat and lon

    data <- data %>%
      mutate(survey_knot = paste(survey, knot, sep = "-")) %>%
      group_by(year,survey_knot) %>%
      mutate(spread_over = n_distinct(new_knot)) %>%
      mutate(dilute_total = 1/spread_over) %>%
      ungroup()

    data <- data %>%
      mutate(rounded_lat = lat,
             rounded_lon = lon,
             knot = new_knot) %>%
      select(-lat, -lon, -recenter_lon, -new_knot,-survey_knot,-spread_over)

    long_data <- data %>%
      gather(variable, value,-dilute_total,-mean_knot_area,-survey, -surveyed_year,-year, -rounded_lat, -rounded_lon,-knot)


    variables <- unique(long_data$variable)

    vars_to_average <- c(variables[
      str_detect(variables,"mean") |
      str_detect(variables, "density") |
        str_detect(variables, "mpa") |
      str_detect(variables, "dist")
    ], "m_below_sea_level","random_var","aggregate_price")

    vars_to_sum <- variables[!variables %in% vars_to_average]

    long_data <- long_data %>%
      rename(mka = mean_knot_area)

    rescaled_data <- long_data %>%
      group_by(survey, year,surveyed_year, rounded_lat, rounded_lon, knot, variable) %>%
      summarise(rescaled_value = rescale_foo(value, variable, dilute_total,mka,vars_to_average, vars_to_sum),
                mean_knot_area = mean(mka))

    rescaled_data <- rescaled_data %>%
      ungroup() %>%
      spread(variable, rescaled_value)

    rescaled_data <- rescaled_data[, orig_names]

    rescaled_data$any_fishing <- rescaled_data$total_hours > 0

        # a <- rescaled_data %>%
    #   select(rounded_lat, rounded_lon, total_hours, year) %>%
    #   mutate(source = "rescaled") %>%
    #   bind_rows(skynet_data %>%  select(rounded_lat, rounded_lon, total_hours, year) %>%
    #               mutate(source = "skynet")) %>%
    #   filter(year == 2012, rounded_lat > 50)
    #
    # a %>%
    #   ggplot(aes(rounded_lon, rounded_lat, fill = total_hours)) +
    #   geom_point(shape = 21) +
    #   facet_wrap(~source) +
    #   scale_fill_viridis()

    rescaled_data$index <- 1:nrow(rescaled_data)
    return(rescaled_data)


}