create_skynet_data <-
  function(gfw_data,
           fish_data,
           fish_knots,
           surveys = c('AIBTS', 'EBSBTS', 'GOABTS', 'WCGBTS', 'WCGHL'),
           ...) {
    # find nearest knot ---------------------------------------------------

    gfw_latlon_coords <- gfw_data %>%
      select(rounded_lon, rounded_lat)
    nearest_knot <-
      RANN::nn2(fish_knots %>% select(-knot), gfw_latlon_coords, k = 1)

    gfw_data$knot <- fish_knots$knot[nearest_knot$nn.idx]

    gfw_data$distance <- nearest_knot$nn.dists

    knot_pairing_plot <- gfw_data %>%
      filter(distance <= quantile(gfw_data$distance, 0.75)) %>%
      select(rounded_lat, rounded_lon, knot) %>%
      unique() %>%
      ggplot() +
      geom_text(
        aes(
          rounded_lon,
          rounded_lat,
          color = factor(knot),
          label = knot
        ),
        show.legend = F,
        alpha = 1,
        size = 1.5
      )


    # merge fishdata ----------------------------------------------------------

    skynet_data <- gfw_data %>%
      left_join(fish_data %>% select(density, year, knot), by = c('year', 'knot')) %>%
      rename(log_density = density)

    # merge additional sata

    covars <- list(...)

    flat_covars <- map_df(covars,
                          process_covars,
                          surveys = surveys,
                          years = unique(skynet_data$year)) %>%
      spread(variable, mean_value)


    skynet_data <- skynet_data %>%
      left_join(flat_covars,
                by = c('rounded_lat' = 'rlat', 'rounded_lon' = 'rlon', 'year'))

    if (nrow(skynet_data) != nrow(gfw_data)){

      stop('joins have gone off the rails')

    }
    return(skynet_data)

  }