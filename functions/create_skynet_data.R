create_skynet_data <-
  function(gfw_data,
           fish_data,
           fish_knots,
           surveys = c('AIBTS', 'EBSBTS', 'GOABTS', 'WCGBTS', 'WCGHL'),
           cap_distance = T,
           distance_quantile = 0.75,
           ...) {
    # find nearest knot ---------------------------------------------------
    gfw_latlon_coords <- gfw_data %>%
      select(rounded_lon, rounded_lat)
    nearest_knot <-
      RANN::nn2(fish_knots %>% select(-knot), gfw_latlon_coords, k = 1)

    gfw_data$knot <- fish_knots$knot[nearest_knot$nn.idx]

    gfw_data$distance <- nearest_knot$nn.dists
    # browser()
    # knot_pairing_plot <- gfw_data %>%
    #   filter(distance <= quantile(gfw_data$distance, 0.75)) %>%
    #   select(rounded_lat, rounded_lon, knot) %>%
    #   unique() %>%
    #   ggplot() +
    #   geom_text(
    #     aes(
    #       rounded_lon,
    #       rounded_lat,
    #       color = factor(knot),
    #       label = knot
    #     ),
    #     show.legend = F,
    #     alpha = 1,
    #     size = 1.5
    #   )


    # merge fishdata ----------------------------------------------------------

    skynet_data <- gfw_data %>%
      left_join(fish_data %>% select(density, year, knot), by = c('year', 'knot')) %>%
      rename(log_density = density) %>%
      mutate(density = exp(log_density))# merge in fishdata

    # merge additional sata

    covars <-
      list(...) #capture all additional data passed to the function


    flat_covars <- map_df(covars,
                          process_covars,
                          surveys = surveys,
                          years = unique(skynet_data$year)) %>%
      spread(variable, mean_value) # collect and spread covariates


    skynet_data <- skynet_data %>%
      ungroup() %>%
      left_join(flat_covars,
                by = c('rounded_lat' = 'rlat', 'rounded_lon' = 'rlon', 'year')) %>%
      mutate(random_var = rnorm(n = nrow(.), mean = 0, sd = 1))#join covariates to data

    if (nrow(skynet_data) != nrow(gfw_data)) {
      stop('joins have gone off the rails')

    }

    skynet_data <-
      purrr::set_names(skynet_data, tolower(colnames(skynet_data)))

    if (cap_distance == T) {
      skynet_data <- skynet_data %>%
        filter(distance <= quantile(skynet_data$distance, 0.75))

    }

    return(skynet_data)

  }