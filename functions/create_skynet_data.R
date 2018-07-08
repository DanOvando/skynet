#' Create Skynet Data
#' \code{create_skynet_data} creates a dataframe stictching together gfw, fishdata, and environmental data
#' @param gfw_data gfw data
#' @param fish_data fishdata data
#' @param fish_knots fishdata knot locations
#' @param surveys the survey being used
#' @param cap_distance logical indicating whether to set maximum distance from knot included
#' @param distance_quantile max distance quantile included in knot distance
#' @param ... space for arbitrary environmental parameters
#'
#' @return a stitched dataframe
#' @export
#'
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

    # merge fishdata ----------------------------------------------------------
    skynet_data <- gfw_data %>%
      left_join(
        fish_data %>% select(
          density,
          lag_density,
          biomass,
          lag_economic_density,
          economic_density,
          year,
          surveyed_year,
          knot,
          mean_knot_area
        ),
        by = c('year', 'knot')
      ) %>%
      mutate(
        log_lag_density = log(lag_density),
        log_lag_economic_density = log(lag_economic_density),
        log_density = log(density),
        log_biomass = log(biomass),
        log_economic_density = log(economic_density)
      )# merge in fishdata

    # merge additional data

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
        filter(distance <= quantile(skynet_data$distance, distance_quantile))

    }

    return(skynet_data)

  }