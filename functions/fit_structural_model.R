fit_structural_model <- function(beta_distance,
                                 beta_size,
                                 beta_depth,
                                 beta_mpant,
                                 beta_mparu,
                                 beta_engine_power,
                                 q,
                                 sigma,
                                 dat,
                                 price,
                                 marginal_profits = 0,
                                 use = 'mle',
                                 log_space = 1) {

  old_sigma <- sigma
  sigma <-  exp(sigma)

  old_q <- q

  q <- exp(q)

  # beta_size <- exp(beta_size)

  # beta_distance <-  exp(beta_distance)

  # beta_size <- exp(beta_size)

  # beta_interaction <- exp(beta_interaction)

  # beta_intercept <- exp(beta_intercept)

  cost <-
    beta_distance * dat$dist_from_port +
    beta_size * dat$mean_vessel_length +
    beta_engine_power * dat$total_engine_power +
    beta_mpant * dat$no_take_mpa +
    beta_mparu * dat$restricted_use_mpa +
    beta_depth * dat$m_below_sea_level


  cost <- exp(cost)

  # estimated_abundance <-
  #   ((cost + marginal_profits)) / (price  * exp(-q * dat$total_hours))

  estimated_abundance <-
    log(cost + marginal_profits) -  log(price  * exp(-q * dat$total_hours))

  # browser()
  if (any(is.finite(estimated_abundance) == F)) {
    browser()
  }


  if (use == 1) {
    out <-
      -sum(stats::dnorm(dat$log_density, estimated_abundance, sigma, log = T))
  } else {
    out <- estimated_abundance
  }
  # print(out)
  return(out)

}