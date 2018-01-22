fit_models <- function(training,
                       test,
                       mtry = 7,
                       best_lm_model,
                       mle_coefs) {
  training <- training %>%
    as_data_frame()

  test <- test %>%
    as_data_frame()

  rf_model <- randomForest(
    density ~ total_engine_hours +
      total_hours +
      large_vessels +
      port_hours +
      port_engine_hours +
      port_numbers +
      shore_hours +
      shore_engine_hours +
      shore_numbers +
      num_vessels +
      mean_vessel_length,
    data = training,
    mtry = 7
  )


  lm_model <-
    lm(best_lm_model ,
       data = training)

  structural_model <-  stats4::mle(
    fit_sffs,
    start = list(
      beta_distance =  mle_coefs['beta_distance'] %>% as.numeric() %>% jitter(2),
      beta_size =  mle_coefs['beta_size'] %>% as.numeric() %>% jitter(2),
      beta_interaction = mle_coefs['beta_interaction'] %>% as.numeric(2) %>% jitter(2),
      sigma =  mle_coefs['sigma'] %>% as.numeric() %>% jitter(2),
      q =  mle_coefs['q'] %>% as.numeric() %>% jitter(2)
    ),
    fixed = list(
      price = 765,
      dat = training,
      marginal_profits = 0,
      use = 1
    ),
    upper = list(sigma = log(100),
                 q = log(1.5e-5))
  )


  performance <- performance %>%
    mutate(
      rf = predict(rf_model, newdata = .),
      lm = predict(lm_model, newdata = .),
      structural = fit_sffs(
        beta_distance = mle_coefs['beta_distance'],
        beta_size = mle_coefs['beta_size'],
        beta_interaction = mle_coefs['beta_interaction'],
        q = mle_coefs['q'],
        sigma = mle_coefs['sigma'],
        dat = .,
        price = 765,
        marginal_profits = 10,
        use = 0
      )
    ) %>%
    gather('model', 'density_hat', rf:structural) %>%
    mutate(se = (density - density_hat) ^ 2) %>%
    group_by(model) %>%
    dplyr::mutate(rmse = sqrt(mean(se))) %>%
    mutate(model = fct_reorder(model, rmse))

  return(performance)

}