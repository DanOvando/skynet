#' Fit Skynet
#' A wrapper function for fitting models to skynet data
#'
#' @param dep_var the dependent variable being considered
#' @param model the model to be used
#' @param training data to be used in training
#' @param testing data to be used in testing
#' @param fitcontrol_method caret options
#' @param fitcontrol_number caret options
#' @param fitcontrol_repeats caret options
#'
#' @return a list of model fits
#' @export
#'
fit_skynet <- function(dep_var,
                       model,
                       training,
                       testing,
                       never_ind_vars,
                       survey_prices,
                       tree_candidate_vars,
                       lm_candidate_vars,
                       structural_vars = c(
                         'log_density',
                         'dist_from_port',
                         'mean_vessel_length',
                         'm_below_sea_level',
                         'restricted_use_mpa',
                         'no_take_mpa',
                         'total_hours',
                         'total_engine_power',
                         'aggregate_price'
                       ),
                       fitcontrol_method = 'repeatedcv',
                       fitcontrol_number = 1,
                       fitcontrol_repeats = 1,
                       tune_model = T) {
  # model_formula <-
  #   paste(dep_var, '~', ind_vars) %>% as.formula() # construct model forumla
  #
  #

  fit_control <- trainControl(method = fitcontrol_method,
                              number = fitcontrol_number,
                              repeats = fitcontrol_repeats)

  independent_data <- training %>%
    as.data.frame() %>%
    select(-matches(paste0(never_ind_vars, collapse = '|'))) %>%
    select(matches(paste0(tree_candidate_vars, collapse = '|')))

  dependent_data <- training %>%
    as.data.frame() %>%
    as.matrix()

  dependent_data <- dependent_data[, dep_var] %>% as.numeric()

  model_formula <-
    paste0(dep_var, '~', paste(colnames(independent_data), collapse = '+')) %>%
    as.formula()

  reg_data <- independent_data %>%
    as_data_frame() %>%
    mutate(!!dep_var := dependent_data)

  if (model == 'cforest') {
    model <- train(
      model_formula,
      data = reg_data,
      method = "cforest",
      trControl = fit_control
    )
    if (tune_model == T) {
      cforest_importance <- varImp(model$finalModel)

      varimp <- cforest_importance %>%
        mutate(varname = rownames(.)) %>%
        arrange(desc(Overall)) %>%
        mutate(varname = str_replace(varname, '(TRUE)|(FALSE)', ''))


      random_imp <-
        varimp$Overall[varimp$varname == 'random_var']

      vars_to_include <- varimp %>%
        filter(Overall > random_imp) %>% {
          .$varname
        } %>% unique()

      reg_data <- reg_data[, c(dep_var, vars_to_include)]

      model_formula <-
        paste0(dep_var, '~', paste(vars_to_include, collapse = '+')) %>%
        as.formula()

      model <- train(
        model_formula,
        data = reg_data,
        method = "cforest",
        trControl = fit_control
      )

    }

    out_training <- training %>%
      as_data_frame() %>%
      add_predictions(model)

    out_testing <- testing %>%
      as_data_frame() %>%
      add_predictions(model)

  } # close cforest

  if (model == 'rf') {
    #fit random forest

    # rf_grid <-  expand.grid(mtry = c(3,5))

    model <- train(
      independent_data,
      dependent_data,
      method = "rf",
      do.trace = 10,
      # na.action = na.omit,
      trControl = fit_control
      # tuneGrid = rf_grid
    )

    if (tune_model == T) {
      impnames <-
        model$finalModel$importance %>% as.data.frame() %>% purrr::pluck(attr_getter("row.names"))
      varimp <-  model$finalModel$importance %>%
        as_data_frame() %>%
        mutate(varname = impnames) %>%
        arrange(desc(IncNodePurity)) %>%
        mutate(varname = str_replace(varname, '(TRUE)|(FALSE)', ''))

      random_imp <-
        varimp$IncNodePurity[varimp$varname == 'random_var']



      vars_to_include <- varimp %>%
        filter(IncNodePurity > random_imp) %>% {
          .$varname
        } %>% unique()

      independent_data <- independent_data[, vars_to_include]


      model <- train(
        independent_data,
        dependent_data,
        method = "rf",
        do.trace = 10,
        na.action = na.omit,
        trControl = fit_control
      )

    }

    out_testing <- testing %>%
      as_data_frame() %>%
      add_predictions(model)

    out_training <- training %>%
      as_data_frame() %>%
      add_predictions(model)

  } # close rf

  if (model == 'gbm') {

    # gbm_grid <-  expand.grid(interaction.depth = c(1, 5),
    #                         n.trees = (1:10)*100,
    #                         shrinkage = 0.01,
    #                         n.minobsinnode = 20)

    model <- train(
      model_formula,
      data = reg_data,
      method = "gbm",
      verbose = T,
      trControl = fit_control#,
      # tuneGrid = gbm_grid
    )

    on.exit(detach('package:plyr'))


    if (tune_model == T) {
      var_importance <- varImp(model$finalModel)

      varimp <- var_importance %>%
        mutate(varname = rownames(.)) %>%
        arrange(desc(Overall)) %>%
        mutate(varname = str_replace(varname, '(TRUE)|(FALSE)', ''))


      random_imp <-
        varimp$Overall[varimp$varname == 'random_var']

      vars_to_include <- varimp %>%
        filter(Overall > random_imp) %>% {
          .$varname
        } %>% unique()
      reg_data <- reg_data[, c(dep_var, vars_to_include)]

      model_formula <-
        paste0(dep_var, '~', paste(vars_to_include, collapse = '+')) %>%
        as.formula()

      model <- train(
        model_formula,
        data = reg_data,
        method = "gbm",
        verbose = T,
        trControl = fit_control
      )

      on.exit(detach('package:plyr'))

    }

    out_testing <- testing %>%
      as_data_frame() %>%
      add_predictions(model)

    out_training <- training %>%
      as_data_frame() %>%
      add_predictions(model)

  }

  if (model == 'structural') {

    independent_data <- training %>%
      as.data.frame() %>%
      select(matches(paste0(structural_vars, collapse = '|')))


    testing <- testing %>%
      as.data.frame() %>%
      select(matches(paste0(structural_vars, collapse = '|')))

    mle_model <- stats4::mle(
      fit_structural_model,
      start = list(
        beta_distance = (2),
        beta_size = 2 ,
        beta_depth = 2 ,
        beta_mpant = 2,
        beta_mparu = 2 ,
        beta_engine_power = 2,
        sigma = log(sd(independent_data$log_density) / 10),
        q = log(1e-4)
      ),
      fixed = list(
        price = independent_data$aggregate_price,
        dat = independent_data,
        marginal_profits = 0,
        use = 1,
        log_space = 1
      ),
      upper = list(
        sigma = log(100),
        q = log(1.5e-3),
        beta_distance = (20),
        beta_size = 20 ,
        beta_depth = 20 ,
        beta_mpant = 20,
        beta_mparu = 20 ,
        beta_engine_power = 20
      )
    )


    # print(mle_model)

    mle_coefs <- mle_model@coef

    for (i in 1:10) {
      mle_model <- stats4::mle(
        fit_structural_model,
        start = list(
          beta_distance =  mle_coefs['beta_distance'] %>% as.numeric() %>% jitter(2),
          beta_size =  mle_coefs['beta_size'] %>% as.numeric() %>% jitter(2),
          beta_depth = mle_coefs['beta_depth'] %>% as.numeric() %>% jitter(2),
          beta_mpant = mle_coefs['beta_mpant'] %>% as.numeric() %>% jitter(2),
          beta_mparu = mle_coefs['beta_mparu'] %>% as.numeric() %>% jitter(2),
          beta_engine_power = mle_coefs['beta_engine_power'] %>% as.numeric() %>% jitter(2),
          sigma =  mle_coefs['sigma'] %>% as.numeric() %>% jitter(2),
          q =  mle_coefs['q'] %>% as.numeric() %>% jitter(2)
        ),
        fixed = list(
          price = independent_data$aggregate_price,
          dat = independent_data,
          marginal_profits = 0,
          use = 1
        ),
        upper = list(
          sigma = log(100),
          q = log(1.5e-4),
          beta_distance = (20),
          beta_size = 20 ,
          beta_depth = 20 ,
          beta_mpant = 20,
          beta_mparu = 20 ,
          beta_engine_power = 20
        )
      )


      mle_coefs <- mle_model@coef
    } # close for loop


    in_sample_prediction <- fit_structural_model(
      beta_distance = mle_coefs['beta_distance'],
      beta_size = mle_coefs['beta_size'],
      beta_depth = mle_coefs['beta_depth'],
      beta_mpant = mle_coefs['beta_mpant'],
      beta_mparu = mle_coefs['beta_mparu'],
      beta_engine_power = mle_coefs['beta_engine_power'],
      q = mle_coefs['q'],
      sigma = mle_coefs['sigma'],
      dat = independent_data,
      price = independent_data$aggregate_price,
      marginal_profits = 0,
      use = 0
    )

    prediction <- fit_structural_model(
      beta_distance = mle_coefs['beta_distance'],
      beta_size = mle_coefs['beta_size'],
      beta_depth = mle_coefs['beta_depth'],
      beta_mpant = mle_coefs['beta_mpant'],
      beta_mparu = mle_coefs['beta_mparu'],
      beta_engine_power = mle_coefs['beta_engine_power'],
      q = mle_coefs['q'],
      sigma = mle_coefs['sigma'],
      dat = testing,
      price = testing$aggregate_price,
      marginal_profits = 0,
      use = 0
    )

    model <-  mle_model

    out_training <- independent_data %>%
      as_data_frame() %>%
      mutate(pred = in_sample_prediction)

    out_testing <- testing %>%
      as_data_frame() %>%
      mutate(pred = prediction)


  } #close structural

  return(list(model = model, test_predictions = out_testing,
         training_predictions = out_training))

}