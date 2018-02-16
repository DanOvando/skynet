#' Fit Skynet
#' A wrapper function for fitting models to skynet data
#'
#' @param dep_var the dependent variable being considered
#' @param model the model to be used
#' @param training data to be used in training
#' @param testing data to be used in testing
#' @param fitcontrol_method caret options
#' @param fitcontrol_number number of k-fold cross validation folds
#' @param fitcontrol_repeats number of times k-fold cross validation is run
#'
#' @return a list of model fits
#' @export
#'
fit_skynet <- function(dep_var,
                       model_name,
                       training,
                       testing,
                       never_ind_vars,
                       survey_prices,
                       tree_candidate_vars,
                       lm_candidate_vars = NA,
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
                       tune_model = T,
                       cores = 4) {
  # model_formula <-
  #   paste(dep_var, '~', ind_vars) %>% as.formula() # construct model forumla
  #
  #


  doMC::registerDoMC(cores = cores)
  # cluster <- makeCluster(cores) # convention to leave 1 core for OS
  # doParallel::registerDoParallel(cluster)

  fit_control <- trainControl(method = fitcontrol_method,
                              number = fitcontrol_number,
                              repeats = fitcontrol_repeats,
                              allowParallel = TRUE)

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

  if (model_name == 'cforest') {
    model <- train(
      model_formula,
      data = reg_data,
      method = "cforest",
      trControl = fit_control,
      preProcess = c("center","scale")
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
        trControl = fit_control,
        preProcess = c("center","scale")

      )

    }

    out_training <- training %>%
      as_data_frame() %>%
      add_predictions(model)

    out_testing <- testing %>%
      as_data_frame() %>%
      add_predictions(model)

  } # close cforest

  if (model_name == 'ranger') {
    #fit random forest

    # rf_grid <-  expand.grid(mtry = c(3,5))
    model <- train(
      independent_data %>% dmap_if(is.logical, as.numeric),
      dependent_data,
      method = "ranger",
      # na.action = na.omit,
      trControl = fit_control,
      preProcess = c("center","scale"),
      importance = "impurity_corrected"
      # tuneGrid = rf_grid
    )

    if (tune_model == T) {
      impnames <-
        model$finalModel$variable.importance %>% as.data.frame() %>% purrr::pluck(attr_getter("row.names"))
      varimp <-  model$finalModel$variable.importance %>%
        as_data_frame() %>%
        mutate(varname = impnames) %>%
        arrange(desc(value)) %>%
        mutate(varname = str_replace(varname, '(TRUE)|(FALSE)', ''))

      random_imp <-
        varimp$value[varimp$varname == 'random_var']


      vars_to_include <- varimp %>%
        filter(value > random_imp) %>% {
          .$varname
        } %>% unique()

      independent_data <- independent_data[, vars_to_include]


      model <- train(
        independent_data %>% dmap_if(is.logical, as.numeric),
        dependent_data,
        method = "ranger",
        # na.action = na.omit,
        trControl = fit_control,
        preProcess = c("center","scale"),
        importance = "impurity_corrected",
        verbose = TRUE
        # tuneGrid = rf_grid
      )


    }
    out_testing <- testing %>%
      as_data_frame() %>%
      add_predictions(model)

    out_training <- training %>%
      as_data_frame() %>%
      add_predictions(model)
  } # close ranger rf


  if (model_name == 'rf') {
    #fit random forest

    # rf_grid <-  expand.grid(mtry = c(3,5))

    model <- train(
      independent_data,
      dependent_data,
      method = "rf",
      do.trace = 10,
      # na.action = na.omit,
      trControl = fit_control,
      preProcess = c("center","scale")

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
        trControl = fit_control,
        importance = T,
        preProcess = c("center","scale")
      )

    }

    out_testing <- testing %>%
      as_data_frame() %>%
      add_predictions(model)

    out_training <- training %>%
      as_data_frame() %>%
      add_predictions(model)

  } # close rf

  if (model_name == 'gbm') {

    # gbm_grid <-  expand.grid(interaction.depth = c(1, 5),
    #                         n.trees = (1:10)*100,
    #                         shrinkage = 0.01,
    #                         n.minobsinnode = 20)

    model <- train(
      model_formula,
      data = reg_data,
      method = "gbm",
      verbose = T,
      trControl = fit_control,
      preProcess = c("center","scale")
      # tuneGrid = gbm_grid
    )
    # on.exit(detach('package:plyr'))


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
        trControl = fit_control,
        preProcess = c("center","scale")
      )

      # on.exit(detach('package:plyr'))

    }
    out_testing <- testing %>%
      as_data_frame() %>%
      add_predictions(model)

    out_training <- training %>%
      as_data_frame() %>%
      add_predictions(model)

  }

  if (model_name == 'structural') {

    compile(here::here('scripts','fit_structural_skynet.cpp'))

    dyn.load(dynlib(here::here("scripts","fit_structural_skynet")))

    independent_data <- training %>%
      as.data.frame() %>%
      select(matches(paste0(structural_vars, collapse = '|')))

    testing_frame <- testing %>%
      as.data.frame() %>%
      select(matches(paste0(structural_vars, collapse = '|')))

    struct_data <- list(
      data = as.matrix(independent_data %>%
                         select(dist_from_port,
                                mean_vessel_length,
                                total_engine_power,
                                no_take_mpa,
                                restricted_use_mpa,
                                m_below_sea_level)),
      log_d = as.numeric(independent_data$log_density),
      effort = as.numeric(independent_data$total_hours),
      price = as.numeric(independent_data$aggregate_price),
      mp = 0,
      n = nrow(independent_data)
    )

    struct_params <- list(
      betas = rep(0, ncol(struct_data$data)),
      log_sigma = log(sd(struct_data$log_d)),
      logit_q = log(.001 / (1 - .001))
    )

    model <- MakeADFun(data=struct_data,parameters=struct_params)

    mle_fit <-
      nlminb(
        model$par,
        objective = model$fn,
        gradient = model$gr,
        control = list("trace" = 100)
      )

    mle_fit_report <- model$report()

    model <-  list(mle_fit = mle_fit,
                   mle_report = mle_fit_report
                  )

  testing_prediction <- predict_structural_model(mle_fit = mle_fit,
                                    data = testing_frame,
                                    mle_vars = colnames(struct_data$data),
                                    mp = 0)

    out_training <- training %>%
      as_data_frame() %>%
      mutate(pred = mle_fit_report$log_d_hat)

    # out_training %>%
    #   ggplot(aes(pred, pred2, color = total_hours > 0)) +
    #   geom_point()

    out_testing <- testing %>%
      as_data_frame() %>%
      mutate(pred = testing_prediction)

    # out_testing %>%
    #   ggplot(aes(log_density, pred, color = total_hours > 0)) +
    #   geom_point()

  } #close structural

  return(list(model = model, test_predictions = out_testing,
         training_predictions = out_training))

}