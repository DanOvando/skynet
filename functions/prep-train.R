#' prep train
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
prep_train <- function(data_subset,
                       dep_var,
                       model_name,
                       training,
                       testing,
                       never_ind_vars,
                       survey_prices,
                       tree_candidate_vars,
                       weight_surveys = T,
                       lm_candidate_vars = NA,
                       structural_vars = c(
                         'log_density',
                         'dist_from_port',
                         "dist_from_shore",
                         'mean_vessel_length',
                         'm_below_sea_level',
                         'restricted_use_mpa',
                         'no_take_mpa',
                         'total_hours',
                         'total_engine_power',
                         "aggregate_price"
                       ),
                       fitcontrol_method = 'repeatedcv',
                       fitcontrol_number = 1,
                       fitcontrol_repeats = 1,
                       tune_model = T,
                       cores = 4,
                       data_sources) {
  doMC::registerDoMC(cores = cores)
  # cluster <- makeCluster(cores) # convention to leave 1 core for OS
  # doParallel::registerDoParallel(cluster)

  fit_control <- trainControl(
    method = fitcontrol_method,
    number = fitcontrol_number,
    repeats = fitcontrol_repeats,
    allowParallel = TRUE
  )
  dat <-
    data_sources$data[[which(data_sources$data_subset == data_subset)]]

  training <- dat %>%
    filter(index %in% training)

  testing <- dat %>%
    filter(index %in% testing)


  independent_data <- training %>%
    # as.data.frame() %>%
    select(-matches(paste0(never_ind_vars, collapse = '|'))) %>%
    select(matches(paste0(tree_candidate_vars, collapse = '|')))

  if (weight_surveys == T) {
    weights <- training %>%
      # as_data_frame() %>%
      group_by(survey) %>%
      count() %>%
      mutate(weight = 1 / n)

    weights <- training %>%
      # as_data_frame() %>%
      select(survey, num_vessels) %>%
      left_join(weights, by = "survey")

    weights <- weights$weight

  } else{
    weights <- rep(1, nrow(training))
  }

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


  train_recipe <- recipes::recipe(model_formula,
                                  data = reg_data) %>%
    step_nzv(all_predictors()) %>%
    step_center(all_predictors()) %>%
    step_scale(all_predictors())

  prep_recipe <- prep(train_recipe, data = reg_data, retain = T)

  if (model_name == 'ranger') {
    #fit random forest
browser()

    rf_grid <-  expand.grid(mtry = c(3,5))



    model <- train(
      independent_data %>% dmap_if(is.logical, as.numeric),
      dependent_data,
      method = "ranger",
      # na.action = na.omit,
      trControl = fit_control,
      tuneGrid = rf_grid,
      preProcess = c("center", "scale"),
      importance = "impurity_corrected",
      weights = weights
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
        preProcess = c("center", "scale"),
        importance = "impurity_corrected",
        verbose = TRUE,
        weights = weights
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

  if (model_name == 'gbm') {
    on.exit(detach('package:plyr'))

    gbm_grid <-  expand.grid(
      interaction.depth = c(3),
      n.trees = c(3000),
      shrinkage = c(0.005),
      n.minobsinnode = c(10)
    )

    model <- caret::train(
      train_recipe,
      data = reg_data,
      method = "gbm",
      verbose = T,
      trControl = fit_control,
      tuneGrid = gbm_grid
    )

  }

  out <- model$bestTune %>% mutate(model = model_name)

  return(out)

}