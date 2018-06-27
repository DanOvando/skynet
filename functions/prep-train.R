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

  tree_candidate_vars <- tree_candidate_vars[tree_candidate_vars != "random_var"]

  doMC::registerDoMC(cores = cores)

  fit_control <- trainControl(
    method = fitcontrol_method,
    number = fitcontrol_number,
    repeats = fitcontrol_repeats,
    savePredictions = "final",
    allowParallel = TRUE
  )
  dat <-
    data_sources$data[[which(data_sources$data_subset == data_subset)]]

  training <- dat %>%
    filter(index %in% training)

  testing <- dat %>%
    filter(index %in% testing)


  independent_data <- training %>%
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


  if (str_detect(model_name,"mars")) {

    dep_var <- glue::glue("log_",dep_var)

  }

  dependent_data <- dependent_data[, dep_var] %>% as.numeric()

  model_formula <-
    paste0(dep_var, '~', paste(colnames(independent_data), collapse = '+')) %>%
    as.formula()

  reg_data <- independent_data %>%
    as_data_frame() %>%
    mutate(!!dep_var := dependent_data) %>%
    dmap_if(is.logical, as.numeric)


  train_recipe <- recipes::recipe(glue::glue("{dep_var} ~ .") %>% as.formula(),
                                  data = reg_data) %>%
    step_nzv(all_predictors()) %>%
    step_BoxCox(all_predictors())

  prep_recipe <- prep(train_recipe, data = reg_data, retain = T)

  if (model_name == 'ranger') {
    #fit random forest

    default <- ncol(independent_data) %>% sqrt() %>% floor()

    rf_grid <-
      expand.grid(
        mtry = c(2, default, ncol(independent_data) - 2),
        splitrule = c("variance", "extratrees"),
        min.node.size = c(5,10,20,50)
      )

    set.seed(42)
    model <- train(
      train_recipe,
      reg_data,
      method = "ranger",
      trControl = fit_control,
      tuneGrid = rf_grid,
      importance = "impurity"
    )

  } # close ranger rf

  if (model_name == 'gbm') {
    on.exit(detach('package:plyr'))

    # gbm_grid <-  expand.grid(
    #   interaction.depth = c(3,5,10),
    #   n.trees = c(3000,6000,8000),
    #   shrinkage = c(0.001,0.005),
    #   n.minobsinnode = c(10,20)
    # )

    gbm_grid <-  expand.grid(
      interaction.depth = c(2,5,8),
      n.trees = c(5000,10000),
      shrinkage = c(0.001),
      n.minobsinnode = c(5,10)
    )

    set.seed(42)
    model <- caret::train(
      train_recipe,
      data = reg_data,
      method = "gbm",
      verbose = T,
      trControl = fit_control,
      tuneGrid = gbm_grid
    )

  }

  if (model_name == 'mars') {

    on.exit(detach('package:earth'))

    default <- ncol(independent_data) %>% sqrt() %>% floor()

    mars_grid <-  expand.grid(nprune = c(ncol(independent_data) + 1, default),
                              # number of terms
                              degree = 1:8)

    set.seed(42)
    model <- caret::train(
      train_recipe,
      data = reg_data,
      method = "earth",
      trControl = fit_control,
      tuneGrid = mars_grid,
      trace = 1
    )

    se_estimate <- sd(model$pred$obs - model$pred$pred)

    model$pred$obs <- exp(model$pred$obs)

    model$pred$pred <- exp(model$pred$pred + se_estimate^2/2)

  }

  if (model_name == 'bagged_mars') {

    mars_grid <-  expand.grid(degree = 1:8)

    set.seed(42)
    model <- caret::train(
      train_recipe,
      data = reg_data,
      method = "bagEarthGCV",
      trControl = fit_control,
      tuneGrid = mars_grid,
      trace = 1,
      B = 50
    )

    se_estimate <- sd(model$pred$obs - model$pred$pred)

    model$pred$obs <- exp(model$pred$obs)

    model$pred$pred <- exp(model$pred$pred + se_estimate^2/2)

  }


  out <- list(tuned_pars = model$bestTune %>% mutate(model = model_name),
              kfold_preds = model$pred, model = model)

  return(out)

}