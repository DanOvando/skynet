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
fit_skynet <- function(data_subset,
                       dep_var,
                       model_name,
                       training,
                       testing,
                       tuned_pars,
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
  # model_formula <-
  #   paste(dep_var, '~', ind_vars) %>% as.formula() # construct model forumla
  #
  #



  doMC::registerDoMC(cores = cores)
  # cluster <- makeCluster(cores) # convention to leave 1 core for OS
  # doParallel::registerDoParallel(cluster)

  dat <-
    data_sources$data[[which(data_sources$data_subset == data_subset)]]

  dat <- dat %>%
    dmap_if(is.logical, as.numeric)

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

  if (model_name == "mars") {
    dep_var <- glue::glue("log_", dep_var)

  }

  dependent_data <- dependent_data[, dep_var] %>% as.numeric()

  model_formula <-
    paste0(dep_var, '~', paste(colnames(independent_data), collapse = '+')) %>%
    as.formula()

  reg_data <- independent_data %>%
    as_data_frame() %>%
    mutate(!!dep_var := dependent_data)


  train_recipe <-
    recipes::recipe(glue::glue("{dep_var} ~ .") %>% as.formula(),
                    data = reg_data) %>%
    step_nzv(all_predictors()) %>%
    step_BoxCox(all_predictors()) %>%
    step_knnimpute(all_predictors())
  # step_center(all_predictors()) %>%
  # step_scale(all_predictors())

  prep_recipe <- prep(train_recipe, data = reg_data, retain = T)


  if (model_name == 'cforest') {
    model <- train(
      model_formula,
      data = reg_data,
      method = "cforest",
      trControl = fit_control,
      preProcess = c("center", "scale"),
      weights = weights
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
        preProcess = c("center", "scale"),
        weights = weights

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

    tuned_pars <- tuned_pars %>%
      pluck(model_name)

    prepped_recipe <-
      recipes::prep(train_recipe, training, retain = T)

    proc_training <- juice(prepped_recipe)

    model <-
      ranger::ranger(
        formula = glue::glue("{dep_var} ~ .") %>% as.formula(),
        data = proc_training,
        mtry = pmin(pmax(1, ncol(proc_training) - 2), tuned_pars$mtry),
        splitrule = tuned_pars$splitrule,
        min.node.size = tuned_pars$min.node.size,
        importance = "impurity_corrected",
        verbose = TRUE
      )

    if (tune_model == T) {
      impnames <-
        model$variable.importance %>% as.data.frame() %>% purrr::pluck(attr_getter("row.names"))

      varimp <-  model$variable.importance %>%
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

      reg_data <- reg_data[, c(dep_var, vars_to_include)]

      train_recipe <-
        recipes::recipe(glue::glue("{dep_var} ~ .") %>% as.formula(),
                        data = reg_data) %>%
        step_nzv(all_predictors()) %>%
        step_BoxCox(all_predictors()) %>%
        step_knnimpute(all_predictors())
      #
      # step_center(all_predictors()) %>%
      # step_scale(all_predictors())

      prepped_recipe <-
        recipes::prep(train_recipe, training, retain = T)

      proc_training <- juice(prepped_recipe)


      model <-
        ranger::ranger(
          formula = glue::glue("{dep_var} ~ .") %>% as.formula(),
          data = proc_training,
          mtry = pmin(pmax(1, ncol(
            proc_training
          ) - 2), tuned_pars$mtry),
          splitrule = tuned_pars$splitrule,
          min.node.size = tuned_pars$min.node.size,
          verbose = TRUE
        )


    }

    proc_testing <- bake(prepped_recipe, newdata = testing)

    training_pred <- model$predictions

    testing_pred <-
      predict(model,
              data =  proc_testing)

    out_testing <- testing %>%
      mutate(pred = testing_pred$predictions)

    out_training <- training %>%
      mutate(pred = training_pred)

  } # close ranger rf


  if (model_name == 'gbm') {
    # on.exit(detach('package:plyr'))

    tuned_pars <- tuned_pars %>%
      pluck(model_name)

    prepped_recipe <-
      recipes::prep(train_recipe, training, retain = T)

    proc_training <- juice(prepped_recipe)

    model <-
      gbm(
        formula = glue::glue("{dep_var} ~ .") %>% as.formula(),
        data = proc_training,
        n.trees = tuned_pars$n.trees,
        interaction.depth = tuned_pars$interaction.depth,
        shrinkage = tuned_pars$shrinkage,
        n.minobsinnode = tuned_pars$n.minobsinnode,
        weights = weights,
        verbose = TRUE,
        keep.data = FALSE
      )


    if (tune_model == T) {
      var_importance <- summary(model)

      varimp <- var_importance %>%
        mutate(varname = rownames(.)) %>%
        arrange(desc(rel.inf)) %>%
        mutate(varname = str_replace(varname, '(TRUE)|(FALSE)', ''))


      random_imp <-
        varimp$rel.inf[varimp$varname == 'random_var']

      vars_to_include <- varimp %>%
        filter(rel.inf > random_imp) %>% {
          .$varname
        } %>% unique()

      reg_data <- reg_data[, c(dep_var, vars_to_include)]

      model_formula <-
        paste0(dep_var, '~', paste(vars_to_include, collapse = '+')) %>%
        as.formula()

      train_recipe <-
        recipes::recipe(glue::glue("{dep_var} ~ .") %>% as.formula(),
                        data = reg_data) %>%
        step_nzv(all_predictors()) %>%
        step_BoxCox(all_predictors()) %>%
        step_knnimpute(all_predictors())
      #
      # step_center(all_predictors()) %>%
      # step_scale(all_predictors())

      prepped_recipe <-
        recipes::prep(train_recipe, training, retain = T)

      proc_training <- juice(prepped_recipe)


      model <-
        gbm(
          formula = glue::glue("{dep_var} ~ .") %>% as.formula(),
          data = proc_training,
          n.trees = tuned_pars$n.trees,
          interaction.depth = tuned_pars$interaction.depth,
          shrinkage = tuned_pars$shrinkage,
          n.minobsinnode = tuned_pars$n.minobsinnode,
          weights = weights,
          verbose = TRUE,
          keep.data = FALSE
        )

    }


    proc_testing <- bake(prepped_recipe, newdata = testing)

    training_pred <-
      predict(model,
              newdata =  proc_training,
              n.trees = model$n.trees)

    testing_pred <-
      predict(model,
              newdata =  proc_testing,
              n.trees = model$n.trees)

    out_testing <- testing %>%
      mutate(pred = testing_pred)

    out_training <- training %>%
      mutate(pred = training_pred)

  }

  if (model_name == "mars") {
    tuned_pars <- tuned_pars %>%
      pluck(model_name)

    prepped_recipe <-
      recipes::prep(train_recipe, training, retain = T)

    proc_training <- juice(prepped_recipe)

    model <-
      earth::earth(
        formula = glue::glue("{dep_var} ~ .") %>% as.formula(),
        data = proc_training,
        nprune = tuned_pars$nprune,
        degree = tuned_pars$degree,
        trace = 1,
        glm = list(family = gaussian)
      )

    proc_testing <- bake(prepped_recipe, newdata = testing)

    se_estimate <- sd(model$residuals)

    training_pred <-
      predict(model,
              newdata =  proc_training)

    testing_pred <-
      predict(model,
              newdata =  proc_testing)

    out_testing <- testing %>%
      mutate(pred = exp(testing_pred + se_estimate ^ 2 / 2) %>% as.numeric())

    out_training <- training %>%
      mutate(pred = exp(training_pred + se_estimate ^ 2 / 2) %>% as.numeric())


  }

  if (model_name == 'structural') {
    compile(here::here('src', 'fit_structural_skynet.cpp'))

    dyn.load(dynlib(here::here("src", "fit_structural_skynet")))

    independent_data <- training %>%
      as.data.frame() %>%
      select(matches(paste0(structural_vars, collapse = '|'))) %>%
      mutate(intercept = 1)

    testing_frame <- testing %>%
      as.data.frame() %>%
      select(matches(paste0(structural_vars, collapse = '|'))) %>%
      mutate(intercept = 1)

    struct_data <- list(
      data = as.matrix(
        independent_data %>%
          select(
            dist_from_port,
            mean_vessel_length,
            total_engine_power,
            no_take_mpa,
            restricted_use_mpa,
            m_below_sea_level,
            intercept,
            dist_from_shore
          )
      ),
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

    model <-
      MakeADFun(data = struct_data,
                parameters = struct_params,
                DLL = "fit_structural_skynet")

    mle_fit <-
      nlminb(
        model$par,
        objective = model$fn,
        gradient = model$gr,
        control = list("trace" = 100)
      )

    mle_fit_report <- model$report()

    model <-  list(mle_fit = mle_fit,
                   mle_report = mle_fit_report)

    testing_prediction <-
      predict_structural_model(
        mle_fit = mle_fit,
        data = testing_frame,
        mle_vars = colnames(struct_data$data),
        mp = 0
      )


    out_training <- training %>%
      as_data_frame() %>%
      mutate(pred = exp(mle_fit_report$log_d_hat + mle_fit_report$sigma ^
                          2 / 2))

    out_testing <- testing %>%
      as_data_frame() %>%
      mutate(pred = exp(testing_prediction + mle_fit_report$sigma ^ 2 /
                          2))

  } #close structural

  if (model_name == 'hours') {
    independent_data <- training %>%
      as.data.frame() %>%
      select(dep_var, structural_vars) %>%
      filter(total_hours > 0)

    testing_frame <- testing %>%
      as.data.frame() %>%
      select(dep_var, structural_vars) %>%
      filter(total_hours > 0)


    model <-
      lm(as.formula(glue::glue("log_{dep_var} ~ log(total_hours)")), data = independent_data %>% filter(total_hours > 0))

    se_estimate <- sd(model$residuals)


    out_training <- training %>%
      as_data_frame() %>%
      mutate(pred = exp(predict(model, newdata = training) + se_estimate ^ 2 / 2))

    out_testing <- testing %>%
      as_data_frame() %>%
      mutate(pred = exp(predict(model, newdata = testing) + se_estimate ^ 2 / 2))

  } # close hours model

  if (model_name == 'engine_power') {
    independent_data <- training %>%
      as.data.frame() %>%
      select(dep_var, structural_vars) %>%
      # select(matches(paste0(structural_vars, collapse = '|'))) %>%
      filter(total_engine_power > 0)

    testing_frame <- testing %>%
      as.data.frame() %>%
      select(dep_var, structural_vars) %>%
      filter(total_engine_power > 0)

    model <-
      lm(as.formula(glue::glue(
        "log_{dep_var} ~ log(total_engine_power)"
      )), data = independent_data)

    se_estimate <- sd(model$residuals)

    out_training <- training %>%
      as_data_frame() %>%
      mutate(pred = exp(predict(model, newdata = training) + se_estimate ^ 2 / 2))

    out_testing <- testing %>%
      as_data_frame() %>%
      mutate(pred = exp(predict(model, newdata = testing) + se_estimate ^ 2 / 2))


  } # close hours model

  return(list(
    test_predictions = out_testing,
    training_predictions = out_training
  ))

}
