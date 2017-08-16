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
                       tree_candidate_vars,
                       lm_candidate_vars,
                       fitcontrol_method = 'repeatedcv',
                       fitcontrol_number = 1,
                       fitcontrol_repeats = 1,
                       tune_model = T) {
  # model_formula <-
  #   paste(dep_var, '~', ind_vars) %>% as.formula() # construct model forumla
  if (model == 'rf') {
    #fit random forest
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

    dependent_data <- dependent_data[,dep_var] %>% as.numeric()

    # check <-  training %>% as.data.frame() %>% na.omit()
    #
    #
    print('hello')
    model <- train(
      independent_data,
      dependent_data,
      method = "rf",
      do.trace = 10,
      na.action = na.omit,
      trControl = fit_control
    )
    # model <- train(
    #   model_formula,
    #   data = training %>% as_data_frame(),
    #   method = "rf",
    #   do.trace = 10,
    #   na.action = na.omit,
    #   trControl = fit_control
    # )

    # model1 <- model

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

      independent_data <- independent_data[,vars_to_include]


      # model_formula <-
      #   paste(dep_var, '~', vars_to_include) %>% as.formula() # construct model forumla
      #
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

  }

  return(list(model = model, test_predictions = out_testing ))

}