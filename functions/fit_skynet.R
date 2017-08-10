#' Fit Skynet
#' A wrapper function for fitting models to skynet data
#'
#' @param dep_var the dependent variable being considered
#' @param ind_vars the independent variabels being considered
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
                       ind_vars,
                       model,
                       training,
                       testing,
                       fitcontrol_method = 'repeatedcv',
                       fitcontrol_number = 1,
                       fitcontrol_repeats = 1) {
  model_formula <-
    paste(dep_var, '~',ind_vars) %>% as.formula() # construct model forumla
  if (model == 'rf') { #fit random forest
    fit_control <- trainControl(method = fitcontrol_method,
                                number = fitcontrol_number,
                                repeats = fitcontrol_repeats)

    # check <-  training %>% as.data.frame() %>% na.omit()
    model <- train(
      model_formula,
      data = training %>% as_data_frame(),
      method = "rf",
      do.trace = 10,
      na.action = na.omit,
      trControl = fit_control
    )

  }

  return(list(model = model))

}