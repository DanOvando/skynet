fit_skynet <- function(dep_var,
                       ind_vars,
                       model,
                       training,
                       testing,
                       fitcontrol_method = 'repeatedcv',
                       fitcontrol_number = 1,
                       fitcontrol_repeats = 1) {
  model_formula <-
    paste(dep_var, '~',ind_vars) %>% as.formula()
  if (model == 'rf') {
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