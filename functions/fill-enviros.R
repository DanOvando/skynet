fill_enviros <- function(variable, data){

  variable[is.nan(variable)] <-  NA

  data <- data %>%
    mutate(thing = variable)

knn_fit <-
  train(
    thing ~ rounded_lon + rounded_lat + year,
    data = data %>% na.omit(),
    method = "knn",
    preProcess = c("center", "scale"),
    tuneLength = 20,
    tuneGrid = data.frame(k = 2:10)
  )

pred <- data %>%
  select(thing,rounded_lon, rounded_lat, year) %>%
  recipe(thing ~ .) %>%
  step_center(all_predictors()) %>%
  step_scale(all_predictors()) %>%
  # step_knnimpute(all_outcomes(), K = knnFit$bestTune %>% as.numeric(),
  #                impute_with =  imp_vars(all_predictors())) %>%
  prep(retain = T) %>%
  juice()

knn_model <- caret::knnreg(thing ~ rounded_lon + rounded_lat + year,
                           data = pred,
                           k = knn_fit$bestTune %>% as.numeric())

pred <- pred %>%
  add_predictions(knn_model) %>%
  mutate(thing = ifelse(is.na(thing), pred, thing))

return(pred$thing)

}

