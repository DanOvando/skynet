
orig_tree <- orig_tree[orig_tree != "random_var"]

# tree_candidate_vars <- sample(orig_tree, length(orig_tree) - 1, replace = FALSE)

test_train <- skynet_models %>%
  filter(
    data_subset == "skynet",
    test_sets == "random",
    model %in% c("gbm")
  )  %>%
  group_by(model, variables) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(candidate_vars = ifelse(
    str_detect(.$data_subset, "delta"),
    list(delta_candidate_vars),
    ifelse(
      variables == "gfw_only",
      list(gfw_only_tree_candidate_vars),
      ifelse(
        variables == "enviro_only",
        list(enviro_only_tree_candidate_vars),
        list(tree_candidate_vars)
      )
    )
  )) %>%   mutate(
    fitted_model = pmap(
      list(
        data_subset = data_subset,
        dep_var = dep_var,
        model_name = model,
        training = train,
        testing = test,
        tree_candidate_vars = candidate_vars
      ),
      prep_train,
      fitcontrol_number = 10,
      fitcontrol_repeats = 2,
      never_ind_vars = never_ind_vars,
      tune_model = T,
      cores = num_cores,
      data_sources = data_sources
    )
  )

kfold_preds <- map(test_train$fitted_model, "kfold_preds") %>%
  set_names(test_train$variables)


kfold_preds_summary <-  kfold_preds %>%
  map_df(~., .id = "variables") %>%
  select(Resample, variables, obs, pred,) %>%
  group_by(Resample, variables) %>%
  summarise(model_rmse = sqrt(mean((obs - pred)^2))) %>%
  rename(id = Resample) %>%
  spread(variables, model_rmse)

rmse_mod <- perf_mod(kfold_preds_summary, seed = 4344, iter = 5000)

rmse_post <- tidy(rmse_mod)

stacked_summary <- kfold_preds_summary %>%
  gather(variables, model_rmse, -id)

ggplot(rmse_post) +
  geom_point(
    data = stacked_summary,
    aes(x = variables, y = model_rmse),
    alpha = .5, col = "blue"
  )


b <- map_dbl(test_train$fitted_model, ~max(.x$model$results$Rsquared))

test_train$r2 <-  b

test_train$fitted_model[[3]]$model %>% plot()


map_df(test_train$fitted_model, ~.x$model$bestTune)


test_train$fitted_model[[3]]$model$finalModel %>% ranger::importance()

test_train %>%
  ggplot(aes(variables,r2)) +
  geom_col()


