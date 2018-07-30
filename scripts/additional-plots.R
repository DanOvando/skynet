data <- skynet_data %>%
  select(any_fishing, dist_from_port,dist_from_shore,
         no_take_mpa,restricted_use_mpa, density, mean_analysed_sst,
         mean_chla )

fish_recipe <- recipe(any_fishing ~ . , data = data) %>%
  step_center(all_predictors(), -no_take_mpa,-restricted_use_mpa) %>%
  step_scale(all_predictors(), -no_take_mpa,-restricted_use_mpa)

prepped_data <- prep(fish_recipe, data = data, retain = T)

prepped_data <- prepped_data %>% juice()

test <- glm(any_fishing ~., data = prepped_data, family = binomial)

summary(test)


# latent spatial variables


a <- candidate_data %>%
  filter(fished_only == TRUE,
         unfished_only == FALSE,
         survey_months_only == TRUE)

data <- a$skynet_data[[1]] %>%
  filter(survey == "wcgbts") %>% mutate(state = ifelse(
    survey == "wcgbts",
    if_else(
      rounded_lat > 46.25,
      "washington",
      if_else(rounded_lat > 42, "oregon", "california")
    ),
    "Alaska"
  ))

oregon <- data %>%
  filter(state == "oregon", total_hours > 0) %>%
  filter(year == max(year))

oregon_cost <- oregon %>%
  select(dist_from_port, dist_from_shore, m_below_sea_level)

cost_recipe <- recipe(oregon_cost) %>%
  step_log(all_numeric())
# step_center(all_numeric()) %>%
# step_scale(all_numeric())
#
oregon_cost <-
  prep(cost_recipe, oregon_cost, retain = T) %>% juice() %>%
  mutate(intercept = 1)


location_year <- oregon %>%
  select(knot, year) %>%
  mutate(knot_year = paste(knot, year,sep = '_')) %>%
  select(knot_year) %>%
  mutate(marker = 1) %>%
  mutate(index = 1:nrow(.)) %>%
  spread(knot_year, marker, fill = 0) %>%
  select(-index)

warmups <- 5000

total_iterations <- 10000

max_treedepth <-  12

adapt_delta <-  0.8

chains <- 1

struct_data <- list(
  cost_data = oregon_cost,
  location_data = location_year,
  log_effort = as.numeric(oregon$total_hours %>% log()),
  price = oregon$aggregate_price,
  mp = 0,
  n = nrow(oregon),
  max_q = .8/max(oregon$total_hours)
)

struct_data$n_betas <- ncol(struct_data$cost_data)

struct_data$n_cpue_betas <- ncol(struct_data$location_data)

# inits <- list(list(cpue_betas = rep(25, struct_data$n_cpue_betas)))


# stan_fit <-
#   stan(
#     file = here::here("src", "fit_effort_model.stan"),
#     data = struct_data,
#     chains = chains,
#     warmup = warmups,
#     iter = total_iterations,
#     cores = 1,
#     refresh = 250,
#     control = list(max_treedepth = max_treedepth,
#                    adapt_delta = adapt_delta)
#   )
#
# fits <- stan_fit %>%
#   tidy()
#
# knots <- fits %>%
#   filter(str_detect(term, "knot"))
#
# effort <- fits %>%
#   filter(str_detect(term, "effort") & !str_detect(term, "log")) %>%
#   mutate(effort = oregon$total_hours)
#
# effort %>%
#   ggplot(aes(effort, estimate)) +
#   geom_point()



# which data subset

best_data_subset <- best_fits %>%
  filter(dep_var == "density", unfished_only == FALSE) %>%
  group_by(data_subset) %>%
  mutate(median_r2 = median(r2_training, na.rm = T)) %>%
  ungroup() %>%
  mutate(data_subset = fct_reorder(data_subset, median_r2)) %>%
  arrange(desc(median_r2))

best_data_subset %>%
  ggplot() +
  geom_point(position = "dodge",aes(data_subset,r2_training), alpha = 0.5) +
  geom_point(aes(data_subset, median_r2), color = "red", size = 2) +
  coord_flip()

best_data_subset <- "skynet"


best_fits <- skynet_models %>%
  filter(model == best_ml_model | model == "structural",
         dep_var == "density",
         variables == "gfw_only") %>%
  mutate(data_subset = if_else(data_subset == "lag_1_skynet_data", "skynet", data_subset))

best_fits <- best_fits %>%
  filter(data_subset == best_data_subset)

best_fits %>%
  group_by(test_set) %>%
  mutate(mean_r2 = mean(r2, na.rm = T)) %>%
  ungroup() %>%
  mutate(test_set = fct_reorder(test_set, mean_r2)) %>%
  gather(r2_source, r2_value, r2, r2_training) %>%
  ggplot(aes(test_set, r2_value)) +
  geom_line() +
  geom_point(aes(color = r2_source), size = 4) +
  coord_flip() +
  theme(axis.title.y = element_blank()) +
  labs(y = bquote(R ^ 2)) +
  scale_color_manual(
    name = "Data Split",
    labels = c("Testing", "Training"),
    values = wes_palette("Zissou1")
  ) +
  facet_wrap(~model)


# spatial performance -----------------------------------------------------

spatial_performance <- best_fits %>%
  filter(test_set == "spatial")  %>% {
    .$test_data[[1]]
  } %>%
  filter(surveyed_year == T) %>%
  group_by(survey) %>%
  mutate(sd_density = sd(density, na.rm = T)) %>%
  group_by(rounded_lat, rounded_lon, survey) %>%
  summarise(resid = mean(pred - density),
            scaled_resid = pmax(-2,pmin(2,resid / mean(sd_density))))

spatial_performance_plot <- spatial_performance %>%
  mutate(recenter_lon = ifelse(rounded_lon < 0, 180 + (180 - abs(rounded_lon)), rounded_lon)) %>%
  dplyr::mutate(geometry = purrr::map2(recenter_lon, rounded_lat, ~ sf::st_point(x = c(.x, .y), dim = 'XY'))) %>%
  ungroup() %>%
  mutate(geometry = sf::st_sfc(geometry, crs =
                                 "+proj=longlat +datum=WGS84 +no_defs")) %>%
  sf::st_sf()

bbox <- sf::st_bbox(spatial_performance_plot)

alaska_bbox <- sf::st_bbox(spatial_performance_plot %>% filter(survey != "wcgbts"))

wc_bbox <- sf::st_bbox(spatial_performance_plot %>% filter(survey == "wcgbts"))


resid_hist <- spatial_performance %>%
  ggplot(aes(scaled_resid)) +
  geom_density(color = "black", fill = "grey", alpha = 0.5) +
  labs(x = "Scaled Residuals", caption = "Residuals divided by standard deviation of observed densities") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"),
        axis.text.x = element_text(size = 8))



alaska_spatial_residual_plot <- ggplot() +
  geom_sf(data = spatial_performance_plot %>% filter(survey != "wcgbts"), aes(color = scaled_resid),size = 0.5, alpha = 0.75) +
  geom_sf(data = pacific_map, shape = 21) +
  coord_sf(xlim = c(alaska_bbox['xmin'], alaska_bbox['xmax']),
           ylim = c(alaska_bbox['ymin']*.95, alaska_bbox['ymax'] * 1.1), expand = F) +
  scale_color_gradient2(low = "tomato", high = "steelblue",midpoint = 0, mid = "grey", name = "Scaled Residuals", guide = guide_colorbar(frame.colour = "black")) +
  map_theme +
  labs(caption = "Alaska") +
  theme(legend.key.height = unit(1,"cm"))


wc_spatial_residual_plot <- ggplot() +
  geom_sf(data = spatial_performance_plot %>% filter(survey == "wcgbts"), aes(color = scaled_resid), size = 0.5, alpha = 0.75) +
  geom_sf(data = pacific_map) +
  coord_sf(xlim = c(wc_bbox['xmin']*.98, wc_bbox['xmax']),
           ylim = c(wc_bbox['ymin'], wc_bbox['ymax'] * 1.05), expand = F) +
  scale_color_gradient2(low = "tomato", high = "steelblue",midpoint = 0, mid = "grey", name = "Residuals", guide = "none") +
  map_theme +
  labs(caption = "West Coast")


spatial_residual_plot <- (
  wc_spatial_residual_plot + alaska_spatial_residual_plot
) / resid_hist  + plot_layout(nrow = 2, ncol = 1, widths = c(2,1),heights = c(2.5,1)) & theme(plot.margin = unit(rep(0.01,4), "points"))


# more spatial ------------------------------------------------------------

train_performance <- best_fits %>%
  filter(test_set == "california", data_subset == "skynet")  %>%{
    .$training_data[[1]]
  } %>%
  mutate(resid = pred - density) %>%
  mutate(split = "Training")

test_performance <- best_fits %>%
  filter(test_set == "california", data_subset == "skynet")  %>% {
    .$test_data[[1]]
  } %>%
  mutate(resid = pred - density) %>%
  mutate(split = "Testing")

spatial_performance <- train_performance %>%
  bind_rows(test_performance) %>%
  filter(surveyed_year == T) %>%
  group_by(survey) %>%
  mutate(scaled_resid = resid / sd(density))

spatial_performance_plot <- spatial_performance %>%
  mutate(recenter_lon = ifelse(rounded_lon < 0, 180 + (180 - abs(rounded_lon)), rounded_lon)) %>%
  dplyr::mutate(geometry = purrr::map2(recenter_lon, rounded_lat, ~ sf::st_point(x = c(.x, .y), dim = 'XY'))) %>%
  ungroup() %>%
  mutate(geometry = sf::st_sfc(geometry, crs =
                                 "+proj=longlat +datum=WGS84 +no_defs")) %>%
  sf::st_sf()

bbox <- sf::st_bbox(spatial_performance_plot)


resid_hist <- spatial_performance %>%
  ggplot(aes(scaled_resid,fill = split)) +
  geom_density(alpha = 0.75) +
  labs(x = "Scaled Residuals", title = "A)") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  theme(legend.position = "right") +
  scale_fill_discrete(name = '')


spatial_residual_plot <- ggplot() +
  geom_sf(data = spatial_performance_plot %>% filter(split == "Testing"), aes(color = scaled_resid), size = 0.5, alpha = 0.75) +
  geom_sf(data = pacific_map) +
  coord_sf(xlim = c(bbox['xmin'], bbox['xmax']),
           ylim = c(bbox['ymin']*.75, bbox['ymax'] * 1.1), expand = F) +
  scale_shape(guide = FALSE) +
  scale_color_gradient2(low = "tomato", high = "steelblue",midpoint = 0, mid = "lightgrey", name = "Residuals") +
  map_theme +
  labs(title = "B)")


resid_hist + spatial_residual_plot + plot_layout(ncol = 2, nrow = 1, widths = c(1,2))



# no idea -----------------------------------------------------------------

resolution_km2 <- 100

sub_skynet_models <- skynet_models %>%
  mutate(
    run_name = glue::glue(
      "{.$dep_var}--{.$model}--{.$weight_surveys}--{.$test_sets}--{.$data_subset}"
    )
  ) %>%
  filter(model %in% c(best_ml_model),
         dep_var == "density",
         data_subset == "skynet",
         test_set %in% c("random","spatial","alaska", "random_west_coast", "west_coast","year_greq_than_2014"),
         variables == "gfw_only")  %>%
  mutate(test_data = map(test_data, ~.x %>% filter(surveyed_year == T)),
         training_data = map(training_data, ~.x %>% filter(surveyed_year == T)))


downscaled_results <-
  cross_df(list(
    run_name = sub_skynet_models$run_name,
    resolution = resolution_km2
  )) %>%
  left_join(
    sub_skynet_models %>% select(
      run_name,
      test_data,
      training_data,
      dep_var,
      model,
      weight_surveys,
      train_set,
      test_set,
      data_subset
    ),
    by = "run_name"
  ) %>%
  arrange(run_name)

downscaled_results <-  downscaled_results %>%
  mutate(
    training_grid  = map2(
      training_data,
      resolution,
      create_grid,
      lon_name = rounded_lon,
      lat_name = rounded_lat
    )
  ) %>%
  mutate(
    training_data = map2(
      training_data,
      training_grid,
      snap_to_grid,
      old_lon_name = rounded_lon,
      old_lat_name = rounded_lat,
      new_lon_name = lon,
      new_lat_name = lat,
      dep_var = "density"
    )
  )

downscaled_results <-  downscaled_results %>%
  mutate(
    testing_grid  = map2(
      test_data,
      resolution,
      create_grid,
      lon_name = rounded_lon,
      lat_name = rounded_lat
    )
  ) %>%
  mutate(
    test_data = map2(
      test_data,
      testing_grid,
      snap_to_grid,
      old_lon_name = rounded_lon,
      old_lat_name = rounded_lat,
      new_lon_name = lon,
      new_lat_name = lat,
      dep_var = "density"
    )
  )

train_performance <- downscaled_results %>%
  select(test_set,train_set, training_data) %>%
  unnest() %>%
  group_by(train_set) %>%
  mutate(rel_pred = agg_pred / max(agg_pred),
         rel_density = agg_mean_density / max(agg_mean_density)) %>%
  mutate(split = "Training")

test_performance <- downscaled_results %>%
  select(test_set,train_set, test_data) %>%
  unnest() %>%
  group_by(test_set) %>%
  mutate(rel_pred = agg_pred / max(agg_pred),
         rel_density = agg_mean_density / max(agg_mean_density)) %>%
  mutate(split = "Testing")

spatial_performance <- train_performance %>%
  bind_rows(test_performance) %>%
  ungroup() %>%
  filter(!str_detect(test_set, "year_")) %>%
  mutate(set_order = as.numeric(factor(test_set))) %>%
  mutate(test_set = fct_reorder(test_set, set_order),
         train_set = fct_reorder(train_set, set_order))

trained_plot <- spatial_performance %>%
  filter(split == "Training") %>%
  ggplot(aes(rel_density, rel_pred)) +
  geom_point() +
  facet_wrap(~ train_set, strip.position = "left", nrow = n_distinct(spatial_performance$test_set), ncol = 1 ) +
  labs(title = "Trained on...", y = "Predicted")  +
  theme(axis.title.x = element_blank())


tested_plot <- spatial_performance %>%
  filter(split == "Testing") %>%
  ggplot(aes(rel_density, rel_pred)) +
  geom_point() +
  facet_wrap(~ test_set, strip.position = "right", nrow = n_distinct(spatial_performance$test_set), ncol = 1 ) +
  labs(title = "Tested on...", x = "Observed") +
  theme(axis.title.y = element_blank())

trained_plot + tested_plot &
  theme(
    panel.spacing = unit(10, "points"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    strip.text = element_text(size = 7),
    plot.margin = unit(rep(.1, 4), "lines")
  )


## Value of information

skynet_models %>%
  filter(data_subset == "skynet", model == best_ml_model,
         test_set == "spatial", dep_var == "density") %>%
  select(variables, r2, r2_training) %>%
  gather(r2_source, r2_value, -variables) %>%
  ggplot(aes(variables, r2_value, fill = r2_source)) +
  geom_col(position = "dodge")

skynet_models %>%
  filter(data_subset == "skynet", model == best_ml_model,
         test_set == "spatial", dep_var == "density") %>%
  select(variables, rmse, rmse_training) %>%
  gather(rmse_source, rmse_value, -variables) %>%
  ggplot(aes(variables, rmse_value, fill = rmse_source)) +
  geom_col(position = "dodge")



## Supplementary predictions