

woah <- skynet_models %>%
  filter(dep_var == "log_density",
         test_set == "west_coast",
         train_set == "not_west_coast",
         model == "ranger",
         data_subset == "skynet")

woah$resolution_test[[1]]

b <- woah$resolution_test[[1]] %>%
  slice(1:2) %>%
  unnest()


b %>%
  ggplot(aes(long, lat, fill = total_predicted)) +
  geom_raster() +
  scale_fill_viridis() +
  facet_grid(year ~ resolution)

b %>%
  filter(survey == 'wcghl') %>%
  ggplot(aes(total_observed, total_predicted, fill = survey)) +
  geom_point(shape = 21) +
  scale_fill_viridis_d() +
  facet_wrap( ~ resolution, scales = "free_y")

