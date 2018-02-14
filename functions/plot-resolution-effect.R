plot_resolution_effect <- function(mod_res) {

res_effect <- mod_res %>%
  mutate(obs = map_dbl(rounded_data, nrow)) %>%
  ggplot(aes(resolution, r2)) +
  geom_point(aes(size = obs), shape = 21, fill = "steelblue") +
  geom_smooth(se = F) +
  scale_y_continuous(limits = c(0, NA),
                     name = expression(R^2)) +
  scale_size_continuous(name = "Observations",range = c(2,7), breaks = seq(0,3000, by = 500)) +
  theme(axis.title.y = element_text(angle = 0)) +
  labs(x = "Lat/Long Resolution", title = "Model: Fish ~ Fishing")

}