plot_skynet <- function(skynet_data, plot_theme, height = 8, width = 10, run_dir,...){

  plot_theme <- theme_ipsum(base_size = 14, axis_title_size = 18)

  theme_set(plot_theme)

  load(file = here::here("data","global_map.Rdata"))

  pacific_map <- global_map %>%
    as("Spatial") %>%
    maptools::nowrapRecenter() %>%
    sf::st_as_sf()

# GIF of abundance/effort over time

naive_trends <- skynet_data %>%
  group_by(rounded_lat, rounded_lon) %>%
  summarise(total_abundance = sum(density),
            total_engine_hours = sum(total_engine_hours),
            survey = unique(survey)[1]) %>%
  ungroup() %>%
  mutate(recenter_lon = ifelse(rounded_lon < 0, 180 + (180 - abs(rounded_lon)), rounded_lon)) # %>%
  # filter(survey == "ebsbts")

map_naive <-  naive_trends %>%
  dplyr::mutate(geometry = purrr::map2(recenter_lon, rounded_lat, ~ sf::st_point(x = c(.x, .y), dim = 'XY'))) %>%
  ungroup() %>%
  mutate(geometry = sf::st_sfc(geometry, crs =
                                 "+proj=longlat +datum=WGS84 +no_defs")) %>%
  sf::st_sf()

bbox <- sf::st_bbox(map_naive)


  a = rbind(pacific_map, pacific_map)

naive_fishing_map <- naive_trends %>%
  ggplot() +
  geom_sf(data = pacific_map, fill = 'grey60') +
  geom_raster(aes(x = recenter_lon, y = rounded_lat, fill = log(total_engine_hours))) +
  # geom_sf(aes(color = survey), size = 1, alpha = 0.5) +
  coord_sf(xlim = c(bbox['xmin'], bbox['xmax']),
           ylim = c(bbox['ymin'], bbox['ymax'])) +
  scale_fill_viridis_c(name = "Fishing")

naive_fish_map <- naive_trends %>%
  ggplot() +
  geom_sf(data = pacific_map, fill = 'grey60') +
  geom_raster(aes(x = recenter_lon, y = rounded_lat, fill = log(total_abundance))) +
  # geom_sf(aes(color = survey), size = 1, alpha = 0.5) +
  coord_sf(xlim = c(bbox['xmin'], bbox['xmax']),
           ylim = c(bbox['ymin'], bbox['ymax'])) +
  scale_fill_viridis_c(name = "Abundance")

# gganimate(naive_map, "output.gif")

# Plot showing naive accuracry of predictive model at different aggregation levels

roundfoo <- function(resolution, skynet_data){

  rounded <- skynet_data %>%
    mutate(lat = round(rounded_lat * (1 / resolution)) / (1 / resolution),
           long = round(rounded_lon * (1 / resolution)) / (1 / resolution)) %>%
    group_by(lat, long) %>%
    summarise(total_fishing = sum(total_engine_hours, na.rm = T),
              total_abundance = sum(density))

}

mod_res <- data_frame(resolution =  seq(0.25,5, by = 0.25)) %>%
  mutate(rounded_data = map(resolution, roundfoo, skynet_data = skynet_data))


mod_res <- mod_res %>%
  mutate(r2 = map_dbl(rounded_data, ~{lm(total_abundance ~ total_fishing, data = .x) %>% modelr::rsquare(.x)}))

res_plot <- mod_res %>%
  unnest() %>%
  ggplot(aes(log(total_fishing), log(total_abundance), fill = factor(resolution), frame = resolution)) +
  geom_point(shape = 21, show.legend = F) +
  geom_smooth(method = "lm", show.legend = F) +
  labs(x = 'log(Fishing)', y = "log(Fish)")

gganimate(res_plot, filename = paste0(run_dir,'naive_resolution_gif.gif'))


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

ggsave(paste0(run_dir,'naive_resolution_r2.pdf'), height = height,
       width = width)


# Plot showing predictive accuracy of model at different aggregation levels


}