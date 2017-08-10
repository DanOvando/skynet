clip_gfw_to_fishdata <- function(knots,gfw_dat) {
u_knots <- knots %>%
  select(approx_long, approx_lat) %>%
  unique()
#
u_gfw <- gfw_dat %>%
  select(rounded_lon, rounded_lat) %>%
  unique()
#
u_knots <-  u_knots %>%
  dplyr::mutate(geometry = purrr::map2(approx_long, approx_lat, ~ sf::st_point(
    x = c(.x, .y), dim = 'XY'
  ))) %>%
  ungroup() %>%
  mutate(geometry = sf::st_sfc(geometry, crs =
                                 "+proj=longlat")) %>%
  sf::st_sf()


u_gfw <-  u_gfw %>%
  dplyr::mutate(geometry = purrr::map2(rounded_lon, rounded_lat, ~ sf::st_point(
    x = c(.x, .y), dim = 'XY'
  ))) %>%
  ungroup() %>%
  mutate(geometry = sf::st_sfc(geometry, crs =
                                 "+proj=longlat")) %>%
  sf::st_sf()


poly_knot <- sf::st_convex_hull(sf::st_union(u_knots))

poly_gfw <- sf::st_convex_hull(u_gfw)

masked_gfw <- sf::st_intersection(poly_knot , poly_gfw)

extract_ll <- function(x) {
  data.frame(rounded_lon = x[1], rounded_lat = x[2])
}

masked_gfw_coords <- as_data_frame(masked_gfw) %>% {
  map_df(.$geometry, extract_ll)
} %>%
  mutate(in_poly = T)

clipped_gfw <- gfw_dat %>%
  left_join(masked_gfw_coords, by = c('rounded_lon', 'rounded_lat')) %>%
  filter(in_poly == T)
}