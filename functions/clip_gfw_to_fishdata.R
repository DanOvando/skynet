#' clip gfw data to the extent of trawl survey knots
#'
#' @param knots a dataframe of knot locations
#' @param gfw_dat a dataframe of GFW locations
#'
#' @return a clipped gfw dataframe
#' @export
#'
clip_gfw_to_fishdata <- function(knots,gfw_dat) {

# get unique locations in knots and gfw
u_knots <- knots %>%
  select(approx_long, approx_lat) %>%
  unique()
#
u_gfw <- gfw_dat %>%
  select(rounded_lon, rounded_lat) %>%
  unique()
# convert to spatial objects
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

# approximate border polygon of both objects
poly_knot <- sf::st_convex_hull(sf::st_union(u_knots))

poly_gfw <- sf::st_convex_hull(u_gfw)

# find intersection of knots and gfw
masked_gfw <- sf::st_intersection(poly_knot , poly_gfw)

extract_ll <- function(x) {
  data.frame(rounded_lon = x[1], rounded_lat = x[2])
}

masked_gfw_coords <- as_data_frame(masked_gfw) %>% {
  map_df(.$geometry, extract_ll)
} %>%
  mutate(in_poly = T)

# filter out gfw position outside of knots polygon
#
clipped_gfw <- gfw_dat %>%
  left_join(masked_gfw_coords, by = c('rounded_lon', 'rounded_lat')) %>%
  filter(in_poly == T)
}