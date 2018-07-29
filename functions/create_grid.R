create_grid <-
  function(data,
           resolution = 100,
           lon_name,
           lat_name,
           units = "km",
           recenter = T) {
    lon_name <- enquo(lon_name)

    lat_name <- enquo(lat_name)

    map_data <-  data %>%
      dplyr::mutate(geometry = purrr::map2(!!lon_name,!!lat_name, ~ sf::st_point(x = c(.x, .y), dim = 'XY'))) %>%
      ungroup() %>%
      mutate(geometry = sf::st_sfc(geometry, crs =
                                     "+proj=longlat +datum=WGS84 +no_defs")) %>%
      sf::st_sf() %>%
      st_transform(crs = 2163) # convert to meters resolution https://epsg.io/2163

    new_grid <-
      st_make_grid(map_data, cellsize = resolution * 1000, what = "centers") %>%
      sf::st_transform(crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
      sf::st_sf()

    new_grid_coords <-  new_grid %>%
      sf::st_coordinates() %>%
      as_data_frame() %>%
      rename(lon = X, lat = Y) %>%
      mutate(recenter_lon = ifelse(lon < 0, 180 + (180 - abs(lon)), lon)) %>%
      mutate(id = 1:nrow(.))

    # check <-  new_grid_coords %>%
    #   dplyr::mutate(geometry = purrr::map2(recenter_lon, lat, ~ sf::st_point(x = c(.x, .y), dim = 'XY'))) %>%
    #   ungroup() %>%
    #   mutate(geometry = sf::st_sfc(geometry, crs =
    #                                  "+proj=longlat +datum=WGS84 +no_defs")) %>%
    #   sf::st_sf()
    #
    #
    # bbox <- sf::st_bbox(check)
    #
    # ggplot() +
    #   geom_sf(data = check) +
    #   geom_sf(data = pacific_map, fill = 'grey60') +
    #   coord_sf(xlim = c(bbox['xmin'], bbox['xmax']),
    #            ylim = c(bbox['ymin'], bbox['ymax']))


  }