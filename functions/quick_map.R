#' quickly map skynet data
#'
#' @param dat
#' @param lat_var the variable that latitude lives in
#' @param lon_var the variable that longitude lives in
#' @param plot_var the variable to plot
#' @param facet_var the variable to facet by
#' @param min_lat bounds on lat
#' @param min_lon bounds on lon
#'
#' @return a ggmap object
#' @export
#'
quick_map <- function(dat,lat_var, lon_var, plot_var, facet_var = 'none',
                      min_lat = -90, min_lon = -160){

  if (facet_var == 'none'){

    dat <- dat %>%
      mutate(!!facet_var := '')

  }

  dat2 <- dat %>%
    rename(lon = !!lon_var,
           lat = !!lat_var,
           var = !!plot_var,
           facet = !!facet_var) %>%
    mutate(lon = pmax(lon, min_lon),
           lat = pmax(lat, min_lat))

  qmplot(
    x = lon,
    y = lat,
    data = dat2,
    color = (var)
  ) +
    scale_color_viridis() +
    facet_wrap( ~ facet) +
    theme_classic() +
    labs(caption = paste0(plot_var)[2])

}