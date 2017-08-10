quick_map <- function(dat,lat_var, lon_var, plot_var, facet_var){

  dat2 <- dat %>%
    rename(lon = !!lon_var,
           lat = !!lat_var,
           var = !!plot_var,
           facet = !!facet_var)

  qmplot(
    x = lon,
    y = lat,
    data = dat2,
    color = log(var)
  ) +
    scale_color_viridis(guide = F) +
    facet_wrap( ~ facet) +
    theme_classic() +
    labs(caption = paste0(plot_var)[2])

}