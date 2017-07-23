query_erddap <- function(desired_data = 'sst',
                         runit = 0.25,
                         min_year,
                         max_year,
                         date_interval = 14,
                         min_lat,
                         max_lat,
                         min_lon,
                         max_lon,
                         space_interval = .25) {

  year_query <- paste0('[(',min_year,'-01-01):',date_interval,':(',max_year,'-12-31)]')

  lat_query <- paste0('[(',min_lat,'):',100*space_interval,':(', max_lat,')]')

  lon_query <- paste0('[(',min_lon,'):',100*space_interval,':(', max_lon,')]')

  query <- paste0(year_query,lat_query,lon_query)

  if (desired_data == 'sst') {
    dat <-
      jsonlite::fromJSON(
        paste0('https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41.json?analysed_sst',query),
        flatten = T
      )
      # paste0('https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41.json?analysed_sst[(2012-01-01):14:(2013-12-29)][(53.00):25:(66.00)][(-179.99):25:(-157.00)]',

  }
  if (desired_data == 'chl-a'){

    # durl <- "http://upwell.pfeg.noaa.gov/erddap/griddap/erdMH1chlamday.json?chlorophyll[(2010-06-01):1:(2010-10-01)][(66):1:(53)][(-180):1:(-157)]"

    dat <-
      jsonlite::fromJSON(
        paste0('http://upwell.pfeg.noaa.gov/erddap/griddap/erdMH1chlamday.json?chlorophyll',query)

,        flatten = T
      )



  }
  datnames <- dat$table$columnNames

  date_cols <- which(dat$table$columnUnits == 'UTC')

  other_cols <- which(dat$table$columnUnits != 'UTC')


  tidy_dat <- dat$table$rows %>%
    as_data_frame() %>%
    set_names(datnames) %>%
    map_at(date_cols, lubridate::as_date) %>%
    as_data_frame() %>%
    map_at(other_cols, as.numeric) %>%
    as_data_frame()

  tidy_dat <- tidy_dat %>%
    mutate(
      year = lubridate::year(time),
      month = lubridate::month(time),
      rlat = round(latitude * (1 / runit)) / (1 / runit),
      rlon = round(longitude * (1 / runit)) / (1 / runit)
    )

  var <- last(dat$table$columnNames)

  var_name <- paste0('mean_', var)

  var_name_units <- paste0(var_name, '_units')

  agg_dat <- tidy_dat %>%
    group_by(year, rlat, rlon) %>%
    summarise(!!var_name := mean(!!!rlang::parse_exprs(var), na.rm = T)) %>%
    mutate(!!var_name_units := last(dat$table$columnUnits))

  return(agg_dat)

}