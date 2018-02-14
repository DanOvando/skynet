round_enviro <- function(tidy_dat,
                         desired_data,
                         runit = 0.25){
  # tidy_dat <- tidy_dat %>%
  #   unnest()

  agg_foo <- function(tidy_dat, desired_data, runit){

tidy_dat <- tidy_dat %>%
  mutate(
    rlat = round(rlat * (1 / runit)) / (1 / runit),
    rlon = round(rlon * (1 / runit)) / (1 / runit)
  )

# var <- last(dat$table$columnNames)


var_name <- colnames(tidy_dat)[!str_detect(colnames(tidy_dat),'units') & str_detect(colnames(tidy_dat),'mean')]

var_name_units <- colnames(tidy_dat)[str_detect(colnames(tidy_dat),'units')]

# var_unit_name <- last(colnames(tidy_dat))

group_vars <- c(quo(year), quo(rlat), quo(rlon))

if (desired_data == 'topo') {
  group_vars <- c(quo(rlat), quo(rlon))

}

if (desired_data == 'wind') {

  tidy_dat <-  tidy_dat %>%
    gather(variable, value, wind_speed, wind_angle)

  group_vars <-
    c(quo(year), quo(rlat), quo(rlon), quo(variable), quo(mean_wind_speed_units))

  agg_dat <- tidy_dat %>%
    group_by(!!!group_vars) %>%
    summarise(value = mean(value, na.rm = T)) %>%
    spread(variable, value) %>%
    ungroup()


} else{

  agg_dat <- tidy_dat %>%
    group_by(!!!group_vars) %>%
    summarise(!!var_name := mean(!!!rlang::parse_exprs(var_name), na.rm = T),
              !!var_name_units := unique(!!!rlang::parse_exprs(var_name_units))[1]) %>%
    ungroup()

}

} # close function

  colnames(tidy_dat) %in% c()

  desired_col <- last(colnames(tidy_dat)[colnames(tidy_dat) != "plots"])

  temp <- map(last(tidy_dat[, desired_col]), agg_foo, desired_data = desired_data, runit = runit)

  tidy_dat[,desired_col] <- list(temp)

  return(tidy_dat)

}
