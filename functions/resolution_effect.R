resolution_effect <- function(data, dep_var){


roundfoo <- function(resolution, dep_var, data){


  rounded <- data %>%
    mutate(lat = round(rounded_lat * (1 / resolution)) / (1 / resolution),
           long = round(rounded_lon * (1 / resolution)) / (1 / resolution)) %>%
    group_by(lat, long, year) %>%
    summarise(nobs = length(pred),
      total_observed = sum(!!rlang::sym(dep_var), na.rm = T),
              total_predicted = sum(pred))

}

mod_res <- data_frame(resolution =  seq(0.25,5, by = 0.25)) %>%
  mutate(rounded_data = map(resolution, roundfoo, data = data,
                             dep_var = dep_var))

mod_res <- mod_res %>%
  mutate(r2 = map_dbl(rounded_data, ~{lm(total_observed ~ total_predicted, data = .x) %>% modelr::rsquare(.x)}))


return(mod_res)

}
