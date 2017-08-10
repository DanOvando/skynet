#' process environmental covariates
#' takes environmental covariates and spreads them for inclusion in skynet data
#' @param covar the data to be dealt with
#' @param surveys survey to include
#' @param years  year
#' @param non_data_vars variables that aren't data
#'
#' @return mean variable by space and time
#' @export
#'
process_covars <- function(covar,
                           surveys,
                           years = 2012:2016,
                           non_data_vars = c('survey', 'year', 'rlat', 'rlon')) {
  covar <- covar %>%
    filter(survey %in% surveys) %>%
    select(which(colnames(.) != 'plots')) %>% #get rid of any plot columns in there
    select(survey, ncol(.)) %>%
    unnest()

  non_data_vars <-
    c(which(colnames(covar) %in% non_data_vars), which(str_detect(colnames(covar), 'units'))) #locations of variables besides data that you want

  covar <- covar %>%
    gather(variable, value,-non_data_vars)


  if (any(str_detect(colnames(covar), 'year')) == F) {
    # repeat time-invariant variables by year if there is no year variable

    names <- c('year', colnames(covar))

    covar <-
      expand.grid(years, list(covar), stringsAsFactors = F) %>%
      as_data_frame() %>%
      unnest() %>%
      set_names(names)

  }

  covar <- covar %>% #gather all actual data
    group_by(year, rlat, rlon, variable) %>%
    summarise(mean_value = mean(value, na.rm = T)) #calculate mean of variable at year and space


}