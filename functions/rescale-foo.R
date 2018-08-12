rescale_foo <- function(value, variable,dilute_total,mean_knot_area, vars_to_average, vars_to_sum){


 # if (str_detect(variable,"density|biomass")){
 #
 #   value = unique(value)
 #
 # }

  if (unique(variable) %in% vars_to_average){

    out <- weighted.mean(value,mean_knot_area, na.rm = T)

  } else{

    if (str_detect(variable,"biomass")){

      out <- sum(unique(value) * unique(dilute_total), na.rm = T)
    } else{

      out <- sum(value, na.rm = T)

    }

  }

  return(out)

}