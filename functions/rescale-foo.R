rescale_foo <- function(value, variable, vars_to_average, vars_to_sum){


 if (str_detect(variable,"density|biomass")){

   value = unique(value)

 }

  if (unique(variable) %in% vars_to_average){

    out <- mean(value, na.rm = T)

  } else{

    out <- sum(value, na.rm = T)
  }

  return(out)

}