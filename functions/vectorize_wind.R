#' vectorize wind
#' calculates wind angle and speed from wind vectors
#' @param x_wind x wind vector
#' @param y_wind y wind vector
#'
#' @return wind angle and speed
#' @export
#'
vectorize_wind <- function(x_wind, y_wind){

  wind_angle <-  180/pi * atan2(y_wind , x_wind)

  wind_speed <- y_wind / sin(atan2(y_wind , x_wind))

  return(list(wind_angle = wind_angle, wind_speed = wind_speed))

}