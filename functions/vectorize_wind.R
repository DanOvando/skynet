vectorize_wind <- function(x_wind, y_wind){
#
#   x_wind <- 10
#
#   y_wind <- -10

  wind_angle <-  180/pi * atan2(y_wind , x_wind)

  wind_speed <- y_wind / sin(atan2(y_wind , x_wind))

  return(list(wind_angle = wind_angle, wind_speed = wind_speed))

}