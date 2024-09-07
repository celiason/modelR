#' Function to determine appropriate center frequencies for simulations
#' Angle varies with frequency. Thus, angles within crystal aren't same as external (specified)
#' range of angles
#' @param thetas angles (in degrees) over which to run Meep simulations
#' @param fcen center frequency of simulations
#' @param df frequency bandwidth of simulations
#'
#' returns wavelengths
calc_fcen <- function(thetas = seq(1,90,by=10), fcen = 2.4, df = 4) {
  theta_max <- max(thetas)/180*pi
  lambda_max <- 1/(fcen-df/2)
  lambda_min <- 1/(fcen+df/2)
  thetas <- thetas/180*pi
  theta_corner <- asin(sin(theta_max*(1/fcen)/lambda_max))
  theta_cen <- thetas
  ifelse(theta_corner <= theta_cen, sin(theta_max) * lambda_min * (2 / sin(theta_cen) - 1), lambda_max)  # lambda_max values
}


# 1 / calc_fcen(thetas=seq(0, 90, by=10), fcen=3.4, df = 4)

