#' Meep angle-resolved spectra interpolation
#'
#' This function reads in MEEP simulation data (data frame Y, with columns for angle of incidence and rows for different wavelengths). Make sure the wavelengths are in increasing order (use Y[order(Y$wl),] if necessary). 'angles0' specifies which angles simulations were done at. 
#'
#' @param Y a list of simulation data, with wavelengths and reflectance values in columns.
#' @param n number of data points to interpolate over
#' @param theta_0 angles over which simulations were performed
#'

MEEPinterp <- function(Y, fcen, theta_0 = seq(0, 90, by=10), n=100) {
  require(akima)
  ks <- sapply(1:length(fcen), function(x) sin(theta_0[x]/180*pi)*fcen[x])
  wl <- do.call(cbind, lapply(Y, "[[", "wl"))/1000
  refl <- do.call(cbind, lapply(Y, "[[", "refl"))
  angles1 <- sapply(1:length(ks), function(x) asin(ks[x]/(1/wl[,x])))
  angles1[is.na(angles1)] <- 0  # remove NaN values
  int <- interp(x=wl, y=angles1, z=refl, duplicate="mean", xo=seq(min(wl), max(wl), length=n), yo=seq(min(angles1), max(angles1), length=n))
  setNames(int, c("wl", "angle", "refl"))
}
