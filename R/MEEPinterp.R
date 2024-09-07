#This function reads in MEEP simulation data (data frame Y, with columns for angle of incidence and rows for different wavelengths). Make sure the wavelengths are in increasing order (use Y[order(Y$wl),] if necessary). 'angles0' specifies which angles simulations were done at.

MEEPinterp <- function(Y, wl, k0=3.4, angles0=seq(0, 90, by=10)) {

ks <- sin(angles0/180*pi)*k0
angles1 <- matrix(data=NA, nrow=length(wl), ncol=length(angles0))

for (i in 1:dim(Y)[2]) angles1[,i] <- asin(ks[i]/(1/(wl/1000)))/pi*180

angles1[is.na(angles1)] <- 0 # this gets rid of NaN values for angle calculations due to the fact that some values will be infinite

sims.interp <- matrix(data=NA, nrow=length(wl), ncol=91)

# this interpolates at numerous angles:
for (i in 1:length(wl)){
  sims.interp[i, ] <- approx(x=angles1[i, ], y=Y[i, ], xout=seq(0,90,by=1))$y
}

sims.interp[is.na(sims.interp)] <- 0
return (list(wl = wl, angle=seq(0,90,by=1), sims=sims.interp))

}
