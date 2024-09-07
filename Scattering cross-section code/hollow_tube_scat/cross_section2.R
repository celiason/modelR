setwd("~/MEEP/Scattering cross-section code/hollow_tube_scat")

a <- 1000e-9 # scaling factor?
rad <- 0.217/2
sx1 <- rad*4 # size of side of "flux box"

scattering <- read.csv("trans_flux_structure_scat.dat", head=F)
scat_pow <- abs(scattering[,2]-scattering[,3]-scattering[,4]+scattering[,5])
freq <- scattering[,1]
lam <- freq^-1
no_structure <- read.csv("incident_flux.dat", head=F) 
incident_pow <- abs(no_structure[,2])/sx1
absorption <- read.csv("trans_flux_structure_abs.dat", head=F)
abs_pow <- absorption[,2]-absorption[,3]-absorption[,4]+absorption[,5]
scat_cross_section <- scat_pow/incident_pow
abs_cross_section <- abs_pow/incident_pow

plot(x=lam, y=a*scat_cross_section, type='l', lwd=3)

# Compare to exact Mie solution
# mie_exact <- read.csv("/Users/chad/Documents/School/PhD/Projects/Polyplectron/data/2013-03-22/cscat_wavelength.csv", head=F)
# plot(mie_exact$V2~c(200:1000), type='l', ylim=c(0, 7e-7))
# lines(a*scat_cross_section~I(1000*lam), lty=2)

# Save data
# te <- scat_cross_section
# tm <- scat_cross_section
# both <- 0.5*(te+tm)
