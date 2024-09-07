# change working directory

cd "/Users/chad/MEEP code/src"  # for meep

# hollow hexagonal or centered rectangular lattice with ellipsoidal rods; ra = amount of air, a1 = lattice constant 1, a2 = lattice 2
meep ra=.57 a1=.535 a2=.614 k1=0.01 k2=0.01 theta=0 mel=4 ker=2.43 res=50 nx=6 no-cortex?=true cor=1.05 ellipsoid.ctl |tee struc0.out
meep ra=.57 a1=.535 a2=.614 k1=0.01 k2=0.01 theta=0 mel=4 ker=2.43 res=50 nx=6 no-cortex?=true cor=1.05 no-struc?=false ellipsoid.ctl |tee struc1.out
grep flux1: struc0.out > flux0.dat
grep flux1: struc1.out > flux1.dat

# hollow hexagonal, with tilt
meep r=.15 r2=.075 nx=10 TE?=false k1=0 k2=0.01 mel=4 ker=2.43 res=30 no-cortex?=true cor=1.05 tri-rods-hollowTILT.ctl |tee struc0.out
meep r=.15 r2=.075 nx=10 TE?=false k1=0 k2=0.01 mel=4 ker=2.43 res=30 no-cortex?=true cor=1.05 no-struc?=false tri-rods-hollowTILT.ctl |tee struc1.out
grep flux1: struc0.out > flux0.dat
grep flux1: struc1.out > flux1.dat

# hollow hexagonal, with tilt (new version 2012-08-11)
# 534,594,0.57; 552,596,0.56
meep ra1=.5 ra=.57 a1=.535 a2=.614 TE?=true k1=0.01 k2=0.1 mel=4 ker=2.43 res=30 tri_hollow_tilt_ellipsoids.ctl |tee struc0.out
meep ra1=.5 ra=.57 a1=.535 a2=.614 TE?=true k1=0.01 k2=0.1 mel=4 ker=2.43 res=30 no-struc?=false tri_hollow_tilt_ellipsoids.ctl |tee struc1.out
grep flux1: struc0.out > flux0.dat
grep flux1: struc1.out > flux1.dat

# solid hexagonal lattice
meep a=.2 ra=.5 k1=0 k2=0.1 theta=0 mel=4 ker=4 res=40 nx=10 no-cortex?=true cor=1.05 tri-rods-solid_incid.ctl |tee struc0.out
meep a=.2 ra=.5 k1=0 k2=0.1 theta=0 mel=4 ker=4 res=40 nx=10 no-cortex?=true cor=1.05 no-struc?=false tri-rods-solid_incid.ctl |tee struc1.out
grep flux1: struc0.out > flux0.dat
grep flux1: struc1.out > flux1.dat

# hollow hexagonal lattice
meep r=.1 ra=0.5 a1=.2 a2=.2 nx=10 TE?=true k1=0.01 k2=0.01 mel=4 ker=2.43 res=50 no-cortex?=true cor=1.5 tri-rods-hollow_incid.ctl |tee struc0.out
meep r=.1 ra=0.5 a1=.2 a2=.2 nx=10 TE?=true k1=0.01 k2=0.01 mel=4 ker=2.43 res=50 no-cortex?=true cor=1.5 no-struc?=false tri-rods-hollow_incid.ctl |tee struc1.out
grep flux1: struc0.out > flux0.dat
grep flux1: struc1.out > flux1.dat


# sims to run:
# a = 200
# d = 100, 125, 150, 175, 200
# r/a = 0 (solid), 0.5 (air)

# hollow multilayer
meep d=.15 ra=0.2 a=.2 nx=10 k1=0.01 k2=0.01 mel=4 ker=2.43 res=50 no-cortex?=true cor=1.5 multilayer.ctl |tee struc0.out
meep d=.15 ra=0.2 a=.2 nx=10 k1=0.01 k2=0.01 mel=4 ker=2.43 res=50 no-cortex?=true cor=1.5 no-struc?=false multilayer.ctl |tee struc1.out
grep flux1: struc0.out > flux0.dat
grep flux1: struc1.out > flux1.dat


#IMPORT DATA FROM MEEP (LOCATED IN /CHAD FOLDER)
setwd("/Users/chad/MEEP\ code/src")  # for R
require(RColorBrewer)
pal <- brewer.pal(11, "Spectral")

struc0 <- read.csv("flux0.dat", head=F)
struc1 <- read.csv("flux1.dat", head=F)
wl <- (1/struc0[,2])*1000
refl <- -1*(struc1[,4]/struc0[,3])
trans <- struc1[,3]/struc0[,3]
#plot(refl~wl, type='l', xlim=c(300, 1000), col=pal[1])
#lines(refl~wl, col="blue")
lines(refl~wl, col="green")

#res <- cbind(wl, r100s = refl)
res <- cbind(res, r200h = refl)


#

setwd("/Users/chad/Documents/School/PhD/Projects/Hollow\ PC")
source("/Users/chad/R/src/Color/gather.spectra-AVANTES_ver2.R")
#witu <- gather.spectra("spectra/2012-04-02", wl=c(300:1000))
#cile <- gather.spectra("spectra/2012-04-06/pol backscat", wl=c(300:1000))
specs <- gather.spectra("/Users/chad/Documents/School/PhD/Projects/Hollow PC/spectra/2012-04-10", wl=c(300:1000))
sims <- read.csv("/Users/chad/Desktop/hollow sims nocortex final.csv", row.names=1)
par(mfrow=c(1,2), mar=c(4,4,2,2))
plot(apply(sims[,2:4], 1, mean)~wl, type='l', xlim=c(300,1000), col="orange3", ylim=c(0,.8))
lines((specs[,4]-25)/250~witu[,1])  # TE refl for turkey
plot(apply(sims[,5:7], 1, mean)~wl, type='l', xlim=c(300,1000), col="orange3")
lines(specs[,3]/105~cile[,1])  # is this TE?
