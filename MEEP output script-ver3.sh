#hollow hexagonal (or centered rectangular) lattice with ellipsoidal rods
meep ra=.5 a1=.58 a2=.62 k1=0 k2=0.05 theta=0 mel=4 ker=2.43 res=30 nx=5 no-cortex?=true cor=1.05 ellipsoid.ctl |tee struc0.out
meep ra=.5 a1=.58 a2=.62 k1=0 k2=0.05 theta=0 mel=4 ker=2.43 res=30 nx=5 no-cortex?=true cor=1.05 no-struc?=false ellipsoid.ctl |tee struc1.out
grep flux1: struc0.out > flux0.dat
grep flux1: struc1.out > flux1.dat

#hollow hex tilt
meep r=.316 ra=.6 nx=10 TE?=true k1=0 k2=0.01 mel=4 ker=2.43 res=30 no-cortex?=true cor=1.05 tri-rods-hollowTILT.ctl |tee struc0.out
meep r=.316 ra=.6 nx=10 TE?=true k1=0 k2=0.01 mel=4 ker=2.43 res=30 no-cortex?=true cor=1.05 no-struc?=false tri-rods-hollowTILT.ctl |tee struc1.out
grep flux1: struc0.out > flux0.dat
grep flux1: struc1.out > flux1.dat

#hollow hex tilt (new version 04-04-2012)
meep r=.316 ra=.52 TE?=true k1=0 k2=0.01 mel=4 ker=2.43 res=30 cor=0 witu_newtilt.ctl |tee struc0.out
meep r=.316 ra=.52 TE?=true k1=0 k2=0.01 mel=4 ker=2.43 res=30 cor=0 no-struc?=false witu_newtilt.ctl |tee struc1.out
grep flux1: struc0.out > flux0.dat
grep flux1: struc1.out > flux1.dat

#solid hex
meep TE?=true r=.06 a1=.15 a2=.15 k1=0.01 k2=0.01 theta=0 mel=4 ker=2.43 res=30 nx=5 no-cortex?=false cor=.22 tri-rods-solid_incid.ctl |tee struc0.out
meep TE?=true r=.06 a1=.15 a2=.15 k1=0.01 k2=0.01 theta=0 mel=4 ker=2.43 res=30 nx=5 no-cortex?=false cor=.22 no-struc?=false tri-rods-solid_incid.ctl |tee struc1.out
grep flux1: struc0.out > flux0.dat
grep flux1: struc1.out > flux1.dat

#1.95 n_eff for melanin-keratin mix

#hollow hex
meep r=.316 ra=.52 nx=6 theta=0 TE?=true k1=0 k2=0.01 mel=4 ker=2.43 res=50 no-cortex?=true cor=1.048 tri-rods-hollow_incid.ctl |tee struc0.out
meep r=.316 ra=.52 nx=6 theta=0 TE?=true k1=0 k2=0.01 mel=4 ker=2.43 res=50 no-cortex?=true cor=1.048 no-struc?=false tri-rods-hollow_incid.ctl |tee struc1.out
grep flux1: struc0.out > flux0.dat
grep flux1: struc1.out > flux1.dat

#hollow disorder
meep res=30 r=.165 ra=0. rot=0 eps1=2.43 eps2=4 k1=0 k2=0 disorder.ctl |tee struc0.out
meep res=30 r=.165 ra=0. rot=0 eps1=2.43 eps2=4 k1=0 k2=0 no-struc?=false disorder.ctl |tee struc1.out
grep flux1: struc0.out > flux0.dat
grep flux1: struc1.out > flux1.dat

require(RColorBrewer)
pal <- brewer.pal(7, "Set1")

#IMPORT DATA FROM MEEP (LOCATED IN /CHAD FOLDER)
struc0 <- read.csv("flux0.dat", head=F)
struc1 <- read.csv("flux1.dat", head=F)
wl <- (1/struc0[,2])*1000
refl <- -1*(struc1[,4]/struc0[,3])
trans <- struc1[,3]/struc0[,3]
abso <- 1-(refl+trans)
lines(refl~wl, type='l', xlim=c(300,1000), col=pal[1])


plot(loess.smooth(y=sim2.avg, x=wl, span=.15, degree=2, family="gaussian", evaluation=500), type='l', xlim=c(300,2000), ylim=c(0,.7), col=pal[1], lwd=2)
lines(loess.smooth(y=sim.avg, x=wl, span=.15, degree=2, family="gaussian", evaluation=500), col=pal[3], lwd=2)


sim2.te <- refl
sim2.tm <- refl
sim2.avg <- apply(cbind(sim2.te, sim2.tm), 1, mean)


#sim.te <- refl
sim.tm <- refl
sim.avg <- apply(cbind(sim.te, sim.tm), 1, mean)

pdf("turkeyrefl.pdf", width=5, height=10)
par(mfrow=c(2,1), mar=c(4,4,2,2))
plot(loess.smooth(y=sim.avg, x=wl, span=.15, degree=2, family="gaussian", evaluation=500), col=pal[3], lwd=2, type='l', xlim=c(300,1000), ylim=c(0,.8), xlab="Wavelength (nm)", ylab="Reflectance (arb. units)", lty=3)
lines((witu[,2]-20)/max(witu[,2])~witu[,1], type='l', col=pal[3], lwd=1)
abline(v=547*2*1.73/c(1:4), col=pal[3])
plot(loess.smooth(y=sim2.avg/.4, x=wl, span=.15, degree=2, family="gaussian", evaluation=500), type='l', xlim=c(300,1000), ylim=c(0,1), col=pal[1], lwd=2, xlab="Wavelength (nm)", ylab="Reflectance (arb. units)", lty=3)
lines((witu2[,18])/max(witu2[,18])~witu2[,1], type='l', lwd=1, col=pal[1])
abline(v=632/3*2*1.73, col=pal[1])
dev.off()
