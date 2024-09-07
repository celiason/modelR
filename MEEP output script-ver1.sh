# MEEP output script

cd "/Users/chad/MEEP code/scripts"
meep r=.1015 a=.203 n1=1 n2=2 k2=.6 res=30 poly3_sm_ver2.ctl |tee struc0.out
meep r=.1015 a=.203 n1=1 n2=2 k2=.6 res=30 no-struc?=false poly3_sm_ver2.ctl |tee struc1.out
grep flux1: struc0.out > flux0.dat
grep flux1: struc1.out > flux1.dat

#IMPORT DATA FROM MEEP (LOCATED IN /CHAD FOLDER)

setwd("/Users/chad/MEEP code/scripts")

struc0 <- read.csv("flux0.dat",head=F)
struc1 <- read.csv("flux1.dat",head=F)
wl <- (1/struc0[,2])*1000
refl <- -1*(struc1[,4]/struc0[,3])
#plot(refl~wl, type='l', xlim=c(250,1000))
lines(refl~wl, col=rainbow(10)[5])

#sims <- cbind(wl=wl, k0.0=refl)
sims <- cbind(sims, k0.6=refl)

require(RColorBrewer)
pal <- brewer.pal(6, "Spectral")
plot(smooth.spline(x=sims[,1], y=log(sims[,2]), spar=.6), type='l', xlim=c(250,850), ylim=c(-4,1), col=pal[1], lwd=3)
lines(smooth.spline(x=sims[,1], y=log(sims[,3]), spar=.6), type='l', lwd=3, col=pal[2])
lines(smooth.spline(x=sims[,1], y=log(sims[,4]), spar=.6), type='l', lwd=3, col=pal[3])
lines(smooth.spline(x=sims[,1], y=log(sims[,5]), spar=.6), type='l', lwd=3, col=pal[4])
lines(smooth.spline(x=sims[,1], y=log(sims[,6]), spar=.6), type='l', lwd=3, col=pal[5])
lines(smooth.spline(x=sims[,1], y=log(sims[,7]), spar=.6), type='l', lwd=3, col=pal[6])


# Simulations for effect of matrix RI

#sims2 <- cbind(wl=wl, n1=refl)
sims2 <- cbind(sims2, n2.0=refl)
pal <- brewer.pal(11, "Spectral")
par(mfrow=c(1,2))
plot(sims2[,2]~sims2[,1], lwd=2, col=pal[1], type='l', xlim=c(300,1000), ylim=c(0,.16))
for (i in seq(4,12,by=2)) lines(sims2[,i]~sims2[,1], lwd=2, col=pal[i]/2)
plot(sims3[,2]~sims3[,1], lwd=2, col=pal[1], type='l', xlim=c(300,1000), ylim=c(0,.16))
for (i in 3:7) lines(sims3[,i]~sims3[,1], lwd=2, col=pal[i-1])


# Simulations for effects of melanin extinction coefficient

#sims <- cbind(wl = wl, k0 = refl)
sims <- cbind(sims, k0.1 = refl)
plot(sims[,2]~sims[,1], type='l', xlim=c(300,1000), col=pal[2], lwd=3)
lines(sims[,3]~sims[,1], lwd=3, col=pal[3])
lines(sims[,4]~sims[,1], lwd=3, col=pal[4])
lines(sims[,5]~sims[,1], lwd=3, col=pal[5])
