# run MIE scattering code

cd "/Users/chad/MEEP code/scripts"
meep theta=45 res=30 mel=4 ker=2.43 k1=0 k2=.05 TE?=true r1=.103 r2=0 mie-hollow.ctl |tee struc0.out
meep theta=45 res=30 mel=4 ker=2.43 k1=0 k2=.05 TE?=true r1=.103 r2=0 no-struc?=false mie-hollow.ctl |tee struc1.out
grep flux1: struc0.out > flux0.dat
grep flux1: struc1.out > flux1.dat

# IMPORT DATA FROM MEEP (LOCATED IN /CHAD FOLDER)

#setwd("/Users/chad/MEEP code/scripts")
struc0 <- read.csv("flux0.dat",head=F)
struc1 <- read.csv("flux1.dat",head=F)
wl <- (1/struc0[,2])*1000
refl <- -1*(struc1[,4]/struc0[,3])
#plot(refl~wl, type='l', xlim=c(300,1000), lwd=2)
lines(refl~wl, col="red")

mie <- refl