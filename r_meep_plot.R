#IMPORT DATA FROM MEEP (LOCATED IN /CHAD FOLDER)
struc0 <- read.csv("flux0.dat",head=F)
struc1 <- read.csv("flux1.dat",head=F)
wl <- (1/struc0[,2])*1000
refl <- -1*(struc1[,4]/struc0[,3])

plot(refl~wl, type='l', xlim=c(300,2000))
lines(refl~wl, col="blue")

