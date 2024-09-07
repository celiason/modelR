d <- seq(40,300,length=25)
mie <- 2*pi*d

neff <- .906*2+(1-.906)

bragg <- 2*neff*(d/2*sqrt(3))  # for a hexagonal lattice, r/a=0.5

require(RColorBrewer)
pal <- brewer.pal(7, "Set1")

par(lwd=2, mfrow=c(1,1), las=1, mar=c(4,4,2,2))
plot(mie~d, type='l', col=pal[1], xlab="diameter (nm)", ylab="hue (nm)", ylim=c(300,700), xlim=c(60,200))
#for (i in 2:6) lines(mie/i~d, col=pal[1])
lines(bragg~d, col=pal[2], type='l', ylim=c(300,700))
#for (i in 2:6) lines(bragg/i~d, col=pal[2])
#legend("topleft", legend=c("mie", "bragg"), lty=1, col=pal[1:2], bty='n')
text(X, c("m=6.28", "m=3.30"))


lm(bragg~d)
lm(mie~d)
X <- locator()
