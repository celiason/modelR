#understanding bandgap diagrams (read Kittel)
library(RColorBrewer)
pal<-brewer.pal(7,"Set1")

par(mfrow=c(2,2),mar=c(2,4,4,2),mex=.75,ps=10)

r1<-10

Y1<-cos(pi*seq(0,.5,.01)/1)^2
Y2<-sin(pi*seq(0,.5,.01)/1)^2
pot<-dexp(seq(0,.5,.01)*pi,rate=r1)

xs<-seq(0,.5,.01)*pi

plot(Y1~xs,type='l',col=pal[1],lwd=3,main="Wave probabilty functions")
lines(Y2~xs,col=pal[1],lwd=3,lty=2)

plot(pot~xs,ylim=c(0,1),type='l',lwd=3,col=pal[3],main="Potential")
points(0,0,pch=16,cex=3)
plot(Y1*pot~xs,col=pal[2],lwd=3,type='l',main="Potential energy of waves")
lines(Y2*pot~xs,col=pal[2],lwd=3,lty=2)