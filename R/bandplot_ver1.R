bandplot <- function(name="tri-hollow", interp.k=8, kpoints=4)

{

library(RColorBrewer)
pal <- brewer.pal(7, "RdYlBu")
lab <- seq(1, dim(tefreqs)[1], by=interp.k+1)
lab <- seq(1,dim(tefreqs)[1], by=65)
#READ OUTPUT OF MPB PHOTONIC BANDS PROGRAM
tefreqs <- read.csv(paste("/Users/chad/",name,".te.dat", sep=""))
tmfreqs <- read.csv(paste("/Users/chad/",name,".tm.dat", sep=""))

#PLOT BAND STRUCTURE DIAGRAM
par(mar=c(4,4,2,2))

plot(tefreqs[,7], type='l', ylab=expression(paste("Frequency (",omega,"a/",2*pi,"c)")), las=1, xlab="", xaxs='i', xaxt='n', ylim=c(0,1))

	axis(side=1, at=lab, labels=c(expression(Gamma), "M", "K", expression(Gamma)))
	lines(tmfreqs[,7], lty=2)
	for(i in 8:14){
		lines(tefreqs[,i])
		lines(tmfreqs[,i], lty=2)
	}

	abline(v=lab, lty=1)

}