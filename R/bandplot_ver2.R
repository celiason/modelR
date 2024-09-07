# Plot photonic band diagrams from MPB output
# Chad M. Eliason, modified 2012-08-15

bandplot <- function(where, kpoints=3) {
  old <- getwd()
  setwd(where)
  require(RColorBrewer)
  pal <- brewer.pal(7, "Set1")
  tefreqs <- read.csv(list.files(pattern=".te.dat"))
  tmfreqs <- read.csv(list.files(pattern=".tm.dat"))
  lab <- seq(1, dim(tefreqs)[1], by=(dim(tefreqs)[1]-1)/kpoints)
  par(mar=c(4,4,2,2))
  plot(tefreqs[,7], type='l', ylab=expression(paste("Frequency (",omega,"a/",2*pi,"c)")), las=1, xlab="", xaxs='i', xaxt='n', ylim=c(0,1.2))
  axis(side=1, at=lab, labels=c(expression(Gamma), "M", "K", expression(Gamma)))
  lines(tmfreqs[,7], lty=2)
  for(i in 8:dim(tefreqs)[2]) {
  	lines(tefreqs[,i])
  	lines(tmfreqs[,i], lty=2)
  }
  abline(v=lab, lty=1)
  setwd(old)
  list(tefreqs, tmfreqs)
}
