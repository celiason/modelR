# GENERATE PHOTONIC BANDGAP PLOTS FROM OUTPUT OF MPB
# MAKE SURE WORKING DIRECTORY CONTAINS *.te.dat and *.tm.dat files
calcbandgap <- function(...) {

  files <- list.files()

  te <- grep("te", files, value=TRUE)
  tm <- grep("tm", files, value=TRUE)

  tefreqs <- read.csv(te)
  tmfreqs <- read.csv(tm)

  plot(tefreqs[,7], ylim=c(0,1), type='l', ylab=expression(paste("Frequency (",omega,"a/",2*pi,"c)")),las=1, xlab="", xaxt='n')

  N <- floor(dim(tefreqs)[1]/3)

  axis(side=1,at=c(1, 1+N, 1+2*N, dim(tefreqs)[1]), labels=c(expression(Gamma), "M", "K", expression(Gamma)))

  lines(tmfreqs[,7], lty=2)

  for(i in 8:14) {
    lines(tefreqs[,i])
    lines(tmfreqs[,i], lty=2)
  }

  legend("bottom", c("TM","TE"), lty=c(2,1), bty='n')

  res <- list(te=tefreqs, tm=tmfreqs)
  
  invisible(res)

}