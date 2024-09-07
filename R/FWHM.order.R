FWHM.order <- function(ft, rnorm)
{
	rnorm = rnorm
	y2 <- ft
	Yi=max(y2)
	Yj=min(y2)
	Xi=which(y2==Yi)
	
	fsthalf=y2[1:Xi]
	sndhalf=y2[Xi:length(y2)]
	
	halfmax=(Yi+Yj)/2
	
	fstHM=which.min(abs(fsthalf-halfmax))
	sndHM=which.min(abs(sndhalf-halfmax))
	
	Xa=rnorm[fstHM]
	Xb=rnorm[Xi+sndHM]

	plot(y2~rnorm, type='l')
	abline(v=rnorm[y2==Yi], col="red")
	abline(h=halfmax, col="red")
	abline(v=Xa, col="red", lty=2)
	abline(v=Xb, col="red", lty=2)
	return(Xb-Xa)
}