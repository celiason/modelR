# Function that runs MEEP from within R and outputs wavelengths and reflectance values 
# Written by Chad M. Eliason, Sept 4, 2018
runMeep <- function(..., file, clean = TRUE, plot = TRUE, img = FALSE, nm=NULL) {
    require(png)
    require(grid)
    if (clean) {
        system('rm output/*.h5')
        system('rm output/*.out')
        system('rm output/*.png')
    }
    expr <- as.list(substitute(list(...)))[-1L]
    class(expr) <- "expression"
    args <- lapply(expr, eval)
    args.name <- names(args)
    args.val <- as.vector(unlist(args))
    if (!is.null(args.val)) {
        run <- paste(args.name, "=", args.val, sep="", collapse=" ")
        run0 <- paste("meep ", run, " ", file, " |tee struc0.out", sep="")
        run1 <- paste("meep ", run, " no-struc?=false ", file, " |tee struc1.out", sep="")
    } else {
        run0 <- paste("meep ", file, " |tee struc0.out", sep="")
        run1 <- paste("meep no-struc?=false ", file, " |tee struc1.out", sep="")
    }
    system(run0)
    system(run1)
    out0 <- system("grep flux1: struc0.out", intern = TRUE)
    out1 <- system("grep flux1: struc1.out", intern = TRUE)
    struc0 <- read.csv(text=out0, head = FALSE)
    struc1 <- read.csv(text=out1, head = FALSE)
    wl <- (1/struc0[,2])*1000
    refl <- -1*(struc1[,4]/struc0[,3])
    trans <- struc1[,3]/struc0[,3]
    par(mfrow=c(1, 2), mgp=c(1.5, .5, 0), ps=8, mex=.75, mar=c(3,3,1,1))
    if (plot) {
        plot(refl~wl, type='l', main=nm, ylim=c(0, max(refl)))
    }
    res <- data.frame(wl, refl, trans, abs=1-(refl+trans))
    if (img) {
        # install.packages("BiocManager")
        # BiocManager::install("rhdf5")
        file <- list.files("output", pattern="eps.*h5", full.name=T)
        x <- rhdf5::h5read(file, name="eps")
        image(x, col=rev(grey.colors(100, start=1, end=0)), xlim=c(1, 0), ylim=c(1, 0), asp=dim(x)[2]/dim(x)[1])
    }
    invisible(res)
}