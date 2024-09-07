# Function that runs MEEP from within R and outputs wavelengths and reflectance values 
# Goal: a flexible function that can be updated for different meep scripts
# Written by Chad M. Eliason, Jan 24, 2013

# NOTE: you must supply an argument before 'file=' for this to work correctly

runMPB <- function(..., file='tri-rods.ctl') {
    # path to mpb
    oldwd <- getwd()
    setwd("~/MPB")
    require(stringr)
    # require(plyr)
    name <- str_extract(file, '^.*?(?=.ctl)')
    # system('rm output/*.h5')  # Clean directory
    expr <- as.list(substitute(list(...)))[-1L]
    class(expr) <- "expression"
    args <- lapply(expr, eval)
    args.name <- names(args)
    args.val <- as.vector(unlist(args))
    run <- paste(args.name, "=", args.val, sep="", collapse=" ")
    run1 <- paste("mpb ", run, " ", file, " >& ", name, ".out", sep="")
    system(run1)
    system(paste("grep tmfreqs ", name, ".out > ", name, ".tm.dat", sep=""))
    system(paste("grep tefreqs ", name, ".out > ", name, ".te.dat", sep=""))
    tedat <- read.csv(dir(pattern="*.te.dat"), head = TRUE)
    tmdat <- read.csv(dir(pattern="*.tm.dat"), head = TRUE)
    # system('h5topng output/eps*.h5')
    # system('open output/eps*.png')
    # res <- ldply(list(te=tedat, tm=tmdat))
    # res
    setwd(oldwd)
    res <- list(te=tedat, tm=tmdat)
    res
}
