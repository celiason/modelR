# simulations of graded layer for Ming

pars <- expand.grid(d=seq(100,500,length=10), cor=seq(30,400,length=10))

# silicon as the substrate
sims <- lapply(1:nrow(pars), function(x) {
  HR(d=pars$d[x],cor=pars$cor[x],nend=3.9,n_mat=1.5,lim=c(300, 800),k_mel=0.1)
  })

# keratin as the substrate
sims2 <- lapply(1:nrow(pars), function(x) {
  HR(d=pars$d[x],cor=pars$cor[x],nend=1.5,n_mat=1.5,lim=c(300, 800),k_mel=0.1)
  })


require(pavo)

wl <- sims[[1]]$wavelength
refls <- sapply(sims, "[[", "model.s")
specs <- as.rspec(cbind(wl, refls))

plot(specs)

# colrs <- spec2rgb(specs)

pars$b3 <- apply(specs, 2, max)[-1]
pars$rgb <- spec2rgb(specs)

p1 <- ggplot(data=pars, aes(x=d, y=cor)) + geom_tile(aes(fill=rgb)) + scale_fill_identity()

p2 <- ggplot(data=pars, aes(x=d, y=cor)) + geom_tile(aes(fill=b3))

require(gridExtra)

grid.arrange(p1, p2)
