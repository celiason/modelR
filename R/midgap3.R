#give band below gap of interest. k.ind is index of wavevector location where gap occurs
midgap <- function(tefreqs, tmfreqs, gap, k.ind) {
  gap = gap+6
  bot <- min(tefreqs[,gap+1][k.ind], tmfreqs[,8][k.ind])
  top <- max(tefreqs[,gap][k.ind], tmfreqs[,7][k.ind])
  bot.te <- min(tefreqs[,gap+1][k.ind])
  top.te <- max(tefreqs[,gap][k.ind])
  bot.tm <- min(tmfreqs[,gap+1][k.ind])
  top.tm <- max(tmfreqs[,gap][k.ind])
  delta.w.te <- bot.te-top.te
  delta.w.tm <- bot.tm-top.tm
  delta.mid.te <- (bot.te+top.te)/2
  delta.mid.tm <- (bot.tm+top.tm)/2
  #list(delta.w.te=delta.w.te, delta.w.tm=delta.w.tm, delta.mid.te=delta.mid.te, delta.mid.tm=delta.mid.tm, midgap.te=delta.w.te/delta.mid.te, midgap.tm=delta.w.tm/delta.mid.tm)
  data.frame(pol=c("TE","TM"), bot=c(top.te,top.tm), top=c(bot.te,bot.tm), freq=c(delta.mid.te, delta.mid.tm), delta.w=c(delta.w.te,delta.w.tm), midgap=c(delta.w.te/delta.mid.te, delta.w.tm/delta.mid.tm))
}
