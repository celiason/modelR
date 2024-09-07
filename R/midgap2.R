midgap <- function(tefreqs, tmfreqs, gap) #give band below gap of interest
{

gap = gap+6
bot <- min(tefreqs[,gap+1][18], tmfreqs[,8][18])
top <- max(tefreqs[,gap][18], tmfreqs[,7][18])

bot.te <- min(tefreqs[,gap+1][18])
top.te <- max(tefreqs[,gap][18])

bot.tm <- min(tmfreqs[,gap+1][18])
top.tm <- max(tmfreqs[,gap][18])

delta.w.te <- bot.te-top.te
delta.w.tm <- bot.tm-top.tm

delta.mid.te <- (bot.te+top.te)/2
delta.mid.tm <- (bot.tm+top.tm)/2

list(delta.w.te=delta.w.te, delta.w.tm=delta.w.tm, delta.mid.te=delta.mid.te, delta.mid.tm=delta.mid.tm, midgap.te=delta.w.te/delta.mid.te, midgap.tm=delta.w.tm/delta.mid.tm)

}