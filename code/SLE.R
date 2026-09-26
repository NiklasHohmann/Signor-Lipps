load('data/fig_3.Rdata')
extln <- split(df_lines, df_lines$panel)
extbar <- split(df_lo, df_lo$panel)
np <- length(extln)
# plot panels
pdf('figs/figure 3.pdf', width=4, height=6)
pannum <- rep(seq(1:(2*np)), time=rep(c(2,1),np))
matrix(pannum, 4, 9, byrow=T)
par(oma=c(2,2,1,1), mar=c(1,1,0.5,0))
layout(matrix(pannum, 4, 9, byrow=T))
coltrue1 <- 'forestgreen' # color for true-extinctions line
coltrue2 <- 'green' # color for true-extinctions background line
borderobs <- 'gray40' # border color for observed extinctions
colobs <- 'gray40' # bar color for observed extinctions
colprob <- 'firebrick1' # fill color for sampling probability
borderprob <- 'firebrick3' # border  color for sampling probability
for (i in 1:np) {
  h <- hist(extbar[[i]]$lo, plot=F, breaks = seq(0,3.2, length.out = 20))$counts
  max(h)/max(extln[[i]]$ext)
  if (i < 10) lnscalex <- max(h)/max(extln[[i]]$ext)
  if (i > 9) lnscalex <- 0.2*max(h)/max(extln[[i]]$ext)
  lnscaley <- length(h)/max(extln[[i]]$t)
  densplot <- density(extln[[i]]$sampling_prob)
  barplot(h, horiz = T, space = 0, col = colobs, 
          border = borderobs, axes = F)
  points(lnscalex*extln[[i]]$ext, lnscaley*extln[[i]]$t,
         type='l', col = coltrue2,lwd = 4, lend=2, xpd=NA)
  points(lnscalex*extln[[i]]$ext, lnscaley*extln[[i]]$t,
         type='l', col = coltrue1,lwd = 2, lend=2, xpd=NA)
  mtext(side=3, line=-0.6, adj=1.75, cex=0.8, LETTERS[i], xpd=NA)
  if (i == 1) mtext(side=2, 'Sudden extinction', cex=0.7, line=0.5, col=coltrue1, xpd=NA)
  if (i == 4) mtext(side=2, 'Pulsed extinction', cex=0.7, line=0.5, col=coltrue1, xpd=NA)
  if (i == 7) mtext(side=2, 'Gradual extinction', cex=0.7, line=0.5, col=coltrue1, xpd=NA)
  if (i == 10) mtext(side=2, 'Constant extinction', cex=0.7, line=0.5, col=coltrue1, xpd=NA)
  if (i == 1) mtext(side=2, adj=0.3, line=-4, bquote("" %->% ""), cex=0.9)
  if (i == 1) mtext(side=2, adj=0.3, line=-5, 'TIME', cex=0.7)
  if (i == 10) mtext(side=1, adj=0.1, line=0.2, 'Ext. rate', cex=0.7, col=coltrue1)
  if (i == 10) mtext(side=1, adj=0.1, line=-0.4, bquote("" %->% ""), cex=0.9, col=coltrue1)
  plot(y = extln[[i]]$t,  x =  extln[[i]]$sampling_prob, xlab='', ylab='', type='n', axes=F, xlim = c(0,8))
  polygon(y = c(0,extln[[i]]$t, 3.2),   x = c(0,extln[[i]]$sampling_prob,0), col = colprob, border= borderprob)
  if (i == 10) mtext(side=1, adj=0.03, line=0.2, 'Prob.', cex=0.7, col=borderprob)
  if (i == 10) mtext(side=1, adj=0.03, line=-0.4, bquote("" %->% ""), cex=0.9, col=borderprob)
}
# add boxes
x_min <- 0.95*grconvertX(0, from = "ndc", to = "user")
x_max <- 0.9*grconvertX(1, from = "ndc", to = "user")
y_min <- 0.6*grconvertY(0, from = "ndc", to = "user")
y_max <- 0.98*grconvertY(1, from = "ndc", to = "user")
x_divide1 <- grconvertX(1.15/3, from = "ndc", to = "user")
x_divide2 <- grconvertX(2.05/3, from = "ndc", to = "user")
y_divide1 <- grconvertY(1.13/4, from = "ndc", to = "user")
y_divide2 <- grconvertY(2.05/4, from = "ndc", to = "user")
y_divide3 <- grconvertY(3.01/4, from = "ndc", to = "user")
mylwd <- 1
mygridcol <- 'black'
segments(x_divide1, y_min, x_divide1, y_max, xpd=NA, col=mygridcol, lwd=mylwd)
segments(x_divide2, y_min, x_divide2, y_max, xpd=NA, col=mygridcol, lwd=mylwd)
segments(x_min, y_divide1, x_max, y_divide1, xpd=NA, col=mygridcol, lwd=mylwd)
segments(x_min, y_divide2, x_max, y_divide2, xpd=NA, col=mygridcol, lwd=mylwd)
segments(x_min, y_divide3, x_max, y_divide3, xpd=NA, col=mygridcol, lwd=mylwd)
rect(xleft = x_min, ybottom = y_min, xright = x_max, ytop = y_max, 
     border = "black", lwd = 1.5, xpd=NA)
dev.off()




