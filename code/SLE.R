# make some fake data
x1 <- hist(c(rnorm(100), rnorm(100, 6)), plot=FALSE)$counts
x2 <- c(6,6,6,6,6,6,10,10,10,6,6,6,6)
x3 <- c(rnorm(10,1), rnorm(30,5), rnorm(10, 10), rnorm(50,20), rnorm(10,25))
x4 <- c(10,10,10,6,6,6,2,2,2,6,6,20,20)
x5 <- c(rnorm(100,1), rnorm(30,5), rnorm(10, 10), rnorm(5,20), rnorm(10,25))
# assemble into lists (each represents one panel)
p1 <- list(obs=x1, true=x2, prob=density(x3))
p2 <- list(obs=x1, true=x4, prob=density(x5))
# make a list of 12 panels
x <- list(p1, p2, p1, p2, p1, p2, p1, p2, p1, p2, p1, p2)
# setup the layout and margins
np <- length(x)
pannum <- rep(seq(1:(2*np)), time=rep(c(3,1),np))
par(oma=c(2,2,1,1), mar=c(1,2,0,0))
layout(matrix(pannum, 4, 12, byrow=T))
# plot panels
for (i in 1:np) {
barplot(x[[i]][[1]], horiz = T, space = 0, 
        col = 'black', axes = F)
  points(x[[i]][[2]], 0:(length(x[[i]][[2]])-1),
        type='l', col = 'skyblue2',lwd = 2, lend=2)
  mtext(side=3, line=-1, adj=1.1, cex=0.8, LETTERS[i], xpd=NA)
plot(x[[i]][[3]]$y, x[[i]][[3]]$x, xlab='', 
     ylab='', type='n', axes=F)
polygon(x[[i]][[3]]$y, x[[i]][[3]]$x,
        col='firebrick1', border='firebrick')
}
# add boxes
x_min <- 0.95*grconvertX(0, from = "ndc", to = "user")
x_max <- 0.8*grconvertX(1, from = "ndc", to = "user")
y_min <- 0.6*grconvertY(0, from = "ndc", to = "user")
y_max <- 0.98*grconvertY(1, from = "ndc", to = "user")
x_divide1 <- grconvertX(1.15/3, from = "ndc", to = "user")
x_divide2 <- grconvertX(2.07/3, from = "ndc", to = "user")
y_divide1 <- grconvertY(1.2/4, from = "ndc", to = "user")
y_divide2 <- grconvertY(2.1/4, from = "ndc", to = "user")
y_divide3 <- grconvertY(3.02/4, from = "ndc", to = "user")
mylwd <- 1.5
mygridcol <- 'darkgray'
segments(x_divide1, y_min, x_divide1, y_max, xpd=NA, col=mygridcol, lwd=mylwd)
segments(x_divide2, y_min, x_divide2, y_max, xpd=NA, col=mygridcol, lwd=mylwd)
segments(x_min, y_divide1, x_max, y_divide1, xpd=NA, col=mygridcol, lwd=mylwd)
segments(x_min, y_divide2, x_max, y_divide2, xpd=NA, col=mygridcol, lwd=mylwd)
segments(x_min, y_divide3, x_max, y_divide3, xpd=NA, col=mygridcol, lwd=mylwd)
rect(xleft = x_min, ybottom = y_min, xright = x_max, ytop = y_max, 
     border = "black", lwd = 1.5, xpd=NA)
