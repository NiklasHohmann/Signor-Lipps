library(StratPal)

set.seed(435821)


col_ext <- "#0072B2"
col_orig <- "#E69F00"
col_sampling <- "#009E73"

avg_rate <- 2
t_min <- 0
t_max <- 3.2
n_ext <- 1000


maximum_sampling <- 10

rescale_function <- function(f, from, to, avg_rate) {
  total <- integrate(f, lower = from, upper = to, subdivisions = 100000)[[1]]
  f_scaled <- function(x) {
    return(f(x) / total * (to - from) * avg_rate)
  }
  return(f_scaled)
}

ext_sudden <- function(x) {
  y <- rep(0.3, length(x))
  y[x > 2.5 & x < 3] <- 10
  return(y)
}

ext_stepwise <- function(x) {
  y <- rep(0.3, length(x))
  y[x > 2.5 & x < 3] <- 10
  y[x > 1.5 & x < 2] <- 10
  return(y)
}

ext_gradual <- approxfun(
  x = c(t_min, 1, 3, 3.05, t_max),
  y = c(1, 1, 10, 1, 1),
  rule = 2
)

ext_constant <- function(x) {
  return(rep(1, length(x)))
}

pres_constant <- function(x) {
  return(rep(1, length(x)))
}
pres_increasing <- approxfun(x = c(t_min, t_max), y = c(1, 5), rule = 2)
pres_decreasing <- approxfun(x = c(t_min, t_max), y = c(5, 1), rule = 2)


# 1. Sudden              Sudden              Perfect record
pres_p1 <- function(x) {
  y <- rep(1, length(x))
  return(y)
}


# 2. Sudden              Stepwise           Nawrot
pres_p2 <- function(x) {
  y <- rep(1, length(x))
  y[x < 2 & x > 1] <- 5
  y[x > 2.8 & x < 3.2] <- 5
  return(y)
}

# 3. Sudden              Gradual             Signor-Lipps
pres_p3 <- function(x) {
  y <- rep(1, length(x))
  return(y)
}

# 4. Stepwise           Sudden      Holland
pres_p4 <- function(x) {
  y <- rep(1, length(x))
  y[x < 1] <- 10
  return(y)
}

#5. Stepwise           Stepwise           Perfect record
pres_p5 <- function(x) {
  y <- rep(1, length(x))
  y[x < 1] <- 5
  y[x > 2.5] <- 5
  return(y)
}

#6. Stepwise           Gradual             Signor-Lipps
pres_p6 <- function(x) {
  y <- rep(1, length(x))
  return(y)
}

# 7 gradual to sudden
pres_p7 <- function(x) {
  y <- rep(1, length(x))
  y[x > 2] <- 5
  return(y)
}

## 8 gradual to stepwise
pres_p8 <- function(x) {
  y <- rep(1, length(x))
  y[x > 2.5] <- 5
  y[x < 2] <- 5
  return(y)
}

## 9 gradual to gradual
pres_p9 <- approxfun(x = c(t_min, t_max), y = c(1, 10), rule = 2)

## 10 constant to sudden
pres_p10 <- function(x) {
  y <- rep(1, length(x))
  y[x > 2.5] <- 10
  return(y)
}

## 11 constant to stepwise
pres_p11 <- function(x) {
  y <- rep(0, length(x))
  y[x < 1] <- 10
  y[x > 2] <- 10
  return(y)
}

# 12 constant to gradual
pres_p12 <- approxfun(x = c(t_min, t_max), y = c(10, 1), rule = 2)

generate_data <- function(
  f_ext, # extinction rate
  f_prob, # sampling probability, unnormalized
  sampling_prob = 2, # normalization of sampling prob, avg. no of fossils per Myr
  t_min = 0,
  t_max = 3.2,
  n_ext = 1000
) {
  # rescale sampling probability
  f_prob_sampling <- rescale_function(
    f = f_prob,
    from = t_min,
    to = t_max,
    avg_rate = sampling_prob
  )
  f_ext2 <- f_ext
  lower_extension <- 10 / f_ext2(t_min)
  upper_extension <- 10 / f_ext2(t_max)

  # get last occurrences
  lo <- c()
  i <- 1
  while (length(lo) < n_ext) {
    # simulate n_ext extinctions
    ext <- StratPal::p3_var_rate(
      x = f_ext2,
      from = t_min - 10 / f_prob_sampling(t_min),
      to = t_max + 10 / f_prob_sampling(t_min),
      f_max = 50,
      n = 1
    )
    occ <- p3_var_rate(
      f_prob_sampling,
      from = t_min - 10 / f_prob_sampling(t_min),
      to = ext,
      f_max = 50
    )
    occ <- occ[occ > t_min & occ < t_max]
    if (length(occ) > 1) {
      # condition on at least one fossil observed
      lo[i] <- max(occ)
      i <- i + 1
    }
  }

  # data for line plots (sampling rate and extinction rate)
  t <- seq(t_min, t_max, by = 0.01)
  df <- data.frame(t = t, ext = f_ext2(t), prob = f_prob_sampling(t))

  return(list("lo" = lo, "lines" = df))
}


li <- list(
  "A" = c("extinction" = ext_sudden, "sampling" = pres_p1),
  "B" = c("extinction" = ext_sudden, "sampling" = pres_p2),
  "C" = c("extinction" = ext_sudden, "sampling" = pres_p3),
  "D" = c("extinction" = ext_stepwise, "sampling" = pres_p4),
  "E" = c("extinction" = ext_stepwise, "sampling" = pres_p5),
  "F" = c("extinction" = ext_stepwise, "sampling" = pres_p6),
  "G" = c("extinction" = ext_gradual, "sampling" = pres_p7),
  "H" = c("extinction" = ext_gradual, "sampling" = pres_p8),
  "I" = c("extinction" = ext_gradual, "sampling" = pres_p9),
  "J" = c("extinction" = ext_constant, "sampling" = pres_p10),
  "K" = c("extinction" = ext_constant, "sampling" = pres_p11),
  "L" = c("extinction" = ext_constant, "sampling" = pres_p12)
)

df_lo <- data.frame()
df_lines <- data.frame()
for (panel in names(li)) {
  sampling_prob <- 2
  if (panel == "A") {
    sampling_prob <- 8
  }
  a <- generate_data(
    f_ext = li[[panel]]$extinction,
    f_prob = li[[panel]]$sampling,
    sampling_prob = sampling_prob
  )
  df_lo <- rbind(df_lo, data.frame(lo = a$lo, panel = rep(panel, length(a$lo))))
  df_lines <- rbind(
    df_lines,
    data.frame(
      t = a$lines$t,
      ext = a$lines$ext,
      sampling_prob = a$lines$prob,
      panel = rep(panel, length(a$lines$t))
    )
  )
}

extln <- split(df_lines, df_lines$panel)
extbar <- split(df_lo, df_lo$panel)
np <- length(extln)
# plot panels
pdf('figs/figure 3.pdf', width = 4, height = 6)
pannum <- rep(seq(1:(2 * np)), time = rep(c(2, 1), np))
matrix(pannum, 4, 9, byrow = T)
par(oma = c(2, 2, 1, 1), mar = c(1, 1, 0.5, 0))
layout(matrix(pannum, 4, 9, byrow = T))
coltrue1 <- 'forestgreen' # color for true-extinctions line
coltrue2 <- 'green' # color for true-extinctions background line
borderobs <- 'gray40' # border color for observed extinctions
colobs <- 'gray40' # bar color for observed extinctions
colprob <- 'firebrick1' # fill color for sampling probability
borderprob <- 'firebrick3' # border  color for sampling probability
for (i in 1:np) {
  h <- hist(
    extbar[[i]]$lo,
    plot = F,
    breaks = seq(0, 3.2, length.out = 20)
  )$counts
  max(h) / max(extln[[i]]$ext)
  if (i < 10) {
    lnscalex <- max(h) / max(extln[[i]]$ext)
  }
  if (i > 9) {
    lnscalex <- 0.2 * max(h) / max(extln[[i]]$ext)
  }
  lnscaley <- length(h) / max(extln[[i]]$t)
  densplot <- density(extln[[i]]$sampling_prob)
  barplot(h, horiz = T, space = 0, col = colobs, border = borderobs, axes = F)
  points(
    lnscalex * extln[[i]]$ext,
    lnscaley * extln[[i]]$t,
    type = 'l',
    col = coltrue2,
    lwd = 4,
    lend = 2,
    xpd = NA
  )
  points(
    lnscalex * extln[[i]]$ext,
    lnscaley * extln[[i]]$t,
    type = 'l',
    col = coltrue1,
    lwd = 2,
    lend = 2,
    xpd = NA
  )
  mtext(side = 3, line = -0.6, adj = 1.75, cex = 0.8, LETTERS[i], xpd = NA)
  if (i == 1) {
    mtext(
      side = 2,
      'Sudden extinction',
      cex = 0.7,
      line = 0.5,
      col = coltrue1,
      xpd = NA
    )
  }
  if (i == 4) {
    mtext(
      side = 2,
      'Pulsed extinction',
      cex = 0.7,
      line = 0.5,
      col = coltrue1,
      xpd = NA
    )
  }
  if (i == 7) {
    mtext(
      side = 2,
      'Gradual extinction',
      cex = 0.7,
      line = 0.5,
      col = coltrue1,
      xpd = NA
    )
  }
  if (i == 10) {
    mtext(
      side = 2,
      'Constant extinction',
      cex = 0.7,
      line = 0.5,
      col = coltrue1,
      xpd = NA
    )
  }
  if (i == 1) {
    mtext(side = 2, adj = 0.3, line = -4, bquote("" %->% ""), cex = 0.9)
  }
  if (i == 1) {
    mtext(side = 2, adj = 0.3, line = -5, 'TIME', cex = 0.7)
  }
  if (i == 10) {
    mtext(
      side = 1,
      adj = 0.1,
      line = 0.2,
      'Ext. rate',
      cex = 0.7,
      col = coltrue1
    )
  }
  if (i == 10) {
    mtext(
      side = 1,
      adj = 0.1,
      line = -0.4,
      bquote("" %->% ""),
      cex = 0.9,
      col = coltrue1
    )
  }
  plot(
    y = extln[[i]]$t,
    x = extln[[i]]$sampling_prob,
    xlab = '',
    ylab = '',
    type = 'n',
    axes = F,
    xlim = c(0, 8)
  )
  polygon(
    y = c(0, extln[[i]]$t, 3.2),
    x = c(0, extln[[i]]$sampling_prob, 0),
    col = colprob,
    border = borderprob
  )
  if (i == 10) {
    mtext(
      side = 1,
      adj = 0.03,
      line = 0.2,
      'Prob.',
      cex = 0.7,
      col = borderprob
    )
  }
  if (i == 10) {
    mtext(
      side = 1,
      adj = 0.03,
      line = -0.4,
      bquote("" %->% ""),
      cex = 0.9,
      col = borderprob
    )
  }
}
# add boxes
x_min <- 0.95 * grconvertX(0, from = "ndc", to = "user")
x_max <- 0.9 * grconvertX(1, from = "ndc", to = "user")
y_min <- 0.6 * grconvertY(0, from = "ndc", to = "user")
y_max <- 0.98 * grconvertY(1, from = "ndc", to = "user")
x_divide1 <- grconvertX(1.15 / 3, from = "ndc", to = "user")
x_divide2 <- grconvertX(2.05 / 3, from = "ndc", to = "user")
y_divide1 <- grconvertY(1.13 / 4, from = "ndc", to = "user")
y_divide2 <- grconvertY(2.05 / 4, from = "ndc", to = "user")
y_divide3 <- grconvertY(3.01 / 4, from = "ndc", to = "user")
mylwd <- 1
mygridcol <- 'black'
segments(
  x_divide1,
  y_min,
  x_divide1,
  y_max,
  xpd = NA,
  col = mygridcol,
  lwd = mylwd
)
segments(
  x_divide2,
  y_min,
  x_divide2,
  y_max,
  xpd = NA,
  col = mygridcol,
  lwd = mylwd
)
segments(
  x_min,
  y_divide1,
  x_max,
  y_divide1,
  xpd = NA,
  col = mygridcol,
  lwd = mylwd
)
segments(
  x_min,
  y_divide2,
  x_max,
  y_divide2,
  xpd = NA,
  col = mygridcol,
  lwd = mylwd
)
segments(
  x_min,
  y_divide3,
  x_max,
  y_divide3,
  xpd = NA,
  col = mygridcol,
  lwd = mylwd
)
rect(
  xleft = x_min,
  ybottom = y_min,
  xright = x_max,
  ytop = y_max,
  border = "black",
  lwd = 1.5,
  xpd = NA
)
dev.off()
