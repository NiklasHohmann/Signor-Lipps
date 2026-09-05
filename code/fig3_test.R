library(StratPal)
library(ggplot2)
library(ggpubr)
library(ggnewscale)


col_ext <- "#0072B2"
col_orig <- "#E69F00"
col_sampling <- "#009E73"

avg_rate <- 2
t_min <- 0
t_max <- 3.2
n_ext = 1000


maximum_sampling = 10

rescale_function <- function(f, from, to, avg_rate) {
  total <- integrate(f,
                     lower = t_min,
                     upper = t_max,
                     subdivisions = 100000)[[1]]
  f_scaled <- function(x) {
    return(f(x) / total * (t_max - t_min) * avg_rate)
  }
  return(f_scaled)
}

plot_ext <- function(f_ext,
                     # extinction rate
                     f_prob,
                     # sampling probability, unnormalized
                     sampling_prob = 2,
                     # normalization of sampling prob, avg. no of fossils per Myr
                     t_min = t_min,
                     t_max = 3.2,
                     n_ext = 1000) {
  # rescale sampling probability
  f_prob_sampling <- rescale_function(
    f = f_prob,
    from = t_min,
    to = t_max,
    avg_rate = sampling_prob
  )
  
  # simulate n_ext extinctions
  ext <- StratPal::p3_var_rate(
    f_ext,
    from = t_min,
    to = t_max,
    f_max = 50,
    n = n_ext
  )
  
  # get last occurrences
  lo <- c()
  i <- 1
  while (length(lo) < n_ext) {
    occ <- p3_var_rate(f_prob_sampling,
                       from = t_min,
                       to = ext[i],
                       f_max = 50)
    if (length(occ) > 1) {
      # condition on at least one fossil observed
      lo[i] <- max(occ)
    } else {
      lo[i] = NA
    }
    
    i <- i + 1
  }
  lo = lo[!is.na(lo)]
  
  # data for line plots (sampling rate and extinction rate)
  t <- seq(t_min, t_max, by = 0.01)
  df_ext <- data.frame(t = t, ext = f_ext(t))
  df_sampling <- data.frame(t = t, prob = f_prob_sampling(t))
  
  df <- data.frame(lo = lo)
  
  p <- df |>
    ggplot(aes(x = lo)) +
    geom_histogram(aes(y = after_stat(ncount)), binwidth = 0.2) +
    geom_line(
      data = df_ext,
      mapping = aes(
        x = t,
        y = ext / max(ext),
        color = "Extinction rate"
      ),
      inherit.aes = FALSE,
      linewidth = 2
    ) +
    geom_line(
      data = df_sampling,
      mapping = aes(
        x = t,
        y = prob / maximum_sampling,
        color = "Relative sampling probability"
      ),
      inherit.aes = FALSE,
      linewidth = 2
    ) +
    scale_color_manual(
      name = NULL,
      values = c(
        "Extinction rate" = "steelblue",
        "Relative sampling probability" = "firebrick"
      ),
      breaks = c("Extinction rate", "Relative sampling probability")
    ) +
    theme(
      legend.position.inside = c(0.5, 0.98),
      legend.position = "inside",
      legend.justification = c(1, 1),
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = 8),
      legend.margin = margin(2, 2, 2, 2),
      legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
    ) +
    labs(x = "Time [Myr]", y = "LO count") +
    scale_y_continuous(
      name = "Count",
      sec.axis = sec_axis( ~ . *  maximum_sampling, name = "Sampling probability [Fossils/Myr]")
    ) +
    coord_flip() +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.text = element_blank()
    )
  return(p)
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

ext_gradual = approxfun(
  x = c(t_min, 1, 3, 3.05, t_max),
  y = c(1, 1, 10, 1, 1),
  rule = 2
)

ext_constant = function(x) {
  return(rep(1, length(x)))
}

pres_constant = function(x) {
  return(rep(1, length(x)))
}
pres_increasing = approxfun(x = c(t_min, t_max),
                            y = c(1, 5),
                            rule = 2)
pres_decreasing = approxfun(x = c(t_min, t_max),
                            y = c(5, 1),
                            rule = 2)



# 1. Sudden              Sudden              Perfect record
pres_p1 = function(x) {
  y = rep(1, length(x))
  return(y)
}
p1 = plot_ext(
  f_ext = ext_sudden,
  f_prob = pres_p1,
  sampling_prob = 8,
  t_min = t_min,
  t_max = t_max,
  n_ext = n_ext
)
p1

# 2. Sudden              Stepwise           Nawrot
pres_p2 = function(x) {
  y = rep(1, length(x))
  y[x < 1.5 & x > 1] = 5
  return(y)
}
p2 = plot_ext(
  f_ext = ext_sudden,
  f_prob = pres_p2,
  sampling_prob = 2,
  t_min = t_min,
  t_max = t_max,
  n_ext = n_ext
)
p2

# 3. Sudden              Gradual             Signor-Lipps
pres_p3 = function(x) {
  y = rep(1, length(x))
  return(y)
}
p3 = plot_ext(
  f_ext = ext_sudden,
  f_prob = pres_p3,
  sampling_prob = 2,
  t_min = t_min,
  t_max = t_max,
  n_ext = n_ext
)
p3

# 4. Stepwise           Sudden      Holland
pres_p4 = function(x) {
  y = rep(1, length(x))
  y[x < 1] = 10
  return(y)
}
p4 = plot_ext(
  f_ext = ext_stepwise,
  f_prob = pres_p4,
  sampling_prob = 2,
  t_min = t_min,
  t_max = t_max,
  n_ext = n_ext
)
p4

#5. Stepwise           Stepwise           Perfect record
pres_p5 = function(x) {
  y = rep(1, length(x))
  y[x < 1] = 5
  y[x > 2.5] = 5
  return(y)
}
p5 = plot_ext(
  f_ext = ext_stepwise,
  f_prob = pres_p5,
  sampling_prob = 2,
  t_min = t_min,
  t_max = t_max,
  n_ext = n_ext
)
p5

#6. Stepwise           Gradual             Signor-Lipps
pres_p6 = function(x) {
  y = rep(1, length(x))
  return(y)
}
p6 = plot_ext(
  f_ext = ext_stepwise,
  f_prob = pres_p6,
  sampling_prob = 1,
  t_min = t_min,
  t_max = t_max,
  n_ext = n_ext
)
p6

# 7 gradual to sudden
pres_p7 = function(x) {
  y = rep(1, length(x))
  y[x > 2] = 5
  return(y)
}
p7 = plot_ext(
  f_ext = ext_gradual,
  f_prob = pres_p7,
  sampling_prob = 2,
  t_min = t_min,
  t_max = t_max,
  n_ext = n_ext
)
p7

## 8 gradual to stepwise
pres_p8 = function(x) {
  y = rep(1, length(x))
  y[x > 2.5] = 5
  y[x < 2] = 5
  return(y)
}
p8 = plot_ext(
  f_ext = ext_gradual,
  f_prob = pres_p8,
  sampling_prob = 4,
  t_min = t_min,
  t_max = t_max,
  n_ext = n_ext
)
p8

## 9 gradual to gradual
pres_p9 = approxfun(x = c(t_min, t_max),
                    y = c(1, 10),
                    rule = 2)
p9 = plot_ext(
  f_ext = ext_gradual,
  f_prob = pres_p9,
  sampling_prob = 2,
  t_min = t_min,
  t_max = t_max,
  n_ext = n_ext
)
p9

## 10 constant to sudden
pres_p10 = function(x) {
  y = rep(1, length(x))
  y[x > 2.5] = 10
  return(y)
}
p10 = plot_ext(
  f_ext = ext_constant,
  f_prob = pres_p10,
  sampling_prob = 2,
  t_min = t_min,
  t_max = t_max,
  n_ext = n_ext
)
p10

## 11 constant to stepwise
pres_p11 = function(x) {
  y = rep(0, length(x))
  y[x < 1] = 5
  y[x > 2] = 10
  return(y)
}
p11 = plot_ext(
  f_ext = ext_constant,
  f_prob = pres_p11,
  sampling_prob = 2,
  t_min = t_min,
  t_max = t_max,
  n_ext = n_ext
)
p11


# 12 constant to gradual
pres_p12 = approxfun(x = c(t_min, t_max),
                     y = c(10, 1),
                     rule = 2)
p12 = plot_ext(
  f_ext = ext_constant,
  f_prob = pres_p12,
  sampling_prob = 1,
  t_min = t_min,
  t_max = t_max,
  n_ext = n_ext
)
p12


p = ggpubr::ggarrange(p1,
                      p2,
                      p3,
                      p4,
                      p5,
                      p6,
                      p7,
                      p8,
                      p9,
                      p10,
                      p11,
                      p12,
                      nrow = 4,
                      ncol = 3)
p
