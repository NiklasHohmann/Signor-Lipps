library(StratPal)
library(ggplot2)
library(ggpubr)


avg_rate <- 3
t_min <- 0
t_max <- 3

f_ext <- function(x) {
  y <- rep(1, length(x))
  y[x > 2.5 & x < 3] <- 10
  return(y)
}

prob1temp <- function(x) {
  return(rep(1, length(x)))
}

prob2temp <- function(x) {
  return(1 * x + 0.1)
}

prob3temp <- function(x) {
  return(-0.3 * x + 1)
}

prob4temp <- function(x) {
  y <- rep(1, length(x))
  y[x > 1 & x < 2] <- 10
  return(y)
}

prob5temp <- function(x) {
  y <- rep(1, length(x))
  y[x > 1 & x < 2.5] <- 0.1
  return(y)
}

rescale_function <- function(f, from, to, avg_rate) {
  total <- integrate(f, lower = t_min, upper = t_max, subdivisions = 100000)[[
    1
  ]]
  f_scaled <- function(x) {
    return(f(x) / total * (t_max - t_min) * avg_rate)
  }
  return(f_scaled)
}

prob1 <- rescale_function(
  f = prob1temp,
  from = t_min,
  to = t_max,
  avg_rate = avg_rate
)
prob2 <- rescale_function(
  f = prob2temp,
  from = t_min,
  to = t_max,
  avg_rate = avg_rate
)
prob3 <- rescale_function(
  f = prob3temp,
  from = t_min,
  to = t_max,
  avg_rate = avg_rate
)
prob4 <- rescale_function(
  f = prob4temp,
  from = t_min,
  to = t_max,
  avg_rate = avg_rate
)
prob5 <- rescale_function(
  f = prob5temp,
  from = t_min,
  to = t_max,
  avg_rate = avg_rate
)

tt <- seq(0, 3.5, by = 0.1)
plot(tt, f_ext(tt))
plot(tt, prob1(tt))
plot(tt, prob2(tt))
plot(tt, prob3(tt))
plot(tt, prob4(tt))


plot_ext <- function(f_ext, f_prob) {
  n <- 1000
  ext <- StratPal::p3_var_rate(
    f_ext,
    from = t_min,
    to = t_max,
    f_max = 50,
    n = n
  )
  lo <- c()
  i <- 1
  while (length(lo) < n) {
    occ <- p3_var_rate(f_prob, from = t_min, to = ext[i], f_max = 50)
    if (length(occ) < 1) {
      next
    }
    lo[i] <- max(occ)
    i <- i + 1
  }

  t <- seq(t_min, t_max, by = 0.01)
  df_ext <- data.frame(t = t, ext = f_ext(t))
  df_sampling <- data.frame(t = t, prob = f_prob(t))

  df <- data.frame(lo = lo)

  p <- df |>
    ggplot(aes(x = lo)) +
    geom_histogram(aes(y = after_stat(ncount)), binwidth = 0.2) +
    geom_line(
      data = df_ext,
      mapping = aes(x = t, y = ext / max(ext), color = "Extinction rate"),
      inherit.aes = FALSE,
      linewidth = 2
    ) +
    geom_line(
      data = df_sampling,
      mapping = aes(
        x = t,
        y = prob / max(prob) * 0.5,
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
    coord_flip()
}

p1 <- plot_ext(f_ext, prob1)
p2 <- plot_ext(f_ext, prob2)
p3 <- plot_ext(f_ext, prob3)
p4 <- plot_ext(f_ext, prob4)
p5 <- plot_ext(f_ext, prob5)

p1
p2
p3
p4
p5

p <- ggpubr::ggarrange(
  p1,
  p2,
  p3,
  p4,
  p5,
  ncol = 3,
  nrow = 2,
  common.legend = TRUE,
  legend = "bottom"
)
p

ggsave(filename = "figs/variable_sampling_probability.png", plot = p)
