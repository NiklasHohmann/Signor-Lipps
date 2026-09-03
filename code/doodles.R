library(StratPal)
library(ggplot2)
library(ggpubr)
library(ggnewscale)

col_ext <- "steelblue"
col_orig <- "darkorange"
col_sampling <- "darkorchid4"


col_ext <- "#0072B2"
col_orig <- "#E69F00"
col_sampling <- "#009E73"

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
  return(pmax(c(-0.3 * x + 1), 0.1))
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


line_plot_sampling <- function(
  f_sampling,
  f_origin,
  f_ext,
  n_taxa,
  t_range,
  names_tax = LETTERS[1:n_taxa],
  sort_by = c("observed extinction", "true extinction"),
  allow_singletons = TRUE
) {
  ifelse(allow_singletons, min_occ <- 1, min_occ <- 2)
  t <- t_range
  stopifnot(n_taxa == length(names_tax))
  ext <- c()
  orig <- c()
  li <- list()
  i <- 1
  while (length(li) < n_taxa) {
    # simulate origination and extinction
    origin <- StratPal::p3_var_rate(
      x = f_orig,
      from = min(t),
      to = max(t),
      n = 1,
      f_max = 50
    )
    e <- StratPal::p3_var_rate(
      f_ext,
      from = origin,
      to = max(t),
      n = 1,
      f_max = 50
    )
    # simulate fossil samples
    x <- StratPal::p3_var_rate(
      x = f_sampling, # sampling frequency
      from = origin,
      to = e,
      f_max = 50
    )
    if (length(x) >= min_occ) {
      # enforce at least 2 fossils, excludes singletons
      li[[i]] <- x # fossil occurrences
      ext[i] <- e # extinctions
      orig[i] <- origin # origins
      i <- i + 1
    }
  }
  # default ordering: by true time of extinction
  if (sort_by == "true extinction") {
    li <- li[order(ext)]
    orig <- orig[order(ext)]
    ext <- ext[order(ext)]
  }

  if (sort_by == "observed extinction") {
    ext_order <- sapply(li, max) |> order()
    orig <- orig[ext_order]
    ext <- ext[ext_order]
    li <- li[ext_order]
  }

  if (sort_by == "true origination") {
    li <- li[order(orig)]
    ext <- ext[order(orig)]
    orig <- orig[order(orig)]
  }

  if (sort_by == "observed origination") {
    orig_order <- sapply(li, min) |> order()
    orig <- orig[orig_order]
    ext <- ext[orig_order]
    li <- li[orig_order]
  }

  true_ext <- ext
  true_origin <- orig
  last_occ <- sapply(li, max)
  first_occ <- sapply(li, min)
  other_occ <- lapply(li, function(x) x[x != max(x) & x != min(x)])

  # first & last occ, ranges etc.
  df_ext_true <- data.frame(taxon = factor(names_tax), t_ext = true_ext)
  df_ext_obs <- data.frame(taxon = factor(names_tax), last_occ = last_occ)
  df_orig_true <- data.frame(taxon = factor(names_tax), t_orig = true_origin)
  df_orig_obs <- data.frame(taxon = factor(names_tax), first_occ = first_occ)
  df_unobs_range_top <- data.frame(
    taxon = factor(rep(names_tax, 2)),
    lim = c(true_ext, last_occ)
  )
  df_unobs_range_bottom <- data.frame(
    taxon = factor(rep(names_tax, 2)),
    lim = c(true_origin, first_occ)
  )
  df_obs_range <- data.frame(
    taxon = factor(rep(names_tax, 2)),
    lim = c(first_occ, last_occ)
  )

  # data frame for occurrences
  names <- c("taxon", "occ")
  df_occ <- data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df_occ) <- names
  for (i in seq_along(li)) {
    df_occ <- rbind(
      df_occ,
      data.frame(
        "taxon" = rep(names_tax[i], length(other_occ[[i]])),
        "occ" = other_occ[[i]]
      )
    )
  }

  # background fill for origination and extinciton rate
  bg_fill_orig <- data.frame(y = seq(min(t), 1, by = 0.01))
  h <- diff(bg_fill_orig$y)[1]
  bg_fill_orig$ymin <- bg_fill_orig$y - h / 2
  bg_fill_orig$ymax <- bg_fill_orig$y + h / 2
  bg_fill_orig$orig <- f_orig(bg_fill_orig$y)

  bg_fill_ext <- data.frame(y = seq(1, max(t), by = 0.01))
  h <- diff(bg_fill_ext$y)[1]
  bg_fill_ext$ymin <- bg_fill_ext$y - h / 2
  bg_fill_ext$ymax <- bg_fill_ext$y + h / 2
  bg_fill_ext$ext <- f_ext(bg_fill_ext$y)

  bg_fill_sampling <- data.frame(y = seq(min(t), max(t), by = 0.01))
  h <- diff(bg_fill_sampling$y)[1]
  bg_fill_sampling$ymin <- bg_fill_sampling$y - h / 2
  bg_fill_sampling$ymax <- bg_fill_sampling$y + h / 2
  bg_fill_sampling$sampling <- f_sampling(bg_fill_sampling$y)

  rect_ext_1 <- data.frame(ymin = 1, ymax = 1.5) # in `lo` units, since you flip
  rect_ext_2 <- data.frame(ymin = 2.5, ymax = 3) # in `lo` units, since you flip
  rect_orig_1 <- data.frame(ymin = 0, ymax = 0.5) # in `lo` units, since you flip

  p <- ggplot(df_ext_true, aes(x = taxon, y = t_ext)) +
    geom_rect(
      data = bg_fill_sampling,
      inherit.aes = FALSE, # sampling rate
      aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf, fill = sampling)
    ) +
    scale_fill_gradient(
      low = "white",
      high = col_sampling,
      name = "Sampling probability",
      guide = guide_colorbar(order = 1),
      limits = c(0, max(bg_fill_sampling$sampling)),
      breaks = range(bg_fill_sampling$sampling),
      labels = c("low", "high")
    ) +
    geom_rect_pattern(
      data = rect_ext_2,
      aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
      pattern = "stripe",
      pattern_angle = 45,
      pattern_density = 0.05,
      pattern_spacing = 0.015,
      pattern_fill = col_ext,
      pattern_colour = NA,
      fill = NA,
      colour = col_ext,
      inherit.aes = FALSE
    ) +
    geom_rect_pattern(
      data = rect_orig_1,
      aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
      pattern = "stripe",
      pattern_angle = 45,
      pattern_density = 0.05,
      pattern_spacing = 0.015,
      pattern_fill = col_orig,
      pattern_colour = NA,
      fill = NA,
      colour = col_orig,
      inherit.aes = FALSE
    ) +
    geom_line(data = df_obs_range, aes(x = taxon, y = lim)) + # observed range
    geom_line(
      data = df_unobs_range_top,
      aes(x = taxon, y = lim),
      linetype = "dashed"
    ) + # unobserved ranges
    geom_line(
      data = df_unobs_range_bottom,
      aes(x = taxon, y = lim),
      linetype = "dashed"
    ) +
    geom_point(data = df_occ, aes(x = taxon, y = occ)) + # fossil samples (not FO & LOs)
    geom_point(data = df_ext_obs, aes(x = taxon, y = last_occ), color = "red") + # Last occ
    geom_point(shape = 4) + # true ext
    geom_point(
      data = df_orig_true,
      aes(x = taxon, y = true_origin),
      shape = 4
    ) + # true orig
    geom_point(
      data = df_orig_obs,
      aes(x = taxon, y = first_occ),
      color = "red"
    ) + # first occ
    scale_y_continuous(expand = expansion(0)) +
    coord_cartesian(ylim = c(min(t), max(t))) +
    labs(x = "Taxon", y = "Time [Myr]")
  return(p)
}

p2 <- line_plot_sampling(
  prob2,
  f_origin,
  f_ext,
  20,
  t_range = c(0, 3.5),
  sort_by = "observed extinction"
)
p2

p3 <- line_plot_sampling(
  prob3,
  f_origin,
  f_ext,
  20,
  t_range = c(0, 3.5),
  sort_by = "observed extinction"
)
p3

p4 <- line_plot_sampling(
  prob4,
  f_origin,
  f_ext,
  20,
  t_range = c(0, 3.5),
  sort_by = "observed extinction"
)
p4

p5 <- line_plot_sampling(
  prob5,
  f_origin,
  f_ext,
  20,
  t_range = c(0, 3.5),
  sort_by = "observed extinction"
)
p5

p <- ggarrange(
  p2,
  p3,
  p4,
  p5,
  ncol = 2,
  nrow = 2,
  common.legend = TRUE,
  legend = "bottom",
  labels = LETTERS[1:4]
)
p
ggsave(filename = "figs/ranges_sampling_prob.png", plot = p)
