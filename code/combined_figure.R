library(ggplot2)
library(ggnewscale)
library(patchwork)

col_ext <- "#0072B2"
col_orig <- "#E69F00"
col_sampling <- "#009E73"

t_max = 3.2
lambdas_used = c(0.25, 1, 2, 8)
times = seq(0, t_max, by = 0.1)

f_orig <- function(x) {
  y <- rep(0.3, length(x))
  y[x > 0 & x < 0.5] <- 10
  return(y)
}

df_true_orig_rate <- data.frame(
  x = c(0, 0, 0.5, 0.5, t_max),
  y = c(0.3, 10, 10, 0.3, 0.3)
)

f_ext <- function(x) {
  y <- rep(0.3, length(x))
  y[x > 2.5 & x < 3] <- 10
  y[x > 1.5 & x < 2] <- 2
  return(y)
}

df_true_ext_rate <- data.frame(
  x = c(0, 1.5, 1.5, 2, 2, 2.5, 2.5, 3, 3, t_max),
  y = c(0.3, 0.3, 2, 2, 0.3, 0.3, 10, 10, 0.3, 0.3)
)

origination_blurring_plot <- function(
  t = times,
  lambdas = lambdas_used
) {
  # convolution
  orig_blur_vector <- function(f, lambda, t) {
    orig_blur_scalar <- function(t) {
      lambda *
        exp(-lambda * t) *
        integrate(
          function(x) {
            f(x) * exp(lambda * x)
          },
          lower = -Inf,
          upper = t,
          rel.tol = 10^-9,
          subdivisions = 100000
        )[[1]]
    }
    return(sapply(t, orig_blur_scalar))
  }
  # calculate observed origination rate
  names <- c("ext", "t", "lambda")
  df <- data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df) <- names
  for (lambda in lambdas) {
    df <- rbind(
      df,
      data.frame(
        orig = orig_blur_vector(f_orig, lambda, t),
        t = t,
        lambda = rep(lambda, length(t))
      )
    )
  }
  # add "true" rate for plotting
  df <- rbind(
    df,
    data.frame(
      orig = df_true_orig_rate$y,
      t = df_true_orig_rate$x,
      lambda = rep(Inf, length(df_true_orig_rate$x))
    )
  )
  df$lambda <- factor(df$lambda, levels = c(lambdas, Inf))

  # black for true rate
  cols <- setNames(
    c(scales::hue_pal()(length(lambdas)), "black"),
    levels(df$lambda)
  )

  p <- df |>
    ggplot(aes(x = t, y = orig, group = lambda, color = lambda)) +
    geom_path(lwd = 1) +
    scale_color_manual(
      values = cols,
      breaks = levels(df$lambda),
      labels = c("Low", "Intermediate", "High", "Very high", "True rate"),
      drop = FALSE
    ) +
    labs(
      x = "Time [Myr]",
      y = "Origination rate [taxa/Myr]",
      color = "Sampling probability",
      title = "Origination rate"
    ) +
    theme(legend.position = "inside", legend.position.inside = c(0.75, 0.75)) +
    coord_flip(xlim = c(-0.05, 3.2), ylim = c(0, 11), expand = FALSE)
  return(p)
}

p <- origination_blurring_plot()
p


extinction_blurring_plot <- function(
  t = times,
  lambdas = lambdas_used
) {
  conv <- function(f, lambda, t) {
    convolution_scalar_fun <- function(t) {
      lambda *
        exp(lambda * t) *
        integrate(
          function(x) {
            f(x) * exp(-lambda * x)
          },
          lower = t,
          upper = Inf,
          rel.tol = 10^-9,
          subdivisions = 100000
        )[[1]]
    }
    convolution_vectorized <- sapply(t, convolution_scalar_fun)
    return(convolution_vectorized)
  }

  names <- c("ext", "t", "lambda")
  df <- data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df) <- names
  for (lambda in lambdas) {
    df <- rbind(
      df,
      data.frame(
        ext = conv(f_ext, lambda, t),
        t = t,
        lambda = rep(lambda, length(t))
      )
    )
  }

  df <- rbind(
    df,
    data.frame(
      ext = df_true_ext_rate$y,
      t = df_true_ext_rate$x,
      lambda = rep(Inf, length(df_true_ext_rate$x))
    )
  )
  df$lambda <- factor(df$lambda, levels = c(lambdas, Inf))

  # black for true rate
  cols <- setNames(
    c(scales::hue_pal()(length(lambdas)), "black"),
    levels(df$lambda)
  )

  p <- df |>
    ggplot(aes(x = t, y = ext, group = lambda, color = lambda)) +
    geom_path(lwd = 1) +
    scale_color_manual(
      values = cols,
      breaks = levels(df$lambda),
      labels = c("Low", "Intermediate", "High", "Very high", "True rate"),
      drop = FALSE
    ) +
    labs(
      x = "Time [Myr]",
      y = "Extinction rate [taxa/Myr]",
      color = "Sampling probability",
      title = "Extinction rate"
    ) +
    theme(legend.position = "inside", legend.position.inside = c(0.75, 0.25)) +
    coord_flip(xlim = c(-0.05, 3.2), ylim = c(0, 11), expand = FALSE)

  return(p)
}
p <- extinction_blurring_plot()
p

line_plot_full <- function(
  lambda,
  f_ext,
  f_orig,
  n_taxa,
  t_range,
  names_tax = LETTERS[1:n_taxa],
  sort_by = "observed extinction",
  allow_singletons = TRUE,
  plot_occ_legend = FALSE,
  plot_range_legend = FALSE
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
      f_orig,
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
    x <- StratPal::p3(
      rate = lambda, # sampling frequency
      from = origin,
      to = e
    )
    if (length(x) >= min_occ) {
      li[[i]] <- x # fossil occurrences
      ext[i] <- e # extinctions
      orig[i] <- origin # origins
      i <- i + 1
    }
  }

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
  df_unobs_range_top <- data.frame(
    taxon = factor(rep(names_tax, 2)),
    lim = c(true_ext, last_occ),
    rtype = rep("uru", length(names_tax) * 2)
  )
  df_unobs_range_bottom <- data.frame(
    taxon = factor(rep(names_tax, 2)),
    lim = c(true_origin, first_occ),
    rtype = rep("lur", length(names_tax) * 2)
  )
  df_obs_range <- data.frame(
    taxon = factor(rep(names_tax, 2)),
    lim = c(first_occ, last_occ),
    rtype = rep("or", length(names_tax) * 2)
  )
  df_ranges <- rbind(df_unobs_range_bottom, df_unobs_range_top, df_obs_range)

  # data frame for occurrences
  df_ext_true <- data.frame(
    taxon = factor(names_tax),
    t = true_ext,
    type = rep("ext", length(names_tax))
  )
  df_ext_obs <- data.frame(
    taxon = factor(names_tax),
    t = last_occ,
    type = rep("LO", length(names_tax))
  )
  df_orig_true <- data.frame(
    taxon = factor(names_tax),
    t = true_origin,
    type = rep("orig", length(names_tax))
  )
  df_orig_obs <- data.frame(
    taxon = factor(names_tax),
    t = first_occ,
    type = rep("FO", length(names_tax))
  )
  names <- c("taxon", "t", "type")
  df_occ <- data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df_occ) <- names
  for (i in seq_along(li)) {
    df_occ <- rbind(
      df_occ,
      data.frame(
        "taxon" = rep(names_tax[i], length(other_occ[[i]])),
        "t" = other_occ[[i]],
        "type" = rep("occ", length(other_occ[[i]]))
      )
    )
  }
  df_occ_all <- rbind(
    df_ext_obs,
    df_ext_true,
    df_orig_obs,
    df_orig_true,
    df_occ
  )

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

  p <- ggplot() +
    geom_rect(
      data = bg_fill_ext,
      inherit.aes = FALSE, # extinction rate background fill
      aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf, fill = ext)
    ) +
    scale_fill_gradient(
      low = "white",
      high = col_ext,
      name = "True\next. rate",
      guide = guide_colorbar(order = 1),
      limits = range(bg_fill_ext$ext),
      breaks = range(bg_fill_ext$ext),
      labels = c("low", "high")
    ) +
    new_scale_fill() +
    geom_rect(
      data = bg_fill_orig,
      inherit.aes = FALSE, # origination rate background fill
      aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf, fill = orig)
    ) +
    scale_fill_gradient(
      low = "white",
      high = col_orig,
      name = "True\norig. rate",
      guide = guide_colorbar(order = 2),
      limits = range(bg_fill_orig$orig),
      breaks = range(bg_fill_orig$orig),
      labels = c("low", "high")
    ) +
    geom_line(
      data = df_ranges,
      aes(x = taxon, y = lim, linetype = rtype),
      show.legend = ifelse(plot_range_legend, NA, FALSE)
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c("lur" = "dashed", "uru" = "dashed", "or" = "solid")
    ) +
    geom_point(
      data = df_occ_all,
      aes(x = taxon, y = t, shape = type, color = type),
      show.legend = ifelse(plot_occ_legend, NA, FALSE)
    ) +
    scale_color_manual(
      values = c(
        ext = "black",
        LO = "red",
        occ = "black",
        FO = "red",
        orig = "black"
      )
    ) +
    scale_shape_manual(
      values = c(ext = 4, LO = 19, occ = 19, FO = 19, orig = 4)
    ) +
    scale_y_continuous(expand = expansion(0)) +
    coord_cartesian(ylim = c(min(t), max(t))) +
    labs(x = "Taxon", y = "Time [Myr]", title = "Range truncation") +
    guides(
      color = guide_legend(position = "bottom", order = 1),
      linetype = guide_legend(position = "bottom", order = 2),
      fill = guide_colourbar(position = "right")
    )
  return(p)
}

p2 <- line_plot_full(
  lambda = 2,
  f_ext = f_ext,
  f_orig = f_orig,
  n_taxa = 20,
  t_range = c(0, 3.2),
  sort_by = "observed extinction",
  allow_singletons = TRUE
)
p2


p1 <- origination_blurring_plot()
p2 <- extinction_blurring_plot()
p3 <- line_plot_full(
  lambda = 2,
  f_ext = f_ext,
  f_orig = f_orig,
  n_taxa = 20,
  t_range = c(0, 3.2),
  sort_by = "observed extinction",
  allow_singletons = TRUE
)
p3
small_legend <- theme(
  legend.key.size = unit(0.35, "cm"),
  legend.text = element_text(size = 6),
  legend.title = element_text(size = 7),
  legend.margin = margin(0, 0, 0, 0),
  legend.box.spacing = unit(1, "pt") # gap between panel and legend
)
p <- p3 /
  (p2 + small_legend | p1 + small_legend) +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(tag_levels = "A")
p
ggsave(filename = "figs/SLE_combo_plot.png", plot = p)

p <- ggpubr::ggarrange(p1, p2, p3, ncol = 2, nrow = 2, labels = LETTERS[1:3])
p
