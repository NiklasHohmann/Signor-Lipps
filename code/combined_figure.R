library(ggplot2)
library(ggnewscale)
library(patchwork)

f_orig <- function(x) {
  y <- rep(0.3, length(x))
  y[x > 0 & x < 0.5] <- 10
  return(y)
}

df_true_orig_rate <- data.frame(
  x = c(0, 0, 0.5, 0.5, 3.5),
  y = c(0.3, 10, 10, 0.3, 0.3)
)

f_ext <- function(x) {
  y <- rep(0.3, length(x))
  y[x > 2.5 & x < 3] <- 10
  y[x > 1 & x < 1.5] <- 2
  return(y)
}

df_true_ext_rate <- data.frame(
  x = c(0, 1, 1, 1.5, 1.5, 2.5, 2.5, 3, 3, 3.5),
  y = c(0.3, 0.3, 2, 2, 0.3, 0.3, 10, 10, 0.3, 0.3)
)

origination_blurring_plot <- function(
  t = seq(0, 3.5, by = 0.1),
  lambdas = 2^seq(-2, 3)
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
      labels = c(as.character(lambdas), "True\nrate"),
      drop = FALSE
    ) +
    ylim(c(0, 10)) +
    labs(
      x = "Time [Myr]",
      y = "Origination rate [taxa/Myr]",
      color = "Fossil sampling\nfrequency [#/Myr]",
      title = "Origination rate"
    ) +
    theme(legend.position = "inside", legend.position.inside = c(0.75, 0.75)) +
    coord_cartesian(xlim = c(-0.05, 3.5), ylim = c(0, 11), expand = FALSE) +
    coord_flip()
  return(p)
}

p <- origination_blurring_plot()
p


extinction_blurring_plot <- function(
  t = seq(0, 3.5, by = 0.1),
  lambdas = 2^seq(-2, 3)
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
        ext = conv(f, lambda, t),
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
      labels = c(as.character(lambdas), "True\nrate"),
      drop = FALSE
    ) +
    ylim(c(0, 10)) +
    labs(
      x = "Time [Myr]",
      y = "Extinction rate [taxa/Myr]",
      color = "Fossil sampling\nfrequency [#/Myr]",
      title = "Extinction rate"
    ) +
    theme(legend.position = "inside", legend.position.inside = c(0.75, 0.25)) +
    coord_cartesian(xlim = c(-0.05, 3.5), ylim = c(0, 11), expand = FALSE) +
    coord_flip()

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

  p <- ggplot(df_ext_true, aes(x = taxon, y = t_ext)) +
    geom_rect(
      data = bg_fill_ext,
      inherit.aes = FALSE, # extinction rate background fill
      aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf, fill = ext)
    ) +
    scale_fill_gradient(
      low = "white",
      high = "steelblue",
      name = "True\next. rate\n[taxa/Myr]",
      guide = guide_colorbar(order = 1)
    ) +
    new_scale_fill() +
    geom_rect(
      data = bg_fill_orig,
      inherit.aes = FALSE, # origination rate background fill
      aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf, fill = orig)
    ) +
    scale_fill_gradient(
      low = "white",
      high = "darkorange",
      name = "True\norig. rate\n[taxa/Myr]",
      guide = guide_colorbar(order = 2)
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
    labs(x = "Taxon", y = "Time [Myr]", title = "Taxon range truncation")
  return(p)
}

p2 <- line_plot_full(
  lambda = 1,
  f_ext = f_ext,
  f_orig = f_orig,
  n_taxa = 20,
  t_range = c(0, 3.5),
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
  t_range = c(0, 3.5)
)

small_legend <- theme(
  legend.key.size = unit(0.35, "cm"),
  legend.text = element_text(size = 6),
  legend.title = element_text(size = 7),
  legend.margin = margin(0, 0, 0, 0),
  legend.box.spacing = unit(1, "pt") # gap between panel and legend
)
p <- (p2 + small_legend | p1 + small_legend) /
  p3 +
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(tag_levels = "A")
ggsave(filename = "figs/SLE_combo_plot.png", plot = p)

p <- ggpubr::ggarrange(p1, p2, p3, ncol = 2, nrow = 2, labels = LETTERS[1:3])
p
