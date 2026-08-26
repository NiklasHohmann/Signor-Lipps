library(ggplot2)

t = seq(0, 3.5, by = 0.1)

f = function(x){
  y = rep(0.3, length(x))
  y[x > 2.5 & x < 3] = 10
  y[x > 1 & x < 1.5] = 2
  return(y)
}

df_true_ext_rate = data.frame(
  x = c(0,1,1,1.5,1.5,2.5,2.5,3,3,3.5),
  y = c(0.3, 0.3,2,2,0.3,0.3,10,10,0.3,0.3)
)

blurring_plot = function(){
  conv = function(f, lambda, t){
    convolution_scalar_fun = function(t){lambda * exp(lambda * t) * integrate(function(x){f(x)* exp(-lambda * x)}, lower = t, upper = Inf, rel.tol = 10^-9,subdivisions = 100000)[[1]]}
    convolution_vectorized = sapply(t, convolution_scalar_fun)
    return(convolution_vectorized)
  }
  
  
  lambdas = 10^seq(log10(0.1), log10(10), length.out = 5)
  lambdas = c(0.25, 0.5, 1, 2, 4, 8)
  t_max = c()
  ext_max = c()
  for (lambda in lambdas){
    ext = conv(f, lambda, t)
    t_max = c(t_max, t[which.max(ext)])
    ext_max = c(ext_max, ext[which.max(ext)])
  }
  
  names = c("ext", "t", "lambda")
  df = data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  for (lambda in lambdas){
    df = rbind(df, data.frame(ext = conv(f, lambda, t), t = t, lambda = rep(lambda, length(t))))
  }
  df$lambda = factor(df$lambda, levels = c(lambdas, Inf))
  
  p = df |>
    ggplot(aes(x = t, y = ext, group = lambda, color = lambda)) +
    geom_path(data = df_true_ext_rate, aes(x = x, y = y), inherit.aes = FALSE, lty = "dotted", lwd = 1) +
    geom_line(lwd = 1) +
    ylim(c(0,10)) +
    labs(x = "Time or stratigraphic position",
         y = "(Observed) Extinction rate", 
         color = "Fossil sampling\nrate [#/unit]") +
    theme(legend.position = "inside",
          legend.position.inside = c(0.2, 0.7)) +
    coord_cartesian(xlim = c(0, 3.5), ylim = c(0, 11) ,expand = FALSE)
  
  return(p)
  
}




line_plot = function(lambda,
                     f_ext,
                     n_taxa,
                     t_range,
                     names_tax = LETTERS[1:n_taxa],
                     sort_by = "observed extinctino"){
  t = t_range
  stopifnot(n_taxa == length(names_tax))
  ext = c()
  li = list()
  i = 1
  while (length(li) < n_taxa){
    e = StratPal::p3_var_rate(f_ext, from = min(t), to = max(t), n = 1, f_max = 50)
    x = StratPal::p3(rate = lambda, from = min(t), to = e)
    if (length(x) > 0){
      li[[i]] = x
      ext[i] = e
      i = i+1
    }
  }
  li = li[order(ext)]
  ext = ext[order(ext)]
  
  sort_by = "observed extinction"
  if (sort_by == "observed extinction"){
    ext_order = sapply(li, max) |> order()
    ext = ext[ext_order]
    li = li[ext_order]
  }
  
  true_ext = ext
  last_occ = sapply(li, max)
  other_occ = lapply(li, function(x) x[x != max(x)])
  
  df_ext_true = data.frame(taxon = factor(names_tax), t_ext = true_ext)
  df_ext_obs = data.frame(taxon = factor(names_tax), last_occ = last_occ)
  
  df_unobs_range = data.frame(taxon = factor(rep(names_tax, 2)),
                              lim = c(true_ext, last_occ))
  df_obs_range = data.frame(taxon = factor(rep(names_tax, 2)),
                            lim = c(rep(min(t), n_taxa), last_occ))
  names = c("taxon", "occ")
  df_occ = data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df_occ) = names
  for (i in seq_along(li)){
    df_occ = rbind(df_occ, data.frame("taxon" = rep(names_tax[i], length(other_occ[[i]])),
                                      "occ" = other_occ[[i]]))
  }
  
  bg_fill = data.frame(y = seq(min(t), max(t), length.out = 400))
  h = diff(bg_fill$y)[1]
  bg_fill$ymin = bg_fill$y - h / 2
  bg_fill$ymax = bg_fill$y + h / 2
  bg_fill$val  = f_ext(bg_fill$y)
  
  p =   ggplot(df_ext_true, aes(x = taxon, y = t_ext)) +
    geom_rect(data = bg_fill, inherit.aes = FALSE,
              aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf, fill = val)) +
    scale_fill_gradient(low = "white", high = "steelblue", name = "Ext. rate") +
    geom_line(data = df_unobs_range, aes(x = taxon, y = lim), linetype = "dashed") +
    geom_line(data = df_obs_range, aes(x = taxon, y = lim)) +
    geom_point(data = df_occ, aes(x = taxon, y = occ)) +
    geom_point(data = df_ext_obs, aes(x = taxon, y = last_occ), color = "red") +
    geom_point(shape = 4) +
    scale_y_continuous(expand = expansion(0)) +
    coord_cartesian(ylim = c(min(t), max(t))) +
    labs(x = "Taxon", y = "Time or stratigrahpic position")
  return(p)
}

p2 = line_plot(lambda = 1,
          f_ext = f,
          n_taxa = 20,
          t_range = c(0,3.5))

p1 = blurring_plot()

p = ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, labels = LETTERS[1:2])

ggsave(filename = "figs/SLE.png",
       plot = p)
p
