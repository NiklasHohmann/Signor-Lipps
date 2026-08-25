library(ggplot2)

t = seq(-2, 4, by = 0.1)

f = function(x){
  y = rep(0.3, length(x))
  y[x > 2 & x < 3] = 3
  y[x > 1 & x < 1.5] = 2
  return(y)
}

plot(t, f(t))

conv = function(f, lambda, t){
  convolution_scalar_fun = function(t){lambda * exp(lambda * t) * integrate(function(x){f(x)* exp(-lambda * x)}, lower = t, upper = Inf, subdivisions = 100000)[[1]]}
  convolution_vectorized = sapply(t, convolution_scalar_fun)
  return(convolution_vectorized)
}
plot(t, conv(f, 2, t))


lambdas = seq(1, 10, by = 1)

names = c("ext", "t", "lambda")
df = data.frame(matrix(nrow = 0, ncol = length(names)))
names(df) = names
for (lambda in lambdas){
  df = rbind(df, data.frame(ext = conv(f, lambda, t), t = t, lambda = rep(lambda, length(t))))
}
df = rbind(df, data.frame(ext = f(t), t = t, lambda = Inf))
df$lambda = factor(df$lambda, levels = c(lambdas, Inf))

p = df |>
  ggplot(aes(x = t, y = ext, group = lambda, color = lambda)) +
  geom_line() +
  ylim(c(0,3)) +
  labs(x = "Time",
       y = "(Observed) Extinction Rate", 
       color = "Fossils per time unit") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.2, 0.6))

ggsave(filename = "figs/SLE_group_level_convolution.png",
       plot = p)




n_taxa = 10
names_tax = LETTERS[1:n_taxa]
lambda = 1
ext = StratPal::p3_var_rate(f, from = min(t), to = max(t), n = n_taxa, f_max = 5) |> sort()
li = list()
for (taxon in seq_along(ext)){
  li[[taxon]] = StratPal::p3(rate = lambda, from = min(t), to = ext[taxon])
}

png(filename = "figs/SLE_lineplot_true_order.png")
plot(NULL, xlim = c(0, 11), ylim = c(min(t), max(t)))
for (i in seq_along(ext)){
  occ = li[[i]]
  lines(x = c(i, i), y = c(min(t), max(occ)))
  lines(x = c(i,i), y = c(max(occ), ext[i]), lty = 3)
  points(x = i, y = ext[i], pch = 19)
  points(x = rep(i, length(li[[i]])), li[[i]])
}
dev.off()

ext_order = sapply(li, max) |> order()
ext_reordered = ext[ext_order]
li_reordered = li[ext_order]

png(filename = "figs/SLE_lineplot_observed_order.png")
plot(NULL, xlim = c(0, 11), ylim = c(min(t), max(t)),
     xlab = "Taxon id",
     ylab = "Time",
     main = "Range truncation")
for (i in 1:n_taxa){
  occ = li_reordered[[i]]
  lines(x = c(i, i), y = c(min(t), max(occ)))
  lines(x = c(i,i), y = c(max(occ), ext_reordered[i]), lty = 3)
  points(x = i, y = ext_reordered[i], pch = 1)
  points(x = i, y = max(occ), pch = 16)
  points(x = rep(i, length(occ)-1), occ[occ != max(occ)], pch = 4)
}
dev.off()
