avg_rate = 3

f_ext = function(x){
  y = rep(1, length(x))
  y[x > 2.5 & x < 3] = 10
  return(y)
}

tt = seq(0, 3.5, by = 0.1)
plot(tt, f_ext(tt))

prob1temp = function(x){
  return(rep(1, length(x)))
}

t_min = 0
t_max = 3
total_1 = integrate(prob1temp, lower = t_min, upper = t_max)[[1]]

prob1 = function(x){
  return(prob1temp(x)/total_1 * (t_max - t_min) * avg_rate)
}

plot(tt, prob1(tt))

prob2temp = function(x){
  return(0.3 * x + 1)
}
total_2 =  integrate(prob2temp, lower = t_min, upper = t_max)[[1]]

prob2 = function(x){
  return(prob2temp(x)/total_2 * (t_max - t_min) * avg_rate)
}

prob3temp = function(x){
  y = rep(1, length(x))
  y[x > 1 & x < 2] = 10
  return(y)
}
total_3 = integrate(prob3temp, lower = t_min, upper = t_max)[[1]]
prob3 = function(x){
  return(prob3temp(x)/total_3 * (t_max - t_min) * avg_rate)
}

plot(tt, prob2(tt))
integrate(prob2, 0, 3)


n = 1000
ext = StratPal::p3_var_rate(f_ext, from = 0, to = 3, n = n)
ext
hist(ext)
library(StratPal)

lo = c()
li = list()
for (i in 1:n){
  lo[i] = p3_var_rate(prob3, from = 0, to = ext[i]) |> max()
}
lo
plot(lo)
hist(lo[!is.null(lo)])
