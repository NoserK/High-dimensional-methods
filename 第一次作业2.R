cauchyRand = function(i){
  x = 2*tan(pi*(runif(i)-0.5))
  return(x)
}
k = rep(1,time = 50000)
b1 = sapply(k,cauchyRand)
hist(b1, nclass = 50, probability = T)
lines(density(b1))

normalRand = function(i){
  z = (-2*log(runif(i)))^0.5*cos(2*pi*runif(i))
}

k = rep(1,time = 50000)
b2 = sapply(k,normalRand)
hist(b2, nclass = 50, probability = T)
lines(density(b2))

f = function(x){   #target distrivutiopn function
  exp(-(x)^(2.5))*(sin(x))^2
}
g1 = function(x){  # instrumental function normal
  1/(2*pi)^0.5*exp(-1/2*x^2)
}

a <- sapply(rep(1,time = 10000), normalRand)
M = 5
b = c()
k = 0
for (i in 1:10000){
  x = a[i]

  if (x > 0 && runif(1) <= f(x)/(M*g1(x))){
    k = k+1
    b = c(b,x)
  }
}
hist(b,nclass = 50,prob = T)
lines(density(b))
k/10000
1/M
(k/10000)*M