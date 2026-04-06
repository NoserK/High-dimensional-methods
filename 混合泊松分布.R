#混合正态分布模型模拟：
a1 = rnorm(10000,0,1)
a2 = rnorm(10000,4,1)
a3 = c()
for (i in 1:10000){
  if (runif(1) <= 0.5){
    a3[i] = a1[i]
  }
  else{
    a3[i] = a2[i]
  }
}

hist(a3, nclass = 50, probability = T)
lines(density(a3))

#混合指数分布模型模拟：
b1 = rexp(10000,1)#1/theta*exp{-1/theta*x}

b2 = rexp(10000,0.125)

b3 = c()
for (i in 1 :10000){
  if (runif(1) <= 0.5){
    b3[i] = b1[i]
  }
  else{
    b3[i] = b2[i]
  }
}
hist(b3, nclass = 50, probability = T)
lines(density(b3))6+

  
#混合泊松分布抽样 
b1 = rpois(10000,1)#1/theta*exp{-1/theta*x}

b2 = rpois(10000,10)

b3 = c()
for (i in 1 :10000){
  if (runif(1) <= 0.5){
    b3[i] = b1[i]
  }
  else{
    b3[i] = b2[i]
  }
}
hist(b3, nclass = 50, probability = T)
lines(density(b3))
  

mydata = b3[1:30]
hist(mydata, probability = T)
#带数据的联合分布
# 1
# (gamma(alpha_1) * gamma(alpha_2) * beta_1 ^ alpha_1 * beta_2 ^alpha_2) * lamda_1 ^ (alpha_1 - 1) * lamda_2 ^ (alpha_2 - 1) * exp(-lamda_1 / beta_1 -lamda_2 / beta_2) 
#
mul = 1
z = rbinom()
for (i in 1 : (length(mydata))){
  mul = mul * (as.numeric(z[i] == 1) * p * lamda_1 ^ mydata[i] * exp(-lamda_1)/gamma(mydata[i] + 1) + as.numeric(z[i] == 0) * (1-p) * lamda_2 ^ mydata[i] * exp(-lamda_2)/gamma(mydata[i] + 1)) 
}
joint = function(p,lamda_1,lamda_2,z,mydata){ #联合分布函数
  mul = 1
  for (i in 1 : (length(mydata))){
    mul = mul * (as.numeric(z[i] == 1) * p * lamda_1 ^ mydata[i] * exp(-lamda_1)/gamma(mydata[i] + 1) + as.numeric(z[i] == 0) * (1-p) * lamda_2 ^ mydata[i] * exp(-lamda_2)/gamma(mydata[i] + 1)) 
  }
  mul * (gamma(alpha_1) * gamma(alpha_2) * beta_1 ^ alpha_1 * beta_2 ^alpha_2) * lamda_1 ^ (alpha_1 - 1) * lamda_2 ^ (alpha_2 - 1) * exp(-lamda_1 / beta_1 -lamda_2 / beta_2)
  
}



#对 p 的满条件抽样

p = rbeta(1,(sum(z)+1), (length(z)-sum(z)+1))

#对 lamda_1 的满条件抽样
alpha = alpha_1 - 1
for( i in 1:length(mydata)){
  alpha = alpha + as.numeric(mydata[i] == 1) * mydata[i]
}

rate = 1/beta_1
for( i in 1:length(mydata)){
  rate = rate + as.numeric(mydata[i] == 1)
}
lamda_1 = rgamma(1,alpha,rate)

#对 lamda_2 的满条件抽样
alpha = alpha_2 - 1
for( i in 1:length(mydata)){
  alpha = alpha + as.numeric(mydata[i] == 0) * mydata[i]
}

rate = 1/beta_1
for( i in 1:length(mydata)){
  rate = rate + as.numeric(mydata[i] == 0)
}
lamda_1 = rgamma(1,alpha,rate)

#对 z 的满条件抽样
