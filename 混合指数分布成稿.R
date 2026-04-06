###正式实验
#混合指数分布抽样 
b1 = rexp(10000,1/50)#1/theta*exp{-1/theta*x}

b2 = rexp(10000,1/10)

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


mydata = b3[1:500]
hist(mydata, probability = T)

lamda_1_ls = c()
lamda_2_ls = c()
alpha_1 = alpha_2 = beta_1 = beta_2 = 1

z = rbinom(500,1,0.5)

for ( ou in 1 : 4000){
  p = rbeta(1,(sum(z)+1), (length(z)-sum(z)+1))
  
  #对 lamda_1 的满条件抽样
  alpha = alpha_1 - 1
  for( i in 1:length(mydata)){
    alpha = alpha + as.numeric(z[i] == 1)
  }
  
  rate = 1/beta_1
  for( i in 1:length(mydata)){
    rate = rate + as.numeric(z[i] == 1) * mydata[i]
  }
  lamda_1 = rgamma(1,alpha,rate)
  
  #对 lamda_2 的满条件抽样
  alpha = alpha_2 - 1
  for( i in 1:length(mydata)){
    alpha = alpha + as.numeric(z[i] == 0)
  }
  
  rate = 1/beta_1
  for( i in 1:length(mydata)){
    rate = rate + as.numeric(z[i] == 0) * mydata[i]
  }
  lamda_2 = rgamma(1,alpha,rate)
  
  #对 z 的满条件抽样
  
  for (i in 1:length(z)){
    prob_z = p * lamda_1 * exp(-lamda_1 * mydata[i] )/(p * lamda_1 * exp(-lamda_1 * mydata[i] ) + (1-p) * lamda_2 * exp(-lamda_2 * mydata[i] ))
    z[i] = rbinom(1,1,prob_z)
  }
  
  lamda_1_ls[ou] = lamda_1
  lamda_2_ls[ou] = lamda_2
}

1/mean(lamda_1_ls[3000:4000])
1/mean(lamda_2_ls[3000:4000])
