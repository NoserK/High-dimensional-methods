library(openxlsx)
library(stargazer)
library(lars)
library(car)
library(magrittr)
library(ridge)
faic<-function(n,x,y,beta,k)
{
  SSeq<-t(y-x%*%beta)%*%(y-x%*%beta)
  aicn<-k*log(SSeq)+2*n
  return(aicn)
}
fbic<-function(n,x,y,beta,k)
{
  SSeq<-t(y-x%*%beta)%*%(y-x%*%beta)
  bicn<-k*log(SSeq)+n/k*log(k)
  return(bicn)
}
p_lambda<-function(theta,lambda=0.69)
{
  p_lambda<-sapply(theta,function(x)
  {
    if(abs(x)<lambda)
    {
      return(lambda^2 - (abs(x) - lambda)^2)
    }
    else
    {
      return(lambda^2)
    }
  }
  )
  return(p_lambda)
}
##定义惩罚项导数
p_lambda_d<-function(theta,a=3.7,lambda=0.69)
{
  if(abs(theta)>lambda)
  {
    if(a*lambda>theta)
    {
      return((a*lambda-theta)/(a-1))
    }
    else
    {
      return(0)
    }
  }
  else
  {
    return(lambda)
  }
}
# ##当beta_j0不等于0，定义惩罚项导数近似
# p_lambda_d_apro<-function(beta_j0,beta_j,a = 3.7, lambda = 2){
#   return(beta_j * p_lambda_d(beta = beta_j0,a = a, lambda = lambda)/abs(beta_j0))
# }
# ##当beta_j0 不等于0，指定近似惩罚项,使用泰勒展开逼近
# p_lambda_apro<-function(beta_j0,beta_j,a = 3.7, lambda = 2){
#   if(abs(beta_j0)< 1e-16){
#     return(0)
#   }else{
#     p_lambda<-p_lambda(theta = beta_j0, lambda = lambda)
#       0.5 * (beta_j^2 - beta_j0^2) * p_lambda_d(theta = beta_j0, a = a, lambda = lambda)/abs(beta_j0)
#   }
# }
#define the log-likelihood function
loglikelihood_SCAD<-function(x,y,b,lambda=0.69)
{
  linear_comb<-as.vector(x%*%b/(nrow(x)*12000000))
  ll<-sum(y*linear_comb)+sum(log(1/(1+exp(linear_comb))))-sum(p_lambda(theta=b))/12000000
  return (ll)
}
data=read.xlsx("E:\\2018-2019\\1.xlsx",
               fill = T, colNames = T) #根据自己的保存位置更改
data=na.omit(data) #去掉缺失值
write.xlsx(data,"E:\\2018-2019\\trimeddata.xlsx")

#summary statistics
ss=summary(data)
stargazer(data) #输出latex文本

#OLS
data2=data[,c(1,2,4,6,8,10,14)] #去掉_24h
AQI=data2$AQI
PM2.5=data2$PM2.5
PM10=data2$PM10
SO2=data2$SO2
NO2=data2$NO2
O3=data2$O3
CO=data2$CO
y<-data2$AQI
x<-matrix(c(data2$PM2.5,data2$PM10,data2$SO2,data2$NO2,data2$O3,data2$CO),nrow = 1442,ncol = 6)
k=1442
xx<-cor(x)
eigen(xx)
kappa(xx,exact = T)
betals<-(solve(t(x)%*%x)%*%t(x))%*%y
yhat=x%*%betals
b0<-betals
b1<-b0
ll_list<-list(ll.old)
lam<-seq(0.01,2,0.01)
m<-matrix(rep(0,200),6,200) 
gcv<-matrix(rep(0,200),1,200) 
aic<-matrix(rep(0,200),1,200)
bic<-matrix(rep(0,200),1,200)
for (i in 1:200)  
{ 
  lambda<-lam[i] 
  ##b.best用于存放历史最大似然值对应系数
  b.best_SCAD<-b0
  # the initial value of loglikelihood
  ll.old<-loglikelihood_SCAD(x=x,y=y,b=b0)
  # initialize the difference between the two steps of theta
  diff<-1  
  #record the number of iterations
  iter<-0
  #set the threshold to stop iterations
  epsi<-1e-6
  #the maximum iterations  
  max_iter<-100000
  #初始化一个列表用于存放每一次迭代的系数结果
  b_history<-list(data.frame(b0))
  #初始化列表用于存放似然值
  #######-------SCAD迭代---------
  while(diff>epsi&iter<max_iter)
  {
    for(j in 1:6)
    {
      if(abs(b0[j,1])<10^-6)
      {
        next()
      }
      else
      {
        #线性部分
        linear_comb=as.vector(x%*%b0)
        #分子
        nominator=sum(y*x[,j]-x[,j])+nrow(x)*b0[j,1]*p_lambda_d(theta = b0[j,1])/abs(b0[j,1])
        #分母,即二阶导部分
        denominator=-sum(x[,j]^2*exp(linear_comb)/(1+exp(linear_comb))^2)+nrow(x)*p_lambda_d(theta = b0[j,1])/abs(b0[j,1])
        #2-(3) :更新b0[j]
        b0[j,1]<-b0[j,1]-nominator/denominator
        #2-(4)
        if(abs(b0[j,1])<10^-6)
        {
          b0[j,1]<-0
        }
        
      }
    }
    #更新似然值
    ll.new<- loglikelihood_SCAD(x=x,y=y,b=b0)
    #若似然值有所增加，则将当前系数保存
    if(ll.new>ll.old)
    {
      #更新系数
      b.best_SCAD<-b0
    }
    
    #求差异
    diff<-abs((ll.new-ll.old)/ll.old)
    ll.old<-ll.new
    iter<-iter+1
    b_history[[iter]]<-data.frame(b0)
    ll_list[[iter]]<-ll.old
  }
  m[,i]=b.best_SCAD
  bet<-b.best_SCAD
  for (l in 1:6)
  { 
    if (abs(b.best_SCAD[l,1]<10^-4)) 
      bet[l,1]<-10^-5
  } 
  p<-x%*%solve((t(x)%*%x)+k*lambda*diag(c(bet)))%*%t(x)
  tr<-sum(diag(p)) 
  t1<-t(y-x%*%m[,i])%*%(y-x%*%m[,i]) 
  gcv[1,i]<-(t(y-x%*%m[,i]))%*%(y-x%*%m[,i])/(1-tr/k)^2
  n=sum(ifelse(m[,i]>10^-5,1,0))
  aic[1,i]<-faic(n,x,y,m[,i],k)
  bic[1,i]<-fbic(n,x,y,m[,i],k)
}
gcvmin<-which.min(gcv)
aicmin<-which.min(aic)
bicmin<-which.min(bic)
betagcv<-m[,gcvmin]
betaaic<-m[,aicmin]
betabic<-m[,bicmin]
betagcv
betabic
betaaic