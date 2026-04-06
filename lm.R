rm(list=ls())
gc()
library(openxlsx)
library(stargazer)
library(lars)
library(car)
library(magrittr)
library(ridge)
library(ncvreg)

data=read.xlsx("E:\\2018-2019\\1.xlsx",
               fill = T, colNames = T) #根据自己的保存位置更改
data=na.omit(data) #去掉缺失值
write.xlsx(data,"E:\\2018-2019\\trimeddata.xlsx")

#data=read.xlsx("C:\\Users\\apple\\Desktop\\1.xlsx",
#               fill = T, colNames = T) #根据自己的保存位置更改
#data=na.omit(data) #去掉缺失值
#write.xlsx(data,"C:\\Users\\apple\\Desktop\\trimeddata.xlsx")

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
model.ols=lm(AQI~.,data=data2)
summary(model.ols)

#influence analysis
influence.measures(model.ols)
op=par(mfrow=c(2,2), mar=0.4+c(4,4,1,1), 
       oma=c(0,0,2,0))
plot(model.ols, 1:4)
par(op)
pickout=c(-1110,-1278,-1467)#异常值

model.ols.correct=lm(AQI~PM2.5+PM10+SO2+NO2+O3+CO,subset=pickout)#去掉异常值

#collinearity test
vif(model.ols.correct)
cor(data2)
plot(data2[,1:7])

#stepwise
stepwise=step(model.ols,direction="backward")
summary(stepwise)

#lasso
x=as.matrix(data2[,2:7])
model.lasso=lars(x,AQI,type="lasso")
coef(model.lasso,mode="norm")
summary(model.lasso)
plot(model.lasso)
CV.lasso=cv.lars(x,AQI,plot.it = T)
best=CV.lasso$index[which.min(CV.lasso$cv)] #交叉验证
coef.lasso=coef.lars(model.lasso, mode='fraction', s=best)

#ridge regression
model.ridge=linearRidge(AQI~.,data=data2)
summary(model.ridge)
plot(model.ridge)
coef.ridge=coef(model.ridge)

##functions
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
falasso<-function(x,y,k)
{
  betals<-solve(t(x)%*%x)%*%t(x)%*%y
  lambda<-seq(0.11,1.60,0.01)
  m2<-matrix(rep(0,150),6,150) 
  gcv2<-matrix(rep(0,150),1,150) 
  aic2<-matrix(rep(0,150),1,150)
  bic2<-matrix(rep(0,150),1,150)
  for (i in 1:150)  
  { 
    lam<-lambda[i] 
    beta0<-(solve(t(x)%*%x)%*%t(x))%*%y
    a<-matrix(rep(0,6),6,1) 
    while(t(beta0-a)%*%(beta0-a)>10^-5)
    { 
      bet<-matrix(rep(0,6),6,1) 
      for(j in 1:6)
      { 
        if(abs(beta0[j,1]<1e-4)) 
        {
          beta0[j,1]<-10^-5
        }
        bet[j,1]<-1/abs(beta0[j,1])
      } 
      beta11<-beta0-solve((t(x)%*%x)+k*lam*diag(c(bet/betals)))%*%((-t(x))%*%(y-x%*%beta0)+k*lam*diag(c(bet/betals))%*%beta0)
      for(g in 1:6)
      { 
        if(abs(beta11[g,1])<1*10^-3)
          beta11[g,1]<-0 
      } 
      a<-beta0;
      beta0<-beta11; 
    }
    m2[,i]<-beta0
    bet<-beta0
    for (l in 1:6)
    { 
      if (abs(beta0[l,1]<10^-4)) 
        bet[l,1]<-10^-5
    } 
    p<-x%*%solve((t(x)%*%x)+k*lam*diag(c(bet)))%*%t(x)
    tr<-sum(diag(p)) 
    t1<-t(y-x%*%beta0)%*%(y-x%*%beta0) 
    gcv2[1,i]<-(t(y-x%*%beta0))%*%(y-x%*%beta0)/(1-tr/k)^2
    n=sum(ifelse(m2[,i]>1e-5,1,0))
    aic2[1,i]<-faic(n,x,y,m2[,i],k)
    bic2[1,i]<-fbic(n,x,y,m2[,i],k)
  }
  gcvmin2<-which.min(gcv2)
  aicmin2<-which.min(aic2)
  bicmin2<-which.min(bic2)
  betagcv2<-m2[,gcvmin2]
  betaaic2<-m2[,aicmin2]
  betabic2<-m2[,bicmin2]
  return(list(betagcv2=betagcv2,betaaic2=betaaic2,betabic2=betabic2,aicmin2=aicmin2,bicmin2=bicmin2,gcvmin2=gcvmin2))
}
fsubsetbic<-function(x,y,n)
{
  q<-1
  s<-combn(c(1:6),1)
  m<-array(0,dim<-c(6,n,1))
  betta1<-matrix(0,6,6)
  BIC1<-c()
  for(i in 1:6){
    m[i, , ]<-cbind(x[,s[1,i]])
    xq<-m[i,,]
    SSeq<-t(y)%*%(diag(n)-xq%*%solve(t(xq)%*%xq)%*%t(xq))%*%y
    BIC1[i]<-log(SSeq/n)+q/n*log(n)
    betta1[i,i]<-solve(t(xq)%*%xq)%*%t(xq)%*%y
  }
  betta1
  BIC1
  q<-2
  s<-combn(c(1:6),2)
  m<-array(0,dim<-c(15,n,2))
  BIC2<-c()
  betta2<-matrix(0,6,15)
  for(i in 1:15){
    m[i, , ]<-cbind(x[,s[1,i]],x[,s[2,i]])
    xq<-m[i,,]
    SSeq<-t(y)%*%(diag(n)-xq%*%solve(t(xq)%*%xq)%*%t(xq))%*%y
    BIC2[i]<-log(SSeq/n)+q/n*log(n)
    p<-solve(t(xq)%*%xq)%*%t(xq)%*%y
    betta2[s[1,i],i]<-p[1,1]
    betta2[s[2,i],i]<-p[2,1]
  }
  BIC2
  betta2
  q<-3
  s<-combn(c(1:6),3)
  m<-array(0,dim<-c(20,n,3))
  betta3<-matrix(0,6,20)
  BIC3<-c()
  for(i in 1:20){
    m[i, , ]<-cbind(x[,s[1,i]],x[,s[2,i]],x[,s[3,i]])
    xq<-m[i,,]
    SSeq<-t(y)%*%(diag(n)-xq%*%solve(t(xq)%*%xq)%*%t(xq))%*%y
    BIC3[i]<-log(SSeq/n)+q/n*log(n)
    p<-solve(t(xq)%*%xq)%*%t(xq)%*%y
    betta3[s[1,i],i]<-p[1,1]
    betta3[s[2,i],i]<-p[2,1]
    betta3[s[3,i],i]<-p[3,1]
  }
  betta3
  BIC3
  q<-4
  s<-combn(c(1:6),4)
  m<-array(0,dim<-c(15,n,4))
  betta4<-matrix(0,6,15)
  BIC4<-c()
  for(i in 1:15){
    m[i, , ]<-cbind(x[,s[1,i]],x[,s[2,i]],x[,s[3,i]],x[,s[4,i]])
    xq<-m[i,,]
    SSeq<-t(y)%*%(diag(n)-xq%*%solve(t(xq)%*%xq)%*%t(xq))%*%y
    BIC4[i]<-log(SSeq/n)+q/n*log(n)
    p<-solve(t(xq)%*%xq)%*%t(xq)%*%y
    betta4[s[1,i],i]<-p[1,1]
    betta4[s[2,i],i]<-p[2,1]
    betta4[s[3,i],i]<-p[3,1]
    betta4[s[4,i],i]<-p[4,1]
  }
  BIC4
  betta4
  q<-5
  BIC5<-c()
  s<-combn(c(1:6),5)
  betta5<-matrix(0,6,6)
  for(i in 1:6){
    xq<-x[,-i]
    SSeq<-t(y)%*%(diag(n)-xq%*%solve(t(xq)%*%xq)%*%t(xq))%*%y
    BIC5[i]<-log(SSeq/n)+q/n*log(n)
    p<-solve(t(xq)%*%xq)%*%t(xq)%*%y
    betta5[s[1,i],i]<-p[1,1]
    betta5[s[2,i],i]<-p[2,1]
    betta5[s[3,i],i]<-p[3,1]
    betta5[s[4,i],i]<-p[4,1]
    betta5[s[5,i],i]<-p[5,1]
  }
  BIC5
  betta5
  q<-6
  xq<-x
  SSeq<-t(y)%*%(diag(n)-xq%*%solve(t(xq)%*%xq)%*%t(xq))%*%y
  BIC6<-log(SSeq/n)+q/n*log(n)
  betta6<-solve(t(xq)%*%xq)%*%t(xq)%*%y
  betta6
  BIC6
  
  betta<-matrix(c(betta1,betta2,betta3,betta4,betta5,betta6),6,63)
  BIC<-c(BIC1,BIC2,BIC3,BIC4,BIC5,BIC6)
  t<-which.min(BIC)
  betabic<-betta[,t]
  return(betabic)
}

##subset and lasso
y<-data2$AQI
x<-matrix(c(data2$PM2.5,data2$PM10,data2$SO2,data2$NO2,data2$O3,data2$CO),nrow = 1442,ncol = 6)
k=1442
xx<-cor(x)
eigen(xx)
kappa(xx,exact = T)
betals<-(solve(t(x)%*%x)%*%t(x))%*%y
yhat=x%*%betals
betasubbic<-fsubsetbic(x,y,k)
betasubbic
betalasso<-falasso(x,y,k)
betalassoaic<-betalasso$betaaic2
betalassobic<-betalasso$betabic2
betalassogcv<-betalasso$betagcv2
betalassoaic
betalassobic
betalassogcv
lambda=0.10+0.01*betalasso$bicmin2
lambda

#SCAD
model.SCAD=ncvreg(x,AQI,penalty="SCAD")
plot(model.SCAD)
cv.SCAD=cv.ncvreg(x,AQI)
plot(cv.SCAD)
lbd=cv.SCAD$lambda.min #best lambda
model.SCAD.correct=ncvreg(x,AQI,penalty="SCAD",lambda=lbd)
model.SCAD.correct$beta
plot(model.SCAD.correct)
