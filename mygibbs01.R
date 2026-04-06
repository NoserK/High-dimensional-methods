gelmanrubin<-function(psi)
{
  #psi[i.j] is the statstic psi(x[i,1:j])
  #for chain in i-th row of X
  psi<-as.matrix(psi)
  n=ncol(psi)
  k=nrow(psi)
  psi.means=rowMeans(psi)
  B=n*var(psi.means)#between variances
  psi.w=apply(psi,1,"var")
  W=mean(psi.w)
  v.hat=W*(n-1)/n+(B/(n*k))
  r.hat=v.hat/W
  return(r.hat)
}
fmydata<-function(i)
{
  b1=rexp(10000,1/50)#1/theta*exp{-1/theta*x}
  b2=rexp(10000,1/10)
  b3=c()
  for(i in 1:10000)
  {
    if(runif(1)<=0.5)
    {
      b3[i]=b1[i]
    }
    else
    {
      b3[i]=b2[i]
    }
  }
  mydata = b3[1:500]
  return(mydata)
}
mysampling<-function(p,N,mydata)
{
  lamda_1_e=c()
  lamda_2_e=c()
  p_e=c()
  alpha_1=alpha_2=10*runif(1)
  beta_1=beta_2=50*runif(1)
  z=rbinom(500,1,p)
  for(ou in 1:N)
  {
    p=rbeta(1,(sum(z)+1),(length(z)-sum(z)+1))
    #对 lamda_1 的满条件抽样
    alpha=alpha_1-1
    for(i in 1:length(mydata))
    {
      alpha=alpha+as.numeric(z[i]==1)
    }
    rate=1/beta_1
    for(i in 1:length(mydata))
    {
      rate=rate+as.numeric(z[i]==1)*mydata[i]
    }
    lamda_1=rgamma(1,alpha,rate)
    #对 lamda_2 的满条件抽样
    alpha=alpha_2-1
    for(i in 1:length(mydata))
    {
      alpha=alpha+as.numeric(z[i]==0)
    }
    rate=1/beta_1
    for(i in 1:length(mydata))
    {
      rate=rate+as.numeric(z[i]==0)*mydata[i]
    }
    lamda_2=rgamma(1,alpha,rate)
    #对 z 的满条件抽样
    for (i in 1:length(z))
    {
      prob_z=p*lamda_1*exp(-lamda_1*mydata[i])/(p*lamda_1*exp(-lamda_1*mydata[i])+(1-p)*lamda_2*exp(-lamda_2*mydata[i]))
      z[i]=rbinom(1,1,prob_z)
    }
    lamda_1_e[ou]=lamda_1
    lamda_2_e[ou]=lamda_2
    p_e[ou]=p
  }
  res<-list(lamda_1_e=lamda_1_e,lamda_2_e=lamda_2_e,p_e=p_e)
  return(res)
}
k<-30
n<-15000
Xlamda1<-matrix(0,nrow = k,ncol = n)
Xlamda2<-matrix(0,nrow = k,ncol = n)
Xp<-matrix(0,nrow = k,ncol = n)
for(i in 1:k)
{
  p<-runif(1)
  mydata<-fmydata(i)
  a<-mysampling(p,n,mydata)
  Xlamda1[i,]<-a$lamda_1_e
  Xlamda2[i,]<-a$lamda_2_e
  Xp[i,]<-a$p_e
}
psilamda1<-t(apply(Xlamda1,1,cumsum))
psilamda2<-t(apply(Xlamda2,1,cumsum))
psip<-t(apply(Xp,1,cumsum))
for(i in 1:nrow(psilamda1))
{
  psilamda1[i,]<-psilamda1[i,]/(1:nrow(psilamda1))
}
print(gelmanrubin(psilamda1))
for(i in 1:nrow(psilamda2))
{
  psilamda2[i,]<-psilamda2[i,]/(1:nrow(psilamda2))
}
print(gelmanrubin(psilamda2))
for(i in 1:nrow(psip))
{
  psip[i,]<-psip[i,]/(1:nrow(psip))
}
print(gelmanrubin(psip))