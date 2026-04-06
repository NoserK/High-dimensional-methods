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
normalchain<-function(sigma,N,X1)
{
  x<-rep(0,N)
  x[1]<-X1
  u<-runif(N)
  for(i in 2:N)
  {
    xt<-x[i-1]
    y<-rnorm(1,xt,sigma)
    r1<-dnorm(y,0,1)*dnorm(xt,y,sigma)
    r2<-dnorm(xt,0,1)*dnorm(y,xt,sigma)
    r<-r1/r2
    if(u[i]<=r)
    {
      x[i]<-y
    }
    else
    {
      x[i]<-xt
    }
  }
  return(x)
}
sigma<-0.2
k<-4
n<-15000
b<-1000
x0<-c(-10,-5,5,10)
X<-matrix(0,nrow = k,ncol = n)
for(i in 1:k)
{
  X[i,]<-normalchain(sigma,n,x0[i])
}
psi<-t(apply(X,1,cumsum))
for(i in 1:nrow(psi))
{
  psi[i,]<-psi[i,]/(1:nrow(psi))
}
print(gelmanrubin(psi))