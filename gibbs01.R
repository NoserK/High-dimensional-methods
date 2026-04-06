rm(list = ls(all = TRUE))
K=2
rmulti = function(p)
{
  K=length(p)
  random_variables=rmultinom(1,size=1,prob=p)
  out=0
    for(k in 1:K)
    {
      if(random_variables[k,1]==1) 
        out=k
    }
  return(out)
}
sample_mean_k=function(x,z)
{
  n_k=c(0,0)
  for(k in 1:K)
  {
    n_k[k]=length(z[z==k])
  }
  x_bar_k=matrix(0,nrow=K,ncol=1)
  x_var_k=matrix(0,nrow=K,ncol=1)
  for(k in 1:K)
  {
      x_bar_k[k]=(1/n_k[k])*sum(x[z==k])
      x_var_k[k]=var(x[z==k])
  }
  out=list(n_k=n_k,x_bar_k=x_bar_k,x_var_k=x_var_k)
  return(out)
}
gibbs_sampling_z=function(x,mu)
{
  d=ncol(x)
  z=rep(0,30)
  for(i in 1:30)
  {
    p=rep(0,length=K)
    for(k in 1:K)
    {
      p[k]=dpois(x[i],mu[k])
    }
    z[i]=rmulti(p=p/sum(p))
  }
  return(z)
}
gibbs_sampling_mu=function(x,z)
{
  sample_mean_k=sample_mean_k(x=x,z=z)
  n_k=sample_mean_k$n_k
  x_bar_k=sample_mean_k$x_bar_k
  x_var_k=sample_mean_k$x_var_k
  mu=matrix(0,nrow=2,ncol=1)
  for(k in 1:K)
  {
    mu[k,]=rgamma(1,(n_k[k]/(n_k[k]+1))*((x_bar_k[k,1])^2/x_var_k[k,1]),(n_k[k]/(n_k[k]+1))*(x_bar_k[k,1]/x_var_k[k,1]))
  }
  mu[is.nan(mu[,])]=0
  return(mu)
}
z_map=function(z)
{
  z_count=matrix(0,nrow=30,ncol=2)
  for(i in 1:30)
  {
    for(k in 1:K)
    {
      z_count[i,k]=length(z[i,z[i,]==k])
    }
  }
  z_map=rep(0,30)
  for(i in 1:30)
  {
    for(k in 1:K)
    {
      if(z_count[i,k]==max(z_count[i,])) 
        z_map[i]=k
    }
  }
  return(z_map)
}
poisson_mixture_model = function(x,iter_max=1000,burn_in=0)
{
  mu=array(0,dim=c(2,iter_max))
  z=array(0,dim=c(30,iter_max))
  cat("number of iteration is:",1,"\n")
  x_random=sample(1:30,size=2,replace=TRUE)
  for(k in 1:K)
  {
    mu[k,1]=x[x_random[k]]
    cat(mu[k,1],"\n")
  }
       z[,1]=gibbs_sampling_z(x=x,mu=mu[,1])
  cat("class is:", z[,1], "\n")
  # gibbs sampling
  for(s in 2:iter_max)
  {
    mu[,s]=gibbs_sampling_mu(x=x,z=z[,s-1])
    z[,s]=gibbs_sampling_z(x=x,mu=mu[,s])
    cat("number of iteration is:",s,"\n")
    cat("class is:",z[,s],"\n")
  }
  z_map = z_map(z=z)
  out=list(z=z[,burn_in:iter_max],z_map=z_map,mu=mu,x=x)
  return(out)
  message("the process has finished")
}
b1=rpois(10000,1)
b2=rpois(10000,10)
b3=matrix(0,nrow=10000,ncol=2)
for(i in 1:10000)
{
  if(runif(1)<=0.5)
  {
    b3[i,1]=b1[i]
    b3[i,2]=1
  }
  else
  {
    b3[i,1]=b2[i]
    b3[i,2]=2
  }
}
mydata=b3[1:30,]
x=as.matrix(mydata,1)
fit_gmm=poisson_mixture_model(x=x)
z_true=b3[1:30,2]
table(fit_gmm$z_map,z_true)