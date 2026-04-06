library(MASS)
##mvrnorm()
##瀹氫箟涓€涓骇鐢熷鍏冩鎬佸垎甯冪殑闅忔満鍚戦噺鍗忔柟宸煩闃?
Simu_Multi_Norm<-function(x_len,sd=1,pho=0.5)
{
  #鍒濆鍖栧崗鏂瑰樊鐭╅樀
  V<-matrix(data =NA,nrow=x_len,ncol=x_len)
  #mean鍙妔d鍒嗗埆涓洪殢鏈哄悜閲弜鐨勫潎鍊煎拰鏂瑰樊
  #瀵瑰崗鏂瑰樊鐭╅樀杩涜璧嬪€紁ho(i,j) = pho^|i-j|
  for(i in 1:x_len)
  { ##閬嶅巻姣忎竴琛?
    for(j in 1:x_len)
    { ##閬嶅巻姣忎竴鍒?
      V[i,j]<-pho^abs(i-j)
    }
  }
  V<-(sd^2)*V
  return(V)
}
##浜х敓妯℃嫙鏁板€艰嚜鍙橀噺X
set.seed(123)
X<-mvrnorm(n=200,mu=rep(0,10),Simu_Multi_Norm(x_len=10,sd =1,pho=0.5))
##浜х敓妯℃嫙鏁板€硷細鍝嶅簲鍙橀噺y
beta<-c(1,2,0,0,3,0,0,0,-2,0)
#alpha<-0
#prob<-exp(alpha X %*% beta)/(1 exp(alpha X %*% beta))
prob<-exp(X%*%beta)/(1+exp(X%*%beta))
y<-rbinom(n=200,size=1,p=prob)
##浜х敓model matrix
mydata<-data.frame(X=X,y=y)
#X<-model.matrix(y~., data = mydata)
##鍖呭惈鎴煩椤圭殑绯绘暟
#b_real<-c(alpha,beta)
b_real<-beta
########----瀹氫箟鎯╃綒椤圭浉鍏崇殑鍑芥暟-----------------
##瀹氫箟鎯╃綒椤?
####杩愯鍙戠幇锛岃嫢lambda璁剧疆涓?2锛屽垯绯绘暟鍏ㄨ鍘嬬缉涓?0.
####鏈▼搴忔牴鎹畆cvreg鐢–V閫夊嚭鏉ョ殑lambda璁剧疆涓€涓緝涓哄悎鐞嗙殑lambda銆?
p_lambda<-function(theta,lambda=0.025)
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
##瀹氫箟鎯╃綒椤瑰鏁?
p_lambda_d<-function(theta,a=3.7,lambda=0.025)
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
# ##褰揵eta_j0涓嶇瓑浜?0锛屽畾涔夋儵缃氶」瀵兼暟杩戜技
# p_lambda_d_apro<-function(beta_j0,beta_j,a = 3.7, lambda = 2){
#   return(beta_j * p_lambda_d(beta = beta_j0,a = a, lambda = lambda)/abs(beta_j0))
# }
# ##褰揵eta_j0 涓嶇瓑浜?0锛屾寚瀹氳繎浼兼儵缃氶」,浣跨敤娉板嫆灞曞紑閫艰繎
# p_lambda_apro<-function(beta_j0,beta_j,a = 3.7, lambda = 2){
#   if(abs(beta_j0)< 1e-16){
#     return(0)
#   }else{
#     p_lambda<-p_lambda(theta = beta_j0, lambda = lambda)
#       0.5 * (beta_j^2 - beta_j0^2) * p_lambda_d(theta = beta_j0, a = a, lambda = lambda)/abs(beta_j0)
#   }
# }
#define the log-likelihood function
loglikelihood_SCAD<-function(X,y,b)
{
  linear_comb<-as.vector(X%*%b)
  ll<-sum(y*linear_comb)+sum(log(1/(1+exp(linear_comb))))-nrow(X)*sum(p_lambda(theta=b))
  return (ll)
}
##鍒濆鍖栫郴鏁?
#b0<-rep(0,length(b_real))

#b0<- b_real rnorm(length(b_real), mean = 0, sd = 0.1)

##灏嗘棤鎯╃綒鏃剁殑浼樺寲缁撴灉浣滀负鍒濆鍊?
b.best_GS<-solve(t(X)%*%X)%*%t(X)%*%y
b0<-b.best_GS
##b1鐢ㄤ簬璁板綍鏇存柊绯绘暟
b1<-b0
##b.best鐢ㄤ簬瀛樻斁鍘嗗彶鏈€澶т技鐒跺€煎搴旂郴鏁?
b.best_SCAD<-b0
# the initial value of loglikelihood
ll.old<-loglikelihood_SCAD(X=X,y=y,b=b0)
# initialize the difference between the two steps of theta
diff<-1  
#record the number of iterations
iter<-0
#set the threshold to stop iterations
epsi<-1e-10
#the maximum iterations  
max_iter<-100000
#鍒濆鍖栦竴涓垪琛ㄧ敤浜庡瓨鏀炬瘡涓€娆¤凯浠ｇ殑绯绘暟缁撴灉
b_history<-list(data.frame(b0))
#鍒濆鍖栧垪琛ㄧ敤浜庡瓨鏀句技鐒跺€?
ll_list<-list(ll.old)
#######-------SCAD杩唬---------
while(diff>epsi&iter<max_iter)
{
  for(j in 1:length(b_real))
  {
    if(abs(b0[j])<1e-06)
    {
      next()
    }
    else
    {
      #绾挎€ч儴鍒?
      linear_comb<-as.vector(X%*%b0)
      #鍒嗗瓙
      nominator<-sum(y*X[,j]-X[,j]*exp(linear_comb)/(1+exp(linear_comb)))+nrow(X)*b0[j]*p_lambda_d(theta = b0[j])/abs(b0[j])
      #鍒嗘瘝,鍗充簩闃跺閮ㄥ垎
      denominator=-sum(X[,j]^2*exp(linear_comb)/(1+exp(linear_comb))^2)+nrow(X)*p_lambda_d(theta = b0[j])/abs(b0[j])
      #2-(3) :鏇存柊b0[j]
      b0[j]<-b0[j]-nominator/denominator
      #2-(4)
      if(abs(b0[j])<1e-06)
      {
        b0[j]<-0
      }
      #       #鏇存柊浼肩劧鍊?
      #       ll.new<- loglikelihood_SCAD(X = X, y = y, b = b0)
      #       
      #       
      #       
      #       #鑻ヤ技鐒跺€兼湁鎵€澧炲姞锛屽垯灏嗗綋鍓嶇郴鏁颁繚瀛?
      #       if(ll.new > ll.old){
      #         #鏇存柊绯绘暟
      #         b.best_SCAD[j]<-b0[j]
      #       }
      #       
      #       #姹傚樊寮?
      #       diff<- abs((ll.new - ll.old)/ll.old)
      #       ll.old <- ll.new
      #       iter<- iter 1
      #       b_history[[iter]]<-data.frame(b0)
      #       ll_list[[iter]]<-ll.old
      #       ##褰撹揪鍒板仠姝㈡潯浠舵椂锛岃烦鍑哄惊鐜?
      #       if(diff < epsi){
      #         break
      #       }
      #       
      
    }
  }
  #鏇存柊浼肩劧鍊?
  ll.new<- loglikelihood_SCAD(X=X,y=y,b=b0)
  #鑻ヤ技鐒跺€兼湁鎵€澧炲姞锛屽垯灏嗗綋鍓嶇郴鏁颁繚瀛?
  if(ll.new>ll.old)
  {
    #鏇存柊绯绘暟
    b.best_SCAD<-b0
  }
  #姹傚樊寮?
  diff<-abs((ll.new-ll.old)/ll.old)
  ll.old<-ll.new
  iter<-iter+1
  b_history[[iter]]<-data.frame(b0)
  ll_list[[iter]]<-ll.old
}

b_hist<-do.call(rbind,b_history)
#b_hist
ll_hist<-do.call(rbind,ll_list)
#ll_hist
#
iter
##
ll.best<-max(ll_hist)
ll.best
##
b.best_SCAD
##瀵规瘮
cbind(coeff_glm,b.best,b.best_SCAD,b_real)
##----------ncvreg楠岃瘉-----------
library(ncvreg)
my_ncvreg<-ncvreg(X,y,family = c("binomial"),penalty = c("SCAD"),lambda = 2)
my_ncvreg$beta
my_ncvreg<-ncvreg(X,y,family = c("binomial"),penalty = c("SCAD"))
summary(my_ncvreg)
my_ncvreg$beta
###鐢╟v鎵炬渶浼樼殑lambda
scad_cv<-cv.ncvreg(X,y,family = c("binomial"),penalty='SCAD')
scad_cv$lambda.min
mySCAD=ncvreg(X,y,family = c("binomial"),penalty='SCAD',lambda=scad_cv$lambda.min)
summary(mySCAD)
ncv_SCAD<-mySCAD$beta[-1]