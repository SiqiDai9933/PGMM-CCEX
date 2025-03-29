
# Input:XT,YT,N,T,r,SW,lambdamax,lambdamin,grids

# NT \times k dimentional XT=(X'_{.1},...,X'_{.T})';
# NT \times 1 dimentional YT=(Y'_{.1},...,Y'_{.T})';
# N: number of cross sections;
# T: number of time periods;
# r: use r_th order spatial lagged term of x to construct IVs. 
# r=1 or 2;
# SW: sptial weighting matrix (must be row-normalized);
# lambdamax: the maximum value of tuning parameter $\lambda$;
# lambdamin: the minimum value of tuning parameter $\lambda$;
# girds: the number of grids in [lambdamin,lambdamax].


# output:
# IC: the minimum value information criterion when \lambda in [lambdamin,lambdamax].
# D_hat: the post-PGMM-CCEX estimator; 
# se: standard error;
# t_value;
# p_value;


PGMM-CCEX<-function(XT,YT,N,T,r,SW,lambdamax,lambdamin,grids){
  
  source("initial_estimation.R")
  source("PGMM_ADMM.R") 
  source("post_estimation.R") 
  options(warn=-1)
  
  phol<-0.1*log(N*T)/(sqrt(N*T)) #hyperparameter $rho_NT$ in IC
  k<-ncol(XT) #The number of explanatory variables excluding spatial lagged term

  out1<-initial_estimation(XT,YT,N,T,k,r,SW)
  QT_W<-out1[,2:(T*(k+1)+1)]
  QT_W<-matrix(QT_W,T*(k+1),T*(k+1))
  RT_W<-out1[,1]
  RT_W<-matrix(RT_W,T*(k+1),1)
  wt<-out1[1:T,(T*(k+1)+2)]
  wt<-matrix(wt,T,1)
  dotD<-out1[,(T*(k+1)+3)] 
  dotD<-matrix(dotD,T*(k+1),1)
  dotDT<-t(matrix(dotD,k+1,T))
  
  
 
  ############## lambda=lambdamin ###################
  log_lambda_grid <- exp(seq(log(lambdamin), log(lambdamax), length.out = grids))
  ICC<-matrix(0,grids,1)
  
  lambda<-lambdamin
  lh<-1
  brea<-0

  
  theta<-PGMM_ADMM(N,T,QT_W,RT_W,wt,k,lambda,dotD)$theta
  thetak<-t(matrix(theta,(k+1),(T-1)))
  bn<-matrix(0,k,1)
  #count the number of breaks for each regressors
  for (p in (i:k)){
    bn[p,]=length(which(thetak[,p]!=0))
  }
  if((sum(bn))==0){
    brea<-0
    ICC[lh,1]=post_estimation(XT,YT,r,N,T,brea,k,SW)$IC
    for (llh in ((lh):grids)){
      ICC[llh,1]=ICC[lh,1]
    }
  }else{
    brea<-which(thetak[,1]!=0)
    ICC[lh,1]=post_estimation(XT,YT,r,N,T,brea,k,SW)$IC
  }
  
  ###################################################################

  for (lambda in log_lambda_grid[-1]){
    lh<-which.min(abs(log_lambda_grid - lambda)) 
    #the lh-th grid/lambda
    theta<-(PGMM_ADMM(N,T,QT_W,RT_W,wt,k,lambda,dotD))$theta
    thetak<-t(matrix(theta,(k+1),(T-1)))
    bn<-matrix(0,k,grids)
    #count the number of breaks for each regressors
    for (p in (i:k)){
      bn[p,lh]=length(which(thetak[,p]!=0))
    }
    if((sum(bn[,lh])==0)){
      brea<-0
      ICC[lh,1]=post_estimation(XT,YT,r,N,T,brea,k,SW)$IC
      for (llh in ((lh):grids)){
        ICC[llh,1]=ICC[lh,1]
      }
    }else{
      if (any(bn[,lh]==bn[,(lh-1)])){
        brea<-which(thetak[,1]!=0)
        ICC[lh,1]=ICC[(lh-1),1]
      }else{
      brea<-which(thetak[,1]!=0)
      ICC[lh,1]=post_estimation(XT,YT,r,N,T,brea,k,SW)$IC}
    }
  }
  
  index0<-which(ICC==min(ICC[ICC!=0]))[1]
  lambda0<-log_lambda_grid[index0]
  
  # we have chosen lambda0 now
  theta<-(PGMM_ADMM(N,T,QT_W,RT_W,wt,k,lambda0,dotD))$theta
  thetak<-t(matrix(theta,(k+1),(T-1)))
  brea<-which(thetak[,1]!=0)
  if (any(brea != 0)) {
    brea <- brea  
  } else {
    brea <- 0  
  }
  
  out<-post_estimation(XT,YT,r,N,T,brea,k,SW)
  
}







