
############################## The second function: PGMM_ADMM() ———— using ADMM algorithm to get PGMM-CCEX estimation ############################

PGMM_ADMM<-function(N,T,QT_W,RT_W,wt,k,lambda,dotD){
  
  I=1000;zeta<-1 #hyperparameters
  
  positive<-function(x){
    if(x>0){
      x
    }
    else{
      0
    }
  }
  
  
  ST<-function(x,y){
    positive(1-y/norm(x,type="2"))%*%x
  }
  
  I=1000
  zeta<-1
  tri<-matrix(0,T,T-1)
  diag(tri[-T,])<- -1
  diag(tri[-1,])<- 1
  tri<-t(tri)
  tri<-kronecker(tri,diag(k+1)) #tri:(T-1,T)
  theta<-matrix(0,(T-1)*(k+1),1)
  
  for (t in (2:T)){
    theta[((t-2)*(k+1)+1):((t-1)*(k+1)),]<-dotD[((t-1)*(k+1)+1):(t*(k+1)),]-dotD[((t-2)*(k+1)+1):((t-1)*(k+1)),]
  }
  
  v<-matrix(0,(T-1)*(k+1),1)
  Zet<-matrix(0,(T-1)*(k+1),1)
  D<-matrix(0,T*(k+1),1) # D=delta
  
  r<-matrix(0,(T-1)*(k+1),1)
  
  #begin iteration
  
  A<-QT_W+zeta*crossprod(tri,tri)
  solveA<-solve(A)
  for (i in (1:I)){
    B<-RT_W+zeta*crossprod(tri,theta-zeta^(-1)*v)
    D<-solveA%*%B
    for (t in (2:T)){
      Zet[((t-2)*(k+1)+1):((t-1)*(k+1)),]<-D[((t-1)*(k+1)+1):(t*(k+1)),]-D[((t-2)*(k+1)+1):((t-1)*(k+1)),]+zeta^(-1)*v[((t-2)*(k+1)+1):((t-1)*(k+1)),]
      theta[((t-2)*(k+1)+1):((t-1)*(k+1)),]<-ST(Zet[((t-2)*(k+1)+1):((t-1)*(k+1)),],((lambda*wt[t])/zeta))
      v[((t-2)*(k+1)+1):((t-1)*(k+1)),]<-v[((t-2)*(k+1)+1):((t-1)*(k+1)),]+zeta*(D[((t-1)*(k+1)+1):(t*(k+1)),]-D[((t-2)*(k+1)+1):((t-1)*(k+1)),]-theta[((t-2)*(k+1)+1):((t-1)*(k+1)),])
    }
    #r<-crossprod(t(tri),D)-theta
    normr<-norm((crossprod(t(tri),D)-theta),type="2")
    if(normr<0.001){
      break;
      print(i)
    }
    out2<-list(D=D,theta=theta,i=i)
  }
  out2
}
