
############################## The first function:initial_estimation() ———— to obatain initial estimatiors \dot{delta}_t ############################

initial_estimation<-function(XT,YT,N,T,k,r,SW){
  
  SWT<-bdiag(replicate(T,SW,simplify=FALSE))
  
  #XT=(X'_{.1},...,X'_{.T})'  #X=((X_1.)',...,(X_N.)')'
  X<-matrix(0,N*T,k)
  Y<-matrix(0,N*T,1)
  X_bar<-matrix(0,T,k) 
  ###Y_bar<-matrix(0,T,1)
  for(i in 1:N){
    select<-rep(0,T)
    for(t in 1:T){
      select[t]<-(t-1)*N+i 
      ##select {i,N+i,...,(T-1)*N+i}th row
      X_bar[t,]<-colMeans(XT[((t-1)*N+1):(t*N),])
    }###Y_bar[t,1]<-mean(Y[select])
    X[((i-1)*T+1):(i*T),]<-XT[select,]
    Y[((i-1)*T+1):(i*T),]<-YT[select,]
  }
  
  # construct defactorisation matrix (M_X_bar)
  I_T<-diag(T)
  I_N<-diag(N)
  M_X_bar<-I_T-tcrossprod(X_bar%*%ginv(crossprod(X_bar,X_bar)),X_bar)
  
  
  # choose (X,W^1*X,..,W^r*X) as Z (IV variables)
  if (r==1){
    Zr<-matrix(0,N*T,k)
    Zr<-kronecker(SW,I_T) %*% X
    #Z1<-Zr[,1:k] #Z2<-Zr[,k+1:2k] #Zr<-Zr[,(r-1)*k+1:rk]
    Z<-cbind(X,Zr)
  } 
  if (r==2){
    Zr1<-matrix(0,N*T,k)
    Zr1<-kronecker(SW,I_T) %*% X
    Zr2<-matrix(0,N*T,k)
    Zr2<-kronecker(SW%*%SW,I_T) %*% X
    #Z1<-Zr[,1:k] #Z2<-Zr[,k+1:2k] #Zr<-Zr[,(r-1)*k+1:rk]
    Z<-cbind(X,Zr1,Zr2)  
  }  
  
  #Z[N*T,(r+1)*k]  #iota= (r+1)*k
  
  # defactor Z
  ZT<-matrix(0,N*T,(r+1)*k) #ZT=(Z_1',...,Z_T')'
  diagZT<-matrix(0,N*T,T*(r+1)*k) #diagZ=diag(Z_1,...,Z_T)
  for (t in 1:T){
    b_t<-matrix(0,T,1) 
    b_t[t,]<-1
    ZT[((t-1)*N+1):(t*N),]<-diagZT[((t-1)*N+1):(t*N),((t-1)*(r+1)*k+1):(t*(r+1)*k)]<-kronecker(I_N,crossprod(b_t,M_X_bar))%*%Z
    #Z_1<-Z[1:N,1:iota]  #Z_2<-Z[N+1:2N,iota+1:2*iota] #Z_T<-Z[(T-1)*N+1:TN,(T-1)*iota+1:T*iota]
  }
  
  # construct X^{\ast}_{.t}=(SW*y_{.t},X_{.t})
  Xast<-cbind(SWT%*%YT,XT) #(N*T,1+k)
  
  # get \dot(\delta_t)_W
  
  # first step: GMM weight=I_N #\dot(\delta_t)_I
  QT_I<-matrix(0,T*(k+1),T*(k+1))
  RT_I<-matrix(0,T*(k+1),1)
  
  dotDT_I<-matrix(0,T*(k+1),1) 
  for (t in 1:T){
    QT_I[((t-1)*(k+1)+1):(t*(k+1)),((t-1)*(k+1)+1):(t*(k+1))]<-matrix(crossprod((1/N)*crossprod(ZT[((t-1)*N+1):(t*N),],Xast[((t-1)*N+1):(t*N),]),(1/N)*crossprod(ZT[((t-1)*N+1):(t*N),],Xast[((t-1)*N+1):(t*N),])),k+1,k+1) ##crossprod(x,)=t(x)%*%x
    RT_I[((t-1)*(k+1)+1):(t*(k+1)),]<-matrix(crossprod((1/N)*crossprod(ZT[((t-1)*N+1):(t*N),],Xast[((t-1)*N+1):(t*N),]),(1/N)*crossprod(ZT[((t-1)*N+1):(t*N),],YT[((t-1)*N+1):(t*N),])),k+1,1)
  }
  dotDT_I<-solve(QT_I)%*%RT_I
  dotDDT_I<-t(matrix(dotDT_I,k+1,T))
  dotDDNT_I<-dotDDT_I[rep(1:nrow(dotDDT_I),each=N),] #dotDDNT:(N*T,k+1)
  
  # constructe time-varing GMM weight W_t #=\tilde(V)_t=(1/N)\sum_i^N{g_it(\dot{\delta_t})g'_it}
  
  #g_it((r+1)*k,1)=z_{it}[(y_{it}-x_{it}'^{\ast}{\delta}_{t}]  
  #z_it:((r+1)*k,1) #xast_it:(k+1,1)
  
  gNT<-matrix(0,(r+1)*k,N*T) #gNT=(g_11,...,g_N1,g_12,...,g_N2,...,g_1T,...,g_NT)
  #W_t((r+1)*k,(r+1)*k)=(1/N)sum_(i=1)^N(g_it%*%g'_it)=sum_(i=1)^N tcrossprod(g_it)/N
  g2NT<-matrix(0,(r+1)*k,(r+1)*k*N*T)
  for (s in 1:(N*T)){
    gNT[,s]<-matrix(ZT[s,])%*%(YT[s,]-tcrossprod(matrix(Xast[s,],1,(k+1)),matrix(dotDDNT_I[s,],1,(k+1)))) #tcrossprod(x,)=x%*%t(x)
    #ZT[s,]=t(z_it):(1,(r+1)*k) #Xast[s,]=t(xast_it):(1,k+1) #dotDDNT[s,]=t({\delta}_{t}):(1,k+1)
    g2NT[,((s-1)*(r+1)*k+1):(s*(r+1)*k)]<-tcrossprod(gNT[,s],gNT[,s]) #g_it%*%g_it'
  }
  
  dim(g2NT)<-c(((r+1)*k),((r+1)*k*N),T) ## split g2NT to T ((r+1)*k,(r+1)*k*N) matrixs
  ##g2NT[,,1]=(g_11g_11',...,g_N1g_N1'):(r*k,r*k*N) #g2NT[,,T]=(g_T1g_T1',...,g_NTg_NT')
  
  diagV<-matrix(0,T*(r+1)*k,T*(r+1)*k)
  diagW<-matrix(0,T*(r+1)*k,T*(r+1)*k)
  for (t in 1:T){
    for (i in 1:N){
      diagV[((t-1)*(r+1)*k+1):(t*(r+1)*k),((t-1)*(r+1)*k+1):(t*(r+1)*k)]<-diagV[((t-1)*(r+1)*k+1):(t*(r+1)*k),((t-1)*(r+1)*k+1):(t*(r+1)*k)]+g2NT[,((i-1)*(r+1)*k+1):(i*(r+1)*k),t]
      #t=1,i=1#diagW[1:(r+1)*k,1:(r+1)*k]<-0+g2NT[,1:(r+1)*k,1]=g_11g_11'
      #t=1,i=2#diagW[1:(r+1)*k,1:(r+1)*k]<-(t=1,i=1)+g2NT[,(r+1)*k+1:2*(r+1)*k,1]=g_11g_11'+g_21g_21'
      #t=1,i=N#diagW[1:(r+1)*k,1:(r+1)*k]<-(t=1,i=1)+(t=1,i=2)+...+g2NT[,(N-1)*(r+1)*k+1:N*(r+1)*k,1]=g_11g_11'+g_21g_21'+...+g_N1g_N1'
      #......
      #t=T,i=1#diagW[(T-1)*(r+1)*k+1:T*(r+1)*k,(T-1)*(r+1)*k+1:T*(r+1)*k]<-0+g2NT[,1:(r+1)*k,T]=g_1Tg_1T'
    }
  }
  ## W_t<-diagW[((t-1)*(r+1)*k)+1:t*(r+1)*k,((t-1)*(r+1)*k)+1:t*(r+1)*k]
  diagV<-diagV/N
  diagW<-solve(diagV)
  #diagW<-diag(T*(r+1)*k)
  # sencond step 
  
  QT_W<-matrix(0,T*(k+1),T*(k+1))
  RT_W<-matrix(0,T*(k+1),1)
  DT_W<-matrix(0,T*(k+1),1)
  for (t in 1:T){
    QT_W[((t-1)*(k+1)+1):(t*(k+1)),((t-1)*(k+1)+1):(t*(k+1))]<-matrix(crossprod((1/N)*crossprod(ZT[((t-1)*N+1):(t*N),],Xast[((t-1)*N+1):(t*N),]),diagW[((t-1)*(r+1)*k+1):(t*(r+1)*k),((t-1)*(r+1)*k+1):(t*(r+1)*k)]%*%((1/N)*crossprod(ZT[((t-1)*N+1):(t*N),],Xast[((t-1)*N+1):(t*N),]))),k+1,k+1) 
    RT_W[((t-1)*(k+1)+1):(t*(k+1)),]<-matrix(crossprod((1/N)*crossprod(ZT[((t-1)*N+1):(t*N),],Xast[((t-1)*N+1):(t*N),]),diagW[((t-1)*(r+1)*k+1):(t*(r+1)*k),((t-1)*(r+1)*k+1):(t*(r+1)*k)]%*%((1/N)*crossprod(ZT[((t-1)*N+1):(t*N),],YT[((t-1)*N+1):(t*N),]))),k+1,1)
  }
  dotDT_W<-ginv(QT_W)%*%RT_W ##dotDT=(\dot{\delta'_1},...,\dot{\delta'_T})':(T*(k+1),1)
  
  # construct data-driven weight \dot(w)
  
  dotDDT_W<-t(matrix(dotDT_W,k+1,T))
  dotw<-matrix(0,T,1)
  for (t in 2:T){
    dotw[t,]<-(norm(dotDDT_W[t,]-dotDDT_W[t-1,],type="2"))^(-2) # kappa=2
  }
  dotw[1,]=1
  wt<-matrix(0,T*(k+1),1)
  wt[1:T,]<-dotw
  out1<-cbind(RT_W,QT_W,wt,dotDT_W)
}
