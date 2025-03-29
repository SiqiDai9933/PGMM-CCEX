
################# The third function:post_estimation() ———— calculate post-PGMM-CCEX estimators and IC under various tuning parameter $\lambda$ #############

post_estimation<-function(XT,YT,r,N,T,brea,k,SW){ 
  
  
  SWT<-bdiag(replicate(T,SW,simplify=FALSE))
  
################# brea: the break dates ######################
# thetak<-t(matrix(theta,(k+1),(T-1)))
# brea<-which(thetak[,1]!=0)
# e.g. brea=2, the structural break occurs in the 3-th period, period 1: t \in {1,2}.
  
  
  if (brea[1]==0){
    breaks<-c(0,T) 
    m_Tilta<-length(breaks)-1
  }else{
    breaks<-c(0,brea,T) # when t=m_Tilta (breaks[m_tilta+1]-1))=T
    m_Tilta<-length(breaks)-1  
  } 
# e.g. brea=c(2,9,15),breaks=(0,2,9,15,T),m_Tilta=4——four period:1:2,3:9,10:15,16:T

  
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
    }
    X[((i-1)*T+1):(i*T),]<-XT[select,]
    Y[((i-1)*T+1):(i*T),]<-YT[select,]
  }
  
  # construct defactorisation matrix (M_X_bar)
  
  I_T<-diag(T)
  M_X_bar<-I_T-tcrossprod(X_bar%*%solve(crossprod(X_bar,X_bar)),X_bar)#CCEX
  I_N<-diag(N)
  
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
  
  ZT<-matrix(0,N*T,(r+1)*k) #ZT=(Z_1',...,Z_T')'
  diagZT<-matrix(0,N*T,T*(r+1)*k) #r*k=iota #diagZ=diag(Z_1,...,Z_T)
  for (t in 1:T){
    b_t<-matrix(0,T,1) 
    b_t[t,]<-1
    ZT[((t-1)*N+1):(t*N),]<-diagZT[((t-1)*N+1):(t*N),((t-1)*(r+1)*k+1):(t*(r+1)*k)]<-kronecker(I_N,crossprod(b_t,M_X_bar))%*%Z
    #Z_1<-Z[1:N,1:iota]  #Z_2<-Z[N+1:2N,iota+1:2*iota] #Z_T<-Z[(T-1)*N+1:TN,(T-1)*iota+1:T*iota]
  }
  
  Xast<-cbind(SWT%*%YT,XT) #(N*T,1+k)
  
  
  QTP_I<-matrix(0,m_Tilta*(k+1),m_Tilta*(k+1))
  RTP_I<-matrix(0,m_Tilta*(k+1),1)
  DTP_I<-matrix(0,m_Tilta*(k+1),1)
  for (m_tilta in (1:m_Tilta)){ 
    #m_tilta in c(1,2,3,4): the m_tilta-th break and final point
    for (t in ((breaks[m_tilta]+1):(breaks[m_tilta+1]))){
      for (tt in ((breaks[m_tilta]+1):(breaks[m_tilta+1]))){
        #t=1:2,3:9,10:15,16:T
        QTP_I[((m_tilta-1)*(k+1)+1):(m_tilta*(k+1)),((m_tilta-1)*(k+1)+1):(m_tilta*(k+1))]<- QTP_I[((m_tilta-1)*(k+1)+1):(m_tilta*(k+1)),((m_tilta-1)*(k+1)+1):(m_tilta*(k+1))]+matrix(crossprod((1/N)*crossprod(ZT[((t-1)*N+1):(t*N),],Xast[((t-1)*N+1):(t*N),]),(1/N)*crossprod(ZT[((tt-1)*N+1):(tt*N),],Xast[((tt-1)*N+1):(tt*N),])),k+1,k+1) ##crossprod(x,)=t(x)%*%x
        RTP_I[((m_tilta-1)*(k+1)+1):(m_tilta*(k+1)),]<-RTP_I[((m_tilta-1)*(k+1)+1):(m_tilta*(k+1)),]+matrix(crossprod((1/N)*crossprod(ZT[((t-1)*N+1):(t*N),],Xast[((t-1)*N+1):(t*N),]),(1/N)*crossprod(ZT[((tt-1)*N+1):(tt*N),],YT[((tt-1)*N+1):(tt*N),])),k+1,1)
      }
    }
  }
  
  
  dotDTP_I<-solve(QTP_I)%*%RTP_I
  
  dotDDMP_I<-t(matrix(dotDTP_I,k+1,m_Tilta)) #dotDDTP:trans (m_tilta,k+1) to (T,k+1)
  re<-matrix(0,(m_Tilta+1),1) #m_Tilta periods + 0
  dotDDTP_I<-matrix(0,T,k+1)
  for (m_tilta in (1:m_Tilta)){ 
    re[(m_tilta+1),]<-breaks[m_tilta+1]-breaks[m_tilta]
    a<-breaks[m_tilta]+1
    b<-breaks[(m_tilta+1)]
    dotDDTP_I[a:b,]<-dotDDMP_I[rep(m_tilta,each=(re[m_tilta+1,])),]
  }
  dotDDNTP_I<-dotDDTP_I[rep(1:nrow(dotDDTP_I),each=N),] #dotDDNTP:trans (N*m_tilta,k+1) to (N*T,k+1)
  
  gNT<-matrix(0,(r+1)*k,N*T) #gNT=(g_11,...,g_N1,g_12,...,g_N2,...,g_1T,...,g_NT)
  #W_t((r+1)*k,(r+1)*k)=(1/N)sum_(i=1)^N(g_it%*%g'_it)=sum_(i=1)^N tcrossprod(g_it)/N
  g2NT<-matrix(0,(r+1)*k,(r+1)*k*N*T)
  for (s in 1:(N*T)){
    gNT[,s]<-matrix(ZT[s,])%*%(YT[s,]-tcrossprod(matrix(Xast[s,],1,(k+1)),matrix(dotDDNTP_I[s,],1,(k+1)))) #tcrossprod(x,)=x%*%t(x)
    #ZT[s,]=t(z_it):(1,(r+1)*k) #Xast[s,]=t(xast_it):(1,k+1) #dotDDNT[s,]=t({\delta}_{t}):(1,k+1)
    g2NT[,((s-1)*(r+1)*k+1):(s*(r+1)*k)]<-tcrossprod(gNT[,s],gNT[,s]) #g_it%*%g_it'
  }
  
  dim(g2NT)<-c(((r+1)*k),((r+1)*k*N),T) ## split g2NT to T ((r+1)*k,(r+1)*k*N) matrixs
  ##g2NT[,,1]=(g_11g_11',...,g_N1g_N1'):(r*k,r*k*N) #g2NT[,,T]=(g_T1g_T1',...,g_NTg_NT')
  
  diagV<-matrix(0,m_Tilta*(r+1)*k,m_Tilta*(r+1)*k) # W_1,..,W_m_Tilda
  #diagVV<-matrix(0,m_Tilta*(r+1)*k,m_Tilta*(r+1)*k)
  for (m_tilta in (1:m_Tilta)){ 
    for (t in ((breaks[m_tilta]+1):(breaks[m_tilta+1]))){ # accumulate residual both in t and i
      for (i in 1:N){
        diagV[((m_tilta-1)*(r+1)*k+1):(m_tilta*(r+1)*k),((m_tilta-1)*(r+1)*k+1):(m_tilta*(r+1)*k)]<-diagV[((m_tilta-1)*(r+1)*k+1):(m_tilta*(r+1)*k),((m_tilta-1)*(r+1)*k+1):(m_tilta*(r+1)*k)]+g2NT[,((i-1)*(r+1)*k+1):(i*(r+1)*k),t]
      }
    }
    #diagVV[((m_tilta-1)*(r+1)*k+1):(m_tilta*(r+1)*k),((m_tilta-1)*(r+1)*k+1):(m_tilta*(r+1)*k)]<-(1/(breaks[m_tilta+1]-breaks[m_tilta]))*diagV[((m_tilta-1)*(r+1)*k+1):(m_tilta*(r+1)*k),((m_tilta-1)*(r+1)*k+1):(m_tilta*(r+1)*k)]
  }
  ## W_m_tilta<-diagW[((m_tilta-1)*(r+1)*k)+1:m_tilta*(r+1)*k,((m_tilta-1)*(r+1)*k)+1:m_tilta*(r+1)*k]
  
  diagW<-matrix(0,T*(r+1)*k,T*(r+1)*k)
  reW<-matrix(0,(m_Tilta+1),1) #trans (m_tilta*(r+1)*k,m_tilta*(r+1)*k) to (T*(r+1)*k,T*(r+1)*k)
  for (m_tilta in (1:m_Tilta)){
    reW[m_tilta+1,]<-breaks[m_tilta+1]-breaks[m_tilta]
    diagW_M_tilta<-solve(diagV/(N*reW[m_tilta+1,]))
    diagW[(((breaks[m_tilta])*(r+1)*k)+1):((breaks[(m_tilta+1)])*(r+1)*k),(((breaks[m_tilta])*(r+1)*k)+1):((breaks[(m_tilta+1)])*(r+1)*k)]<-matrix(bdiag(replicate(reW[m_tilta+1,],diagW_M_tilta[((r+1)*k*(m_tilta-1)+1):((r+1)*k*m_tilta),((r+1)*k*(m_tilta-1)+1):((r+1)*k*m_tilta)],simplify=FALSE)),(reW[(m_tilta+1),]*(r+1)*k),(reW[(m_tilta+1),]*(r+1)*k))
  }
  
  
  # sencond step 
  
  QTP_W<-matrix(0, m_Tilta*(k+1), m_Tilta*(k+1))
  RTP_W<-matrix(0, m_Tilta*(k+1),1)
  DTP_W<-matrix(0, m_Tilta*(k+1),1)
  
  for (m_tilta in (1:m_Tilta)){ #m_tilta in c(1,2,3,4): the m_tilta-th break and final point
    for (t in ((breaks[m_tilta]+1):(breaks[m_tilta+1]))){
      for (tt in ((breaks[m_tilta]+1):(breaks[m_tilta+1]))){
        QTP_W[((m_tilta-1)*(k+1)+1):(m_tilta*(k+1)),((m_tilta-1)*(k+1)+1):(m_tilta*(k+1))]<- QTP_W[((m_tilta-1)*(k+1)+1):(m_tilta*(k+1)),((m_tilta-1)*(k+1)+1):(m_tilta*(k+1))]+matrix(crossprod((1/N)*crossprod(ZT[((t-1)*N+1):(t*N),],Xast[((t-1)*N+1):(t*N),]),diagW[((t-1)*(r+1)*k+1):(t*(r+1)*k),((t-1)*(r+1)*k+1):(t*(r+1)*k)]%*%((1/N)*crossprod(ZT[((tt-1)*N+1):(tt*N),],Xast[((tt-1)*N+1):(tt*N),]))),k+1,k+1)
        RTP_W[((m_tilta-1)*(k+1)+1):(m_tilta*(k+1)),]<-RTP_W[((m_tilta-1)*(k+1)+1):(m_tilta*(k+1)),]+matrix(crossprod((1/N)*crossprod(ZT[((t-1)*N+1):(t*N),],Xast[((t-1)*N+1):(t*N),]),diagW[((t-1)*(r+1)*k+1):(t*(r+1)*k),((t-1)*(r+1)*k+1):(t*(r+1)*k)]%*%((1/N)*crossprod(ZT[((tt-1)*N+1):(tt*N),],YT[((tt-1)*N+1):(tt*N),]))),k+1,1)
      }
    }
  }
  dotDTP_W<-solve(QTP_W)%*%RTP_W
  dotDDMP_W<-t(matrix(dotDTP_W,k+1,m_Tilta)) #post-PGMM-CCEX estimator

########################################## calculate IC and SE #################################################
  
  dotDDTP_W<-matrix(0,T,k+1)
  for (m_tilta in (1:m_Tilta)){
    reW[m_tilta+1,]<-breaks[m_tilta+1]-breaks[m_tilta]
    a=breaks[m_tilta]+1
    b=breaks[(m_tilta+1)]
    dotDDTP_W[a:b,]<-dotDDMP_W[rep(m_tilta,each=(re[m_tilta+1,])),]
  }
  dotDDNTP_W<-dotDDTP_W[rep(1:nrow(dotDDTP_W),each=N),] #dotDDNTP:(N*m_tilta,k+1)
  
  
  gNTP<-matrix(0,(r+1)*k,N*T) #gNT=(g_11,...,g_N1,g_12,...,g_N2,...,g_1T,...,g_NT)
  g2NTP_1<-matrix(0,1,N*T)
  g2NTP_2<-matrix(0,(r+1)*k,(r+1)*k*N*T)
  for (s in (1:(N*T))){
    gNTP[,s]<-matrix(ZT[s,])%*%(YT[s,]-tcrossprod(matrix(Xast[s,],1,(k+1)),matrix(dotDDNTP_W[s,],1,(k+1)))) #tcrossprod(x,)=x%*%t(x)
    g2NTP_1[,s]<-crossprod(matrix(gNTP[,s]))#g_it'%*%g_it
    g2NTP_2[,((s-1)*(r+1)*k+1):(s*(r+1)*k)]<-tcrossprod(gNTP[,s],gNTP[,s]) #g_it%*%g_it'
  }
  
  sigma2<-sum(g2NTP_1)/(N*T)
  IC=log(sigma2)+phol*(k+1)*(m_Tilta)
  
  dim(g2NTP_2)<-c(((r+1)*k),((r+1)*k*N),T) ## split g2NT_2 to T ((r+1)*k,(r+1)*k*N) matrixs
  ##g2NT_2[,,1]=(g_11g_11',...,g_N1g_N1'):(r*k,r*k*N) #g2NT_2[,,T]=(g_T1g_T1',...,g_NTg_NT')
  
  
  diagVP<-matrix(0,m_Tilta*(r+1)*k,m_Tilta*(r+1)*k) # W_1,..,W_m_Tilda
  for (m_tilta in (1:m_Tilta)){ 
    for (t in ((breaks[m_tilta]+1):(breaks[m_tilta+1]))){ # accumulate residual both in t and i
      for (i in 1:N){
        diagVP[((m_tilta-1)*(r+1)*k+1):(m_tilta*(r+1)*k),((m_tilta-1)*(r+1)*k+1):(m_tilta*(r+1)*k)]<-diagVP[((m_tilta-1)*(r+1)*k+1):(m_tilta*(r+1)*k),((m_tilta-1)*(r+1)*k+1):(m_tilta*(r+1)*k)]+g2NTP_2[,((i-1)*(r+1)*k+1):(i*(r+1)*k),t]
      }
    }
  }
  
  
  diagWP<-matrix(0,T*(r+1)*k,T*(r+1)*k)
  reW<-matrix(0,(m_Tilta+1),1) #trans (m_tilta*(r+1)*k,m_tilta*(r+1)*k) to (T*(r+1)*k,T*(r+1)*k)
  for (m_tilta in (1:m_Tilta)){
    reW[m_tilta+1,]<-breaks[m_tilta+1]-breaks[m_tilta]
    diagWP_M_tilta<-solve(diagVP/(N*reW[m_tilta+1,]))
    diagWP[(((breaks[m_tilta])*(r+1)*k)+1):((breaks[(m_tilta+1)])*(r+1)*k),(((breaks[m_tilta])*(r+1)*k)+1):((breaks[(m_tilta+1)])*(r+1)*k)]<-matrix(bdiag(replicate(reW[m_tilta+1,],diagWP_M_tilta[((r+1)*k*(m_tilta-1)+1):((r+1)*k*m_tilta),((r+1)*k*(m_tilta-1)+1):((r+1)*k*m_tilta)],simplify=FALSE)),(reW[(m_tilta+1),]*(r+1)*k),(reW[(m_tilta+1),]*(r+1)*k))
  }
  
  
  DVD<-matrix(0, m_Tilta*(k+1), m_Tilta*(k+1))
  for (m_tilta in (1:m_Tilta)){ #m_tilta in c(1,2,3,4): the m_tilta-th break and final point
    for (t in ((breaks[m_tilta]+1):(breaks[m_tilta+1]))){
      for (tt in ((breaks[m_tilta]+1):(breaks[m_tilta+1]))){
        DVD[((m_tilta-1)*(k+1)+1):(m_tilta*(k+1)),((m_tilta-1)*(k+1)+1):(m_tilta*(k+1))]<- DVD[((m_tilta-1)*(k+1)+1):(m_tilta*(k+1)),
                                                                                               ((m_tilta-1)*(k+1)+1):(m_tilta*(k+1))]+matrix(crossprod((1/N)*crossprod(ZT[((t-1)*N+1):(t*N),],Xast[((t-1)*N+1):(t*N),]),diagWP[((t-1)*(r+1)*k+1):(t*(r+1)*k),
                                                                                                                                                                                                                               ((t-1)*(r+1)*k+1):(t*(r+1)*k)]%*%((1/N)*crossprod(ZT[((tt-1)*N+1):(tt*N),],Xast[((tt-1)*N+1):(tt*N),]))),k+1,k+1)
      }
    }
  }
  
  
  solveWP<-solve(diagWP)
  DWVWD<-matrix(0, m_Tilta*(k+1), m_Tilta*(k+1)) 
  for (m_tilta in (1:m_Tilta)){ #m_tilta in c(1,2,3,4): the m_tilta-th break and final point
    for (t in ((breaks[m_tilta]+1):(breaks[m_tilta+1]))){
      for (tt in ((breaks[m_tilta]+1):(breaks[m_tilta+1]))){
        DWVWD[((m_tilta-1)*(k+1)+1):(m_tilta*(k+1)),((m_tilta-1)*(k+1)+1):(m_tilta*(k+1))]<- DWVWD[((m_tilta-1)*(k+1)+1):(m_tilta*(k+1)),
                                                                                                   ((m_tilta-1)*(k+1)+1):(m_tilta*(k+1))]+matrix(crossprod((1/N)*crossprod(ZT[((t-1)*N+1):(t*N),],Xast[((t-1)*N+1):(t*N),]),
                                                                                                                                                           diagW[((t-1)*(r+1)*k+1):(t*(r+1)*k),((t-1)*(r+1)*k+1):(t*(r+1)*k)]%*%solveWP[((t-1)*(r+1)*k+1):(t*(r+1)*k),((t-1)*(r+1)*k+1):(t*(r+1)*k)]%*%
                                                                                                                                                             diagW[((t-1)*(r+1)*k+1):(t*(r+1)*k),((t-1)*(r+1)*k+1):(t*(r+1)*k)]%*%((1/N)*crossprod(ZT[((tt-1)*N+1):(tt*N),],Xast[((tt-1)*N+1):(tt*N),]))),k+1,k+1)
      }
    }
  }
  
  if(length(reW)==2){
    Lj<-kronecker((reW[-1,]),diag(k+1))
  }else{
    Lj<-kronecker(diag(reW[-1,]),diag(k+1))
  }
  
  
  Omega<-Lj%*%solve(QTP_W)%*%DWVWD%*%solve(t(QTP_W))%*%Lj
  
############################################### t test ################################################################
  
  if (brea[1]==0){
    period<-1
    pl<-T
  }else{
    mhat<-length(brea)
    period<-mhat+1
    pl<-period_length<-diff(breaks)
  }
  p_value<-se<-t_value<-matrix(0,period,(k+1))
  R<-matrix(0,1,(k+1))
  for (t in 1:period){
    for (p in 1:(k+1)){
      R[,p]<-1
      beta<-dotDDMP_W[t,p]
      omega<-Omega[((t-1)*(k+1)+1):(t*(k+1)),((t-1)*(k+1)+1):(t*(k+1))]
      t_value[t,p]=sqrt(pl[t]*N)*beta/(sqrt(R%*%omega%*%t(R)))
      se[t,p]=sqrt(R%*%omega%*%t(R))/sqrt(pl[t]*N)
      p_value[t,p]<- 2 * pt(abs(t_value[t,p]),(pl[t]*N-1),lower.tail = FALSE)
      R<-matrix(0,1,(k+1))
    }
  }
  out3<-list(IC=IC,D_hat=dotDDMP_W,se=se,t_value=t_value,p_value=p_value)
}