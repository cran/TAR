#' @import mvtnorm utils
#' @export
#' 
#' @title
#' Identify the number of regimes and the corresponding thresholds for a Gaussian TAR model.
#' @description 
#' This function identify the number of regimes and the corresponding thresholds for a TAR model with Gaussian noise process.
#' @return 
#' The function returns the identified number of regimes with posterior probabilities and the thresholds with credible intervals.
#' @details
#' The TAR model is given by \deqn{X_t=a_0^{(j)} + \sum_{i=1}^{k_j}a_i^{(j)}X_{t-i}+h^{(j)}e_t} when \eqn{Z_t\in (r_{j-1},r_j]} for som \eqn{j} (\eqn{j=1,\cdots,l}). 
#' the \eqn{\{Z_t\}} is the threshold process, \eqn{l} is the number of regimes, \eqn{k_j} is the autoregressive order in the regime \eqn{j}. \eqn{a_i^{(j)}} with \eqn{i=0,1,\cdots,k_j} denote the autoregressive coefficients, while \eqn{h^{(j)}} denote the variance weights. \eqn{\{e_t\}} is the Gaussian white noise process \eqn{N(0,1)}.
#' @author Hanwen Zhang <hanwenzhang at usantotomas.edu.co>
#' @param Z The threshold series
#' @param X The series of interest
#' @param n.sim Number of iteration for the Gibbs Sampler
#' @param p.burnin Percentage of iterations used for Burn-in
#' @param n.thin Thinnin factor for the Gibbs Sampler
#' 
#' @references 
#' Nieto, F. H. (2005), \emph{Modeling Bivariate Threshold Autoregressive Processes in the Presence of Missing Data}. Communications in Statistics. Theory and Methods, 34; 905-930
#' @seealso \code{\link{LS.norm}}
#' @examples 
#' set.seed(12345678)
#' # Example 1, TAR model with 2 regimes
#' Z<-arima.sim(n=300,list(ar=c(0.5)))
#' l <- 2
#' r <- 0
#' K <- c(2,1)
#' theta <- matrix(c(1,-0.5,0.5,-0.7,-0.3,NA), nrow=l)
#' H <- c(1, 1.5)
#' X <- simu.tar.norm(Z,l,r,K,theta,H)
#' #res <- reg.thr.norm(Z,X)
#' #res$L.est
#' #res$L.prob
#' #res$R.est
#' #res$R.CI
#' 


reg.thr.norm <- function(Z, X, n.sim=500, p.burnin=0.2, n.thin=1){
  set.seed(123456789)
    l_Max=3
    l_Min=2
    k_Max=3
    k_Min=0
 
     N <- length(Z)
    
    # Prios for L
    l_param <- round(mean(c(l_Min,l_Max)))
    DPIL<-dpois(l_Min:l_Max,l_param)/sum(dpois(l_Min:l_Max,l_param))
    
    # Prioris for K
    k_param <- round(mean(c(k_Min,k_Max)))
    DPIK<-dpois(k_Min:k_Max,k_param)/sum(dpois(k_Min:k_Max,k_param))
    
    # Priors for the autoregressive coefficients  
    theta_0j <- array(0,c(l_Max-l_Min+1,l_Max,k_Max+1))
    v_theta <- 0.1
    V_0j <- diag(rep(v_theta,k_Max+1)) 
    
    # Priors for the variance weights
    esp <- 1
    vari<-5
    beta <- esp/vari
    alpha <- esp*beta 
    # Least Square estimation for the parameters
    estima.theta <- function(l,K,R){
      theta.esti <- function(K,j){
        theta.est <- rep(NA,K[j]+1)
        AA<-W_j(l,K,j,R)
        theta.est[1:(K[j]+1)] <- solve(t(AA$W.j)%*%AA$W.j)%*%t(AA$W.j)%*%AA$X.j
        theta.est      
      }
      theta.gorro <- matrix(NA,l,k_Max+1)
      
      for(ll in 1:l){
        theta.gorro[ll,1:(K[ll]+1)] <- theta.esti(K,ll)
      }
      J<-DIZ(l,K,R)$J  
      k<-max(K)
      e<-rep(NA,N)
      if(k>0) {
        e[1:k]<-0
        for(t in (k+1):N){
          jt<-J[t-k]   
          a<-theta.gorro[jt,1:(K[jt]+1)]   
          e[t]<-X[t]-a[1]-sum(a[-1]*c(X[(t-1):(t-K[jt])]))
        }
      }
      
      ##########################
      if(k==0){
        e<-rep(NA,N)
        for(t in 1:N){
          jt<-J[t]   
          a<-theta.gorro[jt,1]   
          e[t]<-X[t]-a
        }
      }
      
      h.est <- rep(NA,l)
      for(ll in 1:l){
        h.est[ll] <- sd(e[which(J==ll)])
      }
      list(theta.est=theta.gorro,h.est=h.est)
      
    }
    # Likelihood function
    LIKHOD<-function(L,K,A,H,R){
      J<-DIZ(L,K,R)$J
      nj<-table(J)
      k<-max(K)
      e<-rep(NA,N)
      if(k>0){
        e[1:k]<-0
        ht<-1
        for(t in (k+1):N){
          jt<-J[t-k]  
          a<-A[jt,1:(K[jt]+1)]  
          e[t]<-(X[t]-a[1]-sum(a[-1]*X[(t-1):(t-K[jt])]))/H[jt]
          ht<-ht*H[jt]
        }
      }
      ##########################
      if(k==0){
        e<-rep(NA,N)
        ht<-1
        for(t in 1:N){
          jt<-J[t]  
          a<-A[jt,1]  
          e[t]<-(X[t]-a)/H[jt]
          ht<-ht*H[jt]
        }
      }
      e.n <- e[(k+1):N]
      f<-(2*pi)^((k-N)/2)*exp(-sum(e.n^2)/2)/ht
      # f.log<--sum(e.n^2)/2-sum(nj*log(H))-(N-k)*log(2*pi)/2
      f.log<--sum(e.n^2)/2
      lista<-list(f=f, f.log=f.log, e=e.n)
      return(lista)
    }
    DIZ<-function(L,K,R){
      l <- L
      J <- rep(l, N)
      r <- c(min(Z)-0.1, R)
      for(j in 1:l){
        J[which(as.double(Z<=r[j+1])-as.double(Z<=r[j])==1)] <- j
      }
      k <- max(K)
      if(k==0){
        J.k <- J
        nj.k <- table(J.k)
      }
      if(k>0){
        J.k <- J[-(1:k)]
        nj.k <- table(J[-(1:k)])
      }
      list(J=J.k,nj=nj.k)
    }
    W_j<-function(L,K,j,R){
      k<-max(K)
      J.k<-DIZ(L,K,R)$J
      nj.k <- DIZ(L,K,R)$nj
      W.j<-matrix(NA,nj.k[j],K[j]+1)
      W.j[,1]<-1
      if(K[j]>0){
        ss<-0
        for(i in 1:(N-k)){
          if(J.k[i]==j){
            ss<-ss+1
            W.j[ss,-1]<-X[(i+k-1):(i+k-K[j])]
          } }       }
      X.j<-matrix(X[which(J.k==j)+k])
      list(W.j=W.j,X.j=X.j)
    }
    # Posterior probabilities for the regime number
    aux_L <- function(L,K,theta,H,R){
      LL.log <- LIKHOD(L,K,theta,H,R)$f.log+log(DPIL[L-1])
      LL <- LIKHOD(L,K,theta,H,R)$f*DPIL[L-1]
      lista<-list(LL=LL,LL.log=LL.log)
      return(lista)
    }
    ### All posible thresholds for L=2
    pos_R2 <- quantile(Z,seq(0.01,0.99,by=0.005))
     ### Eliminte all thresholds that induce too few observations in any regime
    indicador <- c(0)
    
    print("Creating all possible thresholds for 2 regimes...", quote = FALSE)
    n.min <- 40
    pbR2 <- txtProgressBar(min = 0, max = length(pos_R2), style = 3)
    for(iii in 1:length(pos_R2)){
      J2<-rep(NA,N)
      for(ii in 1:N){
        if(Z[ii]<=pos_R2[iii])        {J2[ii]<-1}
        if(Z[ii]>pos_R2[iii])         {J2[ii]<-2  }
      }
      nj <- table(J2)
      if(min(nj)>n.min&length(nj)==2){indicador<-c(indicador,iii)}
      Sys.sleep(0.1)   
      setTxtProgressBar(pbR2, iii)
    }
    close(pbR2)
    indicador<-indicador[-1]
    pos_R2 <- pos_R2[indicador]
    ### All posible thresholds for L=3
    m <- length(pos_R2)
    pos_R3 <- rep(NA,2)
    for(r in 1:(m-1)) {
      pos_R3 <- rbind(pos_R3,cbind(rep(pos_R2[r],m-r),pos_R2[(r+1):m]))}
    pos_R3 <- pos_R3[-1,]
    ### Eliminte all thresholds that induce too few observations in any regime
    indicador <- c(0)
    
    print("Creating all possible thresholds for 3 regimes...", quote = FALSE)
    
    pbR3 <- txtProgressBar(min = 0, max = dim(pos_R3)[1], style = 3)
    for(iii in 1:dim(pos_R3)[1]){
      J3<-rep(NA,N)
      for(ii in 1:N){
        if(Z[ii]<=pos_R3[iii,1])        {J3[ii]<-1}
        if(Z[ii]>pos_R3[iii,1]&Z[ii]<=pos_R3[iii,2])  {J3[ii]<-2}
        if(Z[ii]>pos_R3[iii,2])         {J3[ii]<-3  }
      }
      nj <- table(J3)
      if(min(nj)>n.min&length(nj)==3){indicador<-c(indicador,iii)}
      Sys.sleep(0.1)   
      setTxtProgressBar(pbR3, iii)
    }
    close(pbR3)
    indicador<-indicador[-1]
    pos_R3 <- pos_R3[indicador,]
    ### All posible thresholds for L=4
   # m <- length(pos_R2)
  #  pos_R4 <- rep(NA,2)
  #  R4_aux <- c(NA)
  #  for(mr in 2:(m-1)){
  #    for(r in mr:(m-1)) {
  #      pos_R4 <- rbind(pos_R4,cbind(rep(pos_R2[r],m-r),pos_R2[(r+1):m]))}}
  #  for(mm in 2:(m-1)){
  #    R4_aux <- c(R4_aux,rep(mm-1,(m-mm)*(m-mm+1)/2))}
  #  pos_R4 <- cbind(pos_R2[R4_aux[-1]],pos_R4[-1,])
  #  ### Eliminate all thresholds that induce too few observations in any regime
  #  indicador <- c(0)
  # for(iiii in 1:dim(pos_R4)[1]){
  #    J4<-rep(NA,N)
  #    for(ii in 1:N){
  #      if(Z[ii]<=pos_R4[iiii,1])        {J4[ii]<-1}
  #      if(Z[ii]>pos_R4[iiii,1]&Z[ii]<=pos_R4[iiii,2])  {J4[ii]<-2}
  #      if(Z[ii]>pos_R4[iiii,2]&Z[ii]<=pos_R4[iiii,3])  {J4[ii]<-3}
  #      if(Z[ii]>pos_R4[iiii,3])         {J4[ii]<-4  }
  #    }
  #    nj <- table(J4)
  #    if(min(nj)>n.min&length(nj)==4){indicador<-c(indicador,iiii)}
  #  }
  #  indicador<-indicador[-1]
  #  pos_R4 <- pos_R4[indicador,]
    ### Matrix for all possible thresholds for all regime numbers
    pos_R <- array(10000,c(3,dim(pos_R3)[1],l_Max-1))
    pos_R[1,1:length(pos_R2),1]<-pos_R2
    pos_R[2,1:dim(pos_R3)[1],1:2]<-pos_R3
    # pos_R[3,1:dim(pos_R4)[1],1:3]<-pos_R4
    # Coditional density for thresholds
    aux_R <- function(L,K,theta,H,R){
      k<-max(K)
      J.k<-DIZ(L,K,R)$J
      e.R <-LIKHOD(L,K,theta,H,R)$e
      RR<-(2*pi)^((k-N)/2)*exp(-sum(e.R^2)/2)/(prod(H[J.k]))
      RR.log<--sum(e.R^2)/2-sum(log(H[J.k]))+((k-N)/2)*log(2*pi) 
      lista<-list(RR=RR,RR.log=RR.log)
      return(lista)
    }
    
    # Log posterior probabilities for the autoregressive orders
    aux_k <- function(L,K,theta,H,ki,R){
      KK.log <- LIKHOD(L,K,theta,H,R)$f.log+log(DPIK[ki+1])
      KK <- LIKHOD(L,K,theta,H,R)$f*DPIK[ki+1]
      lista<-list(KK=KK,KK.log=KK.log)
      return(lista)
    }
    # Sample one value for the autoregressive coefficient vector theta_j
    sample_theta <- function(L,K,j,h2,R){
      AA<-W_j(L,K,j,R)
      Sigma <- V_0j[1:(K[j]+1),1:(K[j]+1)]+ t(AA$W.j)%*%AA$W.j/h2[j]
      mu <- solve(Sigma)%*%(t(AA$W.j)%*%AA$X.j/h2[j]+V_0j[1:(K[j]+1),1:(K[j]+1)]%*%theta_0j[L-1,j,1:(K[j]+1)])
      rmvnorm(1,mu,solve(Sigma))
    }
    # Sample one value for the variance weight h2_j
    sample_h2 <- function(L,K,j,theta_j,R){
      k<-max(K)
      J.k<-DIZ(L,K,R)$J
      nj.k<-table(DIZ(L,K,R)$J)
      theta_j<-matrix(theta_j)
      AA<-W_j(L,K,j,R)
      alfa.j <- alpha+nj.k[j]/2
      beta.j <- beta+sum((AA$X.j-AA$W.j%*%theta_j)^2)/2
      1/rgamma(1,shape=alfa.j,scale=1/beta.j)
    }
    # Space to store the values for all the parameters
    res_L <- rep(l_param,n.sim)
    res_R <- array(NA,c(n.sim,l_Max-1,l_Max-1))
    res_R[,1,1] <- sample(pos_R2,n.sim,replace=T)
    for(i in 1:n.sim){
      res_R[i,2,1:2] <- pos_R3[sample(c(1:dim(pos_R3)[1]),1),]
    #  res_R[i,3,1:3] <- pos_R4[sample(c(1:dim(pos_R4)[1]),1),]
    }
    res_K <- array(NA,c(l_Max-l_Min+1,l_Max,n.sim))
    for(rrr in 1:n.sim){
      for(ii in 1:(l_Max-l_Min+1)){
        for(jjj in 1:(ii+1)){
          u<-runif(1)
          priori.K<-c(0,cumsum(DPIK))
          res_K[ii,jjj,rrr] <- min(which(as.double(u<=priori.K)==1))-2
        }}}
    res_theta_j_2 <- array(rnorm(n.sim*(k_Max+1)^(2+1)*2,0,sqrt(1/v_theta) ),c(k_Max+1,k_Max+1,2,n.sim,k_Max+1))
    res_theta_j_3 <- array(rnorm(n.sim*(k_Max+1)^(3+1)*3,0,sqrt(1/v_theta) ),c(k_Max+1,k_Max+1,k_Max+1,3,n.sim,k_Max+1))
   # res_theta_j_4 <- array(rnorm(n.sim*(k_Max+1)^(4+1)*4,0,sqrt(1/v_theta) ),c(k_Max+1,k_Max+1,k_Max+1,k_Max+1,4,n.sim,k_Max+1))
    res_h_j <- array(sqrt(rgamma((l_Max-l_Min+1)*l_Max*n.sim,shape=alpha,scale=1/beta)),c(l_Max-l_Min+1,l_Max,n.sim))
    
    prob_L <- matrix(NA,l_Max-l_Min+1, n.sim)
  
    # ----------------------------
    # Run Gibbs sampler
    # ----------------------------
    
    print("Running Gibbs Sampler...", quote = FALSE)
    
    pbGibbs <- txtProgressBar(min = 0, max = n.sim-1, style = 3)
    
    for(rr in 2:n.sim){
      # Select the regime number
      cond.l <- rep(NA,l_Max-l_Min+1)
      cond.l[1]<-aux_L(2,res_K[1,1:2,rr-1],res_theta_j_2[res_K[1,1,rr-1]+1,res_K[1,2,rr-1]+1,,rr-1,],res_h_j[1,1:2,rr-1],res_R[rr-1,2-1,1])$LL.log
      cond.l[2]<-aux_L(3,res_K[2,1:3,rr-1],res_theta_j_3[res_K[2,1,rr-1]+1,res_K[2,2,rr-1]+1,res_K[2,3,rr-1]+1,,rr-1,],res_h_j[2,1:3,rr-1],res_R[rr-1,3-1,1:2])$LL.log
     # cond.l[3]<-aux_L(4,res_K[3,1:4,rr-1],res_theta_j_4[res_K[3,1,rr-1]+1,res_K[3,2,rr-1]+1,res_K[3,3,rr-1]+1,res_K[3,4,rr-1]+1,,rr-1,],res_h_j[3,1:4,rr-1],res_R[rr-1,4-1,1:3])$LL.log
      cond.ll<-exp(cond.l)
      if(min(cond.l)<= -744){cond.ll<-exp(cond.l-max(cond.l))}
      if(max(cond.l)>700){cond.ll<-exp(cond.l-max(cond.l))}
      temp<-runif(1,0,1)
      aux.l <- c(0,cumsum(cond.ll/sum(cond.ll)))
      prob_L[,rr] <- cond.ll/sum(cond.ll)
      res_L[rr] <- min(which(as.double(temp<=aux.l)==1))
      # Select the thresholds
      cond.R <- rep(NA)
      if(res_L[rr]==2){
        for(rc in 1:length(pos_R2)){
          cond.R[rc] <- aux_R(res_L[rr],res_K[res_L[rr]-1,1:res_L[rr],rr-1],res_theta_j_2[res_K[res_L[rr]-1,1,rr-1]+1,res_K[res_L[rr]-1,2,rr-1]+1,,rr-1,],res_h_j[res_L[rr]-1,1:res_L[rr],rr-1],pos_R2[rc])$RR.log}
        cond.RR<-exp(cond.R)
        if(min(cond.R)<= -744){cond.RR<-exp(cond.R-max(cond.R))}
        if(max(cond.R)>700){cond.RR<-exp(cond.R-max(cond.R))}
        res_R[rr,res_L[rr]-1,1:(res_L[rr]-1)]<-pos_R2[sample(c(1:length(pos_R2)),size=1,prob=cond.RR/(sum(cond.RR)))]
      }
      
      if(res_L[rr]==3){
        for(rc in 1:dim(pos_R3)[1]){
          cond.R[rc] <- aux_R(res_L[rr],res_K[res_L[rr]-1,1:res_L[rr],rr-1],res_theta_j_3[res_K[res_L[rr]-1,1,rr-1]+1,res_K[res_L[rr]-1,2,rr-1]+1,res_K[res_L[rr]-1,3,rr-1]+1,,rr-1,],res_h_j[res_L[rr]-1,1:res_L[rr],rr-1],pos_R3[rc,])$RR.log
          }
        cond.RR<-exp(cond.R)
        if(min(cond.R)<= -744){cond.RR<-exp(cond.R-max(cond.R))}
        if(max(cond.R)>700){cond.RR<-exp(cond.R-max(cond.R))}
        res_R[rr,res_L[rr]-1,1:(res_L[rr]-1)]<-pos_R3[sample(c(1:dim(pos_R3)[1]),size=1,prob=cond.RR/(sum(cond.RR))),]
      }
      
  #    if(res_L[rr]==4){
  #      for(rc in 1:dim(pos_R4)[1]){
  #        cond.R[rc] <- aux_R(res_L[rr],res_K[res_L[rr]-1,1:res_L[rr],rr-1],res_theta_j_4[res_K[res_L[rr]-1,1,rr-1]+1,res_K[res_L[rr]-1,2,rr-1]+1,res_K[res_L[rr]-1,3,rr-1]+1,res_K[res_L[rr]-1,4,rr-1]+1,,rr-1,],res_h_j[res_L[rr]-1,1:res_L[rr],rr-1],pos_R4[rc,])$RR.log}
  #      cond.RR<-exp(cond.R)
  #      if(min(cond.R)<= -744){cond.RR<-exp(cond.R-max(cond.R))}
  #      if(max(cond.R)>700){cond.RR<-exp(cond.R-max(cond.R))}
  #      res_R[rr,res_L[rr]-1,1:(res_L[rr]-1)]<-pos_R4[sample(c(1:dim(pos_R4)[1]),size=1,prob=cond.RR/(sum(cond.RR))),]
  #    }
      
      # Select the autoregressive orders
      for(j in 1:res_L[rr]){
        cond.k <- rep(NA,k_Max-k_Min+1)
        for(ki in k_Min:k_Max){
          K.can <- res_K[res_L[rr]-1,1:res_L[rr],rr-1]
          K.can[j] <- ki
          if(j>1){K.can[1:(j-1)]<-res_K[res_L[rr]-1,1:(j-1),rr]}
          #         if(res_L[rr]==2){
          #           theta_j <- res_theta_j_2[K.can[1]+1,K.can[2]+1,,rr-1,]}
          #         if(res_L[rr]==3){
          #           theta_j <- res_theta_j_3[K.can[1]+1,K.can[2]+1,K.can[3]+1,,rr-1,]}
          #         if(res_L[rr]==4){
          #          theta_j <- res_theta_j_4[K.can[1]+1,K.can[2]+1,K.can[3]+1,K.can[4]+1,,rr-1,]}
          theta_j<-estima.theta(res_L[rr],K.can,res_R[rr,res_L[rr]-1,1:(res_L[rr]-1)])$theta.est
          cond.k[ki+1] <- aux_k(res_L[rr],K.can,theta_j,res_h_j[res_L[rr]-1,1:res_L[rr],rr-1],ki,res_R[rr,res_L[rr]-1,1:(res_L[rr]-1)])$KK.log
        }
        
        cond.kk<-exp(cond.k)
        if(min(cond.k)<= -744){cond.kk<-exp(cond.k-max(cond.k))}
        if(max(cond.k)>700){cond.kk<-exp(cond.k-max(cond.k))}
        
        ## seleccionar un valor para k_j
        temp<-runif(1,0,1)
        aux.k <- c(0,cumsum(cond.kk/sum(cond.kk)))
        res_K[res_L[rr]-1,j,rr] <- min(which(as.double(temp<=aux.k)==1))-2   
        }
      
      # Estimtate autoregressiv coefficients
      K.f <- res_K[res_L[rr]-1,1:res_L[rr],rr]
      if(res_L[rr]==2){
        for(j in 1:(res_L[rr])){
          res_theta_j_2[K.f[1]+1,K.f[2]+1,j,rr,1:(K.f[j]+1)]<-sample_theta(res_L[rr],K.f,j,(res_h_j[res_L[rr]-1,1:res_L[rr],rr-1])^2,res_R[rr-1,res_L[rr]-1,1:(res_L[rr]-1)])
          #res_theta_j[rr,j,1:(K[j]+1)] <- A.real[j,1:(K[j]+1)]
        } }
      if(res_L[rr]==3){
        for(j in 1:(res_L[rr])){
          res_theta_j_3[K.f[1]+1,K.f[2]+1,K.f[3]+1,j,rr,1:(K.f[j]+1)]<-sample_theta(res_L[rr],K.f,j,(res_h_j[res_L[rr]-1,1:res_L[rr],rr-1])^2,res_R[rr-1,res_L[rr]-1,1:(res_L[rr]-1)])
          #res_theta_j[rr,j,1:(K[j]+1)] <- A.real[j,1:(K[j]+1)]
        }}
      if(res_L[rr]==4){
        for(j in 1:(res_L[rr])){
          res_theta_j_4[K.f[1]+1,K.f[2]+1,K.f[3]+1,K.f[4]+1,j,rr,1:(K.f[j]+1)]<-sample_theta(res_L[rr],K.f,j,(res_h_j[res_L[rr]-1,1:res_L[rr],rr-1])^2,res_R[rr-1,res_L[rr]-1,1:(res_L[rr]-1)])
          #res_theta_j[rr,j,1:(K[j]+1)] <- A.real[j,1:(K[j]+1)]
        }}
      # Estimate variance weights
      if(res_L[rr]==2){
        for(jj in 1:res_L[rr]){
          res_h_j[res_L[rr]-1,jj,rr]<-sqrt(sample_h2(res_L[rr],K.f,jj,res_theta_j_2[K.f[1]+1,K.f[2]+1,jj,rr,1:(K.f[jj]+1)],res_R[rr,res_L[rr]-1,1:(res_L[rr]-1)])) 
        }                }
      if(res_L[rr]==3){
        for(jj in 1:res_L[rr]){
          res_h_j[res_L[rr]-1,jj,rr]<-sqrt(sample_h2(res_L[rr],K.f,jj,res_theta_j_3[K.f[1]+1,K.f[2]+1,K.f[3]+1,jj,rr,1:(K.f[jj]+1)],res_R[rr,res_L[rr]-1,1:(res_L[rr]-1)]))
        }
      }
      if(res_L[rr]==4){
        for(jj in 1:res_L[rr]){
          res_h_j[res_L[rr]-1,jj,rr]<-sqrt(sample_h2(res_L[rr],K.f,jj,res_theta_j_4[K.f[1]+1,K.f[2]+1,K.f[3]+1,K.f[4]+1,jj,rr,1:(K.f[jj]+1)],res_R[rr,res_L[rr]-1,1:(res_L[rr]-1)])) 
        }
      }
     # print(res_L[rr])
      Sys.sleep(0.1)
      # update progress bar
      setTxtProgressBar(pbGibbs, rr)
    } # Here ends the Gibbs Sampler
  close(pbGibbs)
    
  iter.final <- seq(p.burnin*n.sim+1, n.sim, by=n.thin)
  # Identify the autoregressive orders
  table(res_L[iter.final])
  L.prob <- matrix(rowMeans(prob_L[,-1]),nrow=1)
  colnames(L.prob) <- c("2 regimes", "3 regimes")
  rownames(L.prob) <- c("Posterior probability")
  L.final <- which(L.prob==max(L.prob))+1
  R.final <- c()
  R.IC <- matrix(NA,2,L.final-1)
  rownames(R.IC) <- c("Limit inferior","Limit superior")
  colnames(R.IC) <- colnames(R.IC, do.NULL = FALSE, prefix = "Regime")
  for(j in 1:(L.final-1)){
    R.final[j] <- mean(res_R[which(res_L[iter.final]==L.final)+p.burnin*n.sim,L.final-1,j])
    R.IC[,j] <- quantile(res_R[which(res_L[iter.final]==L.final)+p.burnin*n.sim,L.final-1,j], c(0.025,0.975))
  }
    impre <- cat("The identified regime number is: ", L.final, ".\n\nA posteriori probability of the number of regimes are:\n")
    print(L.prob)
    cat("\nThe estimated threshold(s) is(are):", R.final)
  indd <- c("1st","2nd","3rd")
  for(j in 1:(L.final-1)){
    cat(".\n\nThe 95% credible interval of the",indd[j],"threshold is:", 
  "\n", "(",R.IC[1,j],",",R.IC[2,j],")")
  }
  results = list(L.est=L.final, L.prob=L.prob, R.est=R.final, R.CI=R.IC, 
                 R.values = res_R[which(res_L[iter.final]==L.final)+p.burnin*n.sim,L.final-1,j], L.values = res_L[iter.final])
  results
}
