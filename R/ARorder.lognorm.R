#' @import mvtnorm utils
#' @export
#' 
#' @title
#' Identify the autoregressive orders for a log-normal TAR model given the number of regimes and thresholds.
#' @description 
#' This function identify the autoregressive orders for a log-normal TAR model given the number of regimes and thresholds.
#' @return 
#' The identified autoregressive orders with posterior probabilities
#' @details
#' The log-normal TAR model is given by \deqn{log X_t=a_0^{(j)} + \sum_{i=1}^{k_j}a_i^{(j)}log X_{t-i}+h^{(j)}e_t} when \eqn{Z_t\in (r_{j-1},r_j]} for some \eqn{j} (\eqn{j=1,\cdots,l}). 
#' the \eqn{\{Z_t\}} is the threshold process, \eqn{l} is the number of regimes, \eqn{k_j} is the autoregressive order in the regime \eqn{j}. \eqn{a_i^{(j)}} with \eqn{i=0,1,\cdots,k_j} denote the autoregressive coefficients, while \eqn{h^{(j)}} denote the variance weights. \eqn{\{e_t\}} is the Gaussian white noise process \eqn{N(0,1)}.
#' @author Hanwen Zhang <hanwenzhang at usantotomas.edu.co>
#' @param Z The threshold series
#' @param X The series of interest
#' @param l The number of regimes.
#' @param r The vector of thresholds for the series \{Z_t\}.
#' @param k_Max The minimum value for each autoregressive order. The default is 3.
#' @param k_Min The maximum value for each autoregressive order. The default is 0.
#' @param n.sim Number of iteration for the Gibbs Sampler
#' @param p.burnin Percentage of iterations used for burn-in
#' @param n.thin Thinnin factor for the Gibbs Sampler
#' 
#' @references 
#' Nieto, F. H. (2005), \emph{Modeling Bivariate Threshold Autoregressive Processes in the Presence of Missing Data}. Communications in Statistics. Theory and Methods, 34; 905-930
#' @seealso \code{\link{simu.tar.lognorm}, \link{ARorder.norm}}
#' @examples 
#' set.seed(12345678)
#' Z<-arima.sim(n=500,list(ar=c(0.5)))
#' l <- 2
#' r <- 0
#' K <- c(2,1)
#' theta <- matrix(c(1,0.5,-0.3,-0.5,-0.7,NA),nrow=l)
#' H <- c(1, 1.3)
#' X <- simu.tar.lognorm(Z,l,r,K,theta,H)
#' #res <- ARorder.lognorm(Z,X,l,r)
#' #res$K.est
#' #res$K.prob

ARorder.lognorm <- function(Z, X, l, r, k_Max=3, k_Min=0, n.sim=500, p.burnin=0.3, n.thin=1){
    if(length(r)!=(l-1))
      stop("The number of thresholds should be l-1")
    N <- length(Z)
    
    # Prioris for the autoregressive orders K
    k_param <- round(mean(c(k_Min,k_Max)))
    DPIK<-dpois(k_Min:k_Max,k_param)/sum(dpois(k_Min:k_Max,k_param))
    
    # Priors for the autoregressive coefficients and variance weights 
    theta_0j <- matrix(0,l,k_Max+1)
    v_theta <- 0.01
    V_0j <- diag(rep(v_theta, k_Max+1)) 
    alpha <- 2
    beta <- 3
    
    # Likelihood function
    LIKHOD<-function(K,A,H){
      J<-DIZ(l,K,r)$J
      nj<-table(J)
      k<-max(K)
      e<-rep(NA,N)
      if(k>0){
        e[1:k]<-0
        ht<-1
        for(t in (k+1):N){
          jt<-J[t-k]   
          a<-A[jt,1:(K[jt]+1)]   
          e[t]<-(log(X[t])-a[1]-sum(a[-1]*log(X[(t-1):(t-K[jt])])))/H[jt]
          ht<-ht*H[jt]*X[t]
        }
      }
      ##########################
      if(k==0){
        e<-rep(NA,N)
        ht<-1
        for(t in 1:N){
          jt<-J[t]   
          a<-A[jt,1]  
          e[t]<-(log(X[t])-a)/H[jt]
          ht<-ht*H[jt]*X[t]
        }
      }
      e.n <- e[(k+1):N]
      f<-(2*pi)^((k-N)/2)*exp(-sum(e.n^2)/2)/ht
      f.log<--sum(e.n^2)/2
      #-log(ht)-(N-k)*log(2*pi)/2
      lista<-list(f=f, f.log=f.log, e=e.n)
      return(lista)
    }

    # Function for regime indicators and regime sizes 
    DIZ<-function(L,K,r){
      J <- rep(l, N)
      r <- c(min(Z)-0.1, r)
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
    
    W_j<-function(K,j){
      k<-max(K)
      J.k<-DIZ(l,K,r)$J
      nj.k <- DIZ(l,K,r)$nj
      W.j<-matrix(NA,nj.k[j],K[j]+1)
      W.j[,1]<-1
      if(K[j]>0){
        ss<-0
        for(i in 1:(N-k)){
          if(J.k[i]==j){
            ss<-ss+1
            W.j[ss,-1]<-log(X[(i+k-1):(i+k-K[j])])
          } }       }
      X.j<-log(matrix(X[which(J.k==j)+k]))
      list(W.j=W.j,X.j=X.j)
    }
    # Space to store the values for the parameters and the autoregressive orders
    res_theta_j <- array(NA, c(rep(k_Max+1,l),l,n.sim,k_Max+1))
    res_h_j<-array(NA, c(rep(k_Max+1,l),l,n.sim))
    for(k in k_Min:k_Max){
      for(kk in k_Min:k_Max){
        J.aux <- c(k,kk)
        LS <- LS.lognorm(Z,X,l,r,c(k,kk))
        for(j in 1:l){
          for(ri in 1:n.sim){
            res_theta_j[k+1, kk+1, j, ri,1:(J.aux[j]+1)] <- LS$theta.est[j,1:(J.aux[j]+1)]
            res_h_j[k+1, kk+1, j, ri]<- LS$h.est[j]
          }
        }
      }
    }
    
    res_K <- matrix(NA,l,n.sim)
    res_K[,1] <- k_Min
    prob_K <- array(NA,c(l,k_Max-k_Min+1,n.sim))
    # log posterior probabilities for the autoregressive orders
    aux_k <- function(K,theta,H,ki){
      LIKHOD(K,theta,H)$f.log+log(DPIK[ki+1])
    }
    # Sample one value for the autoregressive coefficient vector theta_j
    sample_theta <- function(K,j,h2){
      AA<-W_j(K,j)
      Sigma <- V_0j[1:(K[j]+1),1:(K[j]+1)]+ t(AA$W.j)%*%AA$W.j/h2[j]
      mu <- solve(Sigma)%*%(t(AA$W.j)%*%AA$X.j/h2[j]+V_0j[1:(K[j]+1),1:(K[j]+1)]%*%theta_0j[j,1:(K[j]+1)])
      rmvnorm(1,mu,solve(Sigma))
    }
    # Sample one value for the variance weight h2_j
    sample_h2 <- function(K,j,theta_j){
      k<-max(K)
      J.k<-DIZ(l,K,r)$J
      nj.k<-table(DIZ(l,K,r)$J)
      theta_j<-matrix(theta_j)
      AA<-W_j(K,j)
      alfa.j <- alpha+nj.k[j]/2
      beta.j <- beta+sum((AA$X.j-AA$W.j%*%theta_j)^2)/2
      1/rgamma(1, shape=alfa.j, scale=1/beta.j)
    }
    # Progress bar
    pb <- txtProgressBar(min = 0, max = n.sim, style = 3)
    # Run Gibbs sampler
    for(rr in 2:n.sim){
      for(j in 1:l){
        cond.k <- rep(NA,k_Max-k_Min+1)
        for(ki in k_Min:k_Max){
          K.can <- res_K[,rr-1]
         if(j>1){K.can[1:(j-1)]<-res_K[1:(j-1),rr]}
          K.can[j] <- ki
          # AAA<-LS.norm(Z, X, l, r,K.can)
          #  theta_j<-AAA$theta.est
          # h_j<-AAA$h.est
          if(l==2){
            res_theta_j[K.can[1]+1,K.can[2]+1,j,rr,1:(K.can[j]+1)] <- sample_theta(K.can,j,(res_h_j[K.can[1]+1,K.can[2]+1,,rr])^2)
            res_h_j[K.can[1]+1,K.can[2]+1,j,rr] <- sqrt(sample_h2(K.can,j,res_theta_j[K.can[1]+1,K.can[2]+1,j,rr,1:(K.can[j]+1)]))
    #       res_theta_j[K.can[1]+1,K.can[2]+1,j,rr,1:(K.can[j]+1)] <- sample_theta(K.can,j,(res_h_j[K.can[1]+1,K.can[2]+1,,rr])^2)
    #       res_h_j[K.can[1]+1,K.can[2]+1,j,rr] <- sqrt(sample_h2(K.can,j,res_theta_j[K.can[1]+1,K.can[2]+1,j,rr,1:(K.can[j]+1)]))
          }
          if(l==3){
            res_theta_j[K.can[1]+1,K.can[2]+1,K.can[3]+1,j,rr,1:(K.can[j]+1)] <- sample_theta(K.can,j,(res_h_j[K.can[1]+1,K.can[2]+1,K.can[3]+1,,rr])^2)
            res_h_j[K.can[1]+1,K.can[2]+1,K.can[3]+1,j,rr] <- sqrt(sample_h2(K.can,j,res_theta_j[K.can[1]+1,K.can[2]+1,K.can[3]+1,j,rr,1:(K.can[j]+1)]))
          }
          if(l==4){
            res_theta_j[K.can[1]+1,K.can[2]+1,K.can[3]+1,K.can[4]+1,j,rr,1:(K.can[j]+1)] <- sample_theta(K.can,j,(res_h_j[K.can[1]+1,K.can[2]+1,K.can[3]+1,K.can[4]+1,,rr])^2)
            res_h_j[K.can[1]+1,K.can[2]+1,K.can[3]+1,K.can[4]+1,j,rr] <- sqrt(sample_h2(K.can,j,res_theta_j[K.can[1]+1,K.can[2]+1,K.can[3]+1,K.can[4]+1,j,rr,1:(K.can[j]+1)]))
          }
          cond.k[ki+1] <- aux_k(K.can,as.matrix(res_theta_j[K.can[1]+1,K.can[2]+1,,rr,1:(max(K.can)+1)]),res_h_j[K.can[1]+1,K.can[2]+1,,rr],ki)
        }
        cond.kk<-exp(cond.k)
        if(min(cond.k)<= -740){cond.kk<-exp(cond.k-max(cond.k))}
        if(max(cond.k)>700){cond.kk<-exp(cond.k-max(cond.k))}
        prob_K[j,,rr]<-cond.kk/sum(cond.kk)
        temp<-runif(1,0,1)
        cond.kkk <- c(0,cumsum(prob_K[j,,rr]))
        res_K[j,rr] <- min(which(as.double(temp<=cond.kkk)==1))-2 
      }
      
      ## Sample un value for k_j
     # print(rr)
    #  print(res_K[,rr])  
      Sys.sleep(0.1)
      # update progress bar
      setTxtProgressBar(pb, rr)
    } # Here ends the Gibbs Sampler
  close(pb)
    
  iter.final <- seq(p.burnin*n.sim, n.sim, by=n.thin)
  # Identify the autoregressive orders
  K.final <- c()
 # fre.K <- c()
  K.matrix <- matrix(NA, l, k_Max-k_Min+1)
  colnames(K.matrix) <- paste(k_Min:k_Max)
  rownames(K.matrix) <- paste("Regime", 1:l, sep="")
  for(j in 1:l){
    K.matrix[j,] <- rowMeans(prob_K[j,,iter.final])
    K.final[j] <- which(K.matrix[j,]==max(K.matrix[j,]))-1 
  }

  indd <- c("1st","2nd","3rd")
  if(l>3){
    indd <- c(c("1st","2nd","3rd"), paste(4:l,"th",sep=""))
  }
  
  for(j in 1:l){
    impre <- cat("The identified AR order in  the",indd[j],"regime is:" , K.final[j],"\n")
  }
  cat("\nA posteriori probability of the autoregressive orders are:","\n")
  print(K.matrix)
  
  list(K.est=K.final, K.prob=K.matrix, K.res=res_K[,iter.final])
}
