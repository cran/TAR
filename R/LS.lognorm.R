#' @import 
#' @export
#' 
#' @title
#' Estimate a log-normal TAR model using Least Square method given the structural parameters.
#' @description 
#' This function estimate a log-normal TAR model using Least Square method given the structural parameters, i.e. the number of regimes, thresholds and autoregressive orders.
#' @return 
#' The function returns the autoregressive coefficients matrix theta and variance weights H. Rows of the matrix theta represent regimes
#' @details
#' The TAR model is given by \deqn{log X_t=a_0^{(j)} + \sum_{i=1}^{k_j}a_i^{(j)}log X_{t-i}+h^{(j)}e_t} when \eqn{Z_t\in (r_{j-1},r_j]} for som \eqn{j} (\eqn{j=1,\cdots,l}). 
#' the \eqn{\{Z_t\}} is the threshold process, \eqn{l} is the number of regimes, \eqn{k_j} is the autoregressive order in the regime \eqn{j}. \eqn{a_i^{(j)}} with \eqn{i=0,1,\cdots,k_j} denote the autoregressive coefficients, while \eqn{h^{(j)}} denote the variance weights.
#' @author Hanwen Zhang <hanwenzhang at usantotomas.edu.co>
#' @param Z The threshold series
#' @param X The series of interest
#' @param l The number of regimes.
#' @param r The vector of thresholds for the series \eqn{\{Z_t\}}.
#' @param K The vector containing the autoregressive orders of the \eqn{l} regimes.
#' 
#' 
#' @references 
#' Nieto, F. H. (2005), \emph{Modeling Bivariate Threshold Autoregressive Processes in the Presence of Missing Data}. Communications in Statistics. Theory and Methods, 34; 905-930
#' @seealso \code{\link{simu.tar.norm}}
#' @examples 
#' Z<-arima.sim(n=500,list(ar=c(0.5)))
#' l <- 2
#' r <- 0
#' K <- c(2,1)
#' theta <- matrix(c(1,0.5,-0.3,-0.5,-0.7,NA),nrow=l)
#' H <- c(1, 1.3)
#' X <- simu.tar.lognorm(Z,l,r,K,theta,H)
#' ts.plot(X)
#' LS.lognorm(Z,X,l,r,K)

LS.lognorm <- function(Z, X, l, r, K){
  if(length(r)!=(l-1))
    stop("The number of thresholds should be l-1")
  if(length(K)!=l)
    stop("A TAR model with l regimes should have l autoregressive orders")
  N <- length(Z)
  J <- rep(l, N)
  r <- c(min(Z)-0.1, r)
  for(j in 1:l){
    J[which(as.double(Z<=r[j+1])-as.double(Z<=r[j])==1)] <- j
  }
  k <- max(K)
  if(k>0){
    J.k <- J[-(1:k)]
    nj.k <- table(J.k)
  }
  if(k==0){
    J.k <- J
    nj.k <- table(J.k)
  }
  
  theta.est <- matrix(NA,l,k+1)
  
  for(j in 1:l){
    CONT<-nj.k[j]
    Y<-matrix(NA,CONT,1)
    ZZ<-matrix(NA,CONT,(K[j]+1))
    if(K[j]==0){ZZ<-matrix(1,CONT,1)
    IN=0
    for(ix in (k+1):N){
      { if (J.k[ix-k]==j)
      { IN=IN+1
      Y[IN,1]=log(X[ix])
      }
      }
    }
    }
    if(K[j]>0){
      IN=0
      for  (iy in (k+1):N)
      {        if (J.k[iy-k]==j)
      { IN=IN+1
      Y[IN,1]=log(X[iy])
      ZZ[IN,1]=1
      for (ii in (1:K[j])){ZZ[IN,ii+1]=log(X[iy-ii])}
      }
      } }
   theta.est[j,1:(K[j]+1)] <- solve(t(ZZ)%*%ZZ)%*%t(ZZ)%*%Y}
   e<-rep(NA,N)
  if(k>0) {
    e[1:k]<-0
    for(t in (k+1):N){
      jt<-J[t]   ### jt es el r?gimen donde est? el t-?simo dato
      a<-theta.est[jt,1:(K[jt]+1)]   ### coeficientes autoregresivos para el regimen correspondiente al t-?simo dato
      e[t]<-log(X[t])-a[1]-sum(a[-1]*log(X[(t-1):(t-K[jt])]))
    }
  }
  ##########################
  if(k==0){
    e<-rep(NA,N)
    for(t in 1:N){
      jt<-J[t]   ### jt es el r?gimen donde est? el t-?simo dato
      a<-theta.est[jt,1]   ### coeficientes autoregresivos para el regimen correspondiente al t-?simo dato
      e[t]<-log(X[t])-a
    }
  }
  h.est <- rep(NA,l)
  for(ll in 1:l){
    h.est[ll] <- sd(e[which(J==ll)])
  }
  list(theta.est=theta.est,h.est=h.est)
  
}
