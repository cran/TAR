#' @import TSA
#' @export
#' 
#' @title
#' Simulate a series from a TAR model with Gaussian distributed error.
#' @description 
#' This function simulates a serie from a TAR model with Gaussian distributed error given the parameters of the model from a given threshold process \eqn{\{Z_t\}}
#' @return 
#' The time series \eqn{\{X_t\}}.
#' @details
#' The TAR model is given by \deqn{X_t=a_0^{(j)} + \sum_{i=1}^{k_j}a_i^{(j)}X_{t-i}+h^{(j)}e_t} when \eqn{Z_t\in (r_{j-1},r_j]} for som \eqn{j} (\eqn{j=1,\cdots,l}). 
#' the \eqn{\{Z_t\}} is the threshold process, \eqn{l} is the number of regimes, \eqn{k_j} is the autoregressive order in the regime \eqn{j}. \eqn{a_i^{(j)}} with \eqn{i=0,1,\cdots,k_j} denote the autoregressive coefficients, while \eqn{h^{(j)}} denote the variance weights. \eqn{\{e_t\}} is the Gaussian white noise process \eqn{N(0,1)}.
#' @author Hanwen Zhang <hanwenzhang at usantotomas.edu.co>
#' @param Z The threshold series
#' @param l The number of regimes.
#' @param r The vector of thresholds for the series \eqn{\{Z_t\}}.
#' @param K The vector containing the autoregressive orders of the l regimes.
#' @param theta The matrix of autoregressive coefficients of dimension \eqn{l\times\max{K}}. \eqn{j}-th row contains the autoregressive coefficients of regime \eqn{j}.
#' @param H The vector containing the variance weights of the \eqn{l} regimes.
#' 
#' @references 
#' Nieto, F. H. (2005), \emph{Modeling Bivariate Threshold Autoregressive Processes in the Presence of Missing Data}. Communications in Statistics. Theory and Methods, 34; 905-930
#' @seealso \code{\link{simu.tar.norm}}
#' @examples 
#' Z<-arima.sim(n=500,list(ar=c(0.5)))
#' l <- 2
#' r <- 0
#' K <- c(2,1)
#' theta <- matrix(c(1,-0.5,0.5,-0.7,-0.3,NA),nrow=l)
#' H <- c(1, 1.5)
#' X <- simu.tar.norm(Z,l,r,K,theta,H)
#' ts.plot(X)
#' 

simu.tar.norm <- function(Z, l, r, K, theta, H){
  if(length(r)!=(l-1))
    stop("The number of thresholds should be l-1")
  if(length(K)!=l)
    stop("A TAR model with l regimes should have l autoregressive orders")
  if(nrow(theta)!=l)
    stop("The number of rows of the matrix of autoregressive orders theta should be l")
  if(length(H)!=l)
    stop("The number of variance weights should be l")
  N <- length(Z)
  J <- rep(l, N)
  r <- c(min(Z)-0.1, r)
  for(j in 1:l){
    J[which(as.double(Z<=r[j+1])-as.double(Z<=r[j])==1)] <- j
  }
  k <- max(K)
  X <- rep(0,N)
  e <- rnorm(N)
  for(i in (k+1):N){
    j <- J[i]
    X[i] <- sum(c(1,X[i-(1:K[j])])*theta[j,1:(K[j]+1)]) + H[j]*e[i]
  }
  X
}