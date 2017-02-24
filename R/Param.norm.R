#' @import mvtnorm
#' @export
#' 
#' @title
#' Estimate a Gaussian TAR model using Gibbs Sampler given the structural parameters.
#' @description 
#' This function estimate a Gaussian TAR model using Gibbs Sampler given the structural parameters, i.e. the number of regimes, thresholds and autoregressive orders.
#' @return 
#' The function returns the autoregressive coefficients matrix theta and variance weights \eqn{H}. Rows of the matrix theta represent regimes
#' @details
#' The TAR model is given by \deqn{X_t=a_0^{(j)} + \sum_{i=1}^{k_j}a_i^{(j)}X_{t-i}+h^{(j)}e_t} when \eqn{Z_t\in (r_{j-1},r_j]} for som \eqn{j} (\eqn{j=1,\cdots,l}). 
#' the \eqn{\{Z_t\}} is the threshold process, \eqn{l} is the number of regimes, \eqn{k_j} is the autoregressive order in the regime \eqn{j}. \eqn{a_i^{(j)}} with \eqn{i=0,1,\cdots,k_j} denote the autoregressive coefficients, while \eqn{h^{(j)}} denote the variance weights. \eqn{\{e_t\}} is the Gaussian white noise process \eqn{N(0,1)}.
#' @author Hanwen Zhang <hanwenzhang at usantotomas.edu.co>
#' @param Z The threshold series
#' @param X The series of interest
#' @param l The number of regimes.
#' @param r The vector of thresholds for the series \eqn{\{Z_t\}}.
#' @param K The vector containing the autoregressive orders of the \eqn{l} regimes.
#' @param n.sim Number of iteration for the Gibbs Sampler
#' @param p.burnin Percentage of iterations used for burn-in
#' @param n.thin Thinnin factor for the Gibbs Sampler
#' 
#' @references 
#' Nieto, F. H. (2005), \emph{Modeling Bivariate Threshold Autoregressive Processes in the Presence of Missing Data}. Communications in Statistics. Theory and Methods, 34; 905-930
#' @seealso \code{\link{LS.norm}}
#' @examples 
#'  # Example 1, TAR model with 2 regimes
#' Z<-arima.sim(n=500,list(ar=c(0.5)))
#' l <- 2
#' r <- 0
#' K <- c(2,1)
#' theta <- matrix(c(1,-0.5,0.5,-0.7,-0.3,NA), nrow=l)
#' H <- c(1, 1.5)
#' X <- simu.tar.norm(Z,l,r,K,theta,H)
#' # res <- Param.norm(Z,X,l,r,K)
#'
#' # Example 2, TAR model with 3 regimes
#' Z<-arima.sim(n=300, list(ar=c(0.5)))
#' l <- 3
#' r <- c(-0.6, 0.6)
#' K <- c(1, 2, 1)
#' theta <- matrix(c(1,0.5,-0.5,-0.5,0.2,-0.7,NA, 0.5,NA), nrow=l)
#' H <- c(1, 1.5, 2)
#' X <- simu.tar.norm(Z, l, r, K, theta, H)
#' # res <- Param.norm(Z,X,l,r,K)



Param.norm <- function(Z, X, l, r, K, n.sim=500, p.burnin=0.2, n.thin=3){
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
  J.k <- J[-(1:k)]
  nj.k <- table(J[-(1:k)])
  
  
  # Priors 
  theta_0j <- matrix(0,l,k+1)
  V_0j <- diag(rep(0.01,k+1))   
  alpha <- 2
  beta <- 3
  

  
  W_j<-function(j){
    W.j<-matrix(NA,nj.k[j],K[j]+1)
    W.j[,1]<-1
    ss<-0
    for(i in 1:(N-k)){
      if(J.k[i]==j){
        ss <- ss+1
        W.j[ss,-1] <- X[(i+k-1):(i+k-K[j])]
      } }
    X.j<-matrix(X[which(J.k==j)+k])
    list(W.j=W.j,X.j=X.j)
  }
  
  sample_theta <- function(j,h2){
    AA<-W_j(j)
    Sigma <- V_0j[1:(K[j]+1),1:(K[j]+1)]+ t(AA$W.j)%*%AA$W.j/h2[j]       
    mu <- solve(Sigma)%*%(t(AA$W.j)%*%AA$X.j/h2[j]+V_0j[1:(K[j]+1),1:(K[j]+1)]%*%theta_0j[j,1:(K[j]+1)])
    rmvnorm(1,mu,solve(Sigma))
  }
  sample_h2 <- function(j,theta_j){
    theta_j<-matrix(theta_j)
    AA<-W_j(j)
    alfa.j <- alpha+nj.k[j]/2
    beta.j <- beta+sum((AA$X.j-AA$W.j%*%theta_j)^2)/2
    1/rgamma(1,shape=alfa.j,scale=1/beta.j)
  }
  N.sim <- n.sim
  res_theta_j <- array(NA,c(N.sim,l,k+1))
  res_theta_j[1,,] <- 0  ### valores iniciales para iniciar Gibbs
  
  res_h_j <- matrix(NA,N.sim,l)
  res_h_j[1,] <- 1  ### valores iniciales para iniciar Gibbs
  
  
  pbGibbs <- txtProgressBar(min = 0, max = N.sim-1, style = 3)
  
  for(rr in 2:N.sim){
    for(j in 1:l){
      res_theta_j[rr,j,1:(K[j]+1)]<-sample_theta(j,(res_h_j[rr-1,])^2)
      #res_theta_j[rr,j,1:(K[j]+1)] <- A.real[j,1:(K[j]+1)]
    }
    for(jj in 1:l){
      res_h_j[rr,jj]<-sqrt(sample_h2(jj,res_theta_j[rr,jj,1:(K[jj]+1)]))
    }
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pbGibbs, rr)
  } # Here ends Gibbs
  close(pbGibbs)
  
  iter.final <- seq(p.burnin*n.sim, n.sim, by=n.thin)
  
  h.est <- matrix(NA,l,4)
  
  theta.matrix <- vector("list",l)
  for(j in 1:l){
    aux <- matrix(NA,K[j]+1,4)
    aux[1:(K[j]+1),1] <- colMeans(res_theta_j[iter.final,j,1:(K[j]+1)])
    aux[1:(K[j]+1),2] <- apply(res_theta_j[iter.final,j,1:(K[j]+1)],2,sd)
    aux[1:(K[j]+1),3] <- apply(res_theta_j[iter.final,j,1:(K[j]+1)], 2, function(x) quantile(x,0.025))
    aux[1:(K[j]+1),4] <- apply(res_theta_j[iter.final,j,1:(K[j]+1)], 2, function(x) quantile(x,0.975))
    colnames(aux) <- c("Estimate","Std.Err","Limit inferior","Limit superior")
    rownames(aux) <- paste("lag", 0:K[j], sep="")
    theta.matrix[[j]] <- aux
    h.est[j,1] <- mean(res_h_j[iter.final,j])
    h.est[j,2] <- sd(res_h_j[iter.final,j])
    h.est[j,3:4] <- quantile(res_h_j[iter.final,j], c(0.025,0.975))
  }
  colnames(h.est) <- c("Estimate","Std.Err","Limit inferior","Limit superior")
  rownames(h.est) <- paste("Regime", 1:l, sep="")
  
  indd <- c("1st","2nd","3rd")
   if(l>3){
    indd <- c(c("1st","2nd","3rd"), paste(4:l,"th",sep=""))
  }
  
  for(j in 1:l){
    cat("\n\nThe autoregressive coefficients of the",indd[j],"regime are:","\n") 
        print(theta.matrix[[j]])
  }
  cat("\n\nThe variance weights are:", "\n")
  print(h.est)
  
  list(theta.value=res_theta_j[iter.final,,], h.value=res_h_j[iter.final,])
}
