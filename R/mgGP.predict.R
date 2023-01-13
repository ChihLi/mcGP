predict.mgGP <- function(fit, Y, xnew=NULL, sig2.fg=TRUE, array.fg=TRUE){
  
  X <- fit$X
  theta <- fit$theta
  tau2 <- fit$tau2
  q_Z <- fit$q_Z
  nu <- fit$nu
  g <- fit$g
  d <- fit$d
  p <- fit$p
  n <- fit$n
  N <- fit$N
  K <- fit$VI.settings$K
  
  if(is.null(xnew)){
    xnew <- X
  }
  
  n.new <- nrow(xnew)
  
  y.pred <- matrix(0, nrow=N, ncol=n.new)
  if(array.fg) y.pred.array <- y.sig2.array <- array(0, c(N,n.new,K))
  if(sig2.fg) y.sig2 <- matrix(0, nrow=N, ncol=n.new)
  
  tick <- proc.time()
  for(k in 1:K){
    Phi <- matern.kernel(X,theta[k,],nu=nu) + g*diag(1,n)
    Phi.inv <- solve(Phi)
    Phi_Xx <- matern.kernel(X,theta[k,],nu,xnew)
    DD <- t(Phi_Xx)%*%Phi.inv
    if(array.fg) y.pred.array[,,k] <- t(DD%*%t(Y))
    # y.pred <- y.pred + t(t(Phi_Xx)%*%Phi.inv%*%t(Y)) * post.p[,k]
    y.pred <- y.pred + t(DD%*%t(Y)) * q_Z[,k]
    
    if(sig2.fg) {
      y.sig2 <- y.sig2 + t(outer(tau2[k]*(1+g-pmin(1, diag(DD%*%Phi_Xx))), q_Z[,k])) +
        t(DD%*%t(Y))^2 * q_Z[,k]
      if(array.fg) y.sig2.array[,,k] <- t(outer(tau2[k]*(1+g-pmin(1, diag(DD%*%Phi_Xx))), rep(1,N)))
    }
  }
  tock <- proc.time()
  
  if(sig2.fg){
    if(array.fg){
      return(list(mean=y.pred, sig2=y.sig2-(y.pred)^2, 
                  mean.array=y.pred.array, sig2.array=y.sig2.array,
                  time.elapsed=tock-tick))
    }else{
      return(list(mean=y.pred, sig2=y.sig2-(y.pred)^2,
                  time.elapsed=tock-tick))
    }
  }else{
    if(array.fg){
      return(list(mean=y.pred, mean.array=y.pred.array, sig2.array=y.sig2.array,
                  time.elapsed=tock-tick))
    }else{
      return(list(mean=y.pred, time.elapsed=tock-tick))
    }
  }
}
