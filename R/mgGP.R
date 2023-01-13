##### mgGP #####
mgGP <- function(X, Y, S,
                 VI.settings=list(maxit=100, K=10, reltol=eps),
                 priors.para=list(alpha0=0.5,R0=NULL,mu0=NULL,v0=NULL,W0=NULL),
                 GP.settings=list(nu=2.5, g=eps, theta.init=0.1, theta.lower=eps, theta.upper=100),
                 parallel=FALSE, n.cores=detectCores(), trace=FALSE){
  
  ##### VI setting #####
  
  # number of iterations
  if(is.null(VI.settings$maxit)){
    maxit <- 100
  }else{
    maxit <- VI.settings$maxit
  }
  
  # relative convergence tolerance
  if(is.null(VI.settings$reltol)){
    reltol <- 1e-10
  }else{
    reltol <- VI.settings$reltol
  }
  
  # truncated number
  if(is.null(VI.settings$K)){
    K <- 10 
  }else{
    K <- VI.settings$K
  }
  
  d <- ncol(X)
  p <- ncol(S)
  n <- nrow(X)
  N <- nrow(S)
  
  ##### prior setting #####
  if(is.null(priors.para$alpha0)){
    alpha0 <- 0.5
  }else{
    alpha0 <- priors.para$alpha0
  }
  
  if(is.null(priors.para$R0)){
    R0 <- solve(cov(S))
  }else{
    R0 <- priors.para$R0
  }
  
  if(is.null(priors.para$mu0)){
    mu0 <- colMeans(S)
  }else{
    mu0 <- priors.para$mu0
  }
  
  if(is.null(priors.para$v0)){
    v0 <- p
  }else{
    v0 <- priors.para$v0
  }
  
  if(is.null(priors.para$W0)){
    W0 <- R0/p
  }else{
    W0 <- priors.para$W0
  }
  
  ##### GP setting #####
  
  # smoothing parameter for matern kernel
  if(is.null(GP.settings$nu)){
    nu <- 4.5 
  }else{
    nu <- GP.settings$nu
  }
  
  # nugget
  if(is.null(GP.settings$g)){
    g <- eps
  }else{
    g <- GP.settings$g
  }
  
  # theta.init
  if(is.null(GP.settings$theta.init)){
    theta.init <- 0.1
  }else{
    theta.init <- GP.settings$theta.init
  }
  if(length(theta.init)==1) theta.init <- rep(theta.init,d)
  
  # theta.lower
  if(is.null(GP.settings$theta.lower)){
    theta.lower <- eps
  }else{
    theta.lower <- GP.settings$theta.lower
  }
  if(length(theta.lower)==1) theta.lower <- rep(theta.lower,d)
  
  # theta.upper
  if(is.null(GP.settings$theta.upper)){
    theta.upper <- 3
  }else{
    theta.upper <- GP.settings$theta.upper
  }
  if(length(theta.upper)==1) theta.upper <- rep(theta.upper,d)
  
  ##### initialization #####
  if(parallel) {
    cl <- makeCluster(n.cores)
    registerDoParallel(cl)
  }
  tau2.init <- var(c(Y))
  
  theta <- matrix(theta.init, ncol=d, nrow=K)
  tau2 <- rep(tau2.init, K)
  # q_Z <- matrix(1/K,ncol=K,nrow=N)
  q_Z <- matrix(c(1,rep(0,K-1)),ncol=K,nrow=N,byrow=TRUE)
  
  a_u <- rep(1,K)
  b_u <- rep(0,K)
  a_mu <- matrix(0, ncol=p, nrow=K)
  b_mu <- array(0, c(K,p,p))
  a_R <- array(0, c(K,p,p))
  for(k in 1:K) a_R[k,,] <- W0
  b_R <- rep(v0, K)
  log_rho <- matrix(0, ncol=K, nrow=N)
  elbo <- rep(0,maxit)
  i <- 0
  relerr <- Inf
  val <- 0
  max.fg <- FALSE
  
  pb <- txtProgressBar(min = 0, max = maxit, initial = 0, style=3) 
  
  tick <- proc.time()
  ##### E-step: update the variational posterior distribution #####
  while(i < maxit & (relerr < 0 | !max.fg | relerr > reltol * (abs(val) + reltol))){
    i <- i+1
    setTxtProgressBar(pb, i)
    
    for(ii in 1:5){
      tock0 <- proc.time()
      ##### (E1) update q(u_k): Beta(a_u, b_u) #####
      for(k in 1:(K-1)) {
        a_u[k] <- sum(q_Z[,k]) + 1
        b_u[k] <- sum(q_Z[,(k+1):K]) + alpha0
      }
      tock1 <- proc.time()
      if(trace) cat("== update q(u_k):", (tock1 - tock0)[3], "second==\n")
      
      ##### (E2) update q(mu_k): N(a_mu, b_mu) #####
      for(k in 1:K) {
        R_k1 <- a_R[k,,] * b_R[k] # E[R_k]
        R_k2 <- R_k1 %*% t(S) %*% q_Z[,k]
        R_k3 <- sum(q_Z[,k]) * R_k1
        a_mu[k,] <- solve(R0+R_k3, R0%*%mu0+R_k2)
        b_mu[k,,] <- solve(R0+R_k3)
      }
      tock2 <- proc.time()
      if(trace) cat("== update q(mu_k):", (tock2 - tock1)[3], "second==\n")
      
      ##### (E3) update q(R_k): W(a_R, b_R) #####
      for(k in 1:K) {
        u_k1 <- sum(q_Z[,k])
        d_S <- t(S) - a_mu[k,]
        u_k2 <- rep(0, p)
        for(j in 1:N) u_k2 <- u_k2 + q_Z[j,k] * d_S[,j,drop=FALSE] %*% t(d_S[,j,drop=FALSE])
        u_k2 <- u_k2 + sum(q_Z[,k])*b_mu[k,,]
        # u_k2 <- d_S %*% diag(q_Z[,k]) %*% t(d_S) + sum(q_Z[,k])*b_mu[k,,]
        a_R[k,,] <- solve(solve(W0) + u_k2)
        b_R[k] <- v0+u_k1
      }
      tock3 <- proc.time()
      if(trace) cat("== update q(R_k):", (tock3 - tock2)[3], "second==\n")
      
      ##### (E4) update q(z_s): multinomial #####
      E_logu <- digamma(a_u) - digamma(a_u+b_u)
      E_log1_u <- digamma(b_u[1:(K-1)]) - digamma(a_u[1:(K-1)]+b_u[1:(K-1)])
      
      if(FALSE){
        log_rho <- foreach(k = 1:K, .combine = "cbind", .packages = c("CholWishart", "plgp"), .export = "matern.kernel") %dopar% {
          E_logRk <- mvdigamma(b_R[k]/2, p) + p*log(2) + determinant(a_R[k,,],logarithm=TRUE)$modulus
          d_S <- t(S) - a_mu[k,]
          a_ks <- apply(d_S, 2, function(x){
            sum(diag((x %*% t(x)+b_mu[k,,]) %*%  (a_R[k,,]*b_R[k])))
          })
          a_ks <- drop(E_logRk)-a_ks-p*log(2*pi)
          
          Phi <- matern.kernel(X,theta[k,],nu=nu) + g*diag(1,n)
          Phi.chol <- chol(Phi)
          logPhi.det <- determinant(Phi.chol,logarithm=TRUE)$modulus*2
          b_ks <- apply(t(solve(Phi.chol))%*%t(Y), 2, function(x) sum(x^2))
          b_ks <- -n*log(2*pi)-n*log(tau2[k])-drop(logPhi.det) - b_ks/tau2[k]
          
          return(E_logu[k] + sum(E_log1_u[1:(k-1)]) + 
                   1/2*(a_ks+b_ks))
        }
      }else{
        for(k in 1:K) {
          E_logRk <- mvdigamma(b_R[k]/2, p) + p*log(2) + determinant(a_R[k,,],logarithm=TRUE)$modulus
          d_S <- t(S) - a_mu[k,]
          # a_ks <- rep(0, N)
          # for(j in 1:N)  a_ks[j] <- sum(diag((d_S[,j] %*% t(d_S[,j])+b_mu[k,,]) %*%  (a_R[k,,]*b_R[k])))
          
          a_ks <- apply(d_S, 2, function(x){
            sum(diag((x %*% t(x)+b_mu[k,,]) %*%  (a_R[k,,]*b_R[k])))
          })
          
          a_ks <- drop(E_logRk)-a_ks-p*log(2*pi)
          
          Phi <- matern.kernel(X,theta[k,],nu=nu) + g*diag(1,n)
          Phi.chol <- chol(Phi)
          logPhi.det <- determinant(Phi.chol,logarithm=TRUE)$modulus*2
          
          b_ks <- apply(t(solve(Phi.chol))%*%t(Y), 2, function(x) sum(x^2))
          b_ks <- -n*log(2*pi)-n*log(tau2[k])-drop(logPhi.det) - b_ks/tau2[k]
          
          # Phi.inv <- solve(Phi)
          # logPhi.det <- determinant(Phi,logarithm=TRUE)$modulus
          # b_ks <- rep(0, N)
          # for(j in 1:N)  {
          #   chol(Phi) %*% 
          #   b_ks[j] <- -n*log(2*pi)-n*log(tau2[k])-drop(logPhi.det)-t(Y[j,])%*%Phi.inv%*%Y[j,]/tau2[k]
          # }
          
          log_rho[,k] <- E_logu[k] + sum(E_log1_u[1:(k-1)]) + 
            1/2*(a_ks+b_ks)
        }
      }
      
      # q_Z <- t(apply(exp(log_rho), 1, function(x) x/sum(x)))
      q_Z <- t(apply(log_rho, 1, function(x) exp(x-max(x))/sum(exp(x-max(x)))))
      tock4 <- proc.time()
      if(trace) cat("== update q(z_s):", (tock4 - tock3)[3], "second==\n")
    }  
    
    Z <- apply(q_Z,1,which.max)
    
    ##### M-step: maximize L(tau2,theta) #####
    
    # for(ii in 1:5){
    #   ##### (M1) update tau2 ####
    #   for(k in 1:K){
    #     Phi <- matern.kernel(X,theta[k,],nu=nu) + g*diag(1,n)
    #     Phi.inv <- solve(Phi)
    #     numerator <- sum(diag(Y%*%Phi.inv%*%t(Y)) * q_Z[,k])
    #     denominator <- sum(q_Z[,k]) * n
    #     tau2[k] <- numerator/denominator
    #   }
    #   tau2[tau2==0] <- 1e-8
    #   
    #   ##### (M2) update theta ####   
    #   for(k in 1:K){
    #     neg.logl <- function(para){
    #       Phi <- matern.kernel(X, para, nu=nu) + g*diag(1,n)
    #       logPhi.det <- determinant(Phi, logarithm=TRUE)$modulus
    #       neg.ll <- logPhi.det*sum(q_Z[,k])+sum(diag(Y%*%Phi.inv%*%t(Y))*q_Z[,k])/tau2[k]
    #       return(neg.ll)
    #     }
    #     theta.out <- optim(theta.init, neg.logl,
    #                        method="L-BFGS-B", lower=theta.lower, upper=theta.upper)
    #     
    #     theta[k,] <- theta.out$par
    #   }
    # }
    ##### (M) update theta ####
    if(parallel){
      foreach.out <- foreach(k = 1:K, .packages = "plgp", .combine = "rbind", .export = c("matern.kernel","eps")) %dopar% {
        neg.logl <- function(para){
          Phi <- matern.kernel(X, para, nu=nu) + g*diag(1,n)
          Phi.chol <- chol(Phi)
          logPhi.det <- determinant(Phi.chol,logarithm=TRUE)$modulus*2
          yPinvy <- apply(t(solve(Phi.chol))%*%t(Y), 2, function(x) sum(x^2))
          neg.ll <- sum(q_Z[,k])/2*(logPhi.det+n*log(sum(yPinvy * q_Z[,k])+eps))
          return(neg.ll)
        }
        theta.out <- optim(theta.init, neg.logl,
                           method="L-BFGS-B", lower=theta.lower, upper=theta.upper)
        # theta.out <- optim(theta[k,], neg.logl,
        #                    method="L-BFGS-B", lower=theta.lower, upper=theta.upper)
        
        theta.return <- theta.out$par 
        elbo.return <- - theta.out$value + n*sum(q_Z[,k])*log(n*sum(q_Z[,k]))/2  # ELBO: Eq[log p]
        
        Phi <- matern.kernel(X,theta.return,nu=nu) + g*diag(1,n)
        Phi.chol <- chol(Phi)
        yPinvy <- apply(t(solve(Phi.chol))%*%t(Y), 2, function(x) sum(x^2))
        numerator <- sum(yPinvy * q_Z[,k])
        denominator <- sum(q_Z[,k]) * n
        tau2.return <- numerator/denominator
        if(tau2.return<eps) tau2.return <- eps
        
        return(c(theta.return, tau2.return, elbo.return))
      }
      theta <- foreach.out[,1:d,drop=FALSE]
      tau2 <- foreach.out[,d+1]
      elbo[i] <- elbo[i] + sum(foreach.out[,d+2])
    }else{
      for(k in 1:K){
        neg.logl <- function(para){
          Phi <- matern.kernel(X, para, nu=nu) + g*diag(1,n)
          # Phi.inv <- solve(Phi)
          # logPhi.det <- determinant(Phi, logarithm=TRUE)$modulus
          # yPinvy <- apply(Y,1,function(x) t(x)%*% Phi.inv %*% (x))
          Phi.chol <- chol(Phi)
          logPhi.det <- determinant(Phi.chol,logarithm=TRUE)$modulus*2
          yPinvy <- apply(t(solve(Phi.chol))%*%t(Y), 2, function(x) sum(x^2))
          neg.ll <- sum(q_Z[,k])/2*(logPhi.det+n*log(sum(yPinvy * q_Z[,k])+eps))
          return(neg.ll)
        }
        theta.out <- optim(theta.init, neg.logl,
                           method="L-BFGS-B", lower=theta.lower, upper=theta.upper)
        # theta.out <- optim(theta[k,], neg.logl,
        #                    method="L-BFGS-B", lower=theta.lower, upper=theta.upper)
        print(theta.out)
        theta[k,] <- theta.out$par 
        elbo[i] <- elbo[i] - theta.out$value + n*sum(q_Z[,k])*log(n*sum(q_Z[,k]))/2 # ELBO: Eq[log p]
        
        Phi <- matern.kernel(X,theta[k,],nu=nu) + g*diag(1,n)
        Phi.inv <- solve(Phi)
        yPinvy <- apply(Y,1,function(x) t(x)%*% Phi.inv %*% (x))
        numerator <- sum(yPinvy * q_Z[,k])
        denominator <- sum(q_Z[,k]) * n
        tau2[k] <- numerator/denominator
        if(tau2[k]<eps) tau2[k] <- eps
      }
    }
    elbo[i] <- elbo[i] - n*sum(q_Z)*(log(2*pi)+1)/2
    # print(theta)
    # print(tau2)
    
    tock5 <- proc.time()
    if(trace) cat("== M-step: estimating theta and tau2:", (tock5 - tock4)[3], "seconds ==\n")
    
    ##### ELBO ####
    ## -log(q(u))
    entropy.q <- lbeta(a_u,b_u)[1:(K-1)] - (a_u[1:(K-1)]-1)*E_logu[1:(K-1)] - (b_u[1:(K-1)]-1)*E_log1_u
    elbo[i] <- elbo[i] + sum(entropy.q)
    ## -log(q(mu))
    for(k in 1:K){
      entropy.q <- p*log(2*pi)/2+determinant(b_mu[k,,],logarithm=TRUE)$modulus/2+p/2
      elbo[i] <- elbo[i] + entropy.q
    }
    ## -log(q(R_k))
    for(k in 1:K){
      entropy.q <- (p+1)/2*determinant(a_R[k,,],logarithm=TRUE)$modulus + 
        1/2*p*(p+1)*log(2)+lmvgamma(b_R[k]/2, p)-
        (b_R[k]-p-1)/2*mvdigamma(b_R[k]/2, p)+b_R[k]*p/2
      elbo[i] <- elbo[i] + entropy.q
    }
    ## -log(q(z_s))
    # entropy.q <- apply(q_Z,1,function(x) {
    #   x <- pmin(pmax(x,eps),1-eps)
    #   xx <- 0
    #   for(k in 1:K){
    #     xx <- xx + sum(dbinom(0:K,K,x[k])*log(factorial(0:K)))
    #   }
    #   -log(factorial(K))-sum(K*x*log(x))+xx
    # })
    truncate.qZ <- pmin(pmax(q_Z,eps),1-eps)
    entropy.q <- -sum(truncate.qZ*log(truncate.qZ))
    elbo[i] <- elbo[i] + entropy.q
    plot(1:i, elbo[1:i])
    
    if(i > 1) {
      relerr <- (elbo[i] - elbo[i-1])/abs(elbo[i-1])
      val <- elbo[i-1]
      max.fg <- all(elbo[1:(i-1)] < elbo[i])
    }
  }
  
  setTxtProgressBar(pb, maxit)
  close(pb)
  tock <- proc.time()
  if(parallel) stopCluster(cl)
  
  return(list(theta=theta,tau2=tau2,nu=nu,g=g,q_Z=q_Z,X=X,  
              d=d,p=p,n=n,N=N, elbo=elbo, time.elapsed=tock-tick, 
              VI.settings=list(maxit=maxit, K=K),
              priors.para=list(alpha0=alpha0,R0=R0,mu0=mu0,v0=v0,W0=W0),
              GP.settings=list(nu=nu, g=g, theta.init=theta.init, theta.lower=theta.lower, theta.upper=theta.upper)))
}

