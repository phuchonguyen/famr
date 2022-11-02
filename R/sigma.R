update_Sigma_IW <- function(prm, Y, Z, Z_int, binary) {
  Ytilde <- Y - prm$eta
  return(NULL)
}

update_Sigma_HIW <- function(prm, Y, Z, Z_int, binary) {
  t <- ncol(Y)
  n <- nrow(Y)
  v <- 2  # param in prior of Sigma, induces uniform prior on marginals of correlations
  d <- rep(10, t)  # hyperparam in prior of Sigma, large values induces weakly informative priors
  l <- prm$l
  K <- ncol(prm$eta)
  
  # Update Sigma
  # NOTE: variance/covariance might blow up when Y is not centered!
  # NOTE: Sigma might be singular when ...
  df <- n + v + t - 1
  # For MBSP on B model
  eta_int <- get_eta_int(prm$eta, K, Z, Z_int, prm$id)
  q <- ncol(eta_int)
  B0 <- matrix(0, q, t)
  Ytilde <- Y - rep(1, nrow(Y))%*%t(prm$alpha)
  #M <- Ytilde - eta_int%*%prm$B
  #SS <- t(M)%*%M
  #S <- 2*v*diag(1/l) + SS + t(prm$B-B0)%*%diag(1/prm$psi)%*%(prm$B-B0)
  Un <- solve(crossprod(eta_int, eta_int) + diag(1/prm$psi, q, q))
  Mn <- Un %*% (crossprod(eta_int, Ytilde) + B0)
  S <- diag(l, t, t) + crossprod(Ytilde, Ytilde) + crossprod(B0,B0) - t(Mn)%*%solve(Un)%*%Mn # if psi=rep(1,q)
  Sigma <- CholWishart::rInvWishart(1, df = df, Sigma = S)[,,1]
  if (t == 1) {
    Sigma = matrix(Sigma, 1, 1)
  }
  
  if (sum(binary) > 0) {
    d <- 1/sqrt(diag(Sigma))
    d[binary==0] <- 1
    Sigma <- diag(d)%*%Sigma%*%diag(d)
  }
  Sinv <- tryCatch(solve(Sigma), error=function(e) {
    cat('\nwhile sampling Sigma:', message(e), '\n')
    print(Sigma)
    cat('\n Var(Ytilde) = ', apply(Ytilde, 2, var, na.rm=TRUE))
  })
  #Update l
#   for (i in sample(1:t)) {
#     l[i] <- 1/rgamma(1, 0.5*(v + t), v*Sinv[i,i] + 1/d[i]^2)
#   }
  
  #stopifnot(sum(is.na(l))==0)
  
  prm[["l"]] <- l
  prm[["Sigma"]] <- Sigma
  prm[["Sigmainv"]] <- Sinv
  return(prm)
}

update_sigmax_sqinv <- function(prm, Y, X, K) {
  
  p <- ncol(X)
  n <- nrow(X)
  # prior 1/sigma^2 ~ Ga(r/2, rs0/2)
  # since data has been standardized to have variance 1, choose r, s0 s.t. Pr(sigma^2 <= 1) = 0.95
  # the choices below from gpreg lecture by Surya to produce Pr(sigma <= 1/3) = 0.5
  r <- 2.5
  s0 <- 0.084 
  
  sigmax_sqinv <- rep(NA, p)
  alpha <- (n + r)*0.5
  for (j in 1:p) {
    SSRj <- sum((X[,j] - prm$eta%*%prm$Theta[j,])^2)
    betaj <- 0.5*(r*s0 + SSRj)
    sigmax_sqinv[j] <- rgamma(1, alpha, betaj)
    stopifnot(sigmax_sqinv[j] > .Machine$double.eps)
  }
  
  prm[["sigmax_sqinv"]] <- sigmax_sqinv
  return(prm)
}



