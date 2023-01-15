update_Sigma_IW_DL <- function(prm, Y, Z, Z_int, binary) {
  K <- ncol(prm$eta)
  eta_int <- get_eta_int(prm$eta, K, Z, Z_int, prm$id)
  p <- ncol(eta_int)
  n <- nrow(Y)
  q <- ncol(Y)
  Ytilde <- Y - prm$alpha - eta_int%*%prm$B - prm$eta_quad %*% prm$Omega
  YtY <- crossprod(Ytilde, Ytilde)
  df <- n + (q + 3)
  S <- YtY + diag(10, q, q)
  tryCatch(solve(S), error=function(e) {print(S)})
  # Sigma <- CholWishart::rInvWishart(1, df = df, Sigma = S)[,,1]
  Sigmainv <- tryCatch(MCMCpack::rwish(df, solve(S)), error=function(e){print("rwish :(")})
  Sigma <- tryCatch(solve(Sigmainv), error=function(e) {
    cat('\nwhile sampling Sigma:', message(e), '\n')
    print(Sigmainv)
  })
  if (q == 1) {
    Sigma = matrix(Sigma, 1, 1)
  }
  
  if (sum(binary) > 0) {
    d <- 1/sqrt(diag(Sigma))
    d[binary==0] <- 1
    Sigma <- diag(d)%*%Sigma%*%diag(d)
  }
  prm[["Sigma"]] <- Sigma
  prm[["Sigmainv"]] <- Sigmainv
  return(prm)
}

update_Sigma_IW_TPBN <- function(prm, Y, Z, Z_int, binary) {
  K <- ncol(prm$eta)
  eta_int <- get_eta_int(prm$eta, K, Z, Z_int, prm$id)
  p <- ncol(eta_int)
  p_quad <- length(prm$Omega_psi)
  n <- nrow(Y)
  q <- ncol(Y)
  Ytilde <- Y - prm$alpha - eta_int%*%prm$B - prm$eta_quad %*% prm$Omega
  YtY <- crossprod(Ytilde, Ytilde)
  BtB <- t(prm$B) %*% diag(1/prm$psi, p, p) %*% prm$B
  OtO <- t(prm$Omega) %*% diag(1/prm$Omega_psi, p_quad, p_quad) %*% prm$Omega
  df <- n + p + p_quad + (q + 3)
  S <- YtY + BtB + OtO + diag(10, q, q)
  tryCatch(solve(S), error=function(e) {print(S)})
  # Sigma <- CholWishart::rInvWishart(1, df = df, Sigma = S)[,,1]
  Sigmainv <- tryCatch(MCMCpack::rwish(df, solve(S)), error=function(e){print("rwish :(")})
  Sigma <- tryCatch(solve(Sigmainv), error=function(e) {
    cat('\nwhile sampling Sigma:', message(e), '\n')
    print(Sigmainv)
  })
  if (q == 1) {
    Sigma = matrix(Sigma, 1, 1)
  }
  
  if (sum(binary) > 0) {
    d <- 1/sqrt(diag(Sigma))
    d[binary==0] <- 1
    Sigma <- diag(d)%*%Sigma%*%diag(d)
  }
  return(list(Sigma = Sigma, Sigmainv = Sigmainv))
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
  Ytilde <- Y - prm$alpha
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

update_sigmax_sqinv <- function(prm, Y, X, K, s0) {
  
  p <- ncol(X)
  n <- nrow(X)
  # prior 1/sigma^2 ~ Ga(r/2, rs0/2)
  # since data has been standardized to have variance 1, choose r, s0 s.t. Pr(sigma^2 <= 1) = 0.95
  # the choices below from gpreg lecture by Surya to produce Pr(sigma <= 1/3) = 0.5
  r <- 2.5
  # s0 <- 0.084 
  
  sigmax_sqinv <- rep(NA, p)
  alpha <- (n + r)*0.5
  Xtilde <- X - prm$eta%*%t(prm$Theta)
  for (j in 1:p) {
    SSRj <- sum(Xtilde[,j]^2)
    betaj <- 0.5*(r*s0 + SSRj)
    sigmax_sqinv[j] <- rgamma(1, shape=alpha, rate=betaj)
    stopifnot(sigmax_sqinv[j] > .Machine$double.eps)
  }
  return(sigmax_sqinv)
}

