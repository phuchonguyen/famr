update_Sigma_IW_TPBN <- function(prm, Y, Z, Z_int, K, binary) {
  p <- ncol(prm$eta_int)
  p_quad <- ncol(prm$eta_quad)
  n <- nrow(Y)
  q <- ncol(Y)
  Ytilde <- Y - tcrossprod(rep(1, n), prm$alpha) - 
    prm$eta_int%*%prm$B - 
    prm$eta_quad %*% prm$Omega
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

update_sigmax_sqinv <- function(prm, X, K, s0=0.084, r=2.5) {
  p <- ncol(X)
  n <- nrow(X)
  # prior 1/sigma^2 ~ Ga(r/2, rs0/2)
  # since data has been standardized to have variance 1, default choice of r, s0 s.t. Pr(sigma^2 <= 1) = 0.95
  # the choices below from gpreg lecture by Surya to produce Pr(sigma <= 1/3) = 0.5
  sigmax_sqinv <- rep(NA, p)
  alpha <- (n + r)*0.5
  Xtilde <- X - prm$eta %*% t(prm$Theta)
  for (j in 1:p) {
    SSRj <- sum(Xtilde[,j]^2)
    betaj <- 0.5*(r*s0 + SSRj)
    sigmax_sqinv[j] <- rgamma(1, shape=alpha, rate=betaj)
    stopifnot(sigmax_sqinv[j] > .Machine$double.eps)
  }
  return(sigmax_sqinv)
}

# When there is random intercepts
update_Sigma_IW_TPBN_re <- function(prm, Y, Z, Z_int, K, binary) {
  p <- ncol(prm$eta_int)
  p_quad <- 0#ncol(prm$eta_quad)
  n <- nrow(Y)
  q <- ncol(Y)
  Ytilde <- Y - tcrossprod(rep(1, n), prm$alpha) -
    prm$eta_int%*%prm$B -
    prm$eta_quad %*% prm$Omega #-
    # prm$xi[prm$numeric_id, ]
  sinv <- 1/(1/prm$sigmay_sqinv + 1/prm$nu_sqinv)
  YtY <- crossprod(Ytilde, Ytilde) * sinv #prm$sigmay_sqinv
  BtB <- t(prm$B) %*% diag(1/prm$psi, p, p) %*% prm$B
  OtO <- 0#t(prm$Omega) %*% diag(1/prm$Omega_psi, p_quad, p_quad) %*% prm$Omega
  # ItI <- t(prm$xi) %*% prm$xi * prm$nu_sqinv
  df <- n + p + p_quad + (q + 10)
  S <- YtY + BtB + OtO + 15*(diag(1e-5 + apply(Y, 2, var), q, q)) #+ ItI
  tryCatch(solve(S), error=function(e) {print(S)})
  # Sigma <- CholWishart::rInvWishart(1, df = df, Sigma = S)[,,1]
  Sigmainv <- tryCatch(MCMCpack::rwish(df, solve(S)), error=function(e){print("rwish re :(")})
  Sigma <- tryCatch(solve(Sigmainv), error=function(e) {
    cat('\nwhile sampling Sigma re:', message(e), '\n')
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
