#' Sample the factor loading matrix using a Multiplicative Gamma Process prior 
#' 
#' @param prm list of current parameters
#' @param Y n x t matrix of responses
#' @param X n x p matrix of correlated covariates
#' @param K number of latent factor
update_Theta_MGP <- function(prm, Y, X, K) {
  Theta <- update_Theta_MGP_cpp(prm$eta, prm$sigmax_sqinv, 
                                prm$phi, prm$delta, cumprod(prm$delta), 
                                K, ncol(X), X)
  phi <- update_phi_MGP_cpp(Theta, cumprod(prm$delta), K, ncol(X))
  delta <- update_delta_MGP_cpp(prm$delta, cumprod(prm$delta), 
                                Theta, phi, 
                                K, ncol(X), prm$a1, prm$a2)
  a1_out <- update_a1_MGP(prm$a1, delta, prm$a1_eps, prm$a1_n_accepted)
  a2_out <- update_a2_MGP(prm$a2, delta, prm$a2_eps, prm$a2_n_accepted, K)
  return(c(list(Theta = Theta, 
                phi = phi, 
                delta = delta),
           a1_out, a2_out))
}

update_a1_MGP <- function(a1, delta, eps, a1_n_accepted) {
  delta1 <- delta[1]
  # proposal with a truncated normal
  prop <- truncnorm::rtruncnorm(1, a=0, mean=a1, sd=eps)
  # proposal log likelihood
  lprop <- log(delta1^(prop-1) * prop * exp(-prop) / gamma(prop)) +
    log(truncnorm::dtruncnorm(a1, a=0, mean=prop, sd=eps))
  # current log likelihood
  lcur <- log(delta1^(a1-1) * a1 * exp(-a1) / gamma(a1)) +
    log(truncnorm::dtruncnorm(prop, a=0, mean=a1, sd=eps))
  # accept with log probability
  lr = lprop - lcur
  if (log(runif(1)) < lr) {
    a1 <- prop
    a1_n_accepted <- a1_n_accepted + 1
  }
  stopifnot(a1 > 0)
  return(list(a1 = a1, a1_n_accepted = a1_n_accepted))
}

update_a2_MGP <- function(a2, delta, eps, a2_n_accepted, K) {
  delta <- delta[-1]
  # proposal with a truncated normal
  prop = truncnorm::rtruncnorm(1, a=0, mean=a2, sd=eps)
  # proposal log likelihood
  lprop = -(K-1)*log(gamma(prop)) + (prop-1)*sum(log(delta)) - prop + log(prop) +
    log(truncnorm::dtruncnorm(a2, a=0, mean=prop, sd=eps))
  # current log likelihood
  lcur = -(K-1)*log(gamma(a2)) + (a2-1)*sum(log(delta)) - a2 + log(a2) +
    log(truncnorm::dtruncnorm(prop, a=0, mean=a2, sd=eps))
  # accept with log probability
  lr = lprop - lcur
  if (log(runif(1)) < lr) {
    a2 <- prop
    a2_n_accepted <- a2_n_accepted + 1
  }
  stopifnot(a2 > 0)
  return(list(a2 = a2, a2_n_accepted = a2_n_accepted))
}

#' Place Gaussian priors on the elements of the factor loading matrix
update_Theta_normal <- function(prm, Y, X, K) {
  Dinv_j <- diag(1/100, K, K)
  ete <- crossprod(prm$eta, prm$eta)
  p <- ncol(X)
  Theta <- matrix(0, p, K)
  for (j in 1:p) {
    Xj <- X[, j]
    sigmax_sqinv_j <- prm$sigmax_sqinv[j]
    Vj <- solve(Dinv_j + ete / sigmax_sqinv_j)
    mj <- crossprod(prm$eta, Xj) / sigmax_sqinv_j
    Theta[j, ] <- MASS::mvrnorm(1, mu = Vj %*% mj, Sigma = Vj)
  }
  prm[["Theta"]] <- Theta
  return(prm)
}

# DL on the COLS
update_Theta_DL_col <- function(prm, Y, X, K) {
  return(0)
}

# DL on the ROWS
update_Theta_DL <- function(prm, Y, X, K) {
  # DL prior variances
  V0 <- prm$phi^2 * prm$omega * tcrossprod(prm$tau^2, rep(1, K))
  V0[V0 < .Machine$double.eps] <- .Machine$double.eps # for numerical instability
  ete <- crossprod(prm$eta, prm$eta)
  p <- ncol(X)
  Theta <- matrix(0, p, K)
  for (j in 1:p) {
    Dinv_j <- diag(1/V0[j,], K, K)
    Xj <- X[, j]
    sigmax_sqinv_j <- prm$sigmax_sqinv[j]
    Vj <- solve(Dinv_j + ete / sigmax_sqinv_j)
    mj <- crossprod(prm$eta, Xj) / sigmax_sqinv_j
    Theta[j, ] <- MASS::mvrnorm(1, mu = Vj %*% mj, Sigma = Vj)
  }
  prm[["Theta"]] <- Theta
  
  update_tau_DL_cpp(prm$tau, prm$Theta, prm$phi, prm$omega, p, K)
  update_phi_DL_cpp(prm$phi, prm$Theta, p, K)
  update_omega_DL_cpp(prm$omega, prm$phi, prm$Theta, prm$tau)
  return(prm)
}

# DL on the ROWS of Theta
update_Theta_DL_mh <- function(prm, Y, X, K) {
  p <- ncol(X)
  # DL prior variances (stack columns)
  V0 <- as.vector(prm$phi^2 * prm$omega * tcrossprod(prm$tau^2, rep(1, K)))
  V0[V0 < .Machine$double.eps] <- .Machine$double.eps # for numerical instability
  accepted <- update_Theta_mh_cpp(prm$Theta, prm$sigmax_sqinv, X,
                                  V0, prm$Theta_eps)
  update_tau_DL_cpp(prm$tau, prm$Theta, prm$phi, prm$omega, p, K)
  update_phi_DL_cpp(prm$phi, prm$Theta, p, K)
  update_omega_DL_cpp(prm$omega, prm$phi, prm$Theta, prm$tau)
  prm[["Theta_n_accepted"]] <- prm$Theta_n_accepted + accepted
  return(prm)
}

#' Integrate out eta and use Metropolis-Hasting instead
#' DONT USE
update_Theta_normal_mh <- function(prm, Y, X, K, s, adaptiveMWG) {
  # indep. normal with 100 variance priors
  V0 = rep(100, length(prm$Theta))
  update_Theta_mh_cpp(prm$Theta, prm$sigmax_sqinv, X, V0,
                      prm$Theta_n_accepted, prm$Theta_eps,
                      s, adaptiveMWG)
  return(prm)
}