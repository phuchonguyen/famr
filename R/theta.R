#' Sample the factor loading matrix using a Multiplicative Gamma Process prior 
#' 
#' @param prm list of current parameters
#' @param Y n x t matrix of responses
#' @param X n x p matrix of correlated covariates
#' @param K number of latent factor
#' @param v1 scalar or matrix of hyperparams for the MGP
#' @param v2 scalar or matrix for hyperparams for the MGP
update_Theta_MGP <- function(prm, Y, X, K, v1, v2) {
  p <- ncol(X)
  phi <- prm$phi
  delta <- prm$delta
  tau <- cumprod(delta)
  Theta <- update_Theta_MGP_cpp(prm$eta, prm$sigmax_sqinv, phi, delta, tau, K, p, X)
  phi <- update_phi_MGP_cpp(Theta, tau, K, p, v1, v2)
  delta <- update_delta_MGP_cpp(delta, tau, Theta, phi, K, p)
  
  prm[["Theta"]] <- Theta
  prm[["phi"]] <- phi
  prm[["delta"]] <- delta
  
  return(prm)
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


#' Integrate out eta and use Metropolis-Hasting instead
update_Theta_normal_mh <- function(prm, Y, X, K) {
  accepted <- update_Theta_normal_mh_cpp(prm$Theta, prm$sigmax_sqinv, X, prm$Theta_eps)
  prm[["Theta_n_accepted"]] <- prm$Theta_n_accepted + accepted
  return(prm)
}