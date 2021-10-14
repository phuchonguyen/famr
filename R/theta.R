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