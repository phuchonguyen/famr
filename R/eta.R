# Use this when there is NO random intercepts
update_eta_mh <- function(prm, Y, X, Z, Z_int, n, K, p, q, p_z, p_int, 
                          s, adaptiveM, adaptiveMWG, adaptiveMWG_batch) {
  if (is.null(Z_int)) Z_int <- matrix(0, 0, 0)
  Ytilde <- Y - tcrossprod(rep(1, nrow(Y)), prm$alpha)
  if (p_z > 0) {
    Ytilde <- Ytilde - Z %*% prm$B[(nrow(prm$B)-p_z+1) : nrow(prm$B),]
  }
  Omega <- aperm(prm$Omega_tensor@data, c(2, 3, 1)) # for tensor operation in C++
  Omega[is.na(Omega)] <- 0
  update_eta_mh_cpp(prm$eta, prm$B, prm$Theta, Omega,
                    prm$Sigmainv, prm$sigmax_sqinv,
                    Ytilde, X, Z_int,
                    prm$uid, prm$id,
                    K, p, q, n, p_int,
                    prm$eta_n_accepted,
                    prm$eta_eps,
                    prm$eta_A, prm$eta_b, 
                    prm$lpmf, s, 
                    adaptiveM, adaptiveMWG, adaptiveMWG_batch)
  return(prm)
}

# For model with random intercepts
update_eta_mh_re <- function(prm, Y, X, Z, Z_int, time, n, K, p, q, p_z, p_int, 
                          s, adaptiveM, adaptiveMWG, adaptiveMWG_batch, eps_power) {
  if (is.null(Z_int)) Z_int <- matrix(0, 0, 0)
  Ytilde <- Y - tcrossprod(rep(1, nrow(Y)), prm$alpha)
  if (p_z > 0) {
    Ytilde <- Ytilde - Z %*% prm$B[(nrow(prm$B)-p_z+1) : nrow(prm$B),]
  }
  # Omega <- aperm(prm$Omega_tensor@data, c(2, 3, 1)) # for tensor operation in C++
  # Omega[is.na(Omega)] <- 0
  Omega <- array(0, dim=c(K, K, q))
  sinv <- 1/(1/prm$sigmay_sqinv + 0) # TODO nu_sqinv
  update_eta_mh_cpp(prm$eta, prm$Bt, prm$B, prm$Theta, Omega,
                    prm$Sigmainv * sinv,
                    prm$sigmax_sqinv,
                    Ytilde, X, Z_int,
                    prm$uid, prm$id,
                    time-min(time), # re-index to start at 0
                    K, p, q, n, p_int,
                    prm$eta_n_accepted,
                    prm$eta_eps,
                    prm$eta_A, prm$eta_b, 
                    prm$lpmf, s,
                    prm$eta_prop, # for DEBUGGING
                    adaptiveM, adaptiveMWG, adaptiveMWG_batch, eps_power) 
  return(prm)
}