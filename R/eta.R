# Only for when interactions between factors are not considered
# Otherwise, use eta-nuts
update_eta_gibbs <- function(prm, Y, X, Z, Z_int, n, K, p, t, q, q_int) {
  if (is.null(Z)) Z <- matrix(0, 0, 0)
  if (is.null(Z_int)) Z_int <- matrix(0, 0, 0)
  Ytilde <- Y - prm$alpha
  update_eta_gibbs_cpp(prm[["eta"]], prm$B, prm$Theta, 
                       prm$Sigmainv, prm$sigmax_sqinv,
                       Ytilde, X, prm$id, Z, Z_int, 
                       K, p, t, n, q, q_int)
  
  return(prm)
}

update_eta_mh <- function(prm, Y, X, Z, Z_int, n, K, p, t, q, q_int, 
                          s, adaptiveM, adaptiveMWG) {
  if (is.null(Z)) Z <- matrix(0, 0, 0)
  if (is.null(Z_int)) Z_int <- matrix(0, 0, 0)
  Ytilde <- Y - prm$alpha
  Omega <- aperm(prm$Omega_tensor@data, c(2, 3, 1)) # for tensor operation in C++
  Omega[is.na(Omega)] <- 0
  update_eta_mh_cpp(prm$eta, prm$B, prm$Theta, Omega,
                    prm$Sigmainv, prm$sigmax_sqinv,
                    Ytilde, X, Z, Z_int,
                    prm$uid, prm$id,
                    K, p, t, n, q, q_int,
                    prm$eta_n_accepted,
                    prm$eta_eps,
                    prm$eta_A, prm$eta_b, s, 
                    adaptiveM, adaptiveMWG)
  return(prm)
}