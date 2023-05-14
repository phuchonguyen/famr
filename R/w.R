update_W_MGP <- function(prm, Y, H) {
  q <- ncol(Y)
  n <- nrow(Y)
  Ytilde <- Y - 
    tcrossprod(rep(1, n), prm$alpha) - 
    prm$eta_int %*% prm$B - 
    prm$eta_quad %*% prm$Omega
  W <- update_Theta_MGP_cpp(prm$xi[prm$numeric_id, ,drop=F], 
                            prm$sigmay_sqinv, 
                            prm$W_phi, prm$W_delta, cumprod(prm$W_delta), 
                            H, q, Ytilde)
  W_phi <- update_phi_MGP_cpp(W, cumprod(prm$W_delta), H, q)
  W_delta <- update_delta_MGP_cpp(prm$W_delta, cumprod(prm$W_delta), 
                                W, W_phi, 
                                H, q, prm$W_a1, prm$W_a1)
  a1_out <- update_a1_MGP(prm$W_a1, W_delta, prm$a1_eps, prm$W_a1_n_accepted)
  a2_out <- update_a2_MGP(prm$W_a2, W_delta, prm$a2_eps, prm$W_a2_n_accepted, H)
  return(list(W = W, 
              W_phi = W_phi, 
              W_delta = W_delta,
              W_a1 = a1_out$a1,
              W_a2 = a2_out$a2,
              W_a1_n_accepted = a1_out$a1_n_accepted,
              W_a2_n_accepted = a2_out$a2_n_accepted
              ))
}

update_xi_mh <- function(prm, Y, H, n, q, s, adaptiveM, adaptiveMWG, adaptiveMWG_batch) {
  Ytilde <- Y - 
    tcrossprod(rep(1, nrow(Y)), prm$alpha) - 
    prm$eta_int %*% prm$B - 
    prm$eta_quad %*% prm$Omega
  update_xi_mh_cpp(prm$xi, prm$xi_A, prm$xi_b, 
                   prm$xi_n_accepted, prm$xi_eps, prm$xi_lpmf,
                   Ytilde, prm$W, prm$Sigmainv,
                   prm$uid, prm$id, H, n, s, 
                   adaptiveM, adaptiveMWG, adaptiveMWG_batch)
  return(prm)
}
