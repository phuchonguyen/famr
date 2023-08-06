update_Lambda <- function(Y, time, K, L, q, prm) {
  S_inv <- prm$Sigmainv * 1/(1/prm$sigmay_sqinv + 0) # TODO update: marginalize out linear effect B
  Ytilde <- Y - tcrossprod(rep(1, nrow(Y)), prm$alpha) - prm$eta_int %*% prm$B
  update_Lambda_cpp(prm$Bt_Lambda, prm$Bt_U, 
                    Ytilde, prm$eta[prm$numeric_id,,drop=FALSE],
                    S_inv, prm$Bt_phi, cumprod(prm$Bt_delta), 
                    time-min(time), nrow(Y), ncol(Y), K, L)
  prm[['Bt_phi']] <- update_phi_MGP_cpp(prm$Bt_Lambda, 
                                        cumprod(prm$Bt_delta), 
                                        L, q)
  prm[['Bt_delta']] <- update_delta_MGP_cpp(prm$Bt_delta, cumprod(prm$Bt_delta), 
                                            prm$Bt_Lambda, prm$Bt_phi, L, q, 
                                            prm$Bt_a1, prm$Bt_a2)
  # a1_out <- update_a1_MGP(prm$a1, delta, prm$a1_eps, prm$a1_n_accepted)
  # a2_out <- update_a2_MGP(prm$a2, delta, prm$a2_eps, prm$a2_n_accepted, K)
  
  return(prm)
}


# L: number of latent GP
# q: number of outcomes
# TT: number of unique time points
# n: number of unique subjects
update_U <- function(prm, Y, time, K, L, q, TT, n) {
  S0_inv <- chol2inv(chol(diag(1, TT, TT) + matrix(1, TT, TT)*1/prm$nu_sqinv)) # TODO nu_sqinv
  S_inv <- kronecker(prm$Sigmainv, S0_inv)
  stopifnot(all(dim(S_inv), c(TT*q, TT*q)))
  Ytilde <- Y - tcrossprod(rep(1, nrow(Y)), prm$alpha) - prm$eta_int %*% prm$B
  
  update_U_cpp(prm$Bt_U, Ytilde, prm$eta, prm$Bt_Lambda, prm$numeric_id-1,
               time-min(time), S_inv, prm$C_inv,
               L, K, TT, q, n)
  
  return(prm)
}

# Update kappa discrete options
# Equal prior probabilities over the options of kappa
# Returns index of the chosen kappa
update_kappa_discrete <- function(Ci_all, ldetC_all, U, L, K, n_kappa_opts) {
  log_prob <- rep(NA, n_kappa_opts)
  for (i in 1:n_kappa_opts) {
    log_prob[i] <-  llike_kappa(Ci_all[[i]], ldetC_all[[i]], U, L, K)
  }
  # From Kelly Moran's code: https://github.com/kelrenmor/bs3fa/blob/master/R/sample_lind.R
  maxlg <- max(log_prob)
  # https://stats.stackexchange.com/questions/66616/converting-normalizing-very-small-likelihood-values-to-probability
  probls <- exp(log_prob - maxlg)
  probls <- probls/sum(probls)
  which_ls <- sample(x=(1:n_kappa_opts), 1, replace=T, prob=probls)
  return(which_ls)
}