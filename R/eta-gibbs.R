# Only for when interactions between factors are not considered
# Otherwise, use eta-nuts
update_eta_gibbs <- function(prm, Y, X, Z, Z_int, n, K, p, t, q, q_int) {
  if (is.null(Z)) Z <- matrix(0, 0, 0)
  if (is.null(Z_int)) Z_int <- matrix(0, 0, 0)
  prm[["eta"]] <- update_eta_gibbs_cpp(prm$eta, prm$B, prm$Theta, prm$Sigmainv, prm$sigmax_sqinv,
                                       Y, X, K, p, t, n, q, q_int, prm$id, Z, Z_int)
  
  return(prm)
}