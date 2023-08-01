update_random_intercept <- function(prm, Y) {
  Ytilde <- Y - tcrossprod(rep(1, nrow(Y)), prm$alpha) -
    prm$Bt_eta - # time-varying main effect of chemicals latent factor
    prm$eta_int %*% prm$B # linear effect terms: interaction with coveriates, main effects of covariates
  
  update_random_intercept_cpp(prm$xi, Ytilde, prm$id, prm$uid,
                prm$Sigma, prm$sigmay_sqinv, prm$nu_sqinv,
                ncol(Y), length(prm$uid))
  # Update between subject variance
  # TODO: tune priors on a, b
  a <- 3.2 + length(prm$uid)/2
  b <- 1/50 + sum(apply(prm$xi, 1, function(x) t(x) %*% prm$Sigmainv %*% x))/2
  prm$nu_sqinv <- rgamma(1, shape=a, rate=b)
  stopifnot(prm$nu_sqinv > .Machine$double.eps)
  
  return(prm)
}