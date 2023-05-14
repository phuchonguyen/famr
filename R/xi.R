update_random_intercept <- function(prm, Y) {
  Ytilde <- Y - tcrossprod(rep(1, nrow(Y)), prm$alpha) -
    prm$eta_int%*%prm$B -
    prm$eta_quad %*% prm$Omega
  
  update_random_intercept_cpp(prm$xi, Ytilde, prm$id, prm$uid,
                prm$Sigma, prm$sigmay_sqinv, prm$nu_sqinv,
                ncol(Y), length(prm$uid))
  # Update between subject variance
  a <- 0.1 + length(prm$uid)/2
  b <- 0.001 + sum(apply(prm$xi, 1, function(x) t(x) %*% prm$Sigmainv %*% x))/2
  prm$nu_sqinv <- rgamma(1, shape=a, rate=b)
  stopifnot(prm$nu_sqinv > .Machine$double.eps)
  
  return(prm)
}