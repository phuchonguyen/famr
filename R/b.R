# Transform eta for estimating the interaction matrix
# Coefficients are order as : interactions between latent factors and covariates, covariates 
get_eta_int <- function(eta, K, Z, Z_int, id) {
  uid <- unique(id)
  # Duplicate eta for repeated outcome measurements: (Y_11,...Y_1T_1,...Y_nT_n)^T
  eta_dup <- eta[match(id, uid),,drop=F]
  if (is.null(colnames(eta_dup))) {
    colnames(eta_dup) <- paste0("e", 1:K)
  }
  eta_int <- c()
  cnames <- c() #colnames(eta_int)
  if (!is.null(Z_int)) {
    if (is.null(colnames(Z_int))) {
      colnames(Z_int) <- paste0("zint", 1:ncol(Z_int))
    }
    for (h in 1:ncol(Z_int)) {
      for (k in 1:K) {
        eta_int <- cbind(eta_int, Z_int[,h]*eta_dup[,k])
        cnames <- c(cnames, paste0('e',k, colnames(Z_int)[h]))
      }
    }
  }
  if (!is.null(Z)) {
    if (is.null(colnames(Z))) {
      colnames(Z) <- paste0("z", 1:ncol(Z))
    }
    eta_int <- cbind(eta_int, Z)
    cnames <- c(cnames, colnames(Z))
  }
  colnames(eta_int) <- cnames
  return(eta_int)
}

# Use when there is NO random intercept
update_B_TPBN <- function(prm, Y, X, K, Z, Z_int) {
  q <- ncol(Y)
  n <- nrow(Y)
  eta_int <- get_eta_int(prm$eta, K, Z, Z_int, prm$id)
  p <- ncol(eta_int)
  alpha <- update_intercept(Y - eta_int %*% prm$B - prm$eta_quad %*% prm$Omega, prm$Sigma)
  Ytilde <- Y - tcrossprod(rep(1, n), prm$alpha) - prm$eta_quad %*% prm$Omega
  B <- update_B_TPBN_cpp(eta_int, Ytilde, prm$Sigma, prm$psi, p, q)
  psi <- as.vector(update_psi_TPBN_cpp(B, prm$Sigmainv, prm$zeta, p, q))
  zeta <- update_zeta_TPBN_cpp(psi, p, global_shrink=1/(K*sqrt(n*log(n))))
  stopifnot(dim(B) == dim(prm$B))
  stopifnot(sum(is.na(B)) == 0)
  return(list(alpha = alpha,
              B = B,
              psi = psi,
              zeta = zeta,
              eta_int = eta_int))
}

# Use when there is random intercept and GP main effects
# Updates linear regression matrices
update_B_TPBN_re <- function(prm, Y, X, K, Z, Z_int) {
  q <- ncol(Y)
  n <- nrow(Y)  # number of observations total: sum_i T_i
  sinv <- 1/(1/prm$sigmay_sqinv + 1/prm$nu_sqinv)
  eta_int <- get_eta_int(prm$eta, K, Z, Z_int, prm$id)
  p <- ncol(eta_int)
  alpha <- update_intercept(Y - prm$Bt_eta - eta_int %*% prm$B, prm$Sigma/sinv)
  Ytilde <- Y - tcrossprod(rep(1, n), alpha) - prm$Bt_eta
  Ytilde <- Ytilde * sqrt(sinv) #sqrt(prm$sigmay_sqinv)
  eta_int <- eta_int * sqrt(sinv) #sqrt(prm$sigmay_sqinv)
  B <- update_B_TPBN_cpp(eta_int, Ytilde, prm$Sigma, prm$psi, p, q)
  psi <- as.vector(update_psi_TPBN_cpp(B, prm$Sigmainv, prm$zeta, p, q))
  zeta <- update_zeta_TPBN_cpp(psi, p, global_shrink=1/(K*sqrt(n*log(n))))
  stopifnot(dim(B) == dim(prm$B))
  stopifnot(sum(is.na(B)) == 0)
  return(list(alpha = alpha,
              B = B,
              psi = psi,
              zeta = zeta,
              eta_int = eta_int))
}

# Using a flat prior p(a) \prop 1
update_intercept <- function(Y, Sigma) {
  # use flat prior on alpha
  return(MASS::mvrnorm(1, 
                       mu = colMeans(Y, na.rm = TRUE), 
                       Sigma = Sigma / nrow(Y)))
}