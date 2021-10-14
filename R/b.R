# Multivariate Global Local Shrinkage prior
# Paper: https://arxiv.org/pdf/1711.07635.pdf

# Transform eta for estimating the interaction matrix
# Coefficients are order as : latent main effect, latent interactions, interact with covariates, covariates 
get_eta_int <- function(eta, K, Z, Z_int, fa_interact, id) {
  uid <- unique(id)
  eta_dup <- eta[match(id, uid), ]
  if (is.null(colnames(eta_dup))) {
    colnames(eta_dup) <- paste0("e", 1:K)
  }
  if (fa_interact) {
    eta_int <- model.matrix(~. + .^2 - 1, data = as.data.frame(eta_dup))
  } else {
    eta_int <- eta_dup
  }
  cnames <- colnames(eta_int)
  if (!is.null(Z_int)) {
    for (h in 1:ncol(Z_int)) {
      for (k in 1:K) {
        eta_int <- cbind(eta_int, Z_int[,h]*eta_dup[,k])
        cnames <- c(cnames, paste0('e',k, colnames(Z_int)[h]))
      }
    }
  }
  if (!is.null(Z)) {
    eta_int <- cbind(eta_int, Z)
    cnames <- c(cnames, colnames(Z))
  }
  colnames(eta_int) <- cnames
  return(eta_int)
}

# Get matrix (K x K_int) indicating which columns in eta_int contain the kth factor
get_k_I <- function(K_int, K, q_int, fa_interact) {
  k_I <- matrix(0, K, K_int)
  j <- K
  for (k in 1:K) {
    k_I[k,k] <- 1
    h <- k+1
    while(h <= K & fa_interact) {
      j <- j+1
      k_I[k,j] <- 1
      k_I[h,j] <- 1
      h <- h+1
    }
  }
  if (!is.null(q_int)) {
    for (h in 1:q_int) {
      for (k in 1:K) {
        j <- j+1
        k_I[k,j] <- 1
      }
    }
  }
  return(k_I)
}

update_B_TPBN <- function(prm, Y, X, K, Z, Z_int, u, a, tau) {
  eta_int <- get_eta_int(prm$eta, K, Z, Z_int, prm$fa_interact, prm$id)
  out <- update_zeta_nu_TPBN(Y, prm$B, prm$nu, prm$zeta, prm$Sigmainv, u, a, tau)
  nu <- out$nu
  zeta <- out$zeta
  
  B <- core_B_TPBN(Y, eta_int, nu, zeta, prm$Sigma)
  rownames(B) <- colnames(eta_int)
  
  stopifnot(dim(B) == dim(prm$B))
  stopifnot(sum(is.na(B)) == 0)
  
  prm[["B"]] <- B
  prm[["nu"]] <- nu
  prm[["zeta"]] <- zeta
  return(prm)
}

core_B_TPBN <- function(Y, eta, nu, zeta, Sigma) {
  ete <- crossprod(eta, eta)
  Dinv <- diag(1/nu)
  U <- chol2inv(chol(ete + Dinv))
  V <- Sigma
  M <- crossprod(tcrossprod(eta, U), Y)
  B <- MBSP::matrix.normal(M = M, U = U, V = V)
  return(B)
}

update_zeta_nu_TPBN <- function(Y, B, nu, zeta, Sigmainv, u, a, tau) {
  t <- ncol(Y)
  n <- nrow(Y)
  q <- t
  Sinv <- Sigmainv
  for (i in sample(1:length(nu))) {
    # Option 1: B ~ TPN
    zeta[i] <- rgamma(1, a, nu[i] + tau)
    b <- B[i,]
    chi <- max(t(b)%*%Sinv%*%b, .Machine$double.eps) # to prevent chi parameter from collapsing to 0 from MBSP
    psi <- 2*zeta[i]
    lambda <- u-q*0.5
    nu[i] <- .Call(GIGrvg:::"C_rgig", n=1, lambda = lambda, chi = chi, psi = psi)
    #nu[i] <- GIGrvg::rgig(n = 1, lambda = lambda, chi = chi, psi = psi)
  }
  
  stopifnot(sum(is.na(zeta)) == 0)
  stopifnot(sum(is.na(nu)) == 0)
  
  return(list(nu=nu, zeta=zeta))
}