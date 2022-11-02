# Transform eta for estimating the interaction matrix
# Coefficients are order as : latent main effect, interactions between latent factors and covariates, covariates 
get_eta_int <- function(eta, K, Z, Z_int, id) {
  uid <- unique(id)
  eta_dup <- eta[match(id, uid), , drop=F]
  if (is.null(colnames(eta_dup))) {
    colnames(eta_dup) <- paste0("e", 1:K)
  }
  # if (fa_interact) {
  #   eta_int <- model.matrix(~. + .^2 - 1, data = as.data.frame(eta_dup))
  eta_int <- eta_dup
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
get_k_I <- function(K_int, K, q_int) {
  k_I <- matrix(0, K, K_int)
  j <- K
  for (k in 1:K) {
    k_I[k,k] <- 1
    h <- k+1
    # while(h <= K & fa_interact) {
    #   j <- j+1
    #   k_I[k,j] <- 1
    #   k_I[h,j] <- 1
    #   h <- h+1
    # }
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

# TODO: pass u, a (v), tau (r) to RCpp funcntions
update_B <- function(prm, Y, X, K, Z, Z_int, u, a, tau) {
  q <- ncol(Y)
  n <- nrow(Y)
  eta_int <- get_eta_int(prm$eta, K, Z, Z_int, prm$id)
  p <- ncol(eta_int)
  B0 <- matrix(0, p, q)
  
  alpha <- update_intercept(Y - eta_int %*% prm$B, prm$Sigma)
  Ytilde <- Y - tcrossprod(rep(1, n), alpha)
  B <- update_B_TPBN(eta_int, Ytilde, prm$Sigma, B0, prm$psi, p, q)
  # psi <- update_psi_IG(B, B0, prm$Sigmainv, q, p)
  psi <- as.vector(update_psi(B, prm$Sigmainv, B0, prm$zeta, p, q))#rep(1, q)#
  zeta <- update_zeta(psi, p, global_shrink=1/(K*sqrt(n*log(n)))) #rep(1, q)#
  
  #b0 <- rep(0, p)#update_b0_Horseshoe(eta_int, Ytilda, prm$Sigmainv, psi, prm$u0sq, prm$tau0sq, q, p, n)
  #u0sq <- rep(1, p)#update_u0(b0, prm$c0, prm$tau0sq, q, p)
  #c0 <- 1#tryCatch(update_c0(u0sq, q),
                 #error=function(e) {print(u0sq)})
  #tau0sq <- update_tau0(c0, prm$d0, q, p)
  #d0 <- update_d0(tau0sq)
  #tau0sq <- update_tau0_hCauchy(b0, prm$u0sq, prm$d0, q, p)
  #d0 <- update_d0(tau0sq)
  #b0 <- update_b0(B, prm$Sigmainv, psi, prm$tau0sq, q, p)
  #tau0sq <- update_tau0_IG(b0, q, p)
  
  
  # out <- update_zeta_nu_TPBN(Y, prm$B, prm$psi, prm$zeta, prm$Sigmainv, u, a, tau)
  # psi <- out$nu
  # zeta <- out$zeta
  # 
  # B <- core_B_TPBN(Y, eta_int, psi, zeta, prm$Sigma)
  # rownames(B) <- colnames(eta_int)

  stopifnot(dim(B) == dim(prm$B))
  stopifnot(sum(is.na(B)) == 0)
  
  prm[["alpha"]] <- alpha
  prm[["B"]] <- B
  prm[["psi"]] <- psi
  prm[["zeta"]] <- zeta
  # prm[["b0"]] <- b0
  # prm[["u0sq"]] <- u0sq
  # prm[["c0"]] <- c0
  #prm[["tau0sq"]] <- tau0sq
  #prm[["d0"]] <- d0
  
  return(prm)
}

# Using a flat prior p(a) \prop 1
update_intercept <- function(Y, Sigma) {
  a <- MASS::mvrnorm(1, mu = colMeans(Y, na.rm = TRUE), Sigma = Sigma / nrow(Y))
  return(a)
}


# update_b0 <- function(B, Sigmainv, psi, tau0sq, q, p) {
#   b0 <- rep(0, q)
#   isi <- crossprod(rep(1, p), Sigmainv) %*% rep(1, p)
#   is <- crossprod(rep(1, p), Sigmainv)
#   for (i in range(q)) {
#     vi <- 1/(isi/psi[i] + 1/tau0sq)
#     mi <- vi * (is %*% B[i,]) /psi[i]
#     b0[i] <- rnorm(1, mi, sqrt(vi))
#   }
#   
#   return(b0)
# }

update_u0 <- function(b0, c0, tau0_sq, q, p) {
  u0_sq <- rep(0, q)
  for (i in range(q)) {
    chi <- b0[i]^2 / tau0_sq
    lam <- 0
    ps <- 2 * c0[i]
    u0_sq[i] <- .Call(GIGrvg:::"C_rgig", n=1, lambda = lam, chi = chi, psi = ps)
  }
  
  return(u0_sq)
}

update_psi_IG <- function(B, B0, Sigmainv, q, p, a=1, b=1) {
  psi <- rep(NA, q)
  for (i in 1:q) {
    SS <- t(B[i,]-B0[i,]) %*% Sigmainv %*% (B[i,]-B0[i,])
    psi[i] <- 1/rgamma(1, (p+a)/2, (SS + b)/2)
  }
  
  return(psi)
}

# update_tau0 <- function(c0, d0, q, p) {
#   tau0sq <- rgamma(1, 0.5 + 0.5*q, sum(c0, d0))
#   
#   return(tau0sq)
# }
# 
# update_d0 <- function(tau0sq) {
#   d0 <- rgamma(1, 1, 1 + tau0sq)
#   
#   return(d0)
# }

# update_tau0_IG <- function(b0, q, p) {
#   tau0_sq <- 1/rgamma(1, 0.5*(q+1), 0.5*sum(b0^2, 1))
#   
#   return(tau0_sq)
# }

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

# TODO: DEBUG
# Spike and Slab priors, with continuous mixtures or spike-n-slab hierarchical priors
# TODO: Edit notations
# Adding in Spike and Slab priors for each row of B
# TODO: Add SNS prior on each entry in a row, mixed with Normal(0, SigmaY) slab
update_B_HSNS <- function(prm, Y, eta_int, prior="laplace") {
  Y    <- Y
  etay <- eta_int 
  q    <- ncol(Y)
  B    <- prm$B
  psi   <- prm$psi
  zeta <- prm$zeta
  gamma <- prm$gamma
  pi_gamma <- prm$pi_gamma
  Sinv  <- prm$Sigmainv
  mg    <- 1 # TODO: changes if select groups of variables
  
  for (i in sample(1:nrow(B))) {
    ### Probability of gamma_j = 0 ###
    Xg <- etay[,i]
    Sginv <- as.matrix(1/psi[i] + t(Xg)%*%Xg)
    Sg <- 1/Sginv
    Mg <- Sg %*% t(Xg) %*% (Y - etay[,-i] %*% B[-i,])
    lRg <- log(psi[i])*(-q/2) + log(Sg)*(q/2) + 0.5*sum(diag( Sinv %*% t(Mg) %*% Sginv %*% Mg))
    p0 <- pi_gamma/(pi_gamma + (1-pi_gamma)*exp(lRg))

    ### Update B[i,], psi[i], zeta[i] given gamma[i] ###
    if(runif(1) < p0) {
      
      gamma[i] <- 0
      B[i,] <- 0
      if (prior == "laplace") { 
        # B[i,] ~ MLaplace, psi ~ G((q+1) / 2, zeta / 2)
        psi[i] <- rgamma(1, 0.5*(q+1), 0.5*zeta[i])  # MBSGS prior, MLaplace on Bg
        zeta[i] <- rgamma(1, 0.5, 1)    
      } else if (prior == "tpb") {
        # B[i,] ~ TPB, psi ~ G(u, zeta)
        psi[i] <- rgamma(1, 0.5, zeta[i])
        zeta[i] <- rgamma(1, 0.5, 1)
      } else stop(paste("Prior", prior, "not implemented"))
      
    } else {
      
      gamma[i] <- 1
      B[i,] <- MASS::mvrnorm(1, mu = Mg, Sigma = Sg[1,1]*prm$Sigma)
      d <- t(B[i,])%*%Sinv%*%B[i,]
      if (prior == "laplace") {
        psi[i] <- 1/statmod::rinvgauss(1, mean = sqrt(zeta[i]/d), shape = zeta[i])      # MBSGS 
        zeta[i] <- rgamma(1, 0.5, psi[i]*0.5 + 0.5) # zeta ~ G(a, tau) prior
      } else if (prior == "tpb") {
        chi <- max(d, .Machine$double.eps) # to prevent chi parameter from collapsing to 0 from MBSP
        psi <- 2*zeta[i]
        lambda <- 0.5-q*0.5
        psi[i] <- .Call(GIGrvg:::"C_rgig", n=1, lambda = lambda, chi = chi, psi = psi)
        zeta[i] <- rgamma(1, 0.5, psi[i] + 1)
      }
    }
    
  }
  
  ### Update P(gamma_j=0): pi_gamma ###
  pi_gamma <- rbeta(1, 1 + sum(gamma==0), 1 + sum(gamma))
  
  prm[["gamma"]] <- gamma
  prm[["pi_gamma"]] <- pi_gamma
  prm[["B"]]     <- B
  prm[["psi"]]    <- psi
  prm[["zeta"]]  <- zeta
  
  return(prm)
}