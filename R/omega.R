# Interaction tensor
get_eta_quad <- function(eta, id, K) {
  # Duplicate eta for repeated outcome measurements
  eta_dup <- eta[match(id, unique(id)), , drop=F]
  # Fill in pairwise interactions of factors
  # Omega needs to be a lower triangular tensor for O3 to match the following
  # pattern of regressor 11 12 ... 1K 22 ... 2K ... KK
  # eta_quad <- eta_dup * eta_dup[,1]
  # for (k in 2:K) {
  #   eta_quad <- cbind(eta_quad, eta_dup[,k:K] * eta_dup[,k])
  # }
  l <- lapply(1:K, function(k) eta_dup[,k:K] * eta_dup[,k])
  eta_quad <- do.call(cbind, l)
  return(eta_quad)
}

# Omega is an rTensor object, use as.tensor(array)
update_Omega <- function(prm, Y, X, K, Z, Z_int, O3, O3_is_col_na) {
  n <- nrow(Y)
  q <- ncol(Y)
  Omega_psir <- prm$Omega_psir
  Omega_zetar <- prm$Omega_zetar
  Omega_psic <- prm$Omega_psic
  Omega_zetac <- prm$Omega_zetac
  # Get residuals of the main effects
  Ytilde <- Y - prm$eta_int %*% prm$B - prm$alpha
  # Get pairwise interaction regressor
  eta_quad <- get_eta_quad(prm$eta, prm$id, K)
  # Diagonals of O3 column variance
  psi <- sapply(1:K, function(k) Omega_psir[k:K] * Omega_psic[k])
  psi <- unlist(psi) # make a vector
  psi <- sapply(psi, function(x) max(x, 1e-10)) # avoid numerical error
  # Matrix normal update
  Omega <- update_Omega_TPBN_cpp(prm$Omega, eta_quad, Ytilde, prm$Sigma, psi)
  # O3: 3-mode matricization of Omega tensor into a (q x K^2) matrix
  O3[, !O3_is_col_na] <- t(Omega)
  # Fold O3 into 3-tensor Omega
  Omega_tensor <- rTensor::k_fold(O3, m = 1, modes = c(q, K, K))
  
  # Update variances
  update_Omega_psi_cpp(Omega_psir, Omega_psic,
                       Omega_tensor@data, prm$Sigmainv, 
                       Omega_zetar, Omega_zetac, K, q)
  update_Omega_zeta_cpp(Omega_zetar, Omega_psir, K)
  update_Omega_zeta_cpp(Omega_zetac, Omega_psic, K)
  Omega_psi <- unlist(sapply(1:K, function(k) Omega_psir[k:K] * Omega_psic[k]))
  return(list(eta_quad = eta_quad,
              Omega = Omega,
              Omega_tensor = Omega_tensor,
              Omega_psir = Omega_psir,
              Omega_zetar = Omega_zetar,
              Omega_psic = Omega_psic,
              Omega_zetac = Omega_zetac,
              Omega_psi = Omega_psi))
}