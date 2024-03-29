#' Impute values under the limit of detection (LOD) in mixtures
#' @param Ilod is (n x p) indicator matrix of observations below lod
#' @param loda is a vector of p lower limit LOD for variables in X
#' @param lodb is a vector of p upper limit LOD for variables in X
impute_X <- function(prm, X, Ilod, loda, lodb, Xmean) {
  X <- impute_X_lod_cpp(prm$eta, prm$Theta, X, prm$sigmax_sqinv, nrow(X), 
                        Ilod, loda-Xmean, lodb-Xmean)
  
  return(X)
}

impute_Y <- function(prm, X, K, Z, Z_int, Y, O, binary, Yraw, missing_Y) {
  Sigma <- prm$Sigma
  eta_int <- get_eta_int(prm$eta, K, Z, Z_int, prm$fa_interact, prm$id)
  M <- eta_int %*% prm$B
  if (missing_Y) {
    O_cont <- O
    O_cont[, binary == 1] <- 1  # don't impute latent Y
    Y <- impute_Ymis_cpp(Y, M, Sigma, O_cont, nrow(Y), ncol(Y))
    Yraw[, binary == 0] <- Y[, binary == 0]
  }
  
  if (sum(binary) > 0) {
    Y <- impute_Yprobit_cpp(Y, M, Sigma, Yraw, binary, nrow(Y), ncol(Y))
    Yraw[, binary == 1][O[, binary == 1] == 0] <- 1 * (Y[, binary == 1][O[, binary == 1] == 0] > 0)
  }
  
  return(list(Y=Y, Yraw=Yraw))
}