famr <- function(niter, Y, X, Z=NULL, Z_int=NULL, K=2, id=NULL,
                  p_adapt=0.5, epsilon_bar=1, H=0, gamma=0.05, t0=10, kappa=0.75, delta=0.65,
                  v1=NULL, v2=NULL, Ilod=NULL, loda=NULL, lodb=NULL,
                  B_u=1, B_a=0.5, B_tau=NULL,  # Strawderman-Berger prior
                  missing_Y=FALSE, verbose=FALSE, binary=NULL,
                  varimax=F) {
  fa_interact=FALSE #TODO: Remove this feature
  
  if (is.vector(Y)) Y <- matrix(Y, ncol = 1)
  if (nrow(Y) != nrow(X)) stop("Missmatched number of rows in Y and X")
  # id for potential repeated measurements/random intercepts
  if (is.null(id)) id <- 1:nrow(Y)
  X <- X[!duplicated(id), ]   # reduce X to a unique subject per row
  if (is.null(binary)) binary <- rep(0, ncol(Y)) 
  stopifnot(length(binary) == ncol(Y))
  Yraw <- Y
  
  t <- ncol(Y)  # number of outcome times
  n <- length(unique(id)) # number of subjects
  p <- ncol(X)  # number of components in mixtures
  q <- ncol(Z)
  q_int <- ncol(Z_int)
  K <- K  # number of factors
  # Mean-centered X, Y to remove intercepts
  Xmean <- colMeans(X, na.rm = TRUE)
  X <- X - tcrossprod(rep(1, nrow(X)), Xmean)
  Y <- Y - tcrossprod(rep(1, nrow(Y)), colMeans(Y, na.rm = TRUE))
  O <- (!is.na(Y))*1  # indicators of missing values in Y
  if (missing_Y) {
    Y[O == 0] <- 0
  }
  if (is.null(colnames(X))) colnames(X) <- paste0("x", 1:p)
  if (is.null(q)) q <- 0
  if(q > 0 & is.null(colnames(Z))) colnames(Z) <- paste0("z", 1:q)
  if (is.null(q_int)) q_int <- 0
  if (q_int > 0 & is.null(colnames(Z_int))) colnames(Z_int) <- paste0("zi", 1:q_int)
  if (is.null(v1)) v1 <- matrix(1, p, K)
  if (is.null(v2)) v2 <- matrix(1, p, K)
  if (is.null(B_tau)) B_tau <- 1/(t*sqrt(n * log(n))) # as recommended by Bai & Ghosh
  
  # Setting initial conditions
  prm <- list()   
  ### For Regression
  prm[["id"]] <- id
  K_int <- K
  prm[["fa_interact"]] <- fa_interact
  if (!is.null(q)) K_int <- K_int + q
  if (q_int > 0) K_int <- K_int + q_int*K
  prm[["B"]] <- array(rnorm(t * K_int), dim = c(K_int, t))      # main effect & interactions matrix
  prm[["nu"]] <- rep(1/2*1/(t*sqrt(n * log(n))), K_int)          # shrinkage on B
  prm[["zeta"]] <- 1/2*prm$nu                                      # shrinkage on B
  prm[["Sigma"]] <- diag(rgamma(t, 2, 1))           # covariance of y_i
  prm[["Sigmainv"]] <- diag(1/diag(prm$Sigma))
  prm[["l"]] <- rep(1, t)            # hierarchical IW
  ### For Factor model
  prm[["sigmax_sqinv"]] <- rgamma(p, 2.5*0.5, 2.5*0.084*0.5)    # noise variance
  prm[["Theta"]]  <- array(rnorm(p * K, 0, 10), dim = c(p, K))  # loadings
  # MGP Prior
  prm[["phi"]]    <- array(rgamma(p * K, 1), dim = c(p, K))     # local shrinkage for Theta
  prm[["delta"]]  <- rgamma(K, 1)                               # global shrinkage param for factor K
  prm[["eta"]]    <- array(rnorm(n * K), dim = c(n, K))         # random factors
  
  # Storage for outputs
  nwarmup <- round(niter/2)
  sims <- list(B = array(NA, dim=c(niter - nwarmup, K_int, t)),
               Bx = array(NA, dim=c(niter - nwarmup, (p + p*q_int), t)),
               Sigma = array(NA, dim = c(niter - nwarmup, t, t)),
               Theta = array(NA, dim = c(niter - nwarmup, p, K)),
               eta = array(NA, dim = c(niter - nwarmup, n, K)),
               sigmax_sqinv = array(NA, dim = c(niter - nwarmup, p)),
               Y_rep = array(NA, dim=c(niter - nwarmup, dim(Y))) 
  )
  if (varimax) {
    sims <- c(sims,
              list(B_varimax = array(NA, dim=c(niter - nwarmup, K_int, t)),
                   Theta_varimax = array(NA, dim = c(niter - nwarmup, p, K))
              ))
  }
  
  pb <- txtProgressBar(style = 3)
  for (s in 1:niter) {
    setTxtProgressBar(pb, s / niter)
    # Gibbs steps
    prm <- update_sigmax_sqinv(prm, Y, X, K)
    prm <- update_Theta_MGP(prm, Y, X, K, v1, v2)
    prm <- update_B_TPBN(prm, Y, X, K, Z, Z_int, B_u, B_a, B_tau)
    prm <- update_Sigma_HIW(prm, Y, Z, Z_int, binary)
    prm <- update_eta_gibbs(prm, Y, X, Z, Z_int, n, K, p, t, q, q_int)
    
    # Imputation
    if (!is.null(Ilod) & !is.null(loda) & !is.null(lodb)) {
      X <- impute_X(prm, X, Ilod, loda, lodb, Xmean)
    }
    if ((missing_Y & !is.null(O)) | sum(binary) > 0) {
      out_Ymis <- impute_Y(prm, X, K, Z, Z_int, Y, O, binary, Yraw, missing_Y)
      Y <- out_Ymis$Y
      Yraw <- out_Ymis$Yraw
    }
    
    # store samples
    if (s > nwarmup) {
      sims[["sigmax_sqinv"]][s-nwarmup, ] <- prm$sigmax_sqinv
      sims[["Theta"]][s-nwarmup, , ] <- prm$Theta 
      sims[["eta"]][s-nwarmup, , ] <- prm$eta
      
      # Varimax rotation of Theta
      if (varimax) {
        vari_out <- varimax(prm$Theta)
        vari_rot <- vari_out$rotmat
        sims[["Theta_varimax"]][s-nwarmup, , ] <- vari_out$loadings
      }
      
      # Rescaling for probit latent variables
      scale_B <- prm$B
      scale_S <- prm$Sigma
      if (sum(binary) > 0) {
        # Scale Sigma
        d <- sqrt(diag(prm$Sigma))
        temp_s <- d[binary == 1]
        d[binary == 1] <- 1
        R <- cov2cor(prm$Sigma)
        scale_S <- diag(d) %*% R %*% diag(d)
        # Scale B
        scale_B[, binary == 1] <- sweep(prm$B[, binary == 1, drop=FALSE], 2, temp_s, '/')
      }
      sims[["B"]][s-nwarmup, , ] <- scale_B
      sims[["Sigma"]][s-nwarmup, , ] <- scale_S
      
      # NOTE: rotate after scaling. Rotation does not change between-outcome variance
      # Varimax rotation of B
      temp_B_varimax <- array(0, dim=dim(scale_B))
      if (varimax) {
        j <- 1
        while (j < (K + K*q_int)) {
          temp_B_varimax[j:(j+K-1),] <- t(vari_rot) %*% scale_B[j:(j+K-1),]
          j <- j + K
        }
      }
      sims[["B_varimax"]][s-nwarmup, , ] <- temp_B_varimax
      
      # Effects in original predictors
      temp_Bx <- matrix(0, nrow=(p + p*q_int), ncol=t)
      j <- 1
      jx <- 1
      Ainv <- solve( t(prm$Theta) %*% diag(prm$sigmax_sqinv) %*% prm$Theta + diag(1, K))
      Ax <- diag(prm$sigmax_sqinv) %*% prm$Theta %*% Ainv
      while (j < (K + K*q_int)) {
        temp_Bx[jx:(jx+p-1),] <- Ax %*% scale_B[j:(j+K-1),]
        j <- j + K
        jx <- jx + p
      }
      sims[["Bx"]][s-nwarmup, , ] <- temp_Bx
      
      # Fitted values for Y
      eta_int <- get_eta_int(prm$eta, K, Z, Z_int, fa_interact, id)
      Y_rep <- eta_int %*% prm$B
      if (sum(binary) > 0) {
        Y_rep[, binary == 1] <- 1 * (Y_rep[, binary == 1] > 0)
      }
      sims[["Y_rep"]][s-nwarmup, , ] <- Y_rep
    }
  }
  
  return(sims)
}