famr <- function(niter, Y, X, Z=NULL, Z_int=NULL, K=2, id=NULL,
                  p_adapt=0.5, thin=1, epsilon_bar=1, H=0, gamma=0.05, t0=10, kappa=0.75, delta=0.65,
                  v1=NULL, v2=NULL, Ilod=NULL, loda=NULL, lodb=NULL,
                  B_u=1, B_a=0.5, B_tau=NULL,  # Strawderman-Berger prior
                  missing_Y=FALSE, verbose=FALSE, binary=NULL,
                  varimax=F, b0_global_shrinkage=1, 
                 init_eta_eps=1e-2, init_Theta_eps=1e-2,
                 eta=NULL, Theta=NULL, B=NULL, Sigma=NULL,
                 X_pred=NULL, Z_pred=NULL, Z_int_pred=NULL
                 ) {
  
  if (is.vector(Y)) Y <- matrix(Y, ncol = 1)
  if (nrow(Y) != nrow(X)) stop("Missmatched number of rows in Y and X")
  # id for potential repeated measurements/random intercepts
  if (is.null(id)) id <- 1:nrow(Y)
  X <- X[!duplicated(id), ] # reduce X to a unique subject per row (when there are repeated measurements of Y)

  q <- ncol(Y) # number of outcomes
  n <- length(unique(id)) # number of subjects
  p <- ncol(X) # number of components in mixtures
  p_z <- ncol(Z) # number of covariates
  p_int <- ncol(Z_int) # number of covariates interacting with the mixture latent factors
  K <- K # number of factors
  
  # Handle binary outcomes using probit model
  if (is.null(binary)) binary <- rep(0, ncol(Y)) 
  stopifnot(length(binary) == ncol(Y))
  Yraw <- Y
  if (sum(binary) > 0) {
    # Initialize probit latent variables to random normal draws
    Y[,binary==1] <- abs(rnorm(length(Y[,binary==1])))*(-1)^(Y[,binary==1]==0)
  }
  # Mean-centered X to remove intercepts
  Xmean <- colMeans(X, na.rm = TRUE)
  X <- X - tcrossprod(rep(1, nrow(X)), Xmean)
  #Y <- Y - tcrossprod(rep(1, nrow(Y)), colMeans(Y, na.rm = TRUE))
  # If there are missing values in Y
  O <- (!is.na(Y))*1 # indicators of observed Y
  if (missing_Y) { # initialize missing values to be zero
    Y[O == 0] <- 0
    Yraw[O == 0] <- 0
  }
  if (is.null(colnames(X))) colnames(X) <- paste0("x", 1:p)
  if (is.null(p_z)) p_z <- 0
  if (p_z > 0 & is.null(colnames(Z))) colnames(Z) <- paste0("z", 1:p_z)
  if (is.null(p_int)) p_int <- 0
  if (p_int > 0 & is.null(colnames(Z_int))) colnames(Z_int) <- paste0("zi", 1:p_int)
  if (is.null(v1)) v1 <- matrix(1, p, K) # TODO: remove???
  if (is.null(v2)) v2 <- matrix(1, p, K) # TODO: remove???
  # as recommended by Bai & Ghosh
  if (is.null(B_tau)) B_tau <- 1/(q*sqrt(n * log(n))) #TODO: remove??? 
  
  # Parallel Tempering temperatures
  # mc3_temps <- c(1)# seq(1.0, 0.2, -0.2)
  # Parallel Tempering acceptance
  # mc3_accept <- rep(NA, niter) 
  # Setting initial conditions
  prm <- list()
  ### For Regression
  prm[["id"]] <- id
  # K_int is the total number of predictors including 
  # the K factors, K*p_int interactions, and p_z covariates
  K_int <- K
  if (!is.null(p_z)) K_int <- K_int + p_z
  if (p_int > 0) K_int <- K_int + p_int*K
  prm[["alpha"]] <- rep(0, q) # outcome-specific intercept (needed for probit model)
  prm[["B"]] <- array(rnorm(q * K_int), dim = c(K_int, q)) # main effect & interactions matrix
  prm[["psi"]] <- rep(1/2*1/(q*sqrt(n * log(n))), K_int) # shrinkage on B
  prm[["zeta"]] <- 1/2*prm$psi # shrinkage on B
  prm[["is_zero"]] <- rep(1, K_int) # indicator whether the 95% credible region includes zeros vector
  prm[["hdr_thres"]] <- qchisq(0.95, df=q) # source: https://stats.stackexchange.com/questions/354063/computing-highest-density-region-given-multivariate-normal-distribution-with-dim
  # prm[["b0"]] <- rep(0, K_int)
  # prm[["u0sq"]] <- rep(1, K_int)
  # prm[["c0"]] <- rep(1, K_int)
  # prm[["tau0sq"]] <- b0_global_shrinkage
  # prm[["d0"]] <- 1
  prm[["Sigma"]] <- diag(rgamma(q, 2, 1), q, q) # covariance of y_i
  prm[["Sigmainv"]] <- diag(1/diag(prm$Sigma), q, q)
  prm[["l"]] <- apply(Y, 2, var) # hierarchical IW or prior on IS
  prm[["sigmay_sqinv"]] <- rgamma(1, 2.5*0.5, 2.5*0.084*0.5) # gives variance of Y a bit more flexibility
  ### For Factor model
  prm[["sigmax_sqinv"]] <- rgamma(p, 2.5*0.5, 2.5*0.084*0.5)    # noise variance
  prm[["eta"]]    <- array(rnorm(n * K), dim = c(n, K))         # random factors
  prm[["Theta"]]  <- array(rnorm(p * K, 0, 10), dim = c(p, K))  # loadings
  prm[["eta_n_accepted"]] <- 0L # for C++ update by reference
  prm[["eta_eps"]] <- init_eta_eps
  prm[["Theta_n_accepted"]] <- 0L # for C++ update by reference
  prm[["Theta_eps"]] <- init_Theta_eps
  # Dirichlet-Laplace prior on the rows of Theta
  prm[["tau"]] <- rep(0.5, p)
  prm[["phi"]] <- matrix(rep(1/K, p*K), p, K)
  prm[["omega"]] <- statmod::rinvgauss(p*K, (1/(2*K))/as.vector(prm$Theta), rep(1, p*K))
  # MGP Prior
  #prm[["phi"]]    <- array(rgamma(p * K, 1), dim = c(p, K))     # local shrinkage for Theta
  #prm[["delta"]]  <- rgamma(K, 1)                               # global shrinkage param for factor K

  # Storage for outputs
  nwarmup <- round(niter/2)
  sims <- list(B = array(NA, dim=c(niter - nwarmup, K_int, q)),
               Bx = array(NA, dim=c(niter - nwarmup, (p + p*p_int), q)),
               # b0 = array(NA, dim=c(niter - nwarmup, K_int)),
               # b0x = array(NA, dim=c(niter - nwarmup, p)), # main effect only
               alpha = array(NA, dim=c(niter - nwarmup, q)), # outcome-specific intercept
               Sigma = array(NA, dim = c(niter - nwarmup, q, q)),
               Theta = array(NA, dim = c(niter - nwarmup, p, K)),
               eta = array(NA, dim = c(niter - nwarmup, n, K)),
               sigmax_sqinv = array(NA, dim = c(niter - nwarmup, p)),
               Y_rep = array(NA, dim=c(niter - nwarmup, dim(Y))),
               is_zero = array(NA, dim=c(niter - nwarmup, K_int))
  )
  if (!is.null(X_pred)) {
    sims[["Y_pred"]] = array(NA, dim=c(niter - nwarmup, nrow(X_pred), q))
  }
  # TODO: remove varimax???
  if (varimax) {
    sims <- c(sims,
              list(B_varimax = array(NA, dim=c(niter - nwarmup, K_int, q)),
                   Theta_varimax = array(NA, dim = c(niter - nwarmup, p, K))
              ))
  }
  
  pb <- txtProgressBar(style = 3)
  # prm[["Theta"]] <- Theta # TODO FOR DEBUGGING
  # prm[["B"]] <- B # TODO FOR DEBUGGING
  # prm[["Sigma"]] <- Sigma # TODO FOR DEBUGGING
  #prm[["eta"]] <- eta # TODO FOR DEBUGGING
  for (s in 1:niter) {
    setTxtProgressBar(pb, s / niter)
    # "For each gradient update, all chains perform
    # one step of Gibbs sampling after which state swaps
    # are proposed between randomly selected neighbouring
    # chains."
    # " temperature is chosen such that swaps are accepted 
    # between 20 and 60% of the time"
    # TODO: add the temperature parameter
    #       TODO: can I only add it to the update step for eta??? or do i need all???
    # Gibbs steps
    prm <- update_sigmax_sqinv(prm, Y, X, K)
    # prm <- update_Theta_MGP(prm, Y, X, K, v1, v2)
    # prm <- update_Theta_DL(prm, Y, X, K) # TODO implement
    prm <- update_Theta_normal_mh(prm, Y, X, K) # TODO TEST!!!
    prm <- update_B(prm, Y, X, K, Z, Z_int, B_u, B_a, B_tau) # TODO TEST!!!
    prm <- update_Sigma_HIW(prm, Y, Z, Z_int, binary) # TODO implement IW
    prm <- update_eta_mh(prm, Y, X, Z, Z_int, n, K, p, q, p_z, p_int) # TODO TEST!!!
    
    # Parallel Tempering swapping step
    # swap_out <- swap_mc3(prm_list)
    # prm_list <- swap_out[["prm_list"]]
    # mc3_accpet[s] <- swap_out[["accept"]]
    
    # Imputation
    # if (!is.null(Ilod) & !is.null(loda) & !is.null(lodb)) {
    #   X <- impute_X(prm, X, Ilod, loda, lodb, Xmean)
    # }
    # if ((missing_Y & !is.null(O)) | sum(binary) > 0) {
    #   out_Ymis <- impute_Y(prm, X, K, Z, Z_int, Y, O, binary, Yraw, missing_Y)
    #   Y <- out_Ymis$Y
    #   Yraw <- out_Ymis$Yraw
    # }
    
    # store samples
    if (s > nwarmup) {
      sims[["sigmax_sqinv"]][s-nwarmup, ] <- prm$sigmax_sqinv
      sims[["Theta"]][s-nwarmup, , ] <- prm$Theta 
      sims[["eta"]][s-nwarmup, , ] <- prm$eta
      # sims[["b0"]][s-nwarmup, ] <- prm$b0
      sims[["alpha"]][s-nwarmup, ] <- prm$alpha
      sims[["is_zero"]][s-nwarmup, ] <- prm$is_zero
      
      # TODO: Remove varimax???
      # Varimax rotation of Theta
      if (varimax) {
        vari_out <- varimax(prm$Theta)
        vari_rot <- vari_out$rotmat
        sims[["Theta_varimax"]][s-nwarmup, , ] <- vari_out$loadings
      }
      
      # Rescaling for probit latent variables
      scale_B <- prm$B
      scale_S <- prm$Sigma
      if (sum(binary) > 0) { # TODO review this
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
      
      # TODO: Remove varimax???
      # NOTE: rotate after scaling. Rotation does not change between-outcome variance
      # Varimax rotation of B
      if (varimax) {
        temp_B_varimax <- array(0, dim=dim(scale_B))
        j <- 1
        while (j < (K + K*p_int)) {
          temp_B_varimax[j:(j+K-1),] <- t(vari_rot) %*% scale_B[j:(j+K-1),]
          j <- j + K
        }
        sims[["B_varimax"]][s-nwarmup, , ] <- temp_B_varimax
      }
      
      # Effects in original predictors
      temp_Bx <- matrix(0, nrow=(p + p*p_int), ncol=q)
      j <- 1
      jx <- 1
      Ainv <- solve( t(prm$Theta) %*% diag(prm$sigmax_sqinv) %*% prm$Theta + diag(1, K))
      Ax <- diag(prm$sigmax_sqinv) %*% prm$Theta %*% Ainv
      while (j < (K + K*p_int)) {
        temp_Bx[jx:(jx+p-1),] <- Ax %*% scale_B[j:(j+K-1),]
        j <- j + K
        jx <- jx + p
      }
      sims[["Bx"]][s-nwarmup, , ] <- temp_Bx
      #sims[["b0x"]][s-nwarmup, ] <- Ax %*% prm$b0[1:K]
      
      # Fitted values for Y
      eta_int <- get_eta_int(prm$eta, K, Z, Z_int, id)
      Y_rep <- eta_int %*% prm$B + MASS::mvrnorm(nrow(Y), rep(0,q), prm$Sigma)
      if (sum(binary) > 0) {
        Y_rep[, binary == 1] <- 1 * (Y_rep[, binary == 1] > 0)
      }
      sims[["Y_rep"]][s-nwarmup, , ] <- Y_rep
      
      # prediction for new data West 2003
      if (!is.null(X_pred)) {
        X_pred_int <- get_eta_int(X_pred, p, Z_pred, Z_int_pred, 1:nrow(X_pred))
        Y_pred <- X_pred_int[, 1:(K_int-p_z)] %*% temp_Bx
        if (p_z > 0) {
          Y_pred <- Y_pred + Z_pred %*% prm$B[(+1):K_int,]
        }
        sims[["Y_pred"]][s-nwarmup, , ] <- Y_pred
      }
    }
  }
  
  message(paste("\nMetropolis-Hasting acceptance prob. of eta is", 
                round(prm$eta_n_accepted/(niter*n), 3)))
  message(paste("\nMetropolis-Hasting acceptance prob. of Theta is", 
                round(prm$Theta_n_accepted/niter, 3)))
  return(sims)
}