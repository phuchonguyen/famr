famr <- function(niter, Y, X, Z=NULL, Z_int=NULL, K=2, id=NULL, 
                 nthin=1, nwarmup=NULL,
                 Ilod=NULL, loda=NULL, lodb=NULL,
                  missing_Y=FALSE, verbose=FALSE, binary=NULL,
                 init_eta_eps=1e-2, 
                 init_Theta_eps=1e-2, adaptiveM=F, adaptiveMWG=F,
                 init_a1_eps=1, init_a2_eps=1,
                 s0=0.084,
                 eta=NULL, Theta=NULL, B=NULL, Sigma=NULL, sigmax_sqinv=NULL,
                 X_pred=NULL, Z_pred=NULL, Z_int_pred=NULL,
                 include_rep=FALSE, include_interactions=FALSE,
                 random_intercept=FALSE
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
  # as recommended by Bai & Ghosh
  # if (is.null(B_tau)) B_tau <- 1/(q*sqrt(n * log(n))) #TODO: remove??? 
  
  # Parallel Tempering temperatures
  # mc3_temps <- c(1)# seq(1.0, 0.2, -0.2)
  # Parallel Tempering acceptance
  # mc3_accept <- rep(NA, niter) 
  # Setting initial conditions
  prm <- list()
  ### For Regression
  prm[["id"]] <- id
  prm[["uid"]] <- unique(id)
  # K_int is the total number of predictors including 
  # the K factors, K*p_int interactions, and p_z covariates
  K_int <- K
  if (!is.null(p_z)) K_int <- K_int + p_z
  if (p_int > 0) K_int <- K_int + p_int*K
  # ********************************************
  # outcome-specific intercept (needed for probit model)
  # ********************************************
  # prm[["alpha_mu"]] <- rnorm(q, 0, 100) # For random intercept
  # prm[["alpha_v"]] <- rgamma(q, 1, 1)
  prm[["alpha"]] <- MASS::mvrnorm(nrow(Y),
                                  mu=rep(0, q), 
                                  Sigma=diag(100, q, q))
  prm[["B"]] <- array(rnorm(q * K_int, 0, 100), dim = c(K_int, q)) # main effect & interactions matrix
  # ********************************************
  # TPBN prior on B
  # ********************************************
  prm[["psi"]] <- rep(1/2*1/(q*sqrt(n * log(n))), K_int) # shrinkage on B
  prm[["zeta"]] <- 1/2*prm$psi # shrinkage on B
  # ********************************************
  # DL prior on B
  # ********************************************
  # prm[["psi"]] <- matrix(1, K_int, q)
  # prm[["zeta"]] <- matrix(1, K_int, q)
  # prm[["nu"]] <- rep(1, K_int)
  
  # Pairwise interactions between the factors
  if (include_interactions) {
    O_tmp <- array(rnorm(q*K*K, 0, 100), dim = c(q, K, K))
  } else {
    O_tmp <- array(0, dim = c(q, K, K))
  }
  for (j in 1:q) {
    O_tmp[j,,][upper.tri(O_tmp[j,,])] <- NA
  }
  O3 <- rTensor::k_unfold(rTensor::as.tensor(O_tmp), m = 1)
  O3_is_col_na <- (colSums(is.na(O3@data)) == q)
  prm[["Omega_tensor"]] <- rTensor::k_fold(O3, m = 1, modes = c(q, K, K)) # upper triangular part has NAs
  prm[["Omega"]] <- t(O3@data[ ,!O3_is_col_na, drop=F])
  prm[["Omega_psir"]] <- rep(1/2*1/(q*sqrt(n * log(n))), K) # shrinkage on Omega
  prm[["Omega_psic"]] <- rep(1/2*1/(q*sqrt(n * log(n))), K) # shrinkage on Omega
  prm[["Omega_psi"]] <- unlist(sapply(1:K, function(k) prm$Omega_psir[k:K] * prm$Omega_psic[k]))
  prm[["Omega_zetar"]] <- 1/2*prm$Omega_psir
  prm[["Omega_zetac"]] <- 1/2*prm$Omega_psic
  prm[["eta_quad"]] <- matrix(1, nrow(Y), nrow(prm$Omega))
  
  prm[["Sigma"]] <- diag(rgamma(q, 2, 1), q, q) # covariance of y_i
  prm[["Sigmainv"]] <- diag(1/diag(prm$Sigma), q, q)
  prm[["l"]] <- apply(Y, 2, var) # hierarchical IW or prior on IS
  ### For Factor model
  # Vxinv_emp <- solve(cov(X) + diag(1e-10, p, p)) # estimated variance of X
  prm[["sigmax_sqinv"]] <- rgamma(p, 2.5*0.5, 2.5*0.084*0.5)    # noise variance
  prm[["eta"]]    <- array(rnorm(n * K), dim = c(n, K))         # random factors
  prm[["Theta"]]  <- array(rnorm(p * K, 0, 10), dim = c(p, K))  # loadings
  prm[["eta_n_accepted"]] <- rep(0, n) # for C++ update by reference
  prm[["eta_eps"]] <- rep(init_eta_eps, n)
  prm[["eta_A"]] <- array(0, dim = c(K, K, n))
  prm[["eta_b"]] <- array(0, dim = c(K, n))
  prm[["Theta_n_accepted"]] <- 0L # for C++ update by reference
  prm[["Theta_eps"]] <- init_Theta_eps
  # ********************************************
  # Dirichlet-Laplace prior on the columns of Theta
  # ********************************************
  # prm[["tau"]] <- rep(0.5, p)
  # prm[["phi"]] <- matrix(rep(1/p, p*K), p, K)
  # prm[["omega"]] <- matrix(statmod::rinvgauss(p*K, abs((1/(2*p))/as.vector(prm$Theta)), rep(1, p*K)), p, K)
  # ********************************************
  # MGP Prior on Theta
  # ********************************************
  prm[["phi"]]    <- array(rgamma(p * K, 1), dim = c(p, K)) # local shrinkage for Theta
  prm[["delta"]]  <- rgamma(K, 1)                           # global shrinkage param for factor K
  prm[["a1"]] <- 2.5
  prm[["a2"]] <- 3.5
  prm[["a1_n_accepted"]] <- 0L
  prm[["a2_n_accepted"]] <- 0L
  prm[["a1_eps"]] <- init_a1_eps
  prm[["a2_eps"]] <- init_a2_eps
  
  # Storage for outputs
  if (is.null(nwarmup)) {
    nwarmup <- round(niter/2) 
  }
  sims <- list(B = array(NA, dim=c(niter - nwarmup, K_int, q)),
               Bx = array(NA, dim=c(niter - nwarmup, (p + p*p_int), q)),
               # b0 = array(NA, dim=c(niter - nwarmup, K_int)),
               # b0x = array(NA, dim=c(niter - nwarmup, p)), # main effect only
               # alpha = array(NA, dim=c(niter - nwarmup, q)), # outcome-specific intercept
               # Sigma = array(NA, dim = c(niter - nwarmup, q, q)),
               #Theta = array(NA, dim = c(niter - nwarmup, p, K)),
               #eta = array(NA, dim = c(niter - nwarmup, n, K)),
               sigmax_sqinv = array(NA, dim = c(niter - nwarmup, p)),
               NULL
  )
  if (include_rep) {
    sims[["Y_rep"]] = array(NA, dim=c(niter - nwarmup, dim(Y)))
  }
  if (!is.null(X_pred)) {
    sims[["Y_pred"]] = array(NA, dim=c(niter - nwarmup, nrow(X_pred), q))
  }
  if (include_interactions) {
    sims[["Omega"]] = array(NA, dim=c(niter - nwarmup, dim(prm$Omega)))
    sims[["Omega_array_x"]] = array(NA, dim=c(niter - nwarmup, q, p, p))
  }
  
  pb <- txtProgressBar(style = 3)
  # prm[["Theta"]] <- Theta # TODO FOR DEBUGGING
  # prm[["B"]] <- B # TODO FOR DEBUGGING
  # prm[["Sigma"]] <- Sigma # TODO FOR DEBUGGING
  # prm[["Sigmainv"]] <- solve(Sigma)
  # prm[["eta"]] <- eta # TODO FOR DEBUGGING
  # prm[["sigmax_sqinv"]] <- sigmax_sqinv
  n_till_adaptive <- min(3000, round(niter/10)) 
  for (s in 1:niter) {
    setTxtProgressBar(pb, s / niter)
    # MCMC steps
    prm[["sigmax_sqinv"]] <- update_sigmax_sqinv(prm, Y, X, K, s0)
    ## Theta MGP ---------------------------------------------------------------
    tmp <- update_Theta_MGP(prm, Y, X, K)
    prm[names(tmp)] <- tmp
    ## -------------------------------------------------------------------------
    # prm <- update_Theta_DL_mh(prm, Y, X, K) # Not as good as MGP
    # prm <- update_Theta_normal(prm, Y, X, K) # Not as good as MGP
    ## eta using Metroplis-within-Gibbs ----------------------------------------
    prm <- update_eta_mh(prm, Y, X, Z, Z_int, n, K, p, q, p_z, p_int,
                         s, adaptiveM=((s>n_till_adaptive) & adaptiveM),
                         adaptiveMWG=adaptiveMWG)
    ## -------------------------------------------------------------------------
    # prm <- update_eta_gibbs(prm, Y, X, Z, Z_int, n, K, p, q, p_z, p_int)
    ## Outcome regression ---------------------------------_--------------------
    tmp <- update_B_TPBN(prm, Y, X, K, Z, Z_int, random_intercept)
    prm[names(tmp)] <- tmp
    if (include_interactions) {
      tmp <- update_Omega(prm, Y, X, K, Z, Z_int, O3, O3_is_col_na)
      prm[names(tmp)] <- tmp
    }
    tmp <- update_Sigma_IW_TPBN(prm, Y, Z, Z_int, binary)
    prm[names(tmp)] <- tmp
    ## -------------------------------------------------------------------------
    # prm <- update_B_DL(prm, Y, X, K, Z, Z_int) # NOT GOOD AS TPBN
    # prm <- update_Sigma_IW_DL(prm, Y, Z, Z_int, binary) # NOT GOOD AS TPBN
    
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
      # sims[["Theta"]][s-nwarmup, , ] <- prm$Theta
      # sims[["eta"]][s-nwarmup, , ] <- prm$eta
      # sims[["b0"]][s-nwarmup, ] <- prm$b0
      # sims[["alpha"]][s-nwarmup, ] <- prm$alpha
      if (include_interactions) {
        sims[["Omega"]][s-nwarmup,,] <- prm$Omega
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
      # sims[["Sigma"]][s-nwarmup, , ] <- scale_S
      
      # Effects in original predictors
      temp_Bx <- matrix(0, ncol=(p + p*p_int), nrow=q)
      j <- 0
      jx <- 0
      Ainv <- tryCatch(solve(t(prm$Theta) %*% diag(prm$sigmax_sqinv, p, p) %*% prm$Theta + diag(1, K, K)),
                       error = function(e) {
                         print(prm$Theta)
                         print(prm$sigmax_sqinv)
                         })
      Ax <- Ainv %*% t(prm$Theta) %*% diag(prm$sigmax_sqinv, p, p)
      while (j < (K + K*p_int)) {
        temp_Bx[, (jx+1):(jx+p)] <- crossprod(scale_B[(j+1):(j+K),,drop=F], Ax)
        j <- j + K
        jx <- jx + p
      }
      temp_Bx <- t(temp_Bx)
      sims[["Bx"]][s-nwarmup, , ] <- temp_Bx
      if (include_interactions) {
        # zero out NAs
        prm[["Omega_tensor"]]@data[is.na(prm$Omega_tensor@data)] <- 0
        # k-mode tensor multiplication to transform this into coefs in X
        sims[["Omega_array_x"]][s-nwarmup,,,] <- rTensor::ttl(prm$Omega_tensor,
                                                               list(t(Ax), t(Ax)),
                                                               c(2, 3))@data
      }
      #sims[["b0x"]][s-nwarmup, ] <- Ax %*% prm$b0[1:K]
      
      # Fitted values for Y
      if (include_rep) {
        eta_int <- get_eta_int(prm$eta, K, Z, Z_int, id)
        Y_rep <- eta_int %*% prm$B + MASS::mvrnorm(nrow(Y), rep(0,q), prm$Sigma)
        if (sum(binary) > 0) {
          Y_rep[, binary == 1] <- 1 * (Y_rep[, binary == 1] > 0)
        }
        sims[["Y_rep"]][s-nwarmup, , ] <- Y_rep
      }
      
      # prediction for new data West 2003
      if (!is.null(X_pred)) {
        X_pred_int <- get_eta_int(X_pred, p, Z_pred, Z_int_pred, 1:nrow(X_pred))
        tmp_p <- ncol(X_pred_int)
        Y_pred <- X_pred_int[, 1:(tmp_p-p_z)] %*% temp_Bx
        if (p_z > 0) {
          Y_pred <- Y_pred + Z_pred %*% prm$B[(K_int-p_z+1):K_int,]
        }
        sims[["Y_pred"]][s-nwarmup, , ] <- Y_pred
      }
    }
  }
  
  message(paste("Metropolis-Hasting acceptance prob. of eta is", 
                round(sum(prm$eta_n_accepted)/(niter*n), 3)))
  # message(paste("\nMetropolis-Hasting acceptance prob. of Theta is", 
  #               round(prm$Theta_n_accepted/(niter*K), 3)))
  message(paste("\nMetropolis-Hasting acceptance prob. of a1 is",
                round(prm$a1_n_accepted/niter, 3)))
  message(paste("\nMetropolis-Hasting acceptance prob. of a2 is",
                round(prm$a2_n_accepted/niter, 3)))
  return(sims)
}