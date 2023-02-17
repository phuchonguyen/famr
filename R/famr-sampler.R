#' @id vector of length N, unique id of subjects, might contain duplicated id
#' 
famr <- function(niter, Y, X, K=2, Z=NULL, Z_int=NULL, id=NULL, 
                 nthin=1, nwarmup=NULL,
                 adaptiveM=F, adaptiveMWG=F,
                 X_pred=NULL, Z_pred=NULL, Z_int_pred=NULL,
                 include_rep=FALSE, include_interactions=FALSE,
                 init_eta_eps=1e-2, init_a1_eps=1, init_a2_eps=1,
                 Ilod=NULL, loda=NULL, lodb=NULL,
                 s0=0.084, r=2.5, verbose=FALSE, missing_Y=FALSE, binary=NULL,
                 eta=NULL, Theta=NULL, B=NULL, Sigma=NULL, sigmax_sqinv=NULL
                 ) {
  
  if (is.vector(Y)) Y <- matrix(Y, ncol = 1)
  if (nrow(Y) != nrow(X)) stop("Missmatched number of rows in Y and X")
  # when id=NULL, no repeated measurements of Y
  if (is.null(id)) id <- 1:nrow(Y)
  # reduce X to a unique subject per row (when there are repeated measurements of Y)
  # each subjet has ONE UNIQUE measurement of X
  X <- X[!duplicated(id), ]
  q <- ncol(Y) # number of outcomes
  n <- length(unique(id)) # number of subjects, each subject might have T_i measurements of Y
  p <- ncol(X) # number of components in mixtures
  p_z <- ncol(Z) # number of covariates
  p_int <- ncol(Z_int) # number of covariates interacting with the mixture latent factors

  # Mean-centered X to remove intercepts
  X <- scale(X, center = T, scale = F)
  # If there are missing values in Y
  O <- (!is.na(Y))*1 # indicators of observed Y
  if (missing_Y) { # initialize missing values to be zero
    Y[O == 0] <- 0
    Yraw <- Y
  }
  if (is.null(colnames(X))) colnames(X) <- paste0("x", 1:p)
  if (is.null(p_z)) p_z <- 0
  if (p_z > 0 & is.null(colnames(Z))) colnames(Z) <- paste0("z", 1:p_z)
  if (is.null(p_int)) p_int <- 0
  if (p_int > 0 & is.null(colnames(Z_int))) colnames(Z_int) <- paste0("zi", 1:p_int)
  
  # Setting initial conditions
  prm <- list()
  prm[["id"]] <- id
  prm[["uid"]] <- unique(id)
  # K_int is the total number of predictors including 
  # the K factors, K*p_int interactions, and p_z covariates
  K_int <- K + p_z + p_int*K
  # ********************************************
  # outcome-specific intercept (needed for probit model)
  # ********************************************
  prm[["alpha"]] <- rnorm(q)
  # ********************************************
  # TPBN prior on main effects of factors/covariates & interactions with covariates
  # ********************************************
  prm[["B"]] <- array(rnorm(q * K_int, 0, 100), dim = c(K_int, q))
  prm[["psi"]] <- rep(1/2*1/(q*sqrt(n * log(n))), K_int) # shrinkage on rows of B
  prm[["zeta"]] <- 1/2*prm$psi # shrinkage on rows of B
  prm[["Sigma"]] <- diag(rgamma(q, 2, 1), q, q) # covariance of y_i
  prm[["Sigmainv"]] <- diag(1/diag(prm$Sigma), q, q)
  # ********************************************
  # Factor model
  # ********************************************
  prm[["sigmax_sqinv"]] <- rgamma(p, 2.5*0.5, 2.5*0.084*0.5) # noise variance
  prm[["eta"]] <- array(rnorm(n * K), dim = c(n, K)) # random factors
  prm[["eta_n_accepted"]] <- rep(0, n) # MH acceptance
  prm[["eta_eps"]] <- rep(init_eta_eps, n) # MH proposal scales
  prm[["eta_A"]] <- array(0, dim = c(K, K, n)) # adaptive MH 
  prm[["eta_b"]] <- array(0, dim = c(K, n)) # adaptive MH
  prm[["eta_int"]] <- get_eta_int(prm$eta, K, Z, Z_int, id)
  # ********************************************
  # Multiplicative Gamma Process prior on the columns of loadings
  # ********************************************
  prm[["Theta"]] <- array(rnorm(p * K, 0, 10), dim = c(p, K)) # loadings
  prm[["phi"]] <- array(rgamma(p * K, 1), dim = c(p, K)) # local shrinkage for Theta
  prm[["delta"]] <- rgamma(K, 1) # global shrinkage param for factor K
  prm[["a1"]] <- 2.5
  prm[["a2"]] <- 3.5
  prm[["a1_n_accepted"]] <- 0L
  prm[["a2_n_accepted"]] <- 0L
  prm[["a1_eps"]] <- init_a1_eps
  prm[["a2_eps"]] <- init_a2_eps
  # ********************************************
  # TPBN prior on each row/column of the pairwise interactions between the factors
  # stored in a (q K K) tensor
  # ********************************************
  if (include_interactions) {
    O_tmp <- array(rnorm(q*K*K, 0, 100), dim = c(q, K, K))
  } else {
    O_tmp <- array(0, dim = c(q, K, K))
  }
  # Tensor made up of lower-diagonal matrices because the matrices are 
  # symmetric otherwise
  for (j in 1:q) {
    O_tmp[j,,][upper.tri(O_tmp[j,,])] <- NA
  }
  # mode-1 matricization of interaction tensor
  O3 <- rTensor::k_unfold(rTensor::as.tensor(O_tmp), m = 1)
  # convert to rTensor object for tensor operations
  prm[["Omega_tensor"]] <- rTensor::k_fold(O3, m = 1, modes = c(q, K, K)) # upper triangular part has NAs
  prm[["Omega_psir"]] <- rep(1/2*1/(q*sqrt(n * log(n))), K) # shrinkage on Omega
  prm[["Omega_psic"]] <- rep(1/2*1/(q*sqrt(n * log(n))), K) # shrinkage on Omega
  prm[["Omega_psi"]] <- unlist(sapply(1:K, function(k) prm$Omega_psir[k:K] * prm$Omega_psic[k]))
  prm[["Omega_zetar"]] <- 1/2*prm$Omega_psir
  prm[["Omega_zetac"]] <- 1/2*prm$Omega_psic
  # matrix representation of tensor for transformed factors (using model.matrix)
  O3_is_col_na <- (colSums(is.na(O3@data)) == q)
  prm[["Omega"]] <- t(O3@data[ ,!O3_is_col_na, drop=F])
  prm[["eta_quad"]] <- get_eta_quad(prm$eta, prm$id, K) # transformed factors
  
  # Storage for outputs
  if (is.null(nwarmup)) {
    nwarmup <- round(niter/2) 
  }
  sims <- list(B = array(NA, dim=c(niter - nwarmup, K_int, q)),
               Bx = array(NA, dim=c(niter - nwarmup, p, q)),
               # b0 = array(NA, dim=c(niter - nwarmup, K_int)),
               # b0x = array(NA, dim=c(niter - nwarmup, p)), # main effect only
               # alpha = array(NA, dim=c(niter - nwarmup, q)), # outcome-specific intercept
               Sigma = array(NA, dim = c(niter - nwarmup, q, q)),
               #Theta = array(NA, dim = c(niter - nwarmup, p, K)),
               #eta = array(NA, dim = c(niter - nwarmup, n, K)),
               sigmax_sqinv = array(NA, dim = c(niter - nwarmup, p)),
               NULL
  )
  if (p_int > 0) {
    sims[["Bxz"]] = array(NA, dim=c(niter - nwarmup, p_int, p, q))
  }
  if (p_z > 0) {
    sims[["Bz"]] = array(NA, dim=c(niter - nwarmup, p_z, q))
  }
  if (include_rep) {
    sims[["Y_rep"]] = array(NA, dim=c(niter - nwarmup, dim(Y)))
  }
  if (!is.null(X_pred)) {
    sims[["E_Y_pred"]] = array(0, dim=c(nrow(X_pred), q))
    X_pred_int <- get_eta_int(X_pred, p, Z_pred, Z_int_pred, 1:nrow(X_pred))
  }
  if (include_interactions) {
    sims[["Omega_array_x"]] = array(NA, dim=c(niter - nwarmup, q, p, p))
    sims[["psir"]] = array(NA, dim=c(niter - nwarmup, K))
    sims[["psic"]] = array(NA, dim=c(niter - nwarmup, K))
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
    prm[["sigmax_sqinv"]] <- update_sigmax_sqinv(prm, Y, X, K, s0, r)
    ## Theta MGP ---------------------------------------------------------------
    tmp <- update_Theta_MGP(prm, Y, X, K)
    prm[names(tmp)] <- tmp
    ## eta using Adaptive Metroplis-within-Gibbs -------------------------------
    prm <- update_eta_mh(prm, Y, X, Z, Z_int, n, K, p, q, p_z, p_int,
                         s, adaptiveM=((s>n_till_adaptive) & adaptiveM),
                         adaptiveMWG=adaptiveMWG)
    ## Outcome regression ------------------------------------------------------
    tmp <- update_B_TPBN(prm, Y, X, K, Z, Z_int)
    prm[names(tmp)] <- tmp
    ## Factor interactions -----------------------------------------------------
    if (include_interactions) {
      tmp <- update_Omega(prm, Y, X, K, Z, Z_int, O3, O3_is_col_na)
      prm[names(tmp)] <- tmp
    }
    ## Outcome covariance ------------------------------------------------------
    tmp <- update_Sigma_IW_TPBN(prm, Y, Z, Z_int, K, binary)
    prm[names(tmp)] <- tmp
    
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
      # if (include_interactions) {
        # sims[["Omega"]][s-nwarmup,,] <- prm$Omega
      # }
      
      # Rescaling for probit latent variables
      # scale_B <- prm$B
      # scale_S <- prm$Sigma
      # if (sum(binary) > 0) { # TODO review this
      #   # Scale Sigma
      #   d <- sqrt(diag(prm$Sigma))
      #   temp_s <- d[binary == 1]
      #   d[binary == 1] <- 1
      #   R <- cov2cor(prm$Sigma)
      #   scale_S <- diag(d) %*% R %*% diag(d)
      #   # Scale B
      #   scale_B[, binary == 1] <- sweep(prm$B[, binary == 1, drop=FALSE], 2, temp_s, '/')
      # }
      sims[["B"]][s-nwarmup, , ] <- prm$B
      sims[["Sigma"]][s-nwarmup, , ] <- prm$Sigma
      
      # Effects in original predictors
      V <- solve(t(prm$Theta) %*% diag(prm$sigmax_sqinv, p, p) %*% prm$Theta + diag(1, K, K))
      Ax <- V %*% t(prm$Theta) %*% diag(prm$sigmax_sqinv, p, p)
      sims[["Bx"]][s-nwarmup, , ] <- t(crossprod(prm$B[1:K,,drop=F], Ax))
      temp_Bx <- sims[["Bx"]][s-nwarmup, , ]
      if (p_int > 0) {
        for (i in 1:p_int) {
          start_idx <- i*K + 1
          end_idx <- i*K + K
          sims[["Bxz"]][s-nwarmup, i, , ] <- t(crossprod(prm$B[start_idx:end_idx,,drop=F], Ax))
          temp_Bx <- rbind(temp_Bx, sims[["Bxz"]][s-nwarmup, i, , ])
        }
      }
      if (p_z > 0) {
        sims[["Bz"]][s-nwarmup, , ] <- prm$B[(K_int-p_z+1):K_int,]
        temp_Bx <- rbind(temp_Bx, sims[["Bz"]][s-nwarmup, , ])
      }
      stopifnot(dim(temp_Bx) == c(p+p*p_int+p_z, q))
      # temp_Bx <- matrix(0, ncol=p, nrow=q)
      # j <- 0
      # jx <- 0
      # while (j < (K + K*p_int)) {
      #   temp_Bx[, (jx+1):(jx+p)] <- crossprod(prm$B[(j+1):(j+K),,drop=F], Ax)
      #   j <- j + K
      #   jx <- jx + p
      # }
      # temp_Bx <- t(temp_Bx)
      # sims[["Bx"]][s-nwarmup, , ] <- temp_Bx
      
      
      if (include_interactions) {
        sims[["psir"]][s-nwarmup,] <- prm$Omega_psir
        sims[["psic"]][s-nwarmup,] <- prm$Omega_psic
        # zero out NAs
        prm[["Omega_tensor"]]@data[is.na(prm$Omega_tensor@data)] <- 0
        # k-mode tensor multiplication to transform this into coefs in X
        tmp_Omega <- rTensor::ttl(prm[["Omega_tensor"]], 
                                  list(t(Ax), t(Ax)), 
                                  c(2, 3))@data
        sims[["Omega_array_x"]][s-nwarmup,,,] <- abind::abind(
          apply(tmp_Omega,
                1,
                function(x) (x+t(x))/2,
                simplify = F),
          along = -1)
      }
      #sims[["b0x"]][s-nwarmup, ] <- Ax %*% prm$b0[1:K]
      
      # Fitted values for Y
      if (include_rep) {
        Y_rep <- tcrossprod(rep(1, nrow(Y)), prm$alpha) + 
          prm$eta_int %*% prm$B +
          MASS::mvrnorm(nrow(Y), rep(0,q), prm$Sigma)
        if (include_interactions) {
          Y_rep <- Y_rep + prm$eta_quad %*% prm$Omega
        }
        # if (sum(binary) > 0) {
        #   Y_rep[, binary == 1] <- 1 * (Y_rep[, binary == 1] > 0)
        # }
        sims[["Y_rep"]][s-nwarmup, , ] <- Y_rep
      }
      
      # Prediction for new data: returns E(y_i | x_i)
      if (!is.null(X_pred)) {
        Y_pred <- tcrossprod(rep(1, nrow(X_pred)), prm$alpha) + 
          X_pred_int %*% temp_Bx
        if (include_interactions) {
          trOV <- sapply(1:q, function(l) sum(diag(prm[["Omega_tensor"]]@data[l,,] %*% V)))
          Y_pred <- Y_pred + tcrossprod(rep(1, nrow(X_pred)), trOV)
          predict_interactions_cpp(Y_pred, X_pred,
                                   aperm(sims[["Omega_array_x"]][s-nwarmup,,,], c(2, 3, 1)),
                                   nrow(X_pred), q)
        }
        sims[["E_Y_pred"]] <- sims[["E_Y_pred"]] + (Y_pred / (niter-nwarmup))
      }
    }
  }
  message(paste("Metropolis-Hasting acceptance prob. of eta is", 
                round(sum(prm$eta_n_accepted)/(niter*n), 3)))
  message(paste("\nMetropolis-Hasting acceptance prob. of a1 is",
                round(prm$a1_n_accepted/niter, 3)))
  message(paste("\nMetropolis-Hasting acceptance prob. of a2 is",
                round(prm$a2_n_accepted/niter, 3)))
  return(sims)
}