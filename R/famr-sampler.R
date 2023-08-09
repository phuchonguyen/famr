# TODO: bug in update_intercept when Z=NULL and Z_int=NULL
#' @id vector of length N, unique id of subjects, might contain duplicated id
#' 
#' @export
famr <- function(niter, Y, X, K=2, Z=NULL, Z_int=NULL, id=NULL, 
                 nthin=1, nwarmup=NULL,
                 adaptiveM=F, adaptiveMWG=F, adaptiveMWG_batch=50,
                 X_pred=NULL, Z_pred=NULL, Z_int_pred=NULL,
                 include_rep=FALSE, include_interactions=FALSE,
                 init_eta_eps=1e-2, init_a1_eps=1, init_a2_eps=1,
                 eta_eps_power=-0.5,
                 Ilod=NULL, loda=NULL, lodb=NULL,
                 s0=0.084, r=2.5, verbose=FALSE, missing_Y=FALSE, binary=NULL,
                 eta=NULL, Theta=NULL, B=NULL, Sigma=NULL, sigmax_sqinv=NULL,
                 W=NULL, xi=NULL, sigmay_sqinv=NULL,
                 random_intercept=FALSE,
                 time=NULL, Bt=NULL, n_kappa_opts=10, time_pred=NULL,
                 kappa_a=5.03, kappa_b=11.65, kappa_discrete=FALSE, kappa=NULL,
                 L=2, U=NULL, Lambda=NULL
                 ) {
  # stopifnot("K must be larger than 1" = K > 1)
  if (is.vector(Y)) Y <- matrix(Y, ncol = 1)
  if (nrow(Y) != nrow(X)) stop("Missmatched number of rows in Y and X")
  # when id=NULL, no repeated measurements of Y
  if (is.null(id)) id <- 1:nrow(Y) # assumes all subjects observed once
  if (is.null(time)) time <- rep(1, nrow(Y)) # assumes all subjects observed once
  # reduce X to a unique subject per row (when there are repeated measurements of Y)
  # each subjet has ONE UNIQUE measurement of X
  X <- X[!duplicated(id), ]
  q <- ncol(Y) # number of outcomes
  n <- length(unique(id)) # number of subjects, each subject might have T_i measurements of Y
  p <- ncol(X) # number of components in mixtures
  p_z <- ncol(Z) # number of covariates
  p_int <- ncol(Z_int) # number of covariates interacting with the mixture latent factors
  TT <- length(unique(time)) # number of unique time point

  # Mean-centered X to remove intercepts in factor model
  Xmean <- colMeans(X)
  X <- scale(X, center = T, scale = F)
  # Indicator matrix dim=dim(Y) of if there are missing values in Y:
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
  # recode id as numbers from 1 to n
  prm[['numeric_id']] <- as.numeric(factor(prm$id, levels = prm$uid))
  prm[["lpmf"]] <- rep(-99999999, n) # log likelihood
  # K_int is the total number of linear predictors including 
  # the K*p_int interactions, and p_z covariates
  K_int <- p_z + p_int*K
  # ********************************************
  # outcome-specific intercept (needed for probit model)
  # ********************************************
  prm[["alpha"]] <- rnorm(q)
  if (random_intercept) {
    prm[["xi"]] <- array(0, dim=c(n, q))
    prm[["nu_sqinv"]] <- rgamma(1, 2.5*0.5, 2.5*0.084*0.5) # between subject error
    prm[["sigmay_sqinv"]] <- 1 # within subject error
  }
  # ********************************************
  # TPBN prior on main effects of factors/covariates & interactions with covariates
  # ********************************************
  prm[["B"]] <- array(rnorm(q * K_int, 0, 100), dim = c(K_int, q))
  prm[["psi"]] <- rep(1/2*1/(q*sqrt(n * log(n))), K_int) # shrinkage on rows of B
  prm[["zeta"]] <- 1/2*prm$psi # shrinkage on rows of B
  # *********************************************
  # Observation noise in the outcomes
  # *********************************************
  prm[["Sigma"]] <- diag(rgamma(q, 2, 1), q, q) # covariance of y_i
  prm[["Sigmainv"]] <- diag(1/diag(prm$Sigma), q, q)
  # ********************************************
  # GP prior on time-varying main effects of factors
  # ********************************************
  prm[["Bt"]] <- array(rnorm(q * K * TT, 0, 10), dim = c(K, q, TT))
  # +++++++++++++++
  # 'Matrix' GP using Kronecker product and matrix normal distribution
  # Did not work well
  # +++++++++++++++
  prm[["Bt_eta"]] <- array(0, dim=dim(Y))
  # prm[["Bt_zeta"]] <- 1/rgamma(K, 0.5, 1)
  # prm[["Bt_phi"]] <- 1/rgamma(1, 0.5, 1)
  # prm[["Bt_psi"]] <- 1/rgamma(K, 0.5, rate=1/prm$Bt_zeta) # shrinkage on rows of Bt
  # prm[["Bt_tau"]] <- 1/rgamma(1, 0.5, rate=1/prm$Bt_phi) # global shrinkage for Bt
  # +++++++++++++++
  # Using factor model of independent GPs
  # +++++++++++++++
  prm[['Bt_Lambda']] <- matrix(rnorm(q * L), q, L)
  prm[['Bt_phi']] <- array(rgamma(q * L, 1), dim = c(q, L)) # local shrinkage for Theta
  prm[["Bt_delta"]] <- rgamma(L, 1) # global shrinkage param for factor K
  prm[['Bt_a1']] <- 2.1
  prm[['Bt_a2']] <- 3.1
  prm[['Bt_U']] <- array(rnorm(L * K * TT), dim = c(L, K, TT))
  # ***** Length-scale for GP *************************
  # From Kelly Moran's code: https://github.com/kelrenmor/bs3fa/tree/master
  # er_to_kappa <- function(er){sqrt(er/6)} # function to go from effective range to length scale kappa
  # is parameterized as sig^2 exp(-0.5 ||d-d'||^2 / kappa ^2)
  # For sig^2 exp(-phi ||d-d'||^2), back of envelope is effective range is 3/phi
  # phi = 0.5 l^(-2) ----> 3/phi = 6 l^2 (small ranges, i.e. smaller than range of data, cause issues)
  # so l = (effective_range / 6) ^(0.5)
  # er_min <- 1/(TT-1) + 1e-2 # corresponds to minimum effective range
  # er_max <- 1 - 1e-2 # corresponds roughly to eff range spanning all data
  # kappa_opts <- seq(er_to_kappa(er_min), er_to_kappa(er_max), length.out=n_kappa_opts)
  kappa_opts <- seq(1 + 1e-2, TT - 1 - 1e-2, length.out=n_kappa_opts) # TODO
  C_all <- lapply(kappa_opts, function(l) covEQ(1:TT, kappa=l, amplitude=1))
  Ci_all <- lapply(C_all, function(V) chol2inv(chol(V))) # avoid using solve because numerical instability can make Ci not symmetric
  ldetC_all <- lapply(C_all, function(V) log(det(V)))
  prm[["kappa"]] <- ifelse(is.null(kappa), kappa_opts[n_kappa_opts %/% 2], kappa) #lapply(1:K, function(i) kappa_opts[n_kappa_opts %/% 2])
  prm[["C"]] <- covEQ(1:TT, prm$kappa, 1) #lapply(1:K, function(i) C_all[[n_kappa_opts %/% 2]])
  prm[["C_inv"]] <- chol2inv(chol(prm$C)) #lapply(1:K, function(i) Ci_all[[n_kappa_opts %/% 2]])
  prm[["logdetC"]] <- log(det(prm$C))
  prm[["kappa_lpdf"]] <- -99999999
  prm[['kappa_eps']] <- 0.5 # range is exp(exp(2*0.5))=800x or exp(-exp(2*0.5))=1/800x of kappa
  prm[['kappa_n_accepted']] <- 0
  # ********************************************
  # Factor model
  # ********************************************
  prm[["sigmax_sqinv"]] <- rgamma(p, 2.5*0.5, 2.5*0.084*0.5) # noise variance
  prm[["eta"]] <- array(rnorm(n * K, 0, 1), dim = c(n, K)) # random factors
  prm[["eta_prop"]] <- array(0, dim=c(niter, n)) # storing proposal of eta_1 for debugging
  prm[["eta_n_accepted"]] <- rep(0, n) # MH acceptance
  prm[["eta_eps"]] <- rep(init_eta_eps, n) # MH proposal scales
  prm[["eta_A"]] <- array(0, dim = c(K, K, n)) # adaptive MH 
  prm[["eta_b"]] <- array(0, dim = c(K, n)) # adaptive MH
  prm[["eta_int"]] <- get_eta_int(prm$eta, K, Z, Z_int, id)
  # ********************************************
  # Multiplicative Gamma Process prior on the columns of loadings
  # ********************************************
  prm[["Theta"]] <- array(rnorm(p * K, 0, 1), dim = c(p, K)) # loadings
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
  sims <- list(B = array(NA, dim=c(niter - nwarmup, K, q, TT)),
               Bx = array(NA, dim=c(niter - nwarmup, p, q, TT)),
               psi = array(NA, dim=c(niter - nwarmup, K_int)),
               Lambda = array(NA, dim=c(niter - nwarmup, q, L)),
               # Bt_psi = array(NA,  dim=c(niter - nwarmup, K)),
               Bt_tau = array(NA,  dim=c(niter - nwarmup, L)),
               Bt_kappa = array(NA,  dim=c(niter - nwarmup, 1)),
               kappa_eps = array(NA,  dim=c(niter - nwarmup, 1)),
               Sigma = array(NA, dim = c(niter - nwarmup, q, q)),
               # nu_sqinv = matrix(NA, niter - nwarmup, 1),
               Theta = array(NA, dim = c(niter - nwarmup, p, K)),
               Theta_tau = array(NA, dim = c(niter - nwarmup, K)),
               alpha = array(NA, dim=c(niter - nwarmup, q)),
               eta = array(NA, dim = c(niter - nwarmup, n, K)),
               sigmax_sqinv = array(NA, dim = c(niter - nwarmup, p)),
               p_eta_accept = c(),
               eta_eps = c(),
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
    sims[["E_Y_pred"]] = array(0, dim=c(niter - nwarmup, nrow(X_pred), q))
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
  # prm[['xi']] <- xi
  # prm[["sigmax_sqinv"]] <- sigmax_sqinv
  n_till_adaptive <- min(3000, round(niter/10))
  for (s in 1:niter) {
    setTxtProgressBar(pb, s / niter)
    # MCMC steps
    prm[["sigmax_sqinv"]] <- update_sigmax_sqinv(prm, X, K, s0, r)
    ## Theta MGP ---------------------------------------------------------------
    tmp <- update_Theta_MGP(prm, X, K)
    prm[names(tmp)] <- tmp
    ## Outcome updates ---------------------------------------------------------
    ## Include random intercepts
    if (random_intercept) {
      ## eta using Adaptive Metropolis-within-Gibbs -----------------------------
      if (is.null(eta)) {
        # TODO: Update Ytilda to use time-varying main effects B(t)
        prm <- update_eta_mh_re(prm, Y, 
                                X, Z, Z_int, time,
                                n, K, p, q, p_z, p_int,
                                s, adaptiveM = ((s>n_till_adaptive) & adaptiveM),
                                adaptiveMWG = adaptiveMWG,
                                adaptiveMWG_batch = adaptiveMWG_batch,
                                eps_power = eta_eps_power)
      } else {
        prm[['eta']] <- eta # TODO: FOR DEBUGGING
      }
      ## Outcome regression ----------------------------------------------------
      # B_TPBN: update step for linear coefficient effects
      if (is.null(B)) {
        tmp <- tryCatch(update_B_TPBN_re(prm, Y, X, K, Z, Z_int),
                        error=function(e) {
                          cat('\nerror while sampling B TPBN:', message(e), '\n')
                        })
        prm[names(tmp)] <- tmp
      } else {
        prm[['B']] <- B
        prm[['alpha']] <- rep(0, q)
      }
      # B_GP: update step for coefficient functions of time
      if (is.null(Bt)) {
        if (is.null(Lambda)) {
          prm <- update_Lambda(Y, time, K, L, q, prm)
        } else {
          prm[['Bt_Lambda']] <- Lambda
        }
        if (is.null(U)) {
          prm <- update_U(prm, Y, time, K, L, q, TT, n)
        } else {
          prm[['Bt_U']] <- U
        }
        for (kk in 1:K) {
          prm[['Bt']][kk,,] <- prm$Bt_Lambda %*% prm$Bt_U[,kk,]
        }
        if (is.null(kappa)) {
          if (kappa_discrete) {
            #Update length-scale kappa using Bayes factor over a set of plausible values
            which_ls <- update_kappa_discrete(Ci_all, ldetC_all, prm$Bt_U,
                                              L, K, n_kappa_opts)
            prm[['kappa']] <- kappa_opts[which_ls]
            prm[['C']] <- C_all[[which_ls]]
            prm[['C_inv']] <- Ci_all[[which_ls]]
          } else {
            tmp <- update_kappa_cpp(prm$kappa, prm$C_inv, prm$C, prm$logdetC,
                                    prm$kappa_lpdf, 1:TT, prm$Bt_U, kappa_a, kappa_b,
                                    L, K, prm$kappa_eps, s, adaptiveMWG_batch,
                                    prm$kappa_n_accepted)
            prm[names(tmp)] <- tmp
          }
        }
        # +++++++++++++++++++++++++++++++++++
        #
        # Gibbs step for main effect functions
        # tmp <- update_B_GP(prm, Y, time, K, q, TT, n)
        # prm[names(tmp)] <- tmp
        # prm[['Bt_psi']] <- rep(1, K)
        # prm[['Bt_tau']] <- 1
        # tauphi <- update_B_GP_amplitude_cpp(prm$Bt_psi, prm$Bt_zeta,
        #                                     prm$Bt_tau, prm$Bt_phi, K, q, TT,
        #                                     abind::abind(prm$C_inv, along=3),
        #                                     prm$Sigmainv,
        #                                     aperm(prm$Bt, c(2, 3, 1))) # convert B from a (K, q, TT) to a (q, TT, K) array
        # prm[['Bt_tau']] <- tauphi[1]
        # prm[['Bt_phi']] <- tauphi[2]
        # # Update wavelength kappa with NUTS
        # kappa_f <- get_cur_kappa_f(unique(time), EQ_kernel_vec,
        #                            prm$Bt, prm$Sigmainv, q, K)
        # tmp <- mcnuts::nuts_iteration(prm$kappa,
        #                               kappa_f$log_p, kappa_f$gradient,
        #                               kappa_eps,
        #                               warmup = (s <= nwarmup),
        #                               M = c(1), verbose = FALSE,
        #                               t = s,
        #                               epsilon_bar = kappa_eps_bar, H = kappa_H,
        #                               # parameters for adapting kappa_eps
        #                               gamma = 0.05, t0 = 10, kappa = 0.75,
        #                               delta = 0.65, mu_eps = kappa_mu_eps
        # )
        # print(tmp$th)
        # prm[['kappa']] <- tmp$th
        # prm[["C"]] <- make_cov(unique(time), EQ_kernel_vec, kappa=prm$kappa, amplitude=1)
        # prm[["C_inv"]] <- solve(prm$C)
        # # While warmup is true, kappa_eps is adapted, o.w. kappa_eps = kappa_eps_bar
        # kappa_eps <- tmp$epsilon
        # kappa_eps_bar <- tmp$epsilon_bar
        # kappa_H <- tmp$H
      } else {
        # prm[['Bt']] <- Bt
        # prm[['Bt_psi']] <- rep(1, K)
        # prm[['Bt_tau']] <- 1
        # tauphi <- update_B_GP_amplitude_cpp(prm$Bt_psi, prm$Bt_zeta,
        #                                     prm$Bt_tau, prm$Bt_phi, K, q, TT,
        #                                     abind::abind(prm$C_inv, along=3), 
        #                                     prm$Sigmainv,
        #                                     aperm(prm$Bt, c(2, 3, 1))) # convert B from a (K, q, TT) to a (q, TT, K) array
        # prm[['Bt_tau']] <- tauphi[1]
        # prm[['Bt_phi']] <- tauphi[2]
        # which_ls <- lapply(1:K, function(k) {
        #   update_kappa_discrete(Ci_all, ldetC_all, prm$Bt[k,,,drop=F],
        #                         prm$Sigmainv, prm$Bt_psi*prm$Bt_tau,
        #                         q, 1, n_kappa_opts)
        # })
        # prm[['kappa']] <- lapply(which_ls, function(ls) kappa_opts[ls])
        # prm[['C']] <- lapply(which_ls, function(ls) C_all[[ls]])
        # prm[['C_inv']] <- lapply(which_ls, function(ls) Ci_all[[ls]])
      }
      prm[['Bt_eta']] <- sapply(1:nrow(Y), function(i) 
        prm$eta[prm$numeric_id[i],] %*% prm$Bt[,,time[i]]) %>% t()
      if (q == 1) {
        prm[['Bt_eta']] <- matrix(prm[['Bt_eta']], nrow=nrow(Y), ncol=ncol(Y))
      }
      stopifnot(dim(prm[['Bt_eta']]) == dim(Y))
      ## Factor interactions ---------------------------------------------------
      # if (include_interactions) { # TODO
      #   tmp <- update_Omega(prm, Y - prm$ri, X, K, Z, Z_int, O3, O3_is_col_na)
      #   prm[names(tmp)] <- tmp
      # }
      ## Random intercept ------------------------------------------------------
      if (is.null(xi)) {
        # prm <- update_random_intercept(prm, Y)
      } else {
        prm[["xi"]] <- xi # TODO: FOR DEBUGGING
      }
      ## Outcome covariance ----------------------------------------------------
      if (is.null(Sigma)) {
        tmp <- update_Sigma_IW_TPBN_re(prm, Y, Z, Z_int, K, TT, binary, n)
        prm[names(tmp)] <- tmp
      } else {
        prm[['Sigma']] <- Sigma
        prm[['Sigmainv']] <- chol2inv(chol(Sigma))
      }
    ## Model without decomposing the noise variance into between individual variances
    # and within individual variances
    } else {
      ## eta using Adaptive Metroplis-within-Gibbs -----------------------------
      if (is.null(eta)) {
        prm <- update_eta_mh(prm, Y, X, Z, Z_int, n, K, p, q, p_z, p_int,
                             s, adaptiveM = ((s>n_till_adaptive) & adaptiveM),
                             adaptiveMWG = adaptiveMWG,
                             adaptiveMWG_batch=adaptiveMWG_batch)
      } else {
        prm[['eta']] <- eta
      }
      ## Outcome regression ----------------------------------------------------
      tmp <- update_B_TPBN(prm, Y, X, K, Z, Z_int)
      prm[names(tmp)] <- tmp
      ## Factor interactions ---------------------------------------------------
      if (include_interactions) {
        tmp <- update_Omega(prm, Y, X, K, Z, Z_int, O3, O3_is_col_na)
        prm[names(tmp)] <- tmp
      }
      ## Outcome covariance ----------------------------------------------------
      tmp <- update_Sigma_IW_TPBN(prm, Y, Z, Z_int, K, binary)
      prm[names(tmp)] <- tmp
    }
    
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
      # # stop adaptive MH
      # adaptiveMWG = F
      # adaptiveM = F
      # sims[["p_eta_accept"]] <- c(sims[["p_eta_accept"]],
      #                             mean(prm$eta_n_accepted/(s-nwarmup)))
      # store adaptive MH diagnostics
      if ((s %% adaptiveMWG_batch == 0)) {
        sims[["p_eta_accept"]] <- c(sims[["p_eta_accept"]],
                                    mean(prm$eta_n_accepted/adaptiveMWG_batch))
        sims[["eta_eps"]] <- rbind(sims[["eta_eps"]], prm$eta_eps)
      }
      
      sims[["eta"]][s-nwarmup,,] <- prm$eta
      sims[["sigmax_sqinv"]][s-nwarmup, ] <- prm$sigmax_sqinv
      sims[["Theta"]][s-nwarmup, , ] <- prm$Theta
      sims[["Theta_tau"]][s-nwarmup, ] <- cumprod(prm$delta)
      sims[["alpha"]][s-nwarmup, ] <- prm$alpha
      sims[['Bt_kappa']][s-nwarmup, ] <- unlist(prm$kappa)
      sims[['Bt_tau']][s-nwarmup, ] <- cumprod(prm$Bt_delta)
      sims[['kappa_eps']][s-nwarmup, ] <- prm$kappa_eps
      sims[['Lambda']][s-nwarmup,,] <- prm$Bt_Lambda
      
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
      sims[["B"]][s-nwarmup, , ,] <- prm$Bt
      if (p_int + p_z > 0) {
        sims[['psi']][s-nwarmup, ] <- prm$psi
      }
      sims[["Sigma"]][s-nwarmup, , ] <- prm$Sigma
      # if (random_intercept) {
      #   sims[['nu_sqinv']][s-nwarmup, ] <- prm$nu_sqinv
      # }
      
      # transformation to convert latent factor effects to effects in original predictors
      V <- solve(t(prm$Theta) %*% diag(prm$sigmax_sqinv, p, p) %*% prm$Theta + diag(1, K, K))
      Ax <- V %*% t(prm$Theta) %*% diag(prm$sigmax_sqinv, p, p) # K \times p
      # TODO Update so that Bx is a (p * q * TT) array
      for (tt in 1:TT) {
        if (K > 1) {
          sims[["Bx"]][s-nwarmup, , , tt] <- t(crossprod(prm$Bt[,,tt], Ax)) 
        } else {
          sims[["Bx"]][s-nwarmup, , , tt] <- t(prm$Bt[,,tt] %*% Ax)
        }
      }
      # variable to store all linear effects
      temp_Bxz <- c()
      # linear interaction between chemicals and covariates in original predictors
      if (p_int > 0) {
        for (i in 1:p_int) {
          start_idx <- (i-1)*K + 1
          end_idx <- (i-1)*K + K
          sims[["Bxz"]][s-nwarmup, i, , ] <- t(crossprod(prm$B[start_idx:end_idx,,drop=F], Ax))
          temp_Bxz <- rbind(temp_Bxz, sims[["Bxz"]][s-nwarmup, i, , ])
        }
      }
      # Linear main effects of covariates
      if (p_z > 0) {
        sims[["Bz"]][s-nwarmup, , ] <- prm$B[(K_int-p_z+1):K_int,]
        temp_Bxz <- rbind(temp_Bxz, sims[["Bz"]][s-nwarmup, , ])
      }
      # TODO: doesn't work when q = 1
      stopifnot(dim(temp_Bxz) == c(p*p_int+p_z, q))
      
      # Fitted outcome values for observed data
      if (include_rep) {
        Y_rep <- tcrossprod(rep(1, nrow(Y)), prm$alpha) + 
          prm$Bt_eta +
          prm$eta_int %*% prm$B +
          MASS::mvrnorm(nrow(Y), rep(0,q), prm$Sigma)
        if (include_interactions) {
          Y_rep <- Y_rep + prm$eta_quad %*% prm$Omega
        }
        if (random_intercept) {
          Y_rep <- Y_rep + prm$xi[prm$numeric_id, ]
        }
        # if (sum(binary) > 0) {
        #   Y_rep[, binary == 1] <- 1 * (Y_rep[, binary == 1] > 0)
        # }
        sims[["Y_rep"]][s-nwarmup, , ] <- Y_rep
      }
      
      # Predicted outcomes for new data: returns E(y_i | x_i)
      if (!is.null(X_pred)) {
        Y_pred <- tcrossprod(rep(1, nrow(X_pred)), prm$alpha) +
          X_pred_int %*% temp_Bxz +
          sapply(1:nrow(X_pred), function(i) {
          t(X_pred[i,]) %*% sims[["Bx"]][s-nwarmup,,,time_pred[i]]
          }) %>% t()
        if (include_interactions) {
          trOV <- sapply(1:q, function(l) sum(diag(prm[["Omega_tensor"]]@data[l,,] %*% V)))
          Y_pred <- Y_pred + tcrossprod(rep(1, nrow(X_pred)), trOV)
          predict_interactions_cpp(Y_pred, X_pred,
                                   aperm(sims[["Omega_array_x"]][s-nwarmup,,,], c(2, 3, 1)),
                                   nrow(X_pred), q)
        }
        sims[["E_Y_pred"]][s-nwarmup, , ] <- Y_pred
      }
    }
  }
  
  if (s >= adaptiveMWG_batch) {
    message(paste("\nMetropolis-Hasting acceptance prob. of eta is",
                  round(tail(sims[["p_eta_accept"]], 1), 3)))
    sims[["eta_prop"]] <- prm$eta_prop
  }
  message(paste("\nMetropolis-Hasting acceptance prob. of a1 is",
                round(prm$a1_n_accepted/niter, 3)))
  message(paste("\nMetropolis-Hasting acceptance prob. of a2 is",
                round(prm$a2_n_accepted/niter, 3)))
  cat('\nGP length scale kappa =', paste(prm$kappa, sep=" "), '\n')
  return(sims)
}
