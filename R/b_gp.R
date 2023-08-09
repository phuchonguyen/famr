# # Updates the matrix of regression functions B(t)
# update_B_GP <- function(prm, Y, time, K, q, TT, n) {
#   Ytilde <- Y - tcrossprod(rep(1, nrow(Y)), prm$alpha) - prm$Bt_eta - prm$eta_int %*% prm$B
#   # creates a array size (\sum_i T_i, q, K) of \eta_k multiplied by its effect B(t)_k
#   Bt_eta_k <- lapply(1:nrow(Y), function(i) t(prm$eta[prm$numeric_id[i],] * prm$Bt[,,time[i]])) %>%
#     abind::abind(along=-1)
#   sinv <- 1/(1/prm$sigmay_sqinv + 1/prm$nu_sqinv)
#   # M <- array(0, dim=c(q, TT, K))
#   # V <- array(0, dim=c(TT, TT, K))
#   # convert list of k Ck_inv to an (TT, TT, K) array
#   Ci_cube <- abind::abind(prm$C_inv, along=3)
#   # returns a (q, TT, K) array
#   B <- update_B_GP_cpp(Ytilde, time, Bt_eta_k, prm$eta,
#                        Ci_cube, prm$Bt_psi*prm$Bt_tau,
#                        prm$numeric_id, sinv, prm$Sigma,
#                        K, q, TT, n)
#   # B <- array(0, dim=c(q, TT, K))
#   # for (k in 1:K) {
#   #   B[,,k] <- matrixNormal::rmatnorm(s=1, M=M[,,k]%*%V[,,k], U=prm$Sigma, V=V[,,k],
#   #                                    tol=1e-7)
#   # }
#   # psi <- update_psi_GP_cpp(B, prm$Sigmainv, prm$C_inv, prm$Bt_zeta,
#   #                          K, q, TT)
#   # zeta <- update_zeta_GP_cpp(psi, K, global_shrink=0.5)#1/(q*sqrt(n*log(n))))
#   B <- aperm(B, c(3, 1, 2)) # convert B to a (K, q, TT) array
#   return(list(Bt = B
#               # Bt_psi = psi, Bt_zeta = zeta
#               ))
# }

# EQ_kernel <- function(t1, t2, kappa=1, amplitude=1) {
#   # Parametrization: exp{-kappa 0.5(t - t')^2}
#   # amplitude * exp(-(t1 - t2)^2 * kappa / 2)
#   # Parametrization: exp{-(t-t')^2 / 2*kappa^2}
#   amplitude * exp(-0.5 * abs(t1 - t2)^1.9999 / kappa^2) # the 1.9999 is for numerical stability
# }
# 
# EQ_kernel_vec <- Vectorize(EQ_kernel)

# # kernel: vectorized kernel function
# # t: vector of unique time points
# make_cov <- function(t, kernel_vec, kappa, amplitude, eps=1e-5) {
#   outer(t, t, kernel_vec, kappa=kappa, amplitude=amplitude)
#   # diag(c) <- 1.0 # replace the diagonal with 1.0
# }
# 
# matrix_trace <- function(A) {
#   sum(diag(A))
# }

# # Return negative log likelihoods and gradient functions that only takes kappa 
# #   and return values conditional on current observations of B, Si
# get_cur_kappa_f <- function(t, kernel_vec, B, Si, q, K) {
#   A <- ss_kappa(B, Si, K)
#   cur_log_p <- function(kappa) {
#     log_p_kappa(kappa, t, kernel_vec, A, q, K)
#   }
#   
#   cur_gradient <- function(kappa) {
#     gradient_kappa(kappa, t, kernel_vec, A, q, K)
#   }
#   
#   return(list('log_p'=cur_log_p, 'gradient'=cur_gradient))
# }

# helper for calculating the gradient and negative log like
# ss_kappa <- function(B, Si, psi, K) {
#   Reduce('+', lapply(1:K, function(k) t(B[k,,]) %*% Si %*% B[k,,] / psi[k]))
# }

# # Negative log likelihood of the wavelength
# # t: vector of unique time points
# # A: sum of squares returned by ss_kappa
# log_p_kappa <- function(kappa, t, kernel_vec, A, q, K) {
#   if (is.na(kappa) | kappa <= 0)  # returns zero if kappa is negative
#     return(-Inf)
#   C <- make_cov(t, kernel_vec, kappa=kappa, amplitude=1)
#   Ci <- solve(C)
#   # choose this prior specification so that there's at least 90% chance that kappa is below 0.2
#   # kappa \sim Gamma(0.05, 1)
#   log_prior <- dgamma(kappa, 10, 15, log=TRUE) # TODO explain choice of prior specification based on data
#   log_like <- -q*K*0.5*log(det(C)) - 0.5*matrix_trace(Ci %*% A) + log_prior
#   return(log_like)
# }

# gradient_EQ_kernel <- function(t1, t2, kappa=1, amplitude=1) {
#   # Parametrization: exp{-kappa 0.5(t - t')^2}
#   # - EQ_kernel(t1, t2, kappa=kappa, amplitude=amplitude) * (t1 - t2)^2 * 0.5
#   # Parametrization: exp{-(t-t')^2 / 2*kappa^2}
#   EQ_kernel(t1, t2, kappa=kappa, amplitude=amplitude) * (t1 - t2)^2 / kappa^3
# }
# 
# gradient_EQ_vec <- Vectorize(gradient_EQ_kernel)

# # Gradient of the negative log likelihood
# # kappa: wavelength, get gradient for this
# # t: vector of unique time points
# # kernel_vec: vectorized kernel function k(t, t')
# gradient_kappa <- function(kappa, t, kernel_vec, A, q, K) {
#   if (is.na(kappa) | kappa <= 0)  # returns zero if kappa is negative
#     return(0)
#   C <- make_cov(t, kernel_vec, kappa=kappa, amplitude=1)
#   Ci <- solve(C)
#   dCdkappa <- make_cov(t, gradient_EQ_vec, kappa=kappa, amplitude=1)
#   gradient <- -q*K*0.5*matrix_trace(Ci %*% dCdkappa) -
#     0.5*matrix_trace(-Ci %*% A %*% Ci %*% dCdkappa) +
#     (-15 + (10-1)/kappa) # prior
#   cat('\nGradient of ', kappa, ' is ', gradient)
#   return(gradient)
# }
# 
# gradient_kappa_numeric <- function(kappa, t, kernel_vec, A, q, K, d=0.001) {
#   if (is.na(kappa) | kappa <= 0.1)  # returns zero if kappa is negative
#     return(0)
#   lp_upper <- log_p_kappa(kappa+d, t, kernel_vec, A, q, K)
#   lp_lower <- log_p_kappa(kappa-d, t, kernel_vec, A, q, K)
#   return((lp_upper-lp_lower)/(2*0.001))
# }

# Update kappa discrete options
# Equal prior probabilities over the options of kappa
# Returns index of the chosen kappa
# update_kappa_discrete <- function(Ci_all, ldetC_all, B, Si, psi, q, K, n_kappa_opts) {
#   A <- ss_kappa(B, Si, psi, K)
#   log_prob <- rep(NA, n_kappa_opts)
#   for (i in 1:n_kappa_opts) {
#     log_prob[i] <- -q*K*0.5*ldetC_all[[i]] - 0.5*matrix_trace(Ci_all[[i]] %*% A)
#   }
#   # From Kelly Moran's code: https://github.com/kelrenmor/bs3fa/blob/master/R/sample_lind.R
#   maxlg <- max(log_prob)
#   # https://stats.stackexchange.com/questions/66616/converting-normalizing-very-small-likelihood-values-to-probability
#   probls <- exp(log_prob - maxlg)
#   probls <- probls/sum(probls)
#   which_ls <- sample(x=(1:n_kappa_opts), 1, replace=T, prob=probls)
#   return(which_ls)
# }

