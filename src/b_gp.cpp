#include <RcppArmadillo.h>
#include <math.h>
#include "helpers.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;

/* O: (n x TT) array of indicator of existences of Y_it, where n is number 
 *    of subjects, TT is number of unique time points
 * time: \sum_i T_i vector of time points, must be numeric starting from 1 to T
 * id: \sum_i T_i vector of numeric ids for each observation, starting from 1 to n
 * n: number of unique subjects
 * M: (q, T, K)
 * V: (T, T, K)
*/
arma::cube update_B_GP_cpp(const arma::mat& Y, const arma::vec& time,
                     const arma::cube& Bt_eta_k, const arma::mat& eta, 
                     const arma::cube& C_inv, const arma::vec& psi, 
                     const arma::vec& id, 
                     double sinv, const arma::mat& Sigma, 
                     int K, int q, int T, int n) {
  arma::cube B(q, T, K);
  arma::uvec idi;
  arma::mat L;
  arma::mat U;
  arma::mat M_tmp;
  for (int k = 0; k < K; k++) {
    arma::mat M_k = arma::zeros(q, T);
    arma::mat V_k = arma::zeros(T, T);
    for (int i = 0; i < n; i++) {
      // indices of all observations of subject i
      idi = find(id == (i+1)); // id is indexed in R starting from 1 to n
      // Y_{i,-k} * \eta_{ik} at observed time points for subject i 
      M_tmp = (Y.rows(idi) + Bt_eta_k.slice(k).rows(idi)).t() * eta(i, k);
      for (int t = 0; t < idi.n_elem; t++) {
        // get index of time point t = 1,..., T_i
        int time_it = time(idi(t)) - 1; // time is indexed in R starting from 1 to T
        // add contribution of Y_{i,-k, t} * \eta_{ik} to observed time point t of mean
        M_k.col(time_it) += M_tmp.col(t);
        // add contribution of \eta_{ik} to the diagonal indexed by observed time point t of variance
        V_k(time_it, time_it) += eta(i, k);
      }
    }
    M_k *= sinv;
    V_k *= sinv;
    V_k += (C_inv.slice(k) / psi(k)) + 1e-6*arma::eye(T, T); // for numerical stability
    if (!V_k.is_sympd()) {
      cout << "V_k not symmetric positive definite and is " << V_k.is_symmetric() << " symmetric" << std::endl;
    }
    // sample Bk an (q \times T) matrix of regression coefs
    //    of factor k for q outcomes at T time points.
    try {
      L = trimatl(arma::chol(V_k, "lower"));
    } catch(const std::exception &exc) {
      std::cerr << exc.what();
      cout << "\n update_B_GP_cpp chol(Vk)" << std::endl;
    }
    try {
      U = trimatu(arma::chol(Sigma, "upper")); // LL' = Sigma
    } catch(...) {
      cout << "\n update_B_GP_cpp chol(Sigma)" << std::endl;
    }
    arma::mat M = arma::solve(trimatu(L.t()), arma::solve(L, M_k.t()));
    B.slice(k) = (M + arma::solve(trimatu(L.t()), arma::randn<arma::mat>(T, q)) * U).t();
  }
  return B;
}

// Horse shoe prior as in Makalic & Schmidt 2015
// B is a (q x T x K) cube
arma::vec update_B_GP_amplitude_cpp(arma::vec& psi_sq, arma::vec& zeta,
                           double tau_sq, double phi, int K, int q, int T,
                           const arma::cube& Ci, const arma::mat& Si,
                           const arma::cube& B) {
  arma::mat CBSB(T, T);
  double tau_rate = 0;
  for (int k=0; k < K; k++) {
    // psi_k \sim IG()
    CBSB = Ci.slice(k) * B.slice(k).t() * Si * B.slice(k);
    double rate = 1.0/zeta(k) + 0.5*arma::trace(CBSB)/tau_sq;
    psi_sq(k) = 1.0 / arma::randg<double>(distr_param(0.5*(q*T+1), 1.0/rate));
    // zeta_k \sim IG(1, 1+1/psi_k)
    zeta(k) = 1.0 / arma::randg<double>(distr_param(1.0, psi_sq(k)/(1+psi_sq(k))));
    tau_rate += arma::trace(CBSB)/psi_sq(k);
  }
  // tau_sq
  tau_rate *= 0.5;
  tau_rate += 1.0 / phi;
  tau_sq = 1.0 / arma::randg<double>(distr_param(0.5*(K*q*T+1), 1.0/tau_rate));
  // phi \sim IG(1, 1+1/tau_sq)
  phi = 1.0 / arma::randg<double>(distr_param(1.0, tau_sq/(1+tau_sq)));

  arma::vec tauphi(2);
  tauphi(0) = tau_sq;
  tauphi(1) = phi;
  return tauphi;
}

// q: number of outcomes; T: number of time points, B_k: (q, T, K) array
arma::vec update_psi_GP_cpp(const arma::cube& B, const arma::mat& Sigmainv,
                            const arma::mat& C_inv, const arma::vec& zeta,
                            int K, int q, int T) {
  arma::vec psi(K);
  double u = 0.5;
  double lam = u - 0.5 * q * T;
  for (int k = 0; k < K; k++) {
    double chi = arma::trace(C_inv * B.slice(k).t() * Sigmainv * B.slice(k));
    double rgig_psi = 2 * zeta(k);
    psi(k) = rgig_cpp(lam, rgig_psi, chi);
  }
  return psi;
}

arma::vec update_zeta_GP_cpp(arma::vec psi, int K, double global_shrink=1) {
  arma::vec zeta(K);
  double v = 0.5;
  for (int i = 0; i < K; i++) {
    double rate = psi(i) + global_shrink;
    zeta(i) = arma::randg<double>(arma::distr_param(v, 1/rate));
  }

  return zeta;
}
