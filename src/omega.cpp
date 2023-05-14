#include <RcppArmadillo.h>
#include <math.h>
#include "helpers.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat update_Omega_TPBN_cpp(const arma::mat& Omega, const arma::mat& X, 
                           const arma::mat& Y, const arma::mat& Sigma, 
                           const arma::vec& psi) {
  arma::mat Psii = arma::diagmat(1 / psi);
  arma::mat LL = X.t() * X + Psii;
  arma::mat Lu = trimatl(arma::chol(LL, "lower"));
  arma::mat Rv = trimatu(arma::chol(Sigma, "upper"));
  arma::mat M = arma::solve(trimatu(Lu.t()), arma::solve(Lu, X.t()*Y));
  return M + arma::solve(trimatu(Lu.t()), arma::randn<arma::mat>(Omega.n_rows, Omega.n_cols)) * Rv;
}

//[[Rcpp::export]]
void update_Omega_psi_cpp(arma::vec& psir, arma::vec& psic, 
                          const arma::cube& Omega, const arma::mat& Sigmainv, 
                          const arma::vec& zetar, const arma::vec& zetac, 
                          int K, int q) {
  arma::mat O;
  arma::vec s;
  arma::vec o;
  // update row-wise variances
  for (int r = 0; r < K; r++) {
    // parameters for the conditional
    double lam = 0.5 - 0.5 * q * (r + 1); // parameter for the conditional, u=0.5
    double rgig_psi = 2 * zetar(r);
    double chi = 0;
    for (int c = 0; c <= r; c++) {
      o = Omega.slice(c).col(r); // Omega is (q x K rows x K cols)
      chi += arma::as_scalar((o.t() * Sigmainv * o) / psic(c));
    }
    // sample new psi(i)
    psir(r) = std::max(rgig_cpp(lam, rgig_psi, chi), datum::eps);
  }
  
  // update column-wise variances
  for (int c = 0; c < K; c++) {
    // parameters for the conditional
    double chi = 0;
    for (int r = c; r < K; r++) {
      o = Omega.slice(c).col(r); // Omega is (q x K rows x K cols) = (row x col x slice)
      chi += arma::as_scalar((o.t() * Sigmainv * o) / psir(r));
    }
    double lam = 0.5 - 0.5 * q * (K - c); // parameter for the conditional, u=0.5
    double rgig_psi = 2 * zetac(c);
    // sample new psi(i)
    psic(c) = std::max(rgig_cpp(lam, rgig_psi, chi), datum::eps);
  }
}

// [[Rcpp::export]]
void update_Omega_zeta_cpp(arma::vec& zeta, const arma::vec& psi, 
                           int K, double global_shrink=1) {
  double v = 0.5;
  for (int i = 0; i < K; i++) {
    double rate = psi(i) + global_shrink;
    zeta(i) = arma::randg<double>(arma::distr_param(v, 1/rate));
  }
}

// [[Rcpp::export]]
void predict_interactions_cpp(arma::mat& Y, const arma::mat& X, 
                              const arma::cube& Omega, int n, int q) {
  arma::rowvec b(q);
  for (int i=0; i<n; i++) {
    for (int j=0; j<q; j++) {
      b(j) = arma::as_scalar(X.row(i) * Omega.slice(j) * X.row(i).t());
    }
    Y.row(i) += b;
  }
}