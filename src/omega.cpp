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
arma::uvec is_finite_rows(const arma::mat& A) {
// returns a vector of indices of rows with at least one nonfinite (NA, inf) element
  arma::uvec res(A.n_rows, fill::zeros);
  for (unsigned int i = 0; i < A.n_rows; i++) {
    if (A.row(i).is_finite()) {
      res(i) = 1;
    }
  }
  return res;
}

//[[Rcpp::export]]
void update_Omega_psi_cpp(arma::vec& psi, const arma::cube& Omega, 
                          const arma::mat& Sigmainv, const arma::vec& psi2, 
                          const arma::vec& zeta, int K, int q) {
  arma::mat O;
  arma::vec s;
  arma::uvec finite(K);
  double lam = 0.5 - 0.5 * K * q; // parameter for the conditional
  for (int i = 0; i < K; i++) {
    // get matrix Omega[i,,]
    O = Omega.row(i); // K x q matrix
    if (q == 1 && K > 1) { // in this case, Omega.row(i) returns a 1xK row vector
      O = arma::mat(O.t()); // turn into K x 1 matrix
    }
    // remove rows with only NAs
    finite = is_finite_rows(O);
    O.shed_rows(find(finite == 0)); // TODO find rows with only NA values
    s = psi2.elem(find(finite == 1));
    // parameters for the conditional
    double rgig_psi = 2 * zeta(i);
    double chi = arma::trace(Sigmainv * O.t() * (O.each_col() / s));
    // sample new psi(i)
    psi(i) = std::max(rgig_cpp(lam, rgig_psi, chi), datum::eps);
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