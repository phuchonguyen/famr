#include <RcppArmadillo.h>
#include <math.h>
#include "helpers.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;

// *******************************************************
// Random intercept
// *******************************************************
// [[Rcpp::export]]
void update_intercept_cpp(arma::mat& alpha, arma::vec& mu, arma::vec& v,
                          const arma::mat& Y, const arma::uvec& id, const arma::uvec& uid,
                          const arma::mat& Sigma, const arma::mat& Sigmainv,
                          double v0, double s0, bool random_intercept, 
                          int n, int N, int q) {
  if (random_intercept) {
    mu = arma::mvnrnd(arma::mean(Y, 0).t(), (Sigma + arma::diagmat(v)) / N);
    arma::vec Yib(q);
    arma::mat C(q, q);
    arma::uvec idi;
    arma::mat Vinv = diagmat(1/v);
    arma::vec m0 = mu / v;
    for (int i = 0; i < n; i++) {
      idi = find(id == uid(i));
      Yib = arma::sum(Y.rows(idi), 0).t();
      C = (Sigmainv*idi.n_elem + Vinv).i();
      alpha.rows(idi) = ones(idi.n_elem) * arma::mvnrnd(C * (Sigmainv*Yib + m0), C).t();
    }
    double vn = (v0 + n)*0.5;
    arma::mat ualpha = alpha.rows(arma::find_unique(id));
    arma::rowvec ssr = arma::sum(arma::square(ualpha.each_row() - mu.t()), 0);
    for (int j = 0; j < q; j++) {
      v(j) = arma::randg(distr_param(vn, 2/(ssr(j) + v0*s0)));
    }
  } else {
    mu = arma::mvnrnd(arma::mean(Y, 0).t(), Sigma / N);
    v.fill(0);
    alpha = ones(N) * mu.t();
  }
}

// *******************************************************
// Dirichlet Laplace prior
// *******************************************************
// [[Rcpp::export]]
void update_B_DL_cpp(arma::mat& B, const arma::mat& Y, 
                     const arma::mat& X, const arma::mat& Sigma_inv, 
                     const arma::mat& V, int n, int q, int p) {
  arma::mat eps(p, q, fill::value(datum::eps));
  arma::mat V_eps = arma::max(V, eps); // add stability to later inversion
  arma::mat Vj_inv(q, q);
  arma::mat E = Y - X*B;
  arma::mat Ej(E);
  arma::mat Cj(q, q);
  arma::vec mj(q);
  for (int j = 0; j < p; j++) {
    Vj_inv = diagmat(1/V_eps.row(j)); // tau_j^2(psi_{j1}^2zeta_{j1},...,)
    Ej = E + X.col(j)*B.row(j);
    Cj = inv_sympd(Sigma_inv * as_scalar(accu(pow(X.col(j), 2))) + Vj_inv);
    mj = sum((Ej.each_col() % X.col(j)) * Sigma_inv, 0).t(); // CHECK
    B.row(j) = mvnrnd(Cj * mj, Cj).t();
  }
  B = arma::max(B, eps); // prevent numerical rounding to zero
}

// [[Rcpp::export]]
void update_nu_DL_cpp(arma::vec& nu, const arma::mat& B, 
                       const arma::mat& psi, const arma::mat& zeta,
                       int p, int q) {
  arma::vec a = sum(abs(B) / psi, 1);
  for (int j = 0; j < p; j++) {
    nu(j) = std::max(rgig_cpp(1 - q, 1, 2 * a(j)), datum::eps);
  }
}

// [[Rcpp::export]]
void update_psi_DL_cpp(arma::mat& psi, const arma::mat& B, int p, int q) {
  double a = 0.1; // prior on lambda parameter
  arma::mat A = abs(B);
  arma::rowvec t(q);
  for (int j = 0; j < p; j++) {
    for (int k = 0; k < q; k++) {
      t(k) = std::max(rgig_cpp(a - 1, 1, 2 * A(j, k)), datum::eps); // prevent numerical rounding to zero
    }
    psi.row(j) = t / accu(t);
  }
}

// [[Rcpp::export]]
void update_zeta_DL_cpp(arma::mat& zeta, const arma::mat& psi, 
                        const arma::mat& B, const arma::vec& nu) {
  zeta = (psi.each_col() % nu) / abs(B);
  zeta.transform(rig_cpp);
}

// *******************************************************
// Three parameter beta normal prior conditional on Sigma
// *******************************************************

// [[Rcpp::export]]
arma::vec update_zeta_TPBN_cpp(arma::vec psi, int q, double global_shrink=1) {
  arma::vec zeta(q);
  double v = 0.5;
  for (int i = 0; i < q; i++) {
    double rate = psi(i) + global_shrink;
    zeta(i) = arma::randg<double>(arma::distr_param(v, 1/rate));
  }
  
  return zeta;
}

// q: number of predictors, p: number of outcomes; B: q x p matrix
// [[Rcpp::export]]
arma::vec update_psi_TPBN_cpp(arma::mat B, arma::mat Sigmainv,
                         arma::vec zeta, int q, int p) {
  arma::vec psi(q);
  double u = 0.5;
  double lam = u - 0.5 * p;
  for (int i = 0; i < q; i++) {
    double chi = arma::as_scalar(B.row(i)*Sigmainv*(B.row(i).t()));
    double rgig_psi = 2 * zeta(i); 
    psi(i) = rgig_cpp(lam, rgig_psi, chi);
  }
  return psi;
}

// Paper: https://arxiv.org/pdf/1711.07635.pdf
// B: p x q matrix
// [[Rcpp::export]]
arma::mat update_B_TPBN_cpp(arma::mat X, arma::mat Y, arma::mat Sigma, 
                        arma::vec psi,
                        int p, int q) {
  arma::mat Psii = arma::diagmat(1 / psi);
  arma::mat LL = (X.t() * X) + Psii;
  arma::mat Lu = trimatl(arma::chol(LL, "lower"));
  arma::mat Rv = trimatu(arma::chol(Sigma, "upper"));
  arma::mat M = arma::solve(trimatu(Lu.t()), arma::solve(Lu, X.t()*Y));
  arma::mat B = M + arma::solve(trimatu(Lu.t()), arma::randn<arma::mat>(p, q)) * Rv;
  return B;
}

// // Returns u0_square
// // [[Rcpp::export]]
// arma::vec update_u0_hCauchy(arma::vec b0, arma::vec c0, double tau0_sq, int q, int p) {
//   arma::vec u0_sq(q);
//   for (int i = 0; i < q; i++) {
//     double chi = b0(i) * b0(i) / tau0_sq;
//     double lam = 0;
//     double ps = 2 * c0(i); 
//     u0_sq(i) = rgig_cpp(lam, chi, ps);
//   }
//   
//   return u0_sq;
// }
// 
// // Returns tau0_square
// // [[Rcpp::export]]
// double update_tau0_hCauchy(arma::vec b0, arma::vec u0_sq, double d0, int q, int p) {
//   double chi = arma::sum(b0 % b0 / u0_sq);
//   double lam = 0.5 * (1 - q);
//   double ps = 2 * d0; 
//   double tau0_sq = rgig_cpp(lam, chi, ps);
//   
//   return tau0_sq;
// }
// 
// // [[Rcpp::export]]
// arma::vec update_c0(arma::vec u0_sq, int q) {
//   arma::vec c0(q);
//   for (int i = 0; i < q; i++) {
//     double rate = u0_sq(i) + 1;
//     c0(i) = arma::randg<double>(arma::distr_param(0.5, 1/rate));
//   }
//   
//   return c0;
// }
// 
// // [[Rcpp::export]]
// double update_d0(double tau0_sq) {
//   double rate = tau0_sq + 1;
//   return arma::randg<double>(arma::distr_param(0.5, 1/rate));
// }
// 
// // [[Rcpp::export]]
// arma::vec update_b0_Horseshoe(arma::mat X, arma::mat Y, arma::mat Sigmainv,
//                               arma::vec psi, arma::vec u0_sq, double tau0_sq,
//                               int q, int p, int n) {
//   arma::mat XtX = X.t() * X;
//   arma::mat U0i = arma::inv(arma::diagmat(1 / psi) + XtX);
//   arma::mat XU = X.t() - XtX * U0i * X.t();
//   arma::mat XUX = XtX - XtX * U0i * XtX;
//   arma::mat Un = arma::inv(XUX + arma::diagmat(1 / u0_sq) / tau0_sq);
//   double vn = 1 / (arma::as_scalar(ones(p).t() * Sigmainv * ones(p)) + 1);
//   arma::vec mn = Un * XU * Y * Sigmainv * ones(p) * vn;
//   arma::vec b0 = arma::mvnrnd(mn, Un * vn);
//   
//   return b0;
// }