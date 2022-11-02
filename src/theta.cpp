#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;

// [[Rcpp::export]]
arma::mat update_Theta_MGP_cpp(arma::mat eta, arma::vec sigmax_sqinv, 
                               arma::mat phi, arma::vec delta, arma::vec tau, 
                               int K, int p, arma::mat X) {
  arma::mat Theta(p, K);
  arma::mat ete = eta.t() * eta;
  for (int j = 0; j < p; ++j) {
    arma::mat Dj = arma::diagmat(phi.row(j).t() % tau);
    arma::mat L = trimatl(arma::chol(Dj + sigmax_sqinv(j) * ete, "lower"));
    arma::vec m =  arma::solve(trimatu(L.t()), arma::solve(L, eta.t() * X.col(j) * sigmax_sqinv(j)));
    Theta.row(j) = arma::solve(L, randn<arma::vec>(K)).t() + m.t();
  }
  return Theta;
}

// [[Rcpp::export]]
arma::mat update_phi_MGP_cpp(arma::mat Theta, arma::vec tau, int K, int p,
                             arma::mat v1, arma::mat v2) {
  arma::mat phi(p, K);
  for (int j = 0; j < p; j++) {
    // Update phi, given prior phi_jh ~ Ga(v1/2, v2/2)
    for (int h = 0; h < K; h++) {
      double shape = (v1(j, h) + 1) * 0.5;
      double rate = (v2(j, h) + tau(h) * pow(Theta(j, h), 2)) * 0.5;
      phi(j, h) = arma::randg(distr_param(shape, 1/rate));
    }
  }
  
  return phi;
}

// [[Rcpp::export]]
arma::vec update_delta_MGP_cpp(arma::vec delta, arma::vec tau, arma::mat Theta, 
                               arma::mat phi, int K, int p) {
  int a1 = 2; // hyperparams for prior as recommended by Dunson 2011
  int a2 = 3;
  arma::mat pT2 = phi % square(Theta);
  double shape = a1 + p * K * 0.5;
  double rate = 1 + 0.5 * sum(tau.t() % sum(pT2, 0)) / delta(0);
  delta(0) = arma::randg( distr_param(shape, 1/rate));
  tau = arma::cumprod(delta);
  
  for (int h = 1; h < K; h++) {
    shape = a2 + p * (K - h) * 0.5;
    arma::vec tau_hk = tau.subvec(h, K-1);
    rate = 1 + 0.5 * sum(tau_hk.t() % sum(pT2.cols(h, K-1), 0)) / delta(h);
    delta(h) = arma::randg( distr_param(shape, 1/rate));
    tau = arma::cumprod(delta);
    arma::vec r = zeros<vec>(K);
  }
  
  return delta;
}
