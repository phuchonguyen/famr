#include <RcppArmadillo.h>
#include <math.h>
#include "helpers.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;


// For MGP Gibbs sampler
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
arma::mat update_phi_MGP_cpp(arma::mat Theta, arma::vec tau, int K, int p) {
  arma::mat phi(p, K);
  double v = 3.0; // from paper
  for (int j = 0; j < p; j++) {
    // Update phi, given prior phi_jh ~ Ga(v/2, v/2)
    for (int h = 0; h < K; h++) {
      double shape = (v + 1) * 0.5;
      double rate = (v + tau(h) * pow(Theta(j, h), 2)) * 0.5;
      phi(j, h) = arma::randg(distr_param(shape, 1/rate));
    }
  }
  
  return phi;
}

// [[Rcpp::export]]
arma::vec update_delta_MGP_cpp(arma::vec delta, arma::vec tau, arma::mat Theta, 
                               arma::mat phi, int K, int p, 
                               double a1=2.0, double a2=3.0) {
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

// For DL priors on the rows of Theta
// [[Rcpp::export]]
void update_tau_DL_cpp(arma::vec& tau, const arma::mat& Theta, 
                       const arma::mat& phi, const arma::mat& omega,
                       int p, int K) {
  arma::vec a = sum(abs(Theta) / phi, 1);
  for (int j = 0; j < p; j++) {
    tau(j) = std::max(rgig_cpp(1 - K, 1, 2 * a(j)), datum::eps); 
  }
}

// [[Rcpp::export]]
void update_phi_DL_cpp(arma::mat& phi, const arma::mat& Theta, int p, int K) {
  double a = 0.5; // prior on lambda parameter
  arma::mat A = abs(Theta);
  arma::rowvec t(K);
  for (int j = 0; j < p; j++) {
    for (int k = 0; k < K; k++) {
      t(k) = std::max(rgig_cpp(a - 1, 1, 2 * A(j, k)), datum::eps); // prevent numerical rounding to zero
    }
    phi.row(j) = t / accu(t); // division by zero here...
  }
}

// [[Rcpp::export]]
void update_omega_DL_cpp(arma::mat& omega, const arma::mat& phi, 
                         const arma::mat& Theta, const arma::vec& tau) {
  omega = (phi.each_col() % tau) / abs(Theta);
  omega.transform(rig_cpp);
}

// Integrate out eta and use Metropolis-Hasting to sample Theta with Gaussian priors
// DONT USE

double llike_T(const arma::mat& T, const arma::vec& sigmax_sqinv,
               const arma::mat& X, int p, int n) {
  double l = 0.0;
  arma::mat V = inv_sympd(T*T.t() + arma::diagmat(sigmax_sqinv));
  for (int i = 0; i < n; i++) {
    l += as_scalar(X.row(i) * V * X.row(i).t());
  }
  return -0.5 * l;
}

double lprior_T(const arma::mat& T, const arma::vec& V0) {
  arma::vec vT = arma::vectorise(T); // stack columns of T
  return -0.5 * arma::accu(pow(vT, 2) / V0); // V0 is indep. prior variances of vec(Theta)
}

// [[Rcpp::export]]
void update_Theta_mh_cpp(arma::mat& Theta, const arma::vec& sigmax_sqinv,
                         const arma::mat& X, const arma::vec& V0, 
                         arma::mat& n_accepted, arma::mat& eps,
                         int s, bool adaptiveMWG = true) {
  int p = X.n_cols;
  int n = X.n_rows;
  int K = Theta.n_cols;
  double lp_cur;
  double lp_prop;
  double prop;
  double cur;
  // update each element of Theta
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < K; j++) {
      // get proposal
      if (adaptiveMWG) {
        if (s % 50 == 0) { // update scaling parameter every 50 iterations
          double d = min(0.01, std::pow((int)s/50, -0.5));
          if (n_accepted(i,j) / s > 0.44) {eps(i,j) += d;} // 0.44 is theoretically optimal acceptance rate
          else {eps(i,j) -= d;}
        }
        prop = Theta(i,j) + std::exp(2*eps(i,j)) * randn<double>();
      } else {
        // proposal from normal with scaling parameter eps
        prop = Theta(i,j) + eps(i,j) * arma::randn<double>();
      }
      cur = Theta(i,j);
      // log likelihood
      lp_cur = llike_T(Theta, sigmax_sqinv, X, p, n) + lprior_T(Theta, V0);
      Theta(i,j) = prop;
      lp_prop = llike_T(Theta, sigmax_sqinv, X, p, n) + lprior_T(Theta, V0);
      // accept with log probability
      if (log(randu<double>()) < (lp_prop-lp_cur)) {
        n_accepted(i,j)++;
      } else {
        Theta(i,j) = cur;
      }
    }
  }
}
