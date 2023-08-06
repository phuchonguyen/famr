#include <RcppArmadillo.h>
#include <math.h>
#include "helpers.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;

// U: L x K x T cube
// Y: nT x q matrix
// eta: nT x K matrix
// Lambda: q x L matrix
// time: nT vector of index of time
// [[Rcpp::export]]
void update_Lambda_cpp(arma::mat& Lambda, const arma::cube& U, 
                       const arma::mat& Y, const arma::mat& eta,
                       const arma::mat& S_inv, const arma::mat& phi, const arma::vec& tau, 
                       const arma::vec time, int n, int q, int K, int L) {
  arma::vec M(q);
  arma::mat C(q, q);
  arma::mat X(L, n);
  arma::rowvec m(q);
  for (int i = 0; i < n; i++) {
    X.col(i) = U.slice(time(i)) * eta.row(i).t();
  }
  for (int l = 0; l < L; l++) {
    arma::mat S0_inv = arma::diagmat(1 / phi.col(l)) / tau(l);
    double s = accu(pow(X.row(l), 2));
    // multiply each row of Y, or Y(i) with X(l, i)
    // return sum of elements in each column of Y
    m = sum(Y.each_col() % X.row(l).t(), 0);
    C = (S_inv * s + diagmat(S0_inv.col(l))).i();
    M = C * S_inv * m.t();
    Lambda.col(l) = arma::mvnrnd(M, C);
  }
}

// U: L x K x T cube
// Y: nT x q
// eta: n x K
// Lambda: q x L
// id: nT vector
// time: nT vector
// S_inv: q x q positive definite precision matrix
// C_inv: T x T kernel matrix
// n: number of unique subjects
// q: number of outcomes
// K: number of latent factor of exposures
// L: number of latent factor of GPs
// T: number of unique time points
// [[Rcpp::export]]
void update_U_cpp(arma::cube& U, const arma::mat& Y, const arma::mat& eta, 
                  const arma::mat& Lambda, const arma::vec& id,
                  const arma::vec& time, const arma::mat& S_inv,
                  const arma::mat& C_inv, int L, int K, int T, int q, int n) {
  arma::mat y(T*q, n); // vectorized y minus expected y as function of Lambda, U, eta
  arma::mat o(T, n, fill::zeros);
  for (int i = 0; i < n; i++) {
    arma::uvec id_i = find(id == i);
    // fill y in with zero at unobserved time points
    arma::mat Y_i(q, T, fill::zeros);
    for (int ii = 0; ii < id_i.n_elem; ii++) {
      // create an indicator vector of time indices when outcomes are observed for subject i
      o.col(i)(time(id_i(ii))) = 1; // return time index of observation i
      Y_i.col(time(id_i(ii))) = Y.row(id_i(ii)).t();
    }
    // create vectorize y
    y.col(i) = vectorise(Y_i);
  }
  
  arma::field<arma::cube> D(n, K); // Tq x T x L x n x K array
  arma::field<arma::mat> V(n, K); // Tq x L x n x K array
  arma::cube DD(T*q, T, L, fill::zeros); // an element in D
  arma::mat VV(T*q, L); // an element in V
  arma::vec u_lk(T);
  // Create D_i^(lk) matrix for each subject
  for (int k = 0; k < K; k++) {
    for (int l = 0; l < L; l++) {
      // create a diagonal matrix where each diagonal element is Lambda.col(l)
      DD.slice(l) = kron(eye(T, T), Lambda.col(l)); // Tq x T matrix
    }
    for (int i = 0; i < n; i++) {
      // arma::vec o_i(T, fill::zeros);
      // arma::uvec id_i = find(id == i);
      // // fill y in with zero at unobserved time points
      // arma::mat Y1_i(q, T, fill::zeros);
      // for (int ii = 0; ii < id_i.n_elem; ii++) {
      //   // create an indicator vector of time indices when outcomes are observed for subject i
      //   o_i(time(id_i(ii))) = 1; // return time index of observation i
      //   Y1_i.col(time(id_i(ii))) = Y.row(id_i(ii)).t();
      // }
      // // create vectorize y
      // y1.col(i) = vectorise(Y1_i);
      // cout << "\nO_i: " << all(o_i == o.col(i)) << " k=" << k << " i=" << i << endl;
      // cout << "\nvec(Y0_i)==vec(Y1_i): " << all(y.col(i) == y1.col(i)) << endl;
      // if (!all(y.col(i) == y1.col(i))) {
      //   cout << "\nvec(Y0_i): " << y.col(i) << endl;
      //   cout << "\nvec(Y1_i): " << y1.col(i) << endl;
      // }
      
      // creates the D_i^(lk) matrix as in paper
      D(i, k) = DD.each_slice() * diagmat(o.col(i)) * eta(i, k);
      
      // create elements that sum to the regression contribution of factor k
      for (int l = 0; l < L; l++) {
        u_lk = U.tube(l, k);
        VV.col(l) = D(i, k).slice(l) * u_lk;
      }
      V(i, k) = VV;
      // subject regression contribution of factor k
      y.col(i) -= sum(VV, 1);
    }
  }
  
  arma::mat tmp;
  arma::mat DSD(T, T);
  arma::mat DSD_inv(T, T);
  arma::vec DSy(T);
  for (int k = 0; k < K; k++) {
    for (int l = 0; l < L; l++) {
      // initialize posterior covariance and mean for u_lk
      DSD = C_inv;
      DSy = zeros(T);
      for (int i = 0; i < n; i++) {
        tmp = D(i, k).slice(l).t() * S_inv;
        // add to posterior covariance
        DSD += tmp * D(i, k).slice(l);
        // y is residual after regressing on K factors, 
        // adding removes contribution of u_lk regressor
        DSy += tmp * (y.col(i) + V(i, k).col(l));
      }
      DSD_inv = DSD.i();
      U.tube(l, k) = arma::mvnrnd(DSD_inv * DSy, DSD_inv);
    }
  }
}

// Use 1-D Exponential quadratic kernel to make a covariance matrix
// [[Rcpp::export]]
arma::mat covEQ(arma::vec t, double kappa, double amplitude) {
  int D = t.n_elem;
  arma::mat C(D, D);
  for (int i = 0; i < D; i++) {
    for (int j = 0; j < D; j++) {
      C(i, j) = exp(-0.5 * pow( abs(t(i) - t(j)), 1.9999) / pow(kappa, 2)); // for stability
      C(j, i) = C(i, j);
      if (i == j) {C(i, j) = 1.0;}
    }
  }
  return amplitude * C;
}

// [[Rcpp::export]]
double llike_kappa(const arma::mat& Ci, double logdetC, const arma::cube& U,
                   int L, int K) {
  double s = 0;
  for (int l = 0; l < L; l++) {
    for (int k = 0; k < K; k++) {
      arma::vec u_lk = U.tube(l, k);
      s += as_scalar(u_lk.t() * Ci * u_lk);
    }
  }
  double l = - 0.5 * K * L * logdetC - 0.5 * s;
  return l;
}


// time: all unique time points
// [[Rcpp::export]]
Rcpp::List update_kappa_cpp(double kappa , arma::mat& Ci, arma::mat& C, 
                      double logdetC, double lpdf, arma::vec time,
                      const arma::cube& U, double a, double b, int L, int K,
                      double eps, int s, int batch_size, int n_accepted) {
  
  // update scaling parameter every batch_size iterations
  if ((s - 1) % batch_size == 0) { // minus 1 because R index starts at 1
    double d = min(0.05, std::pow((int)s/batch_size, -0.5));
    // 0.44 is theoretically optimal acceptance rate for the latest 50 interations
    if ((n_accepted / batch_size) > 0.44) {
      eps += d;
    } else {
      eps -= d;
    }
    // reset acceptances counter
    n_accepted = 0;
  }
  // sample proposal on log scale gives better mixing and makes prop >= 0
  double prop = exp( log(kappa) + std::exp(2*eps) * randn<double>());
  arma::mat C_prop = covEQ(time, prop, 1);
  arma::mat Ci_prop = C_prop.i();
  double logdetC_prop = log(det(C_prop));
  double l_prop = llike_kappa(Ci_prop, logdetC_prop, U, L, K);
  l_prop += ldinvgam(prop, a, b);
  // log acceptance probability
  double logr = l_prop - lpdf;
  if (log(randu<double>()) < logr) {
    kappa = prop;
    Ci = Ci_prop;
    logdetC = logdetC_prop;
    C = C_prop;
    n_accepted += 1;
    lpdf = l_prop;
  }
  
  return Rcpp::List::create(Rcpp::Named("kappa") = kappa,
                            Rcpp::Named("C") = C,
                            Rcpp::Named("C_inv") = Ci,
                            Rcpp::Named("logdetC") = logdetC,
                            Rcpp::Named("kappa_lpdf") = lpdf,
                            Rcpp::Named("kappa_eps") = eps,
                            Rcpp::Named("kappa_n_accepted") = n_accepted
                            );
}

