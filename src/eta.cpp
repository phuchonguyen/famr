#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;

// a (t, K) matrix B^T + \sum_l z_{il}C^T_l
void get_Aij(arma::mat& Aij, const arma::mat& B, int K, int t, int q_int = 0, 
             const arma::rowvec& z_int = zeros<arma::rowvec>(0)) {
  Aij = B.submat(0, 0, K - 1, t - 1);
  int s;
  int e;
  for (int l = 1; l < q_int; l++) {
    s = K * l;
    e = K * l + K - 1;
    Aij += B.submat(s, 0, e, t - 1) * z_int(l - 1);
  }
}

void get_Si_mi(arma::mat& Si, arma::vec& mi,
               const arma::mat& B, const arma::mat& Theta,
               const arma::mat& Sigmay_inv, const arma::vec& sigmax_sqinv,
               const arma::mat& Y, const arma::rowvec& xi, 
               const arma::mat& Z, const arma::mat& Z_int,
               const arma::uvec& ids, int t, int K, int q, int q_int) {
  arma::mat Aij(K, t);
  arma::rowvec yij(t);
  arma::mat L(K, K);
  arma::mat Theta_sigmax = Theta.each_col() % sigmax_sqinv;
  Si.zeros();
  mi.zeros();
  // loop through repeated outcome measurements for subject i
  for (int j = 0; j < (int) ids.n_elem; j++) {
    // coefficient matrix when factor interacts with covariates
    if (q_int > 0) {
      get_Aij(Aij, B, K, t, q_int, Z_int.row(ids(j)));
    } else {
      get_Aij(Aij, B, K, t);
    }
    yij = Y.row(ids(j));
    // outcome minus effect of covariates
    if (q > 0) {
      yij -= Z.row(ids(j)) * B.submat(B.n_rows - q, 0, B.n_rows - 1, t - 1);
    }
    // contribution of Y to likelihood
    Si += Aij * Sigmay_inv * Aij.t();
    mi += Aij * Sigmay_inv * yij.t();
  }
  // contribution of X to likelihood
  Si += Theta_sigmax.t() * Theta;
  mi += Theta_sigmax.t() * xi.t();
}

// Metropolis Hasting sampler for eta
// Can include interaction between factors

// log prior of latent factors for one subject i
// eta_i \sim N(0, I_K)
double lprior(const arma::rowvec etai) {
  return as_scalar(-0.5*etai*etai.t());
}

// log likelihood of latent factors for one subject i
double llike(const arma::rowvec& etai,
             const arma::mat& B, const arma::mat& Theta,
             const arma::mat& Sigmay_inv, const arma::vec& sigmax_sqinv,
             const arma::mat& Y, const arma::rowvec& xi, 
             const arma::mat& Z, const arma::mat& Z_int,
             const arma::uvec& ids, int t, int K, int q, int q_int) {
  arma::mat Si(K, K);
  arma::vec mi(K);
  get_Si_mi(Si, mi, B, Theta, Sigmay_inv, sigmax_sqinv,
            Y, xi, Z, Z_int, ids, t, K, q, q_int);
  return as_scalar(-0.5*etai*Si*etai.t() + etai*mi);
}

// [[Rcpp::export]]
void update_eta_mh_cpp(arma::mat& eta, const arma::mat& B, const arma::mat& Theta, 
                       const arma::mat& Sigmay_inv, const arma::vec& sigmax_sqinv,
                       const arma::mat& Y, const arma::mat& X, const arma::vec& uid, 
                       const arma::mat& Z, const arma::mat& Z_int,
                       int K, int p, int t, int n, int q, int q_int, 
                       Rcpp::IntegerVector n_accepted, double eps = 1e-3) {
  arma::rowvec prop(K);
  double logr;
  double llike_prop;
  double llike_cur;
  arma::uvec ids;
  for (int i = 0; i < n; i++) {
    // ids of all observations of subject i
    ids = find(uid == i);
    // random walk proposal with scaling parameter eps
    prop = eps * randn<rowvec>(K);
    // log likelihood of the proposal
    llike_prop = llike(prop, B, Theta, Sigmay_inv, sigmax_sqinv,
                       Y, X.row(i), Z, Z_int, ids, t, K, q, q_int);
    // log likelihood of the current state
    llike_cur = llike(eta.row(i), B, Theta, Sigmay_inv, sigmax_sqinv,
                      Y, X.row(i), Z, Z_int, ids, t, K, q, q_int);
    // log acceptance probability
    logr = llike_prop + lprior(prop) - llike_cur - lprior(eta.row(i));
    if (log(randu<double>()) < logr) {
      eta.row(i) = prop;
      n_accepted = n_accepted + 1;
    }
  }
}

// Gibbs sampler for eta
// Assuming no interaction between factors

// TODO TEST!!!
// [[Rcpp::export]]
void update_eta_gibbs_cpp(arma::mat& eta, const arma::mat& B, const arma::mat& Theta, 
                          const arma::mat& Sigmay_inv, const arma::vec& sigmax_sqinv,
                          const arma::mat& Y, const arma::mat& X, const arma::vec& uid, 
                          const arma::mat& Z, const arma::mat& Z_int,
                          int K, int p, int t,
                          int n, int q, int q_int) {
  arma::mat Si(K, K);
  arma::vec mi(K);
  arma::uvec ids;
  arma::mat L(K, K);
  for (int i = 0; i < n; i++) {
    // ids of all observations of subject i
    ids = find(uid == i);
    // get posterior mean and variance
    get_Si_mi(Si, mi, B, Theta, Sigmay_inv, sigmax_sqinv,
              Y, X.row(i), Z, Z_int, ids, t, K, q, q_int);
    // add contribution from the prior
    Si += eye(K, K);
    // use Cholesky decomposition for efficiency
    L = trimatl(chol(Si, "lower"));
    eta.row(i) = solve(L, randn<vec>(K)).t() + solve(trimatu(L.t()), solve(L, mi)).t();
  }
}


// posterior of eta_i given x_i
arma::mat predict_eta(const arma::mat& X, const arma::mat& Theta, 
                      const arma::vec& sigmax_sqinv, int n, int K) {
  arma::mat eta(n, K);
  arma::mat U = Theta.each_col() % sigmax_sqinv;
  arma::mat L = trimatl(arma::chol(U.t() * Theta + arma::eye(K, K), "lower"));
  arma::mat V = arma::solve(trimatu(L.t()), arma::solve(L, arma::eye(K, K)));
  for (int i = 0; i < n; i++) {
    eta.row(i) = X.row(i)*U*V;
  }
  return eta;
}
