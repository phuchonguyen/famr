#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;

// Assuming no interaction between factors

arma::mat get_Aij(arma::mat B, int K, int t, int q_int, arma::rowvec z_int = arma::zeros(0)) {
  arma::mat Aij = B.submat(0, 0, K - 1, t - 1);
  int s;
  int e;
  for (int l = 1; l < q_int; l++) {
    s = K * l;
    e = K * l + K - 1;
    Aij += B.submat(s, 0, e, t - 1) * z_int(l);
  }
  return Aij.t();
}

// [[Rcpp::export]]
arma::mat update_eta_gibbs_cpp(arma::mat eta, arma::mat B, arma::mat Theta, arma::mat Sigmay_inv, arma::vec sigmax_sqinv,
                               arma::mat Y, arma::mat X, int K, int p, int t, int n, int q, int q_int, 
                               arma::vec uid, arma::mat Z, arma::mat Z_int) {
  arma::mat Si(K, K);
  arma::vec mi(K);
  arma::mat Aij(t, K);
  arma::rowvec yij(t);
  arma::mat L(K, K);
  arma::mat Thetat_sigmax = Theta.t() * arma::diagmat(sigmax_sqinv);
  for (int i = 0; i < n; i++) {
    Si.zeros();
    mi.zeros();
    arma::uvec ids = find(uid == i);
    for (int j = 0; j < (int) ids.n_elem; j++) {
      Aij = get_Aij(B, K, t, q_int, Z_int.row(ids(j)));
      yij = Y.row(ids(j));
      if (q > 0) {
        yij -= Z.row(ids(j)) * B.submat(B.n_rows - q, 0, B.n_rows - 1, t - 1);
      }
      Si += Aij.t() * Sigmay_inv * Aij;
      mi += Aij.t() * Sigmay_inv * yij.t();
    }
    Si += arma::eye(K, K) + Thetat_sigmax * Theta;
    L = trimatl(arma::chol(Si, "lower"));
    mi += Thetat_sigmax * X.row(i).t();
    eta.row(i) = arma::solve(L, randn<arma::vec>(K)).t() + arma::solve(trimatu(L.t()), arma::solve(L, mi)).t();
  }
  
  return eta;
}