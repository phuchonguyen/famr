#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace arma;

// Like the indexing notation in base R
// [[Rcpp::export]]
arma::mat submat_cpp(arma::mat M, arma::uvec idx, arma::uvec idy) {
  int nx = (int) idx.n_elem, ny = (int) idy.n_elem;
  arma::mat G(nx, ny);
  for (int i = 0; i < ny; i++) {
    for (int j = 0; j < nx; j++) {
      G(j, i) = M(idx(j), idy(i));
    }
  }
  return G;
}

// [[Rcpp::export]]
arma::mat impute_X_lod_cpp(arma::mat eta, arma::mat Theta, arma::mat X, arma::vec sigmax_sqinv,
                           int n, arma::mat Ilod, arma::vec llod, arma::vec ulod) {
  arma::uvec pos;
  arma::rowvec m;
  double ui, ai, bi;
  for (int i = 0; i < n; i++) {
    pos = find(Ilod.row(i) == 1); 
    if (!pos.is_empty()) {
      m = eta.row(i) * Theta.t();
      for (int j = 0; j < (int) pos.n_elem; j++) {
        ai = normcdf(llod(pos(j)), m(pos(j)), sqrt(1/sigmax_sqinv(pos(j))));
        bi = normcdf(ulod(pos(j)), m(pos(j)), sqrt(1/sigmax_sqinv(pos(j))));
        ui = randu<double>() * (bi - ai) + ai;
        X(i, pos(j)) = R::qnorm(ui, m(pos(j)), sqrt(1/sigmax_sqinv(pos(j))), true, false);
      }
    }
  }
  
  return X;
}

// [[Rcpp::export]]
arma::mat impute_Ymis_cpp(arma::mat Y, arma::mat M, arma::mat Sigma, 
                          arma::mat O, int n, int q) {
  arma::uvec mis_id = find(sum(O, 1) > 0);
  for (int j = 0; j < (int) mis_id.n_elem; j++) {
    int i = mis_id(j);
    arma::uvec mis = find(O.row(i) == 1);
    arma::uvec obs = find(O.row(i) == 0);
    arma::rowvec mi = M.row(i);
    if ((int) mis.n_elem == q) { // all values missing 
      Y.row(i) = arma::mvnrnd(mi.t(), Sigma).t();
    } else if ((int) mis.n_elem > 0) {
      arma::rowvec yi = Y.row(i);
      arma::mat Saa_inv = arma::inv(submat_cpp(Sigma, obs, obs));
      arma::mat Sba = submat_cpp(Sigma, mis, obs);
      arma::mat Sigma_b = submat_cpp(Sigma, mis, mis) - Sba * Saa_inv * Sba.t();
      arma::rowvec m_b = mi.elem(mis) + (yi.elem(obs) - mi.elem(obs)).t() * Saa_inv * Sba.t();
      if (Sigma_b.n_rows == 1) {
        yi(mis(0)) = randn<double>() * sqrt(Sigma_b.eval()(0, 0)) + m_b.eval()(0, 0);
      } else {
        yi.elem(mis) = arma::mvnrnd(m_b.t(), Sigma_b).t();
      }
      Y.row(i) = yi;
    }
  }
  
  return Y;
}

// [[Rcpp::export]]
arma::mat impute_Yprobit_cpp(arma::mat Y, arma::mat M, arma::mat Sigma, 
                             arma::mat Yraw, arma::vec binary, int n, int t) {
  arma::uvec mis = find(binary == 1); 
  arma::uvec ids = regspace<uvec>(0, t - 1);
  arma::uvec misj(1);
  for (int j = 0; j < sum(binary); j++) {
    // Conditional mean and covariance
    arma::uvec obs = find(ids != mis(j));
    misj.fill(mis(j));
    arma::mat Saa_inv = arma::inv(submat_cpp(Sigma, obs, obs));
    arma::mat Sba = submat_cpp(Sigma, misj, obs);
    double sigma_mis = as_scalar(submat_cpp(Sigma, misj, misj) - Sba * Saa_inv * Sba.t());
    arma::vec M_mis = M.col(mis(j)) + (Y.cols(obs) - M.cols(obs)) * Saa_inv * Sba.t();
    for (int i = 0; i < n; i++) {
      // Truncated normal
      double l = (Yraw(i, mis(j)) == 1) ? 0.5 : 0;
      double u = (Yraw(i, mis(j)) == 1) ? 1 : 0.5;
      double p = randu<double>() * (u - l) + l;
      Y(i, mis(j)) = R::qnorm(p, M_mis(i), sqrt(sigma_mis), true, false);
    }
  }
  
  return Y;
}