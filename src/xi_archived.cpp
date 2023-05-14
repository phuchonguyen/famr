#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;

// Metropolis Hasting sampler for xi

// log prior of latent factors for one subject i
// xi_i \sim N(0, I_H)
double lp_xi(const arma::rowvec& xii) {
  return as_scalar(-0.5*xii*xii.t());
}

// log likelihood of latent factors for one subject i
// Yi : (T_i \times q) matrix of repeated observations of outcome of subject i
// xii : H-vector of latent variable
// W : (q \times H) factor loading matrix for xi
// Sigma_inv : (q \times q) outcome idiosyncratic noise inversed covariance matrix
double ll_xi(const arma::mat& Yi, const arma::rowvec& xii, 
             const arma::mat& W, const arma::mat& Sigma_inv) {
  int Ti = Yi.n_rows; // number of repeated observations of subject i
  arma::rowvec si = arma::sum(Yi, 0); // sum_t Y_it
  return as_scalar(-0.5*Ti*xii*W.t()*Sigma_inv*W*xii.t() + xii*W.t()*Sigma_inv*si.t());
}

// [[Rcpp::export]]
void update_xi_mh_cpp(arma::mat& xi, 
                      arma::cube& A, arma::mat& b,
                      arma::vec& n_accepted, arma::vec& eps, arma::vec& lpmf,
                      const arma::mat& Y, const arma::mat& W, const arma::mat& Sigma_inv,
                      const arma::vec& uid, const arma::vec& id,
                      int H, int n, int s,
                      bool adaptiveM = true, bool adaptiveMWG = false,
                      int batch_size=50
                       ) {
  arma::rowvec prop(H);
  double logr;
  double llike_prop;
  double llike_cur;
  arma::uvec idi;
  arma::vec M(H, fill::zeros); // for Adaptive Metropolis
  arma::mat C(H, H);
  for (int i = 0; i < n; i++) {
    // Adaptive Metropolis scaling proposal
    if ((s > 1) & adaptiveM) {
      A.slice(i) += xi.row(i).t() * xi.row(i);
      b.col(i) += xi.row(i).t();
      C = (std::pow(2.38, 2) / H) * 
        (A.slice(i) - (b.col(i) * b.col(i).t()) / s) / (s-1) +
        (1e-10 * eye(H, H));
      prop = xi.row(i) + arma::mvnrnd(M, C).t();
    } else if (adaptiveMWG) {
      if ((s-1) % batch_size == 0) { // update scaling parameter every 50 iterations
        double d = min(0.05, std::pow((int)s/batch_size, -0.5));
        // 0.44 is theoretically optimal acceptance rate
        if ((n_accepted(i) / batch_size) > 0.44) {eps(i) += d;} 
        else {eps(i) -= d;}
        // reset acceptance counter
        n_accepted(i) = 0;
      }
      prop = xi.row(i) + std::exp(2*eps(i)) * randn<rowvec>(H);
    } else {
      // Random walk proposal with scaling parameter eps
      prop = xi.row(i) + std::exp(2*eps(i)) * randn<rowvec>(H);
    }
    // indices of all observations of subject i
    arma::mat Yi = Y.rows(find(id == uid(i)));
    // log likelihood of the proposal
    llike_prop = ll_xi(Yi, prop, W, Sigma_inv);
    // log likelihood of the current state
    llike_cur = lpmf(i);
    // log acceptance probability
    logr = llike_prop + lp_xi(prop) - llike_cur - lp_xi(xi.row(i));
    if (log(randu<double>()) < logr) {
      xi.row(i) = prop;
      n_accepted(i) += 1;
      lpmf(i) = llike_prop;
    }
  }
}