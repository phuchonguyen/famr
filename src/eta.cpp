#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;

// a (t, K) matrix B^T + \sum_l z_{il}C^T_l
void get_Aij(arma::mat& Aij, const arma::cube Bt, const arma::mat& B, 
             int t, int K, int q, int p_int = 0, 
             const arma::rowvec& z_int = zeros<arma::rowvec>(0)) {
  Aij = Bt.slice(t);
  int s;
  int e;
  for (int l = 0; l < p_int; l++) {
    s = K * l;
    e = K * (l + 1) - 1;
    Aij += B.submat(s, 0, e, q - 1) * z_int(l);
  }
}

void get_Si_mi(arma::mat& Si, arma::mat& Ri, arma::vec& mi, arma::vec& ti,
               const arma::cube Bt, const arma::mat& B, const arma::mat& Theta,
               const arma::mat& Sigmay_inv, const arma::vec& sigmax_sqinv,
               const arma::mat& Y, const arma::rowvec& xi, 
               const arma::mat& Z_int, const arma::uvec& idi, const arma::uvec& time,
               int q, int K, int p_int) {
  arma::mat Aij(K, q);
  arma::rowvec yij(q);
  arma::mat Theta_sigmax = Theta.each_col() % sigmax_sqinv;
  // loop through repeated outcome measurements for subject i
  for (int j = 0; j < (int) idi.n_elem; j++) {
    // coefficient matrix when factor interacts with covariates
    if (p_int > 0) {
      get_Aij(Aij, Bt, B, time(idi(j)), K, q, p_int, Z_int.row(idi(j)));
    } else {
      get_Aij(Aij, Bt, B, time(idi(j)), K, q);
    }
    yij = Y.row(idi(j));
    // contribution of Y to likelihood
    Si += Aij * Sigmay_inv * Aij.t();
    Ri += Aij.t();
    mi += Aij * Sigmay_inv * yij.t();
    ti += yij.t();
  }
  // contribution of X to likelihood
  Si += Theta_sigmax.t() * Theta;
  mi += Theta_sigmax.t() * xi.t();
}

// Metropolis Hasting sampler for eta
// Can include interaction between factors

// log prior of latent factors for one subject i
// eta_i \sim N(0, I_K)
double lprior(const arma::rowvec etai, int K) {
  return as_scalar(-0.5*log(2*M_PI)*K -0.5*log(K) -0.5*etai*etai.t());
}

// (q x K x K) Omega mode-2 eta_i^T mode-3 eta_i^T
arma::vec get_Omega_tilde(const arma::cube& Omega, const arma::rowvec& etai, int q) {
  arma::vec o(q);
  for (int i = 0; i < q; i++) {
    o(i) = as_scalar(etai * Omega.slice(i) * etai.t());
  }
  return(o);
}

// log likelihood of latent factors for one subject i
double llike(const arma::rowvec& etai, const arma::cube Bt,
             const arma::mat& B, const arma::mat& Theta, const arma::cube& Omega,
             const arma::mat& Sigmay_inv, const arma::vec& sigmax_sqinv,
             const arma::mat& Y, const arma::rowvec& xi, 
             const arma::mat& Z_int, const arma::uvec& idi, const arma::uvec& time,
             int q, int K, int p_int, int p) {
  arma::mat Si(K, K, fill::zeros);
  arma::mat Ri(q, K, fill::zeros);
  arma::vec mi(K, fill::zeros);
  arma::vec ti(q, fill::zeros);
  get_Si_mi(Si, Ri, mi, ti, Bt, B, Theta, Sigmay_inv, sigmax_sqinv,
            Y, xi, Z_int, idi, time, q, K, p_int);
  arma::vec Oi = get_Omega_tilde(Omega, etai, q);
  arma::rowvec OS = Oi.t()*Sigmay_inv;
  return as_scalar(-0.5*log(2*M_PI)*(q*idi.n_elem + p) + 
                   0.5*log(sum(sigmax_sqinv)) + 
                   idi.n_elem*0.5*log(det(Sigmay_inv)) +
                   -0.5*etai*Si*etai.t() + etai*mi +
                   OS*ti - OS*Ri*etai.t() - 0.5*OS*Oi);
}

// Bt: K x q x T array
// [[Rcpp::export]]
void update_eta_mh_cpp(arma::mat& eta, const arma::cube Bt,
                       const arma::mat& B, const arma::mat& Theta,
                       const arma::cube& Omega,
                       const arma::mat& Sigmay_inv, const arma::vec& sigmax_sqinv,
                       const arma::mat& Y, const arma::mat& X,
                       const arma::mat& Z_int,
                       const arma::vec& uid, const arma::vec& id,
                       const arma::uvec& time,
                       int K, int p, int q, int n, int p_int,
                       arma::vec& n_accepted, arma::vec& eps,
                       arma::cube& A, arma::mat& b, 
                       arma::vec& lpmf, int s,
                       arma::mat& eta_prop,
                       bool adaptiveM = true, bool adaptiveMWG = false,
                       int batch_size = 50, double eps_power = -0.5
                       ) {
  arma::rowvec prop(K);
  double logr;
  double llike_prop;
  double llike_cur;
  arma::uvec idi;
  arma::vec M(K, fill::zeros); // for Adaptive Metropolis
  arma::mat C(K, K);
  for (int i = 0; i < n; i++) {
    // indices of all observations of subject i
    idi = find(id == uid(i));
    // Adaptive Metropolis scaling proposal
    if ((s > 1) & adaptiveM) {
      A.slice(i) += eta.row(i).t() * eta.row(i);
      b.col(i) += eta.row(i).t();
      C = (std::pow(2.38, 2) / K) * 
        (A.slice(i) - (b.col(i) * b.col(i).t()) / s) / (s-1) +
        (1e-10 * eye(K, K));
      prop = eta.row(i) + arma::mvnrnd(M, C).t();
    } else if (adaptiveMWG) {
      // update scaling parameter every 50 iterations
      if ((s-1) % batch_size == 0) {
        double d = min(0.05, std::pow((int)s/batch_size, eps_power));
        // 0.44 is theoretically optimal acceptance rate for the latest 50 interations
        if ((n_accepted(i) / batch_size) > 0.44) {
          eps(i) += d;
        } else if ((n_accepted(i) / batch_size) < 0.23){
          eps(i) -= d;
        }
        // reset acceptances counter
        n_accepted(i) = 0;
      }
      prop = eta.row(i) + std::exp(2*eps(i)) * randn<rowvec>(K);
    } else {
      // Random walk proposal with scaling parameter eps
      prop = eta.row(i) + std::exp(2*eps(i)) * randn<rowvec>(K);
    }
    // for DEBUGGING
    eta_prop(s-1,i) = prop(0);
    // log likelihood of the proposal
    llike_prop = llike(prop, Bt, B, Theta, Omega, Sigmay_inv, sigmax_sqinv,
                       Y, X.row(i), Z_int, idi, time, q, K, p_int, p);
    // log likelihood of the current state
    // llike_cur = llike(eta.row(i), B, Theta, Omega, Sigmay_inv, sigmax_sqinv,
    //                   Y, X.row(i), Z_int, idi, time, q, K, p_int, p);
    llike_cur = lpmf(i);
    // log acceptance probability
    logr = llike_prop + lprior(prop, K) - llike_cur - lprior(eta.row(i), K);
    if (log(randu<double>()) < logr) {
      eta.row(i) = prop;
      n_accepted(i) += 1;
      lpmf(i) = llike_prop;
    }
  }
}