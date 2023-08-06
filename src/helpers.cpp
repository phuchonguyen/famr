#include <RcppArmadillo.h>
#include <math.h>
#include "helpers.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;

// From infinitefactor
// [[Rcpp::export]]
double rgig_cpp(double lam, double psi, double chi){
  double omega = std::max(std::pow(psi*chi, 0.5), 1e-7); // so that omega*omega > 0
  double alpha = std::pow(chi/psi, 0.5);
  bool flip = false;
  long mode;
  if(lam<0) {lam *= -1; flip = true;}
  if(lam>=1) mode = (std::pow((lam-1)*(lam-1) + omega*omega, 0.5)+(lam-1))/omega;
  else mode = omega/(std::pow((lam-1)*(lam-1) + omega*omega, 0.5)+(lam-1));
  double lm = std::log(mode);
  double cons = 0.5 * (lam-1) * lm - 0.25*omega*(mode + 1/mode);
  long maxp = ((lam+1)+std::pow((lam+1)*(lam+1)+omega*omega,0.5))/omega;
  long eq = 0.5*(lam+1)*std::log(maxp) - 0.25*omega*(maxp + 1/maxp) - cons;
  double ext = std::exp(eq);
  float U, L, prop;
  do{
    U = ext * randu<float>();
    L = randu<float>();
    prop = U/L;
  } while (
      ((std::log(L)) > 
         (0.5*(lam-1)*std::log(prop) - 
         0.25*omega*(prop + 1/prop) - 
         cons)));
  if(flip) return alpha/prop;
  else return prop/alpha;
}

// [[Rcpp::export]]
double rig_cpp(double mu) {
  double y = randn<double>();
  y *= y;
  double mu2 = std::pow(mu, 2);
  double quad = 4 * mu * y + mu2 * std::pow(y, 2);
  double x = mu + y * mu2 / 2 - mu / 2  * std::pow(quad, 0.5);
  double u = randu<double>();
  if(u < (mu / (x + mu))) return x;
  else return mu2 / x;
}


// Evaluate the density of an inverse gamma at x
// a : shape, > 0
// b : scale, > 0
// log : return the log
// [[Rcpp::export]]
double ldinvgam(double x, double a, double b) {
  if (a < 0.0) {std::cerr << "\nshape a must be non-negative\n";}
  if (b < 0.0) {std::cerr << "\nscale b must be non-negative\n";}
  double l = -(a - 1) * log(x) - b / x;
  return l;
}

// // [[Rcpp::export]]
// arma::mat lpmf_cpp(const arma::mat& Y, const arma::mat& X, const arma::mat& Z, 
//                    const arma::mat& Z_int, int n,
//                    const arma::mat& B, const arma::mat& Theta, 
//                    const arma::mat& Omega, const arma::mat& eta_int, 
//                    const arma::mat& eta_quad,
//                    const arma::mat& Sigmainv, const arma::vec& sigmax_sqinv,
//                    bool include_interactions = false) {
//   arma::vec ll(n);
//   arma::vec my = eta_int * B;
//   if (include_interactions) {
//     my += eta_quad * 1;
//   }
//   arma::vec mx = eta_int *
//   for (int i=0; i<n; i++) {
//     my = 1;
//   }
// }

