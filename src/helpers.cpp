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