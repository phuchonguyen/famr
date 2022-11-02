#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;

// From infinitefactor
// [[Rcpp::export]]
double rgig_cpp(double lam, double chi, double psi){
  double omega = std::pow(psi*chi, 0.5);
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
arma::vec update_zeta(arma::vec psi, int q, double global_shrink=1) {
  arma::vec zeta(q);
  double v = 0.5;
  for (int i = 0; i < q; i++) {
    double rate = psi(i) + global_shrink;
    zeta(i) = arma::randg<double>(arma::distr_param(v, 1/rate));
  }
  
  return zeta;
}

// q: number of predictors, p: number of outcomes; B: q x p matrix
// [[Rcpp::export]]
arma::vec update_psi(arma::mat B, arma::mat Sigmainv, arma::mat B0, 
                         arma::vec zeta, int q, int p) {
  arma::vec psi(q);
  double u = 0.5;
  double lam = u - 0.5 * p;
  for (int i = 0; i < q; i++) {
    double chi = arma::as_scalar((B.row(i)-B0.row(i))*Sigmainv*(B.row(i)-B0.row(i)).t());
    double rgig_psi = 2 * zeta(i); 
    psi(i) = rgig_cpp(lam, chi, rgig_psi);
  }
  return psi;
}

// Paper: https://arxiv.org/pdf/1711.07635.pdf
// B: q x p matrix
// TODO: change p to q, q to p to match with notation in Overleaf
// [[Rcpp::export]]
arma::mat update_B_TPBN(arma::mat X, arma::mat Y, arma::mat Sigma, 
                        arma::mat B0, arma::vec psi, int q, int p) {
  arma::mat Psii = arma::diagmat(1 / psi);
  arma::mat Lu = trimatl(arma::chol(X.t() * X + Psii, "lower"));
  arma::mat Rv = trimatu(arma::chol(Sigma, "upper"));
  arma::mat M = arma::solve(trimatu(Lu.t()), arma::solve(Lu, X.t()*Y + Psii*B0));
  arma::mat B = M + arma::solve(trimatu(Lu.t()), arma::randn<arma::mat>(q, p)) * Rv;
  
  return B;
}

// // Returns u0_square
// // [[Rcpp::export]]
// arma::vec update_u0_hCauchy(arma::vec b0, arma::vec c0, double tau0_sq, int q, int p) {
//   arma::vec u0_sq(q);
//   for (int i = 0; i < q; i++) {
//     double chi = b0(i) * b0(i) / tau0_sq;
//     double lam = 0;
//     double ps = 2 * c0(i); 
//     u0_sq(i) = rgig_cpp(lam, chi, ps);
//   }
//   
//   return u0_sq;
// }
// 
// // Returns tau0_square
// // [[Rcpp::export]]
// double update_tau0_hCauchy(arma::vec b0, arma::vec u0_sq, double d0, int q, int p) {
//   double chi = arma::sum(b0 % b0 / u0_sq);
//   double lam = 0.5 * (1 - q);
//   double ps = 2 * d0; 
//   double tau0_sq = rgig_cpp(lam, chi, ps);
//   
//   return tau0_sq;
// }
// 
// // [[Rcpp::export]]
// arma::vec update_c0(arma::vec u0_sq, int q) {
//   arma::vec c0(q);
//   for (int i = 0; i < q; i++) {
//     double rate = u0_sq(i) + 1;
//     c0(i) = arma::randg<double>(arma::distr_param(0.5, 1/rate));
//   }
//   
//   return c0;
// }
// 
// // [[Rcpp::export]]
// double update_d0(double tau0_sq) {
//   double rate = tau0_sq + 1;
//   return arma::randg<double>(arma::distr_param(0.5, 1/rate));
// }
// 
// // [[Rcpp::export]]
// arma::vec update_b0_Horseshoe(arma::mat X, arma::mat Y, arma::mat Sigmainv,
//                               arma::vec psi, arma::vec u0_sq, double tau0_sq,
//                               int q, int p, int n) {
//   arma::mat XtX = X.t() * X;
//   arma::mat U0i = arma::inv(arma::diagmat(1 / psi) + XtX);
//   arma::mat XU = X.t() - XtX * U0i * X.t();
//   arma::mat XUX = XtX - XtX * U0i * XtX;
//   arma::mat Un = arma::inv(XUX + arma::diagmat(1 / u0_sq) / tau0_sq);
//   double vn = 1 / (arma::as_scalar(ones(p).t() * Sigmainv * ones(p)) + 1);
//   arma::vec mn = Un * XU * Y * Sigmainv * ones(p) * vn;
//   arma::vec b0 = arma::mvnrnd(mn, Un * vn);
//   
//   return b0;
// }