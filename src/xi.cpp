#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;

void update_random_intercept_cpp(arma::mat& xi, const arma::mat& Y,
                                 const arma::vec& id, const arma::vec& uid,
                                 const arma::mat& Sigma,
                                 const double sigmay_sqinv, const double nu_sqinv,
                                 int q, int n) {
  arma::mat Yi;
  arma::vec yi;
  for (int i = 0; i < n; i++) {
    Yi = Y.rows(find(id == uid(i)));
    yi = arma::sum(Yi, 0).t();
    int ti = Yi.n_rows;
    double si = 1/(nu_sqinv + ti * sigmay_sqinv);
    xi.row(i) = arma::mvnrnd(si * sigmay_sqinv * yi, si * Sigma).t();
  }
}
