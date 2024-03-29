% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/theta.R
\name{update_Theta_MGP}
\alias{update_Theta_MGP}
\title{Sample the factor loading matrix using a Multiplicative Gamma Process prior}
\usage{
update_Theta_MGP(prm, Y, X, K, v1, v2)
}
\arguments{
\item{prm}{list of current parameters}

\item{Y}{n x t matrix of responses}

\item{X}{n x p matrix of correlated covariates}

\item{K}{number of latent factor}

\item{v1}{scalar or matrix of hyperparams for the MGP}

\item{v2}{scalar or matrix for hyperparams for the MGP}
}
\description{
Sample the factor loading matrix using a Multiplicative Gamma Process prior
}
