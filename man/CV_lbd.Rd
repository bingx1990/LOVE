% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CV.R
\name{CV_lbd}
\alias{CV_lbd}
\title{Cross validation to select \eqn{\lambda}}
\usage{
CV_lbd(X, lbdGrids, AI, pureVec, diagonal)
}
\arguments{
\item{X}{A \eqn{n} by \eqn{p} data matrix.}

\item{lbdGrids}{A vector of numerical constants.}

\item{AI}{A \eqn{p} by \eqn{K} matrix.}

\item{pureVec}{The estimated set of pure variables.}

\item{diagonal}{Logical. If TRUE, the covariance matrix of \eqn{Z} is
diagonal; else FALSE.}
}
\value{
The selected \eqn{\lambda}.
}
\description{
Cross-validation to select \eqn{\lambda} for estimating the precision
  matrix of \eqn{Z}. Split the data into two parts. Estimating \eqn{Cov(Z)} on two datasets.
  Then, for each value in \code{lbdGrids}, calculate \eqn{Omega} on the first dataset
  and calculate the loss on the second dataset. Choose the value which minimizes
   \deqn{<Cov(Z), \Omega> - log(det(\Omega)).}
}
