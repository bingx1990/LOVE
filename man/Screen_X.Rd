% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PreScreen.R
\name{Screen_X}
\alias{Screen_X}
\title{Screen features that are pure noise}
\usage{
Screen_X(X, thresh_grid = NULL, nfolds = 10, nthresh = 50, max_prop = 0.5)
}
\arguments{
\item{X}{A \eqn{n} by \eqn{p} data matrix.}

\item{thresh_grid}{A numeric vector of thresholds. Default is \code{NULL}.}

\item{nfolds}{The number of folds. Default is 10.}

\item{nthresh}{Integer. The length of \code{thresh_grid} when \code{thresh_grid} is
\code{NULL}.}

\item{max_prop}{A numeric value between [0, 1] specifying the maximal
proportional of pure noise features. Default is 0.5 meaning that at most
50\% of features are pure noise.}
}
\value{
When only one value is provided in \code{thresh_grid}, \code{Screen_X}
  returns a vector of indices that are detected as pure noise. When either
  \code{thresh_grid} is \code{NULL} or multiple values are provided in
  \code{thresh_grid}, \code{Screen_X} returns a list including
  \itemize{
    \item\code{foldid} The indices of observations used for cv.
    \item\code{thresh_min} The value of \code{thresh_grid} that has the
      minimum cv error.
    \item\code{thresh_1se} The largest value of \code{thresh_grid} such that
     the errors are wihtin one standard error of the minimum cv error.
    \item\code{thresh_grid} The used \code{thresh_grid} sequence.
    \item\code{cv_mean} The averaged cv errors.
    \item\code{cv_sd} The standard errors of cv errors.
    \item\code{noise_ind} a vector of indices that are detected as pure noise
      by using \code{thresh_min}.
  }
}
\description{
\code{Screen_X} finds features that are close to pure noise via
  k-fold cross-validation.
}
\seealso{
\code{ \link{KfoldCV_delta}}.
}
