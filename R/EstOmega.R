################################################################################
######                                                                   #######
######     Code to estimate the precision matrix of Z via solving LPs    #######
######                                                                   #######
################################################################################
# library(linprog)


#' Estimation of the precision matrix
#'
#' \code{estOmega} estimates the inverse matrix of \code{C} via regularizations.
#'
#' @param lbd A numeric constant.
#' @param C A \eqn{K} by \eqn{K} matrix.
#'
#' @return A \eqn{K} by \eqn{K} matrix.
#' @details For a given matrix \code{C}, \code{estOmega} finds its inverse
#'   by solving the following linear program
#'   \deqn{\min_{\Omega} ||\Omega||_{\infty,1}}
#'   subject to
#'   \deqn{||C \Omega - I|| \le \lambda.}
#'   The above LP is solved by decoupling into several linear programs, each of
#'   which solves one row of \eqn{\Omega}.
#'
#' @export

estOmega <- function(lbd, C) {
  K <- nrow(C)
  omega <- matrix(0, K, K)
  for (i in 1:K) {
    omega[,i] <- solve_row(i, C, lbd)
  }
  return(omega)
}


#' Estimate each row by solving a LP
#'
#' @param col_ind An integer.
#' @inheritParams estOmega
#'
#' @return A vector of length \eqn{K}.
#' @noRd

solve_row <- function(col_ind, C_hat, lbd) {
  K <- nrow(C_hat)
  cvec <- c(1, rep(0, 2*K))
  Amat <- -cvec
  Amat <- rbind(Amat, c(-1, rep(1, 2*K)))
  tmp_constr <- C_hat %x% t(c(1,-1))
  Amat <- rbind(Amat, cbind(-1 * lbd, rbind(tmp_constr, -tmp_constr)))
  tmp_vec <- rep(0, K)
  tmp_vec[col_ind] <- 1
  bvec <- c(0, 0, tmp_vec, -tmp_vec)

  lpResult <- linprog::solveLP(cvec, bvec, Amat, lpSolve = T)$solution
  while (length(lpResult) == 0) {
    cat("The penalty lambda =", lbd, "is too small and increased by 0.01...\n")
    lbd <- lbd + 0.01
    Amat[-(1:2), 1] <- lbd
    lpResult <- linprog::solveLP(cvec, bvec, Amat, lpSolve = T)$solution[-1]
  }
  ind <- seq(2, 2*K, 2)
  return(lpResult[ind] - lpResult[ind + 1])
}



