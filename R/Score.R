# library(linprog)


#' Function to calculate the score matrix
#'
#' @param R The correlation matrix.
#' @inheritParams KfoldCV_delta
#'
#' @return A list including: \itemize{
#'   \item \code{score} A score matrix.
#'   \item \code{moments} A matrix of the crossproduct of \code{R}.
#' }
#'
#' @noRd

Score_mat <- function(R, q = 2, exact = F) {
  p <- nrow(R)
  score_mat <- matrix(NA, p, p)
  M <- crossprod(R)

  if (q == 2) {
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        V_ij <- M[c(i,j), c(i,j)] - crossprod(R[c(i,j),c(i,j)])
        if (V_ij[1,1] == 0 | V_ij[2,2] == 0)
          score_mat[i,j] = 0
        else {
          score_ij <- min(V_ij[1,1], V_ij[2,2]) / (p - 2) * (1 - V_ij[1,2] ** 2 / V_ij[1,1] / V_ij[2,2])
          score_mat[i,j] = sqrt(abs(score_ij))
        }
      }
    }
  } else if (q == Inf) {
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        score_mat[i,j] = min(LP_Score(R[-c(i,j), c(i,j)], 1, exact),
                             LP_Score(R[-c(i,j), c(i,j)], 2, exact))
      }
    }
  }
  return(list(score = score_mat, moments = M))
}




#' Function to calculate the \code{Inf} score
#'
#' @param R_ij A \eqn{(p-2)} by \eqn{2} matrix.
#' @param ind Either \code{1} or \code{2}.
#' @inheritParams Score_mat
#'
#' @return A numeric value.
#' @noRd

LP_Score <- function(R_ij, ind, exact = F) {
  if (exact) {
    cvec <- c(1, rep(0, 2))
    Amat <- c(0, 1, 1)     #  |v1| <= 1
    constr_mat <- R_ij[,ind] %x% t(c(1, -1))   #  (p-2) by 2 constraint matrix
    Amat <- rbind(Amat, cbind(-1, rbind(constr_mat, -constr_mat)))   # t >= || v1 R_i + R_j ||
    bvec <- c(1, -R_ij[,-ind], R_ij[,-ind])
    lpResult <- linprog::solveLP(cvec, bvec, Amat, lpSolve = T)
    return(lpResult$opt)
  } else {
    v_grid <- seq(-1, 1, length.out = 100)
    score <- min(vapply(v_grid, function(x, R_ij, ind) {
      max(abs(x * R_ij[,ind] + R_ij[,-ind]))
    }, numeric(1), R_ij = R_ij, ind = ind))
    return(score)
  }
}



