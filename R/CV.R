#####     Functions to select tuning parameters via cross validation


### Functions to select delta for homogenous pure loadings


#' @title Cross validation to select \eqn{\delta}
#'
#' @description Cross validation for choosing the tuning parameter \eqn{\delta}.
#'   For each value of \code{deltaGrids}, first split the data into two parts
#'   and calculate \eqn{I}, \eqn{A_I} and \eqn{Cov(Z)}. Then calculate the fit
#'   \eqn{A_I Cov(Z) A_I'} to find the value which minimizes the loss criterion
#'   \deqn{||\Sigma - A_I Cov(Z) A_I'||_{F-off}/(|I|(|I|-1))}.
#'
#' @inheritParams LOVE
#' @param se_est The vector of the standard deviations of \eqn{p} features.
#' @param deltaGrids A vector of numerical constants.
#'
#' @return A numeric constant. The selected optimal \eqn{\delta}.

CV_Delta <- function(X, deltaGrids, diagonal, se_est, merge) {
  n <- nrow(X); p <- ncol(X)
  sampInd <- sample(n, floor(n / 2))
  X1 <- X[sampInd, ]
  X2 <- X[-sampInd, ]
  Sigma1 <- crossprod(X1) / nrow(X1);
  diag(Sigma1) <- 0
  Sigma2 <- crossprod(X2) / nrow(X2)

  result_Ms <- FindRowMax(abs(Sigma1))
  Ms <- result_Ms$M
  arg_Ms <- result_Ms$arg_M

  loss <- c()
  for (i in 1:length(deltaGrids)) {
    resultFitted <- CalFittedSigma(Sigma1, deltaGrids[i], Ms, arg_Ms, se_est,
                                   diagonal, merge)
    fittedValue <- resultFitted$fitted
    estPureVec <- resultFitted$pureVec
    if (is.null(dim(fittedValue)) && fittedValue == -1)
      loss[i] <- Inf
    else {
      denom <- length(estPureVec) * (length(estPureVec) - 1)
      # loss[i] <- 2 * offSum(Sigma2[estPureVec, estPureVec], fittedValue, 1) / denom
      loss[i] <- 2 * offSum(Sigma2[estPureVec, estPureVec] - fittedValue, se_est[estPureVec]) / denom
    }
  }
  # cat(loss)
  return(deltaGrids[which.min(loss)])
}




#' Calculate the fitted value \eqn{A_I Cov(Z) A_I'}.
#'
#' @inheritParams LOVE
#' @inheritParams EstAI
#' @inheritParams FindPureNode
#'
#' @return A list including: \itemize{
#'    \item \code{pureVec} A vector of the indices of the estimated pure variables.
#'    \item \code{fitted} The fitted value \eqn{A_I Cov(Z) A_I'}.
#' }
#' @noRd

CalFittedSigma <- function(Sigma, delta, Ms, arg_Ms, se_est, diagonal, merge) {
  resultPureNode <- FindPureNode(abs(Sigma), delta, Ms, arg_Ms, se_est, merge)

  estPureIndices <- resultPureNode$pureInd
  # lapply(estPureIndices, function(x) cat(x, "\n"))

  if (singleton(estPureIndices))
    return(list(pureVec = NULL, fitted = -1))

  estSignPureIndices <- FindSignPureNode(estPureIndices, Sigma)
  AI <- RecoverAI(estSignPureIndices, length(se_est))
  C <- EstC(Sigma, AI, diagonal)

  if (length(estPureIndices) == 1)
    fitted <- -1
  else {
    subAI <- AI[resultPureNode$pureVec, ]
    fitted <- subAI %*% C %*% t(subAI)
  }
  return(list(pureVec = resultPureNode$pureVec, fitted = fitted))
}








###   Functions to select delta for heterogeneous pure loadings


#' @title Cross-validation for selecting \eqn{delta}
#'
#' @description \code{KfoldCV_delta} uses \code{nfolds} cross-validation to
#'   select \eqn{delta}.
#'
#' @inheritParams LOVE
#' @param ndelta Integer. The length of the grid of \code{delta}.
#' @param q Either \code{2} or \code{Inf} to specify the type of score.
#' @param exact Logical. Only active for compute the \code{Inf} score.
#'   If TRUE, compute the \code{Inf} score exactly via solving a linear program.
#'   Otherwise, use approximation to compute \code{Inf} score.
#' @param nfolds The number of folds. Default is 10.
#' @param max_pure A numeric value between (0, 1] specifying the maximal
#'   proportion of pure variables. Default is NULL. When not specified,
#'   \code{max_pure} = 1 if \eqn{n > p}, \code{max_pure} = 0.8 otherwise.
#'
#' @return A list of objects including: \itemize{
#'   \item\code{foldid} The indices of observations used for cv.
#'   \item\code{delta_min} The value of \code{delta} that has the minimum cv
#'     error.
#'   \item\code{delta_1se} The smallest and largest value of \code{delta} such
#'     that the cv errors are wihtin one standard error of the minimum cv error.
#'   \item\code{delta} The used \code{delta} sequence.
#'   \item\code{cv_mean} The averaged cv errors.
#'   \item\code{cv_sd} The standard errors of cv errors.
#'   \item\code{est_pure} A list including: \itemize{
#'     \item\code{K} The cardinality of parallel rows.
#'     \item\code{I} The index of parallel rows.
#'     \item\code{I_part} The partition of parallel rows.
#'   }
#'   \item\code{score} The score matrix.
#'   \item\code{moments} The crossproduct matrix \eqn{R'R}.
#' }
#'
#' @export


KfoldCV_delta <- function(X, delta = NULL, ndelta = 50, q = 2, exact = FALSE,
                          nfolds = 10, max_pure = NULL) {
  n_total <- nrow(X)
  p_total <- ncol(X)
  R <- cor(X)

  score_res <- Score_mat(R, q, exact)
  score_mat <- score_res$score
  moments_mat <- score_res$moments

  if (is.null(delta)) {
    if (is.null(max_pure))
      max_pure <- ifelse(n_total > p_total, 1, 0.8)

    delta_max <- quantile(apply(score_mat[-p_total,], 1, min, na.rm = T),
                          probs = max_pure)
    delta <- seq(delta_max, min(score_mat, na.rm = T), length.out = ndelta)
  }

  indicesPerGroup = extract(sample(1:n_total), partition(n_total, nfolds))

  loss <- matrix(NA, nfolds, length(delta))
  for (i in 1:nfolds) {
    valid_ind <- indicesPerGroup[[i]]
    trainX <- X[-valid_ind,, drop = F]
    validX <- X[valid_ind,, drop = F]
    R1 <- cor(trainX);  R2 <- cor(validX)

    score_res <- Score_mat(R1, q, exact)
    score_mat_1 <- score_res$score
    moments_1 <- score_res$moments


    for (j in 1:length(delta)) {
      delta_j <- delta[j]
      if (j == 1) {
        pure_res <- Est_Pure(score_mat_1, delta_j)
        I <- pure_res$I
        I_part <- pure_res$I_part
      } else {
        pure_res <- Est_Pure(score_mat_1[pre_I, pre_I], delta_j)
        I <- pre_I[pure_res$I]
        I_part <- lapply(pure_res$I_part, function(x, pre_I) {pre_I[as.numeric(x)]}, pre_I = pre_I)
      }
      pre_I <- I

      if (length(I_part) == 0)
        break
      else {
        result <- Est_BI_C(moments_1, R1, I_part, I)
        B_hat <- result$B
        C_hat <- result$C
        B_hat_left_inv <- result$B_left_inv
        tmp_R1 <- R1
        tmp_R1[I,I] <- B_hat[I,,drop = F] %*% tcrossprod(C_hat, B_hat[I,, drop = F])
        if (length(I) != p_total) {
          tmp <- B_hat_left_inv %*% tmp_R1[I,-I]
          tmp_prime <- try(solve(C_hat, tmp), silent = F)
          if (class(tmp_prime)[1] == "try-error")
            tmp_prime <- ginv(C_hat) %*% tmp
          tmp_R1[-I, -I] <- crossprod(tmp, tmp_prime)
        }
        loss[i,j] <- offSum(tmp_R1 - R2, 1) / p_total / (p_total - 1)
      }
    }
  }
  cv_mean <- apply(loss, 2, mean)
  cv_sd <- apply(loss, 2, sd)

  ind_min <- which.min(cv_mean)
  delta_min <- delta[ind_min]
  # delta_1se <- delta[min(which(cv_mean <= (cv_mean[ind_min] + cv_sd[ind_min])))]
  delta_1se <- delta[range(which(cv_mean <= (cv_mean[ind_min] + cv_sd[ind_min])))]
  return(list(foldid = indicesPerGroup,
              delta_min = delta_min,
              delta_1se = delta_1se,
              delta = delta,
              cv_mean = cv_mean,
              cv_sd = cv_sd,
              est_pure = Est_Pure(score_mat, delta_min),
              score = score_mat,
              moments = moments_mat))
}





Compute_Loss <- function(delta, score_mat_1, moments_1, R1, R2) {
  pure_res <- Est_Pure(score_mat_1, delta)
  I_part <- pure_res$I_part

  if (length(I_part) == 0)
    loss <- NA
  else {
    I <- pure_res$I
    result <- Est_BI_C(moments_1, R1, I_part, I)
    B_hat <- result$B
    C_hat <- result$C
    B_hat_left_inv <- result$B_left_inv
    R1[I,I] <- B_hat[I,,drop = F] %*% tcrossprod(C_hat, B_hat[I,, drop = F])
    if (length(I) != nrow(R1)) {
      # J <- setdiff(1:nrow(R1), I)
      # R1[J, J] <- R1[J, I] %*% ginv(R1[I,I]) %*% R1[I, J]

      tmp <- B_hat_left_inv %*% R1[I,-I]
      tmp_prime <- try(solve(C_hat, tmp), silent = T)
      if (class(tmp_prime)[1] == "try-error")
        tmp_prime <- ginv(C_hat) %*% tmp
      R1[-I, -I] <- crossprod(tmp, tmp_prime)
    }
    loss <- Mean_Off_Diag(R1, R2)
  }
  loss
}






### Functions to select lambda for estimating the precision matrix


#' @title Cross validation to select \eqn{\lambda}
#'
#' @description Cross-validation to select \eqn{\lambda} for estimating the precision
#'   matrix of \eqn{Z}. Split the data into two parts. Estimating \eqn{Cov(Z)} on two datasets.
#'   Then, for each value in \code{lbdGrids}, calculate \eqn{Omega} on the first dataset
#'   and calculate the loss on the second dataset. Choose the value which minimizes
#'    \deqn{<Cov(Z), \Omega> - log(det(\Omega)).}
#'
#' @inheritParams LOVE
#' @param lbdGrids A vector of numerical constants.
#' @param AI A \eqn{p} by \eqn{K} matrix.
#' @param pureVec The estimated set of pure variables.
#'
#' @return The selected \eqn{\lambda}.

CV_lbd <- function(X, lbdGrids, AI, pureVec, diagonal) {
  sampInd <- sample(nrow(X), floor(nrow(X) / 2))
  X1 <- X[sampInd, ]
  X2 <- X[-sampInd, ]
  Sigma1 <- crossprod(X1) / nrow(X1)
  Sigma2 <- crossprod(X2) / nrow(X2)
  C1 <- EstC(Sigma1, AI, diagonal)
  C2 <- EstC(Sigma2, AI, diagonal)
  loss <- c()
  for (i in 1:length(lbdGrids)) {
    Omega <- estOmega(lbdGrids[i], C1)
    det_Omega <- det(Omega)
    loss[i] <- ifelse(det_Omega <= 0, Inf, sum(Omega * C2) - log(det_Omega))
  }
  return(lbdGrids[which.min(loss)])
}
