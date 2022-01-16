###     Pre-screening


#' @title Screen features that are pure noise
#'
#' @description \code{Screen_X} finds features that are close to pure noise via
#'   k-fold cross-validation.
#'
#' @inheritParams LOVE
#' @param thresh_grid A numeric vector of thresholds. Default is \code{NULL}.
#' @param nthresh Integer. The length of \code{thresh_grid} when \code{thresh_grid} is
#'   \code{NULL}.
#' @param max_prop A numeric value between [0, 1] specifying the maximal
#'   proportional of pure noise features. Default is 0.5 meaning that at most
#'   50\% of features are pure noise.
#'
#' @return When only one value is provided in \code{thresh_grid}, \code{Screen_X}
#'   returns a vector of indices that are detected as pure noise. When either
#'   \code{thresh_grid} is \code{NULL} or multiple values are provided in
#'   \code{thresh_grid}, \code{Screen_X} returns a list including
#'   \itemize{
#'     \item\code{foldid} The indices of observations used for cv.
#'     \item\code{thresh_min} The value of \code{thresh_grid} that has the
#'       minimum cv error.
#'     \item\code{thresh_1se} The largest value of \code{thresh_grid} such that
#'      the errors are wihtin one standard error of the minimum cv error.
#'     \item\code{thresh_grid} The used \code{thresh_grid} sequence.
#'     \item\code{cv_mean} The averaged cv errors.
#'     \item\code{cv_sd} The standard errors of cv errors.
#'     \item\code{noise_ind} a vector of indices that are detected as pure noise
#'       by using \code{thresh_min}.
#'   }
#' @export
#' @seealso \code{ \link{KfoldCV_delta}}.


Screen_X <- function(X, thresh_grid = NULL, nfolds = 10, nthresh = 50,
                     max_prop = 0.5) {
  n_total <- nrow(X)
  R <- cor(X)
  diag(R) <- 0

  row_scale <- rowSums(R ** 2)

  if (is.null(thresh_grid) || length(thresh_grid) > 1) {

    if (is.null(thresh_grid)) {
      thresh_range <- quantile(row_scale, c(0, max_prop))
      thresh_grid <- seq(thresh_range[1], thresh_range[2], length.out = nthresh)
    }

    indicesPerGroup = extract(sample(1:n_total), partition(n_total, nfolds))

    loss <- matrix(NA, nfolds, nthresh)
    for (i in 1:nfolds) {
      valid_ind <- indicesPerGroup[[i]]
      trainX <- X[-valid_ind,, drop = F]
      validX <- X[valid_ind,, drop = F]

      R_train <- cor(trainX)
      R_valid <- cor(validX)
      diag(R_train) <- diag(R_valid) <- 0

      loss_i <- sapply(thresh_grid, function(x, row_scale, off_R1, off_R2) {
        noise_ind <- which(row_scale < x)
        pred_R <- off_R1
        pred_R[noise_ind,] = 0
        pred_R[,noise_ind] = 0
        mean((off_R2 - pred_R) ** 2)
      }, row_scale = row_scale, off_R1 = R_train, off_R2 = R_valid)

      loss[i,] <- loss_i
    }

    cv_mean <- apply(loss, 2, mean)
    cv_sd <- apply(loss, 2, sd)

    ind_min <- which.min(cv_mean)
    thresh_min <- thresh_grid[ind_min]
    thresh_1se <- thresh_grid[max(which(cv_mean <= (cv_mean[ind_min] + cv_sd[ind_min])))]
    noise_ind = which(row_scale < thresh_min)

    return(list(foldid = indicesPerGroup,
                thresh_min = thresh_min,
                thresh_1se = thresh_1se,
                thresh_grid = thresh_grid,
                cv_mean = cv_mean,
                cv_sd = cv_sd,
                noise_ind = noise_ind))

  } else if (length(thresh_grid) == 1) {
    return(which(row_scale < thresh_grid))
  }
}





# Screen_X <- function(X, max_prop = 0.5) {
#   p <- ncol(X)
#   n <- nrow(X)
#
#   samp_ind <- sample(1:n, floor(n / 2))
#   X1 <- X[samp_ind,]
#   X2 <- X[-samp_ind,]
#
#   R1 <- cor(X1)
#   off_R1 <- R1
#   diag(off_R1) <- 0
#
#   row_scale <- rowSums(off_R1 ** 2)
#   thresh_range <- quantile(row_scale, c(0, max_prop))
#   thresh_grid <- seq(thresh_range[1], thresh_range[2], length.out = 50)
#
#   R2 <- cor(X2)
#   off_R2 <- R2
#   diag(off_R2) <- 0
#
#   loss <- sapply(thresh_grid, function(x, row_scale, off_R1, off_R2) {
#     noise_ind <- which(row_scale < x)
#     pred_R <- off_R1
#     pred_R[noise_ind, ] = 0
#     pred_R[,noise_ind] = 0
#     mean((off_R2 - pred_R) ** 2)
#   }, row_scale = row_scale, off_R1 = off_R1, off_R2 = off_R2)
#
#   thresh <- thresh_grid[which.min(loss)]
#   which(row_scale < thresh)
# }

