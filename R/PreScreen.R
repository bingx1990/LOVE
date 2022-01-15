###     Pre-screening


Screen_X <- function(X, max_prop = 0.5) {
  p <- ncol(X)
  n <- nrow(X)

  samp_ind <- sample(1:n, floor(n / 2))
  X1 <- X[samp_ind,]
  X2 <- X[-samp_ind,]

  R1 <- cor(X1)
  off_R1 <- R1
  diag(off_R1) <- 0

  row_scale <- rowSums(off_R1 ** 2)
  thresh_range <- quantile(row_scale, c(0, max_prop))
  thresh_grid <- seq(thresh_range[1], thresh_range[2], length.out = 50)

  # sapply(thresh_grid, function(x, row_scale, X1, X2) {
  #   noise_ind <- which(row_scale < x)
  #   pred_X <- X1
  #   pred_X[,noise_ind] = 0
  #   Mse(X2, pred_X)
  # }, row_scale = row_scale, X1 = X1, X2 = X2)


  R2 <- cor(X2)
  off_R2 <- R2
  diag(off_R2) <- 0

  loss <- sapply(thresh_grid, function(x, row_scale, off_R1, off_R2) {
    noise_ind <- which(row_scale < x)
    pred_R <- off_R1
    pred_R[noise_ind, ] = 0
    pred_R[,noise_ind] = 0
    mean((off_R2 - pred_R) ** 2)
  }, row_scale = row_scale, off_R1 = off_R1, off_R2 = off_R2)

  thresh <- thresh_grid[which.min(loss)]
  which(row_scale < thresh)
}
