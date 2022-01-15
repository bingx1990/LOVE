#' Function to estimate parallel rows
#'
#' @param score_mat The score matrix.
#' @param delta A numeric constant.
#'
#' @return A list including: \itemize{
#'   \item\code{K} The cardinality of parallel rows.
#'   \item\code{I} The index of parallel rows.
#'   \item\code{I_part} The partition of parallel rows.
#' }
#'
#' @noRd

Est_Pure <- function(score_mat, delta) {
  score_flag <- which(score_mat <= delta, arr.ind = T)
  score_graph <- igraph::graph_from_data_frame(score_flag, directed = F)

  G <- igraph::groups(igraph::components(score_graph))
  return(list(K = length(G),
              I = as.numeric(unlist(G)),
              I_part = lapply(G, as.numeric)))
}


#' Function to estimate the submatrix \eqn{B_I} and \eqn{Corr(Z)}
#'
#' @param M The product matrix returned by \code{moments} of \code{\link{Score_mat}}.
#' @inheritParams Score_mat
#' @param I_part The partition of parallel rows.
#' @param I The set of parallel rows.
#' @noRd


Est_BI_C <- function(M, R, I_part, I) {
  K <- length(I_part)
  B_square <- matrix(0, nrow(R), K)
  signs <- Gamma <- rep(0, nrow(R))
  for (k in 1:K) {
    I_k <- I_part[[k]]
    for (ell in 1:length(I_k)) {
      i <- I_k[ell]

      j <- ifelse(ell < length(I_k), I_k[ell + 1], I_k[1])
      cross_Vij <- M[c(i,j), c(i,j)] - crossprod(R[c(i,j), c(i,j)])
      B_square[i, k] <- abs(R[i,j]) * sqrt(cross_Vij[1,1] / cross_Vij[2,2])

      Gamma[i] = R[i,i] - B_square[i, k]
      signs[i] = ifelse(ell == 1, 1, sign(R[i, I_k[[1]]]))
    }
  }
  B <- sqrt(B_square)
  BI <- signs[I] * B[I,,drop = F]
  B[I,] = BI
  cross_BI <- crossprod(BI)
  B_left_inv <- solve(cross_BI, t(BI))

  C_hat <- B_left_inv %*% tcrossprod(R[I,I] - diag(Gamma[I]), B_left_inv)
  diag(C_hat) <- 1


  return(list(B = B,
              C = C_hat,
              B_left_inv = B_left_inv,
              Gamma = Gamma))
}






# \code{Re_Est_Pure} re-estimates the pure variables from the selected parallel
# rows

Re_Est_Pure <- function(X, Sigma, M, I_part, Gamma) {
  L_hat <- sapply(I_part, function(x, M) {
    x <- as.numeric(x)
    row_norms <- vapply(x, function(i, M) {
      crossprod(M[i,-i])
    }, numeric(1), M = M)
    x[which.max(row_norms)]
  }, M = M)

  Gamma_LL <- Gamma[L_hat]
  K_est <- Est_K(X, L_hat, Gamma_LL)

  if (K_est < length(I_part) & K_est >= 1) { # Re-select the pure variables
    I_part_tilde <- Post_Est_Pure(Sigma, Gamma_LL, L_hat, I_part, K_est)
  } else
    I_part_tilde <- I_part

  return(I_part_tilde)
}



Post_Est_Pure <- function(Sigma, Gamma_LL, L_hat, I_part, K_tilde) {
  D_Sigma <- diag(Sigma)
  Sigma_E_LL <- Gamma_LL * D_Sigma[L_hat]
  Theta_LL <- Sigma[L_hat, L_hat, drop = F] - diag(Sigma_E_LL, length(L_hat), length(L_hat))

  D_Theta_vec <- diag(Theta_LL)
  D_Theta_LL <- diag(D_Theta_vec, length(L_hat), length(L_hat))

  L_tilde <- which.max(D_Theta_vec)
  L_tilde_comp <- setdiff(1:length(I_part), L_tilde)

  if (K_tilde > 1) {
    for (k in 2:K_tilde) {
      Theta_LL_schur <- D_Theta_LL[-L_tilde, -L_tilde, drop = F] - Theta_LL[-L_tilde, L_tilde, drop = F] %*%
        solve(Theta_LL[L_tilde, L_tilde, drop = F], Theta_LL[L_tilde, -L_tilde, drop = F])
      i_k <- L_tilde_comp[which.max(diag(Theta_LL_schur))]
      L_tilde <- c(L_tilde, i_k)
      L_tilde_comp <- setdiff(L_tilde_comp, i_k)
    }
  }

  lapply(L_tilde, function(x, I_part) {I_part[[x]]}, I_part = I_part)
}




#' @title Function to estimate the number of latent factors
#'
#' @description Function to (re-)estimate the number of latent factors from a set
#'   of representative parallel rows. This should only be used after a given set of
#'   representative parallel rows are selected.
#'
#' @inheritParams KfoldCV_delta
#' @param L_hat A vector of he representative indices of parallel rows.
#' @param Gamma_LL A vector of numeric values.
#'
#' @return An integer. The estimated \eqn{K}.
#' @noRd


Est_K <- function(X, L_hat, Gamma_LL) {

  n <- nrow(X)
  K_hat <- length(L_hat)
  n_ind <- sample(1:n, floor(n / 2))
  X1 <- X[n_ind,]
  X2 <- X[-n_ind,]

  R1 <- cor(X1)
  R2 <- cor(X2)

  Gamma_LL_mat <- diag(Gamma_LL, K_hat, K_hat)
  M1 <- R1[L_hat, L_hat, drop = F] - Gamma_LL_mat
  M2 <- R2[L_hat, L_hat, drop = F] - Gamma_LL_mat

  res_eig <- eigen(M1)

  K_tilde <- which.min(sapply(1:K_hat, function(x, res_eig, M2) {
    U <- res_eig$vectors[,1:x, drop = F]
    M_tilde <- U %*% diag(res_eig$values)[1:x, 1:x, drop = F] %*% t(U)
    sum((M_tilde - M2) ** 2)
  }, res_eig = res_eig, M2 = M2))

  return(K_tilde)
}




