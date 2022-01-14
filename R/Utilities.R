#####      This script contains some helper functions

#' Recover clusters based on a given \eqn{p} by \eqn{K} loading matrix.
#'
#' @param A A \eqn{p} by \eqn{K} matrix.
#'
#' @return A list of group indices with sign sub-partition.
#' @noRd

recoverGroup <- function(A) {
  Group <- list()
  for (i in 1:ncol(A)) {
    column <- A[,i]
    posInd <- which(column > 0)
    negInd <- which(column < 0)
    Group[[i]] <- list(pos = posInd, neg = negInd)
  }
  return(Group)
}


#' Function to check if there exists any element in the given list that has length
#' equal to 1
#'
#' @param estPureIndices A list of indices of the estimated pure variables.
#'
#' @return Logical. If exists at least one, return TRUE; otherwise return FALSE
#' @noRd

singleton <- function(estPureIndices) {
  if (length(estPureIndices) == 0)
    return(T)
  else
    ifelse(sum(sapply(estPureIndices, FUN = function(x) {length(x)}) == 1) > 0, T, F)
}


#' @title Function to hard-threshold
#'
#' @description Threshold the estimated \code{A} based on the given \code{mu}.
#'   If \code{scale} is TRUE, then normalize each row of \code{A} such that the
#'   \eqn{\ell-1} norm of each row is no larger than 1.
#'
#' @inheritParams recoverGroup
#' @param mu A numeric value.
#' @param scale Logical. Normalize the row-wise \eqn{\ell-1} norm if TRUE.
#'
#' @return A matrix with the same dimension as \code{A}.
#' @noRd

threshA <- function(A, mu, scale = FALSE) {
  scaledA <- A
  for (i in 1:nrow(A)) {
    colInd <- abs(A[i, ]) <= mu
    scaledA[i,colInd] = 0
    if (scale && sum(abs(scaledA[i, ])) > 1)
      scaledA[i, ] <- scaledA[i, ] / sum(abs(scaledA[i, ]))
  }
  return(scaledA)
}



#' Calculate the weighted sum of squares of the upper off-diagonal elements
#'
#' @param M A given symmetric matrix
#' @param weights A vector of length equal to the number of rows of \code{M}.
#'
#' @return A numeric value.
#' @noRd

offSum <- function(M, weights) {
  tmp <- M / weights
  tmp <- t(t(tmp) / weights)
  return(sum((tmp[row(tmp) <= (col(tmp) - 1)])^2))
}




partition = function(totalNumb, numbGroup) {
  # This function returns a vector of numbers for each group given the total number of
  # observations and the total number of groups
  # eg: divide 9 obs into 4 groups should give c(3:2:2:2)
  # args: {@code totalNumb} "integer", # of obs
  # return: a vector of length equal to {@code numbGroup}

  remainder = totalNumb %% numbGroup # get the remainder
  numbPerGroup = totalNumb %/% numbGroup
  rep(numbPerGroup,numbGroup) + c(rep(1,remainder),rep(0,numbGroup-remainder))
}


extract = function(preVec, indices) {
  # Extract the indices from the previous vector and return as a list with length equal to length of {@code indices}
  # args: preVec: a vector
  # args: indices: contains the length of each group to be extracted
  # return: list[length=length(indices)]

  newVec = vector("list",length(indices))
  newVec[[1]] = preVec[1:indices[1]]
  for (i in 2:length(indices)) {
    tmpIndices = sum(indices[1:(i-1)]) + 1
    newVec[[i]] = preVec[tmpIndices:(tmpIndices+indices[i]-1)]
  }
  return(newVec)
}







# Thresh <- function(M, thresh) {
#   if (thresh == 0)
#     return(M)
#   else {
#     M_copy <- M
#     M_copy[abs(M) <= thresh] = 0
#     return(M_copy)
#   }
# }
#
#
# Mean_Off_Diag <- function(M, N) {
#   # Calculate the sum of squares of the upper off-diagonal elements of two matrices
#   # require: M and N have the same dimensions
#   tmp <- M-N
#   diag(tmp) <- 0
#   return(mean(tmp ** 2))
# }
