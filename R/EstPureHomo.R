####  This script contains functions related with the estimation of the submatrix
####  of A corresponding to the pure variables.


#' @title Estimate the submatrix of \eqn{A} corresponding to the pure variables.
#'
#' @description Function to calculate the fitted submatrix \eqn{A_I} and to estimate the set and
#' partition of the pure variables.
#'
#' @param Sigma A \eqn{p} by \eqn{p} matrix.
#' @param optDelta A numerical value.
#' @inheritParams CV_Delta
#' @inheritParams LOVE
#'
#' @return A list including: \itemize{
#'    \item \code{AI} The estimated \eqn{p} by \eqn{K} submatrix \eqn{A_I} of \eqn{A}.
#'    \item \code{pureVec} The vector of the indices of estimated pure variables.
#'    \item \code{pureSignInd} The list of the indices of estimated pure variables.
#' }
#' @noRd


EstAI <- function(Sigma, optDelta, se_est, merge) {
  off_Sigma <- abs(Sigma)
  diag(off_Sigma) <- 0
  result_Ms <- FindRowMax(off_Sigma)
  Ms <- result_Ms$M
  arg_Ms <- result_Ms$arg_M

  resultPure <- FindPureNode(off_Sigma, optDelta, Ms, arg_Ms, se_est, merge)
  estPureIndices <- resultPure$pureInd
  estPureVec <- resultPure$pureVec

  estSignPureIndices <- FindSignPureNode(estPureIndices, Sigma)
  AI <- RecoverAI(estSignPureIndices, nrow(off_Sigma))

  return(list(AI = AI, pureVec = estPureVec, pureSignInd = estSignPureIndices))
}


#' Estimate the covariance matrix of \eqn{Z}.
#'
#' @inheritParams EstAI
#' @inheritParams LOVE
#' @inheritParams CV_lbd
#'
#' @return A \eqn{K} by \eqn{K} matrix.
#' @noRd

EstC <- function(Sigma, AI, diagonal) {
  K <- ncol(AI)
  C <- diag(0, K, K)
  for (i in 1:K) {
    groupi <- which(AI[ ,i] != 0)
    sigmai <- as.matrix(abs(Sigma[groupi,groupi]))
    tmpEntry <- sum(sigmai) - sum(diag(sigmai))
    C[i,i] <- tmpEntry / (length(groupi) * (length(groupi) - 1))
    if (!diagonal && i < K) {
      for (j in (i+1):K) {
        groupj <- which(AI[ ,j]!=0)
        # adjust the sign for each row
        sigmaij <- AI[groupi,i] * as.matrix(Sigma[groupi, groupj])
        sigmaij <- t(AI[groupj, j] * t(sigmaij))
        C[i,j] <- C[j,i] <- sum(sigmaij) / (length(groupi) * length(groupj))
      }
    }
  }
  return(C)
}



#' Function to calculate the maximal absolute value for each row of the given matrix.
#'
#' @inheritParams EstAI
#'
#' @return A numerical vector of \eqn{p} elements.
#' @noRd

FindRowMax <- function(Sigma) {
  p <- nrow(Sigma)
  M <- arg_M <- rep(0, p)
  for (i in 1:p) {
    row_i <- Sigma[i,]
    arg_M[i] <- which.max(row_i)
    M[i] <- row_i[arg_M[i]]
  }
  return(list(arg_M = arg_M, M = M))
}


#' Function to estimate the pure variables for a given \code{delta}.
#'
#' @param off_Sigma A \eqn{p} by \eqn{p} matrix.
#' @param delta A numeric constant.
#' @param Ms A vector containing the largest absolute values of each row of \code{off_Sigma}.
#' @param arg_Ms A vector of the indices of the largest absolute values of each
#'   row of \code{off_Sigma}.
#' @inheritParams EstAI
#' @inheritParams LOVE
#'
#' @return A list of two objects including: \itemize{
#'  \item \code{pureInd} A list of the estimated indices of the pure variables.
#'  \item \code{pureVec} A vector of the estimated indices of the pure variables.
#' }
#' @noRd

FindPureNode = function(off_Sigma, delta, Ms, arg_Ms, se_est, merge) {
  G <- list()
  for (i in 1:nrow(off_Sigma)) {
    row_i <- off_Sigma[i,]
    Si <- FindRowMaxInd(i, Ms[i], arg_Ms[i], row_i, delta, se_est)
    if (length(Si) != 0) {
      pureFlag <- TestPure(row_i, i, Si, Ms, arg_Ms, delta, se_est)
      if (pureFlag) {
        if (merge)
          G <- Merge(G, c(Si, i))
        else
          G <- Merge_union(G, c(Si, i))
      }
    }
  }
  return(list(pureInd = G, pureVec = unlist(G)))
}



#' @title Find indices of a given row.
#'
#' @description Function to calculate indices of the ith row such that the absolute
#' values of the corresponding entries are within \eqn{2\times}\code{delta} difference
#' from the given value \code{M}.
#'
#' @param i An integer.
#' @param M The maximal value of the given row.
#' @param arg_M The index of the maximal value.
#' @param vector A numeric vector.
#' @param delta A numerical constant.
#' @inheritParams EstAI
#'
#' @return A vector of indices.
#' @noRd

FindRowMaxInd <- function(i, M, arg_M, vector, delta, se_est) {
  lbd <- delta * se_est[i] * se_est[arg_M] + delta * se_est[i] * se_est
  indices <- which(M <= lbd + vector)
  return(indices)
}



#' @title Test pure variable.
#'
#' @description Function to check if a given row corresponds to a pure variable.
#'
#' @param Sigma_row A given row of \code{Sigma}.
#' @param rowInd An integer index.
#' @param Si A vector of indices.
#' @inheritParams FindPureNode
#' @inheritParams EstAI
#'
#' @return Logical. TRUE or FALSE.
#' @noRd

TestPure <- function(Sigma_row, rowInd, Si, Ms, arg_Ms, delta, se_est) {
  for (i in 1:length(Si)) {
    j <- Si[i]
    delta_j <- (se_est[rowInd] + se_est[arg_Ms[j]]) * se_est[j] * delta
    if (abs(Sigma_row[j] - Ms[j]) > delta_j)
      return(FALSE)
  }
  return(TRUE)
}



#' @title Estimate the signs of pure variables
#'
#' @description Function to estimate the sign sub-partition of the pure variables.
#' If one group has no pure variables with negative sign, then an empty list is
#' inserted in that position.
#'
#' @param pureList A list of indices of pure variables.
#' @inheritParams EstAI
#'
#' @return A list of sign sub-partition of indices.
#' @noRd

FindSignPureNode <- function(pureList, Sigma) {
  signPureList <- list()
  for (i in 1:length(pureList)) {
    purei <- pureList[[i]]
    # purei <- sort(pureList[[i]])   ### For simulation purpose only.
    if (length(purei) != 1) {
      firstPure <- purei[1]
      pos <- firstPure
      neg <- c()
      for (j in 2:length(purei)) {
        purej <- purei[j]
        if (Sigma[firstPure, purej] < 0)
          neg <- c(neg, purej)
        else
          pos <- c(pos, purej)
      }
      if (length(neg) == 0)
        neg <- list()
      signPureList[[i]] <- list(pos = pos, neg = neg)
    } else
      signPureList[[i]] <- list(pos = purei, neg = list())
  }
  return(signPureList)
}



#' Function to merge pure variables via "intersection".
#'
#' @param groupList An existing list of indices of pure variables.
#' @param groupVec A new vector of indices of pure variables.
#'
#' @return A list of indices of pure variables.
#' @noRd

Merge <- function(groupList, groupVec) {
  # merge the new group with the previous ones which have common nodes
  if (length(groupList) != 0) {
    for (i in 1:length(groupList)) {
      common_nodes <- intersect(groupList[[i]], groupVec)
      if (length(common_nodes) != 0) {
        groupList[[i]] <- common_nodes
        return(groupList)
      }
    }
  }
  groupList <- append(groupList, list(groupVec))
  return(groupList)
}

#' @describeIn Merge Merge pure variables via "union".
#' @noRd

Merge_union <- function(groupList, groupVec) {
  # merge the new group with the previous ones which have common nodes
  if (length(groupList) != 0) {
    common_groups <- sapply(groupList, FUN = function(x, y) {
      length(intersect(x, y))
    }, y = groupVec)
    common_inds <- which(common_groups > 0)
    if (length(common_inds) > 0){
      new_group <- unlist(lapply(common_inds,
                                        FUN = function(x, y){y[[x]]}, y = groupList))
      remain_group <- lapply(which(common_groups == 0),
                             FUN = function(x, y){y[[x]]}, y = groupList)
      groupList <- append(remain_group, list(union(groupVec, new_group)))
      return(groupList)
    }
  }
  groupList <- append(groupList, list(groupVec))
  return(groupList)
}



#' Function to return the estimated submatrix \eqn{A_I} from the partition of
#' pure variables.
#'
#' @param estGroupList A list of group indices of the pure variables.
#' @param p An integer.
#'
#' @return A \eqn{p} by \eqn{K} matrix.
#' @noRd

RecoverAI <- function(estGroupList, p) {
  K <- length(estGroupList)
  A <- matrix(0, p, K)
  for (i in 1:K) {
    groupi <- estGroupList[[i]]
    A[groupi[[1]],i] <- 1
    groupi2 <- groupi[[2]]
    if (length(groupi2) != 0)
      A[groupi2, i] <- -1
  }
  return(A)
}






