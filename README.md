
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LOVE

<!-- badges: start -->
<!-- badges: end -->

LOVE performs overlapping clustering of features under a structured
latent factor model.

## Installation

You can install the released version of LOVE from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("LOVE")
```

And the development version from [GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("bingx1990/LOVE")
```

## Example

This is a basic example which shows you how to use two main functions of
LoveER: LOVE and ER. We start by generating a synthetic data set.

``` r
p <- 6
n <- 50
K <- 2
A <- rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, 1), c(1/3, 2/3), c(1/2, -1/2))
Z <- matrix(rnorm(n * K, sd = 2), n, K)
E <- matrix(rnorm(n * p), n, p)
X <- Z %*% t(A) + E
```

The following code calls the LOVE function to perform overlapping
clustering of the columns of the matrix.

``` r
# library(LOVE)
# # basic example code
# res_LOVE <- LOVE(X)
# res_LOVE <- LOVE(X, pure_homo = T)
```
