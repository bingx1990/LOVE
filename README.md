
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LOVE

<!-- badges: start -->
<!-- badges: end -->

LOVE performs overlapping clustering of feature variables under a
structured latent factor model.

## Installation

<!-- the released version of LOVE from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ``` r -->
<!-- install.packages("LOVE") -->
<!-- ``` -->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("bingx1990/LOVE")
```

## Example

This is a basic example which shows you how to use the main function of
LOVE. We start by generating a synthetic data set.

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
library(LOVE)

res_LOVE <- LOVE(X, pure_homo = FALSE)
res_LOVE <- LOVE(X, pure_homo = TRUE, delta = seq(0.1, 1.1 ,0.1))
```

## Practical pre-screening

In practice, we recommend a pre-screening procedure before calling the
function. The function detects the features that are close to pure
noise. The following example demonstrates the usage of .

``` r
aug_A <- rbind(A, c(0, 0))
aug_p <- nrow(aug_A)
E <- matrix(rnorm(n * aug_p), n, aug_p)
X <- Z %*% t(aug_A) + E

noise_ind <- Screen_X(X)$noise_ind
cat("Features with indices in", noise_ind, "are detected as pure noise.\n")
#> Features with indices in 7 are detected as pure noise.
feature_ind <- setdiff(1:aug_p, noise_ind)

res_LOVE <- LOVE(X[,feature_ind,drop = F], pure_homo = FALSE)
```
