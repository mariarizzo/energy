## use the Rcpp exported function U_center or D_center
## the utilities in this file are provided for reference and historical reasons

Dcenter <- function(x) {
  ## x is a dist object or data matrix
  if (!inherits(x, "dist")) x <- dist(x)
  d <- as.matrix(x)
  n <- nrow(d)
  m <- rowSums(d)
  M <- sum(m) / n^2
  m <- m / n
  a <- sweep(d, 1, m)
  b <- sweep(a, 2, m)
  B <- b + M
}

Ucenter <- function(x) {
  ## x is a dist object or data matrix
  if (!inherits(x, "dist")) x <- dist(x)
  d <- as.matrix(x)
  n <- nrow(d)
  m <- rowSums(d)
  M <- sum(m) / ((n-1)*(n-2))
  m <- m / (n-2)
  a <- sweep(d, 1, m)
  b <- sweep(a, 2, m)
  B <- b + M
  diag(B) <- 0
  B
}

