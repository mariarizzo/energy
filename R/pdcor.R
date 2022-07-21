## pdcor.R
##
##

pdcor <- function(x, y, z) {
  x <- .arg2dist.matrix(x)
  y <- .arg2dist.matrix(y)
  z <- .arg2dist.matrix(z)
  partial_dcor(x, y, z)["pdcor"]
}

pdcov <- function(x, y, z) {
  x <- .arg2dist.matrix(x)
  y <- .arg2dist.matrix(y)
  z <- .arg2dist.matrix(z)
  partial_dcov(x, y, z)
}

