## dcovu.R
## unbiased dcov^2 and bias-corrected dcor^2
##


bcdcor <- function(x, y) {
  ## compute bias corrected distance correlation
  dcorU(x, y)
}

dcovU <-
  function(x, y) {
    ## unbiased dcov^2
    if (!inherits(x, "dist")) x <- dist(x)
    if (!inherits(y, "dist")) y <- dist(y)    
    x <- as.matrix(x)
    y <- as.matrix(y)
    n <- nrow(x)
    m <- nrow(y)
    if (n != m) stop("sample sizes must agree")
    if (! (all(is.finite(c(x, y)))))
      stop("data contains missing or infinite values")

    estimates <- dcovU_stats(x, y) #RcppExports
    return (estimates[1])
  }

dcorU <-
function(x, y) {
  ## unbiased dcov^2
  if (!inherits(x, "dist")) x <- dist(x)
  if (!inherits(y, "dist")) y <- dist(y)
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  m <- nrow(y)
  if (n != m) stop("sample sizes must agree")
  if (! (all(is.finite(c(x, y)))))
    stop("data contains missing or infinite values")

  estimates <- dcovU_stats(x, y) #RcppExports
  return (estimates[2])
}
