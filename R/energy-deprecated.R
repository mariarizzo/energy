## deprecated functions in energy package


DCOR <-
  function(x, y, index=1.0) {
    # distance covariance and correlation statistics
    # alternate method, implemented in R without .C call
    # this method is usually slower than the C version
    
    
    .Deprecated(new = "dcor", package = "energy",
                msg = "DCOR is deprecated, replaced by dcor or dcov")
    
    if (!inherits(x, "dist")) x <- dist(x)
    if (!inherits(y, "dist")) y <- dist(y)
    x <- as.matrix(x)
    y <- as.matrix(y)
    n <- nrow(x)
    m <- nrow(y)
    if (n != m) stop("Sample sizes must agree")
    if (! (all(is.finite(c(x, y)))))
      stop("Data contains missing or infinite values")
    if (index < 0 || index > 2) {
      warning("index must be in [0,2), using default index=1")
      index=1.0}
    
    stat <- 0
    dims <- c(n, ncol(x), ncol(y))
    
    Akl <- function(x) {
      d <- as.matrix(x)^index
      m <- rowMeans(d)
      M <- mean(d)
      a <- sweep(d, 1, m)
      b <- sweep(a, 2, m)
      return(b + M)
    }
    
    A <- Akl(x)
    B <- Akl(y)
    dCov <- sqrt(mean(A * B))
    dVarX <- sqrt(mean(A * A))
    dVarY <- sqrt(mean(B * B))
    V <- sqrt(dVarX * dVarY)
    if (V > 0)
      dCor <- dCov / V else dCor <- 0
    return(list(dCov=dCov, dCor=dCor, dVarX=dVarX, dVarY=dVarY))
  }

