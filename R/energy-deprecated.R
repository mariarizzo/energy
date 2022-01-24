## deprecated functions in energy package

dcor.ttest <- function(x, y, distance=FALSE) {
  # x and y are observed samples or distance
  # distance arg is checked in bcdcor
  .Deprecated(new = "dcorT.test", package = "energy",
              msg = "dcort.ttest is deprecated, replaced by dcorT.test")
  if (distance == TRUE) {
    x <- as.dist(x)
    y <- as.dist(y)
  }
  return(dcorT.test(x, y))
}


dcor.t <- function(x, y, distance=FALSE) {
  # computes the t statistic for corrected high-dim dCor
  # should be approximately student T
  # distance arg is checked in bcdcor
  .Deprecated(new = "dcorT", package = "energy",
              msg = "dcort.t is deprecated, replaced by dcorT")
  if (distance == TRUE) {
    x <- as.dist(x)
    y <- as.dist(y)
  }
  return(dcorT(x, y))
}


DCOR <-
  function(x, y, index=1.0) {
    ## deprecated from energy 1.7-9 
    # distance covariance and correlation statistics
    # originally: alternate method, implemented in R
    
    .Deprecated(new = "dcor", package = "energy",
                msg = "DCOR is deprecated, replaced by dcor")
    aa <- dcov(x, x, index=1.0)
    bb <- dcov(y, y, index=1.0)
    ab <- dcov(x, y, index=1.0)
    if (aa*ab > 0.0)
      r <- c(ab/sqrt(aa*ab), ab, aa, bb) else r <- 0.0
    return(r)
  }

