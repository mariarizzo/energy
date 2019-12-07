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

