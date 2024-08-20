## defunct functions from the energy package

dcor.ttest <- function(x, y, distance=FALSE) {
  .Defunct(new = "dcorT.test", package = "energy",
              msg = "dcort.ttest replaced by dcorT.test")
}

dcor.t <- function(x, y, distance=FALSE) {
  .Deprecated(new = "dcorT", package = "energy",
              msg = "dcor.t replaced by dcorT")
}


