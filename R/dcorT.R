### dcorT.R
### implementation of the distance correlation t-test

Astar <- function(d) {
  ## d is a distance matrix or distance object
  ## modified or corrected doubly centered distance matrices
  ## denoted A* (or B*) in JMVA t-test paper (2013)
  d <- as.matrix(d)
  n <- nrow(d)
  if (n != ncol(d)) stop("Argument d should be distance")
  m <- rowMeans(d)
  M <- mean(d)
  a <- sweep(d, 1, m)
  b <- sweep(a, 2, m)
  A <- b + M  #same as plain A
  #correction to get A^*
  A <- A - d/n
  diag(A) <- m - M
  (n / (n-1)) * A
}

BCDCOR <- function(x, y, distance=FALSE) {
  ## compute bias corrected distance correlation
  ## attempt to check if distance flag is valid
  if (distance==FALSE) {
    if (class(x)=="dist" || class(y)=="dist")
      stop("distance==FALSE but argument is a dist object")
    x <- as.matrix(dist(x))
    y <- as.matrix(dist(y))
  } else {
    x <- as.matrix(x)
    y <- as.matrix(y)
    if (distance == TRUE)
      if (!isSymmetric(x) || !isSymmetric(y))
        stop("distance==TRUE but matrices non-symmetric")
  }
  n <- NROW(x)
  AA <- Astar(x)
  BB <- Astar(y)
  XY <- sum(AA*BB) - (n/(n-2)) * sum(diag(AA*BB))
  XX <- sum(AA*AA) - (n/(n-2)) * sum(diag(AA*AA))
  YY <- sum(BB*BB) - (n/(n-2)) * sum(diag(BB*BB))
  list(bcR=XY / sqrt(XX*YY), XY=XY/n^2, XX=XX/n^2, YY=YY/n^2, n=n)
}

dcor.t <- function(x, y, distance=FALSE) {
  # computes the t statistic for corrected high-dim dCor
  # should be approximately student T
  # distance arg is checked in bcdcor
  r <- BCDCOR(x, y, distance)
  Cn <- r$bcR
  n <- r$n
  M <- n*(n-3)/2
  sqrt(M-1) * Cn / sqrt(1-Cn^2)
}

dcor.ttest <- function(x, y, distance=FALSE) {
  # x and y are observed samples or distance
  # distance arg is checked in bcdcor
  dname <- paste(deparse(substitute(x)),"and",
           deparse(substitute(y)))
  stats <- BCDCOR(x, y, distance)
  bcR <- stats$bcR
  n <- stats$n
  M <- n * (n-3) / 2
  df <- M - 1
  names(df) <- "df"
  tstat <-  sqrt(M-1) * bcR / sqrt(1-bcR^2)
  names(tstat) <- "T"
  estimate <- bcR
  names(estimate) <- "Bias corrected dcor"
  pval <- 1 - pt(tstat, df=df)
  method <- "dcor t-test of independence"
  rval <- list(statistic = tstat, parameter = df, p.value = pval,
          estimate=estimate, method=method, data.name=dname)
  class(rval) <- "htest"
  return(rval)
}
