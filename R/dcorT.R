### dcorT.R
### implementation of the distance correlation t-test
### for high dimension

Astar <- function(d) {
  ## d is a distance matrix or distance object
  ## modified or corrected doubly centered distance matrices
  ## denoted A* (or B*) in JMVA t-test paper (2013)
  if (inherits(d, "dist"))
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

BCDCOR <- function(x, y) {
  ## compute bias corrected distance correlation
  ## internal function not in NAMESPACE (external: use bcdcor) 
  ## revised version from v. 1.7-7 
  if (!inherits(x, "dist")) {
    x <- as.matrix(dist(x))
  } else {
    x <- as.matrix(x)
  }
  if (!inherits(y, "dist")) {
    y <- as.matrix(dist(y))
  } else {
    y <- as.matrix(y)
  }
  
  n <- NROW(x)
  AA <- Astar(x)
  BB <- Astar(y)
  XY <- sum(AA*BB) - (n/(n-2)) * sum(diag(AA*BB))
  XX <- sum(AA*AA) - (n/(n-2)) * sum(diag(AA*AA))
  YY <- sum(BB*BB) - (n/(n-2)) * sum(diag(BB*BB))
  list(bcR=XY / sqrt(XX*YY), XY=XY/n^2, XX=XX/n^2, YY=YY/n^2, n=n)
}



dcorT <- function(x, y) {
  # computes the t statistic for corrected high-dim dCor
  # should be approximately student T
  # x and y are observed samples or distance objects
  r <- BCDCOR(x, y)
  Cn <- r$bcR
  n <- r$n
  M <- n*(n-3)/2
  sqrt(M-1) * Cn / sqrt(1-Cn^2)
}

dcorT.test <- function(x, y) {
  # x and y are observed samples or distance objects
  dname <- paste(deparse(substitute(x)),"and",
                 deparse(substitute(y)))
  stats <- BCDCOR(x, y)
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
  method <- "dcor t-test of independence for high dimension"
  rval <- list(statistic = tstat, parameter = df, p.value = pval,
               estimate=estimate, method=method, data.name=dname)
  class(rval) <- "htest"
  return(rval)
}

