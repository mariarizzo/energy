dcor2d<- function(x, y, type = c("V", "U")) {
  ## computes dcor^2 or bias-corrected dcor^2 by O(n log n) algorithm
  ## bivariate data only: (x,y) in R^2
  ## should be faster than direct calc. for big n
  type <- match.arg(type)
  ## argument checking in dcov2d
  stat <- dcov2d(x, y, type, all.stats=TRUE)
  dvarX <- stat[2]
  dvarY <- stat[3]
  R2 <- 0.0
  if (abs(dvarX*dvarY > 10*.Machine$double.eps))
    R2 <- stat[1] / sqrt(dvarX*dvarY)
  return (R2)
}

dcov2d<- function(x, y, type=c("V", "U"), all.stats=FALSE) {
  ## O(n log n) computation of dcovU or dcov^2 (V^2) for (x, y) in R^2 only
  type <- match.arg(type)
  if (!is.vector(x) || !is.vector(y)) {
    if (NCOL(x) > 1 || NCOL(y) > 1)
      stop("this method is only for univariate x and y")
  }
  x <- as.vector(x)
  y <- as.vector(y)
  n <- length(x)
  if (n != length(y))
    stop("sample sizes must agree")
  
  Sums <- .dcovSums2d(x, y, all.sums=all.stats)
  if (type =="V") {
    d1 <- n^2
    d2 <- n^3
    d3 <- n^4
  } else {
    d1 <- n * (n - 3)
    d2 <- d1 * (n - 2)
    d3 <- d2 * (n - 1)
  }
  dCov2d <- Sums$S1/d1 - 2*Sums$S2/d2 + Sums$S3/d3
  if (all.stats) {
    dvarX <- Sums$S1a/d1 - 2*Sums$S2a/d2 + Sums$S3a/d3
    dvarY <- Sums$S1b/d1 - 2*Sums$S2b/d2 + Sums$S3b/d3
  }
  rval <- ifelse(type=="V", c(V=dCov2d), c(U=dCov2d))
  if (all.stats)
    rval <- c(rval, dvarX=dvarX, dvarY=dvarY) 
  return (rval)
}


.dcovSums2d <- function(x, y, all.sums = FALSE) {
  ## compute the sums S1, S2, S3 of distances for dcov^2
  ## dCov^2 <- S1/d1 - 2 * S2/d2 + S3/d3  
  ## denominators differ for U-statistic, V-statisic
  ## if all.sums==TRUE, also return sums for dVar and kernel 
  if (is.matrix(x) || is.matrix(y)) {
    if (ncol(x) > 1 || ncol(y) > 1)
      stop("Found multivariate (x,y) in .dcovSums2d, expecting bivariate")
  }
  n <- length(x)
  SRx <- sortrank(x)
  SRy <- sortrank(y)
  ## compute the rowSums of the distance matrices
  a. <- .rowSumsDist1(x, SRx)
  b. <- .rowSumsDist1(y, SRy)
  S2 <- sum(a. * b.)
  a.. <- sum(a.)
  b.. <- sum(b.)
  S3 <- sum(a.) * sum(b.)
  
  ## also need order and rank for y[order(x)] in gamma1()
  x1 <- SRx$x
  y1 <- y[SRx$ix]
  SRy1 <- sortrank(y1)
  ones <- rep(1, n)
  g_1 <-  .gamma1(x1=x1, y1=y1, z1=ones,  SRx=SRx, SRy1=SRy1)
  g_x <-  .gamma1(x1=x1, y1=y1, z1=x1,    SRx=SRx, SRy1=SRy1)
  g_y <-  .gamma1(x1=x1, y1=y1, z1=y1,    SRx=SRx, SRy1=SRy1)
  g_xy <- .gamma1(x1=x1, y1=y1, z1=x1*y1, SRx=SRx, SRy1=SRy1)
  S1 <- sum(x * y * g_1 + g_xy - x * g_y - y * g_x)

  L <- list(S1=S1, S2=S2, S3=S3, 
            S1a=NA, S1b=NA, S2a=NA, S2b=NA, S3a=NA, S3b=NA,
            rowsumsA=NA, rowsumsB=NA, sumA=NA, sumB=NA)
  if (all.sums) {
    L$S1a <- 2 * n * (n-1) * var(x)
    L$S1b <- 2 * n * (n-1) * var(y)
    L$S2a <- sum(a.^2)
    L$S2b <- sum(b.^2)
    L$S3a <- a..^2
    L$S3b <- b..^2
    L$rowsumsA <- a.
    L$rowsumsB <- b.
    L$sumA <- a..
    L$sumB <- b..
  }
  return (L);
}

.dvarU2 <- function(x, SRx = NULL) {
  ## O(n log n) computation of dvarU for univariate x only
  ## this is an internal function that will do a stand-alone dVar calc.
  ## but it is not faster than dcovU2(x, x) unless we supply
  ## the precomputed sort + rank results in SRx
  n <- length(x)
  ## compute the rowSums of the distance matrices
  if (is.null(SRx))
    SRx <- sortrank(x)
  a. <- .rowSumsDist1(x, SRx)
  S2 <- sum(a. * a.)
  S3 <- sum(a.)^2

  ## also need order and rank for y[order(x)] in gamma1()
  x1 <- SRx$x
  x2 <- x1
  SRx1 <- sortrank(x1)
  ones <- rep(1, n)
  g_1 <-  .gamma1(x1=x1, y1=x2, z1=ones,  SRx, SRx1)
  g_x <-  .gamma1(x1=x1, y1=x2, z1=x1,    SRx, SRx1)
  g_xx <- .gamma1(x1=x1, y1=x2, z1=x1*x2, SRx, SRx1)
  S1 <- sum(x^2 * g_1 + g_xx - 2 * x * g_x)
  d1 <- n * (n - 3)
  d2 <- d1 * (n - 2)
  d3 <- d2 * (n - 1)
  dVar <- S1/d1 - 2 * S2/d2 + S3/d3
  return(dVar)
}

.gamma1 <- function(x1, y1, z1, SRx, SRy1) {
  # computes the terms of the sum (ab) in dcovU
  # original sample (x_i, y_i, z_i)
  # triples (x1_i, y1_i, z1_i) are sorted by ix=order(x)
  # SRx is the result of sortrank(x), original order
  # SRy1 is the result of sortrank(y1), y1=y[order(x)]
  # pre-compute SRx, SRy1 to avoid repeated sort and rank
  #
  n <- length(x1)
  ix <- SRx$ix      #order(x)
  rankx <- SRx$r    #ranks of original sample x

  ## ranks and order vector for this permutation of sample y1
  iy1 <- SRy1$ix       #order(y1)
  ranky1 <- SRy1$r     #rank(y1)

  ## the partial sums in the formula g_1
  psumsy1 <- (cumsum(as.numeric(z1[iy1])) - z1[iy1])[ranky1]
  psumsx1 <- cumsum(as.numeric(z1)) - z1

  gamma1 <- Btree_sum(y=ranky1, z=z1)   #y1 replaced by rank(y1)
  g <- sum(z1) - z1 - 2 * psumsx1 - 2 * psumsy1 + 4 * gamma1
  g <- g[rankx]
}

.rowSumsDist1 <- function(x, Sx = NULL) {
  ## for univariate samples, equivalent to rowSums(as.matrix(dist(x)))
  ## but much faster
  ## Sx is a sortrank object usually pre-computed here
  ## x is the data vector, Sx$x is sort(x)
  if (is.null(Sx))
  Sx <- sortrank(x)
  n <- length(x)
  r <- Sx$r  #ranks
  z <- Sx$x  #ordered sample x
  psums1 <- (cumsum(as.numeric(z)) - z)[r]
  (2*(r-1)-n)*x + sum(x) - 2*psums1
}
