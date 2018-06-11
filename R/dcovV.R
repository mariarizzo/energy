dcovV.test <- function(x, y, max.eigen = 250) {
  # distance covariance test for multivariate independence
  # estimate eigenvalues of limit distribution for p-value
  n <- NROW(x)
  if (n < 50) {
    ## sample size too small for the limit distr method
    return(dcov.test(x, y, R=499))
  }
  
  method <- "dCov V test of independence: estimated eigenvalues"
  dataname <- paste(n, deparse(substitute(x)), "and", deparse(substitute(y)))
  rval <- .dcovV(x, y, max.eigen)
  nV <- rval$nV
  estimate = sqrt(rval$nV/n)
  statistic <- nV
  names(estimate) <- "dCov"
  names(statistic) <- "nV^2"
  e <- list(
    statistic = statistic,
    method = method,
    estimate = estimate,
    estimates = rval$lambda,
    p.value = rval$p,
    n = n,
    data.name = dataname)
    class(e) <- "htest"
    return(e)
}


dcovU.test <- function(x, y, max.eigen = 250) {
  # distance covariance test for multivariate independence
  # estimate eigenvalues of limit distribution for p-value
  n <- NROW(x)
  if (n < 50) {
    ## sample size too small for the limit distr method
    return(dcov.test(x, y, R=499))
  }
  
  method <- "dCov U test of independence: estimated eigenvalues"
  dataname <- paste(n, deparse(substitute(x)), "and", deparse(substitute(y)))
  nU <- n * dcovU(x, y)
  rval <- .dcovU_prob(nU, x, y, max.eigen)
  statistic <- nU
  estimate <- nU / n
  names(estimate) <- "U"
  names(statistic) <- "n U"
  e <- list(
    statistic = statistic,
    method = method,
    estimate = estimate, 
    estimates = rval$lambda,
    p.value = rval$p,
    n = n,
    data.name = dataname)
  class(e) <- "htest"
  return(e)
}


.dcovV <- function(x, y, max.eigen = 250) {
  # compute dcov V-statistics and
  # estimate eigenvalues of limit distribution for p-value
  n <- NROW(x)
  if (n < 50) return (NA)
  x1 <- x
  y1 <- y
  if (n > max.eigen) {
    i <- sample(1:n, replace=FALSE, size=max.eigen)
    n <- max.eigen
    if (is.vector(x))
      x1 <- x[i] else x1 <- x[i,]
      if (is.vector(y))
        y1 <- y[i] else y1 <- y[i,]
  }
  
  Dx <- as.matrix(dist(x1))
  Dy <- as.matrix(dist(y1))
  A <- D_center(Dx)
  B <- D_center(Dy)
  H <- A * B / n

  lambda <- svd(H, nu=0, nv=0)$d
  
  if (NROW(x) > n) {
    V <- n * dcov(x, y)^2    #compute V using complete sample
    } else {
      V <- sum(H)            #n V^2 = n dcov^2
    }
  
  p <- NA
  if (requireNamespace("CompQuadForm", quietly=TRUE)) {
    p <- CompQuadForm::imhof(V, lambda)$Qq
  } else {
    warning("package CompQuadForm required to compute probability")
  }
  return (list(nV=V, pval=p, lambda=lambda))
}

.dcovU_prob <- function(U, x, y, max.eigen = 250) {
  # compute Pr(U > q) for dcov U-statistic using
  # estimated eigenvalues of V
  n <- NROW(x)
  x1 <- x
  y1 <- y
  if (n < 50) return (NA)
  if (n > max.eigen) {
    i <- sample(1:n, replace=FALSE, size=max.eigen)
    n <- max.eigen
    if (is.vector(x))
      x1 <- x1[i] else x1 <- x1[i,]
    if (is.vector(y))
      y1 <- y1[i] else y1 <- y1[i,]
    }
  Dx <- as.matrix(dist(x1))
  Dy <- as.matrix(dist(y1))
  A <- D_center(Dx)
  B <- D_center(Dy)
  H <- A * B / n
  lambda <- svd(H, nu=0, nv=0)$d
  qU <- U + sum(rev(lambda))
  p <- NA
  if (requireNamespace("CompQuadForm", quietly=TRUE)) {
    p <- CompQuadForm::imhof(qU, lambda)$Qq
  } else {
    warning("package CompQuadForm required to compute probability")
  }
  return (list(p = p, lambda = lambda))
}

