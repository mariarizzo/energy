dcovV.test <- function(x, y, max.eigen = 250) {
  # distance covariance test for multivariate independence
  # estimate eigenvalues of limit distribution for p-value
  n <- NROW(x)
  if (is.vector(x))
    x1 <- matrix(x, ncol=1) else x1 <- x
    if (is.vector(y))
      y1 <- matrix(y, ncol=1) else y1 <- y
  stats <- dcov_UV(x1, y1)
  nV <- n * stats["dcovV"]
    
  method <- "dCov V test of independence (use limit)"
  dataname <- paste(n, deparse(substitute(x)), "and", deparse(substitute(y)))
  rval <- .dcovV(nV, x, y, max.eigen)
  estimate = stats["dcovV"]
  statistic <- nV
  names(estimate) <- "dCov V statistic"
  names(statistic) <- "n V^2"
  e <- list(
    statistic = statistic,
    method = method,
    estimate = estimate,
    eigenvalues = rval$lambda,
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
  if (is.vector(x))
    x1 <- matrix(x, ncol=1) else x1 <- x
  if (is.vector(y))
    y1 <- matrix(y, ncol=1) else y1 <- y
  stats <- dcov_UV(x1, y1)
  nV <- n * stats["dcovV"]   #need V for eigenvalues
  nU <- n * stats["dcovU"]
  
  method <- "dCov U test of independence (use limit)"
  dataname <- paste(n, deparse(substitute(x)), "and", deparse(substitute(y)))
  statistic <- nU
  rval <- .dcovU_prob(nU, x, y, max.eigen)
  estimate <- stats["dcovU"]
  names(estimate) <- "dcov U statistic"
  names(statistic) <- "n U"
  e <- list(
    statistic = statistic,
    method = method,
    estimate = estimate, 
    eigenvalues = rval$lambda,
    p.value = rval$p,
    n = n,
    data.name = dataname)
  class(e) <- "htest"
  return(e)
}


.dcovV <- function(V, x, y, max.eigen = 250) {
  # estimate eigenvalues of dcov V-statistic limit distribution for p-value
  n <- NROW(x)
  x1 <- x
  y1 <- y
  if (n > max.eigen) {
    i <- sample(1:n, replace=FALSE, size=max.eigen)
    n <- max.eigen
    if (is.vector(x))
      x1 <- matrix(x1[i], ncol=1) else x1 <- x1[i,]
    if (is.vector(y))
      y1 <- matrix(y1[i], ncol=1) else y1 <- y1[i,]
  }
  
  Dx <- as.matrix(dist(x1))
  Dy <- as.matrix(dist(y1))
  A <- D_center(Dx)
  B <- D_center(Dy)
  H <- A * B / n
  lambda <- svd(H, nu=0, nv=0)$d
  ## informal check for convergence
  err <- mean(tail(lambda, 3))
  if (err > .Machine$double.eps^.2)
    warning("sample size too small for limit method; use permutation test")
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
      x1 <- x[i] else x1 <- x1[i,]
    if (is.vector(y))
      y1 <- y1[i] else y1 <- y1[i,]
    }
  Dx <- as.matrix(dist(x1))
  Dy <- as.matrix(dist(y1))
  A <- D_center(Dx)
  B <- D_center(Dy)
  H <- A * B / n
  lambda <- svd(H, nu=0, nv=0)$d
  err <- mean(tail(lambda, 3))
  if (err > .Machine$double.eps^.2)
    warning("sample size too small for limit method; use permutation test")
  qU <- U + sum(rev(lambda))
  p <- NA
  if (requireNamespace("CompQuadForm", quietly=TRUE)) {
    ## informal check for convergence
    p <- CompQuadForm::imhof(qU, lambda)$Qq
  } else {
    warning("package CompQuadForm required to compute probability")
  }
  return (list(p = p, lambda = lambda))
}

