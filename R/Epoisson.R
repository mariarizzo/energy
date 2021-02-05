poissonGOF.test <-
function(x, R) {
  # parametric bootstrap energy tests of Poisson distribution
  # poisson.e is the energy GOF statistic
  # poisson.m is the mean distance statistic
  method <- c("Poisson M-test", "Poisson E-test")
  stat <- .poisson.stats

  if (any(!is.integer(x)) || any(x < 0)) {
      warning("sample must be non-negative integers")
      return(NULL)
    }
  n <- length(x)
  lambda <- mean(x)
  if (missing(R)) {
    R <- 0
    warning("Specify R > 0 replicates for MC test")
  }
    
  bootobj <- boot::boot(x, statistic = stat, R = R, 
                   sim = "parametric",
                   ran.gen = function(x, y) {rpois(n, lambda)})

  if (R > 0) {
    p1 <- 1 - mean(bootobj$t[,1] < bootobj$t0[1])
    p2 <- 1 - mean(bootobj$t[,2] < bootobj$t0[2])
    p <- c(p1, p2)
    } else {
      p <- c(NA, NA)
    }
  
  # a data frame, not an htest object  
  # comparable to broom::tidy on an htest object
  RVAL <- data.frame(estimate=lambda, statistic=bootobj$t0,
                       p.value=p, method=method)
  return(RVAL)
}

poisson.mtest <- function(x, R) {
  # for backward compatibility
  poissonM.test(x, R)
}

poissonM.test <- 
  function(x, R) {
    # parametric bootstrap mean distance test of Poisson distribution
    if (any(!is.integer(x)) || any(x < 0)) {
      warning("sample must be non-negative integers")
      return(NULL)
    }
    n <- length(x)
    lambda <- mean(x)
    bootobj <- boot::boot(x, statistic = poisson.m, R = R, sim = "parametric",
                          ran.gen = function(x, y) {rpois(n, lambda)})
    if (R > 0)
      p <- 1 - mean(bootobj$t < bootobj$t0) else p <- NA
    names(bootobj$t0) <- "test statistic"
    names(lambda) <- "mean"
    e <- list(
      method = paste("Poisson M-test", sep = ""),
      statistic = bootobj$t0,
      p.value = p,
      data.name = paste("sample size ", n, ", replicates ", R, sep=""),
      estimate = lambda)
    class(e) <- "htest"
    e
  }


poissonE.test <-
  function(x, R) {
    # parametric bootstrap mean distance test of Poisson distribution
    if (any(!is.integer(x)) || any(x < 0)) {
      warning("sample must be non-negative integers")
      return(NULL)
    }
    n <- length(x)
    lambda <- mean(x)
    bootobj <- boot::boot(x, statistic = poisson.e, R = R, sim = "parametric",
                          ran.gen = function(x, y) {rpois(n, lambda)})
    if (R > 0)
      p <- 1 - mean(bootobj$t < bootobj$t0) else p <- NA
    names(bootobj$t0) <- "test statistic"
    names(lambda) <- "mean"
    e <- list(
      method = paste("Poisson E-test", sep = ""),
      statistic = bootobj$t0,
      p.value = p,
      data.name = paste("sample size ", n, ", replicates ", R, sep=""),
      estimate = lambda)
    class(e) <- "htest"
    e
  }

poisson.m <-
function(x) {
    # mean distance statistic for Poissonity
    if (any(!is.integer(x)) || any(x < 0)) {
      warning("sample must be non-negative integers")
      return(NULL)
    }
    n <- length(x)
    stat <- 0
    e <- .C("poisMstat",
            x = as.integer(x),
            nx = as.integer(n),
            stat = as.double(stat),
            PACKAGE = "energy")$stat
    e
}

poisson.e <-
function(x) {
  # energy GOF statistic for Poissonity
  if (any(!is.integer(x)) || any(x < 0)) {
    warning("sample must be non-negative integers")
    return(NULL)
  }
  lambda <- mean(x)
  n <- length(x)
  
  ## E|y-X| for X Poisson(lambda) (vectorized)
  Px <- ppois(x, lambda)
  Px1 <- ppois(x-1, lambda)
  meanvec <- 2*x*Px - 2*lambda*Px1 + lambda - x
  
  ## second mean E|X-X'|
  a <- 2 * lambda
  EXX <- a * exp(-a) * (besselI(a, 0) + besselI(a, 1))
  
  ## third mean = sum_{i,j} |x_i - x_j| / n^2
  K <- seq(1 - n, n - 1, 2)
  y <- sort(x)
  meanxx <- 2 * sum(K * y) / n^2
  
  e <- n * (2 * mean(meanvec) - EXX - meanxx)
}

.poisson.stats <- function(x, lambda) {
  c(poisson.m(x), poisson.e(x))
}
