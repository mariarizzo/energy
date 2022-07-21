poisson.tests <-
function(x, R, test="all") {
  # parametric bootstrap tests of Poisson distribution
  # poisson.e is the energy GOF statistic
  # poisson.m is the mean distance statistic
  # (not related to the test stats::poisson.test)
  if (!is.integer(x) || any(x < 0)) {
    warning("sample must be non-negative integers")
    return(NULL)
  }
  test <- tolower(test)
  poisson.stats <- function(x) {
    c(poisson.m(x), poisson.e(x))
  }
  stat <- switch(test,
                 "m" = poisson.m,
                 "e" = poisson.e,
                 poisson.stats)
  
  method <- switch(test, 
                   m=c("M-CvM","M-AD"), 
                   e="Energy", 
                   c("M-CvM","M-AD","Energy"))
  method <- paste(method, " test", sep="")
  n <- length(x)
  lambda <- mean(x)
  if (missing(R) || is.null(R)) {
    R <- 0
    message("Specify R > 0 replicates for MC test")
  }

  bootobj <- boot::boot(x, statistic = stat, R = R, 
                   sim = "parametric",
                   ran.gen = function(x, y) {rpois(n, lambda)})
  
  N <- length(bootobj$t0)
  p <- rep(NA, times=N)
  if (R > 0) {
    for (i in 1:N) {
    p[i] <- 1 - mean(bootobj$t[,i] < bootobj$t0[i])
    }
  }
  
  # a data frame, not an htest object  
  # comparable to broom::tidy on an htest object
  RVAL <- data.frame(estimate=lambda, statistic=bootobj$t0,
                       p.value=p, method=method)
  return(RVAL)
}

poisson.mtest <-
function(x, R=NULL) {
  if (is.null(R)) R <- 0
  rval <- poisson.tests(x, R, test="M")
  DNAME <- paste(deparse1(substitute(x)), "replicates: ", R)
  stat <- rval$statistic[1]
  names(stat) <- "M-CvM"
    e <- list(
      method = paste("Poisson M-test", sep = ""),
      statistic = stat,
      p.value = rval$p.value[1],
      data.name = DNAME,
      estimate = rval$estimate[1])
    class(e) <- "htest"
    e
}

poisson.etest <- 
function(x, R=NULL) {
  if (is.null(R)) R <- 0
    rval <- poisson.tests(x, R, test="E")
    DNAME <- paste(deparse1(substitute(x)), "replicates: ", R)
    stat <- rval$statistic
    names(stat) <- "E"
    e <- list(
      method = paste("Poisson E-test", sep = ""),
      statistic = stat,
      p.value = rval$p.value,
      data.name = paste("replicates: ", R,  sep=""),
      estimate = rval$estimate)
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
  stats <- .poisMstat(x)
  names(stats) <- c("M-CvM", "M-AD")
  return(stats)
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
  stat <- n * (2 * mean(meanvec) - EXX - meanxx)
  names(stat) <- "E"
  return(stat)
}


