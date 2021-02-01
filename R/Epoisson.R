poissonGOF.test <-
function(x, R, method=c("M","E")) {
    # parametric bootstrap energy tests of Poisson distribution
    # poisson.e is the energy GOF statistic
    # poisson.m is the mean distance statistic
    n <- length(x)
    lambda <- mean(x)
    method <- match.arg(method)
    if (method == "M") {
      bootobj <- boot::boot(x, statistic = poisson.m, R = R, 
                   sim = "parametric",
                   ran.gen = function(x, y) {rpois(n, lambda)})
      test <- "Mean distance test of Poisson distribution"
    } else {
      bootobj <- boot::boot(x, statistic = poisson.e, R = R, 
                            sim = "parametric",
                            ran.gen = function(x, y) {rpois(n, lambda)})
      test <- "Energy test of Poisson distribution"
    }
    if (R > 0)
      p <- 1 - mean(bootobj$t < bootobj$t0) else p <- NA
    names(bootobj$t0) <- "test statistic"
    names(lambda) <- "mean"
    e <- list(
        method = test,
        statistic = bootobj$t0,
        p.value = p,
        data.name = paste("sample size ", n, ", replicates ", R, sep=""),
        estimate = lambda)
    class(e) <- "htest"
    e
}

poisson.mtest <- function(x, R) {
  return(poissonGOF.test(x, R=R, method="M"))
}

poisson.m <-
function(x) {
    # mean distance statistic for Poissonity
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

