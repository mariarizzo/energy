mvnorm.test <- function(x, R) {
  # parametric bootstrap E-test for multivariate normality
  if (missing(R)) {
    method = "Energy test of multivariate normality: (Specify R > 0 for MC test)"
    R <- 0
  } else {
    method = "Energy test of multivariate normality: estimated parameters"
  }
  
  if (is.vector(x) || NCOL(x)==1) {
    n <- NROW(x)
    d <- 1
    bootobj <- boot::boot(x, statistic = normal.e, R = R, sim = "parametric",
                          ran.gen = function(x, y) {
                            return(rnorm(n))
                          })
  } else {
    n <- nrow(x)
    d <- ncol(x)
    bootobj <- boot::boot(x, statistic = mvnorm.e, R = R, sim = "parametric",
                          ran.gen = function(x, y) {
                            return(matrix(rnorm(n * d), nrow = n, ncol = d))
                          })
  }
  if (R > 0)
    p <- 1 - mean(bootobj$t < bootobj$t0) else p <- NA

  names(bootobj$t0) <- "E-statistic"
  e <- list(statistic = bootobj$t0, p.value = p,
            method = method,
            data.name = paste("x, sample size ", n, ", dimension ", d, ", replicates ",
                              R, sep = ""))
  class(e) <- "htest"
  e
}

mvnorm.etest <- function(x, R) {
  return(mvnorm.test(x, R))
}

mvnorm.e <- function(x) {
  # E-statistic for multivariate normality
  if (is.vector(x) || NCOL(x)==1)
    return(normal.e(x))
  n <- nrow(x)
  d <- ncol(x)
  if (n < 2)
    return(normal.e(x))
  # subtract column means and compute S^(-1/2)
  z <- scale(x, scale = FALSE)
  ev <- eigen(var(x), symmetric = TRUE)
  P <- ev$vectors
  lambda <- ev$values
  D <- diag(d)
  diag(D) <- 1 / sqrt(lambda)
  y <- z %*% (P %*% D %*% t(P))
  if (any(!is.finite(y)))
    return(NA)
  return(mvnEstat(y))
}

normal.e <- function(x) {
  ## Case 4: unknown parameters
  x <- as.vector(x)
  n <- length(x)
  y <- (x - mean(x)) / sd(x)
  y <- sort(y)
  K <- seq(1 - n, n - 1, 2)
  if (y[1] == y[n])
    return(NA)
  return(2 * (sum(2 * y * pnorm(y) + 2 * dnorm(y)) -
                n/sqrt(pi) - mean(K * y)))
}


normal.test <- function(x, method=c("mc", "limit"), R) {
  ## implements the test for for d=1 
  ## Case 4: composite hypothesis
  method <- match.arg(method)
  estimate <- c(mean(x), sd(x))
  names(estimate) <- c("mean", "sd")
  
  if (method == "mc") {
    ## Monte Carlo approach
    if (missing(R)) R <- 0
    e <- energy::mvnorm.etest(x, R=R)
    e$method <- "Energy test of normality"
    e$method <- ifelse(R > 0,
      paste0(e$method,": estimated parameters"),
        paste0(e$method, "  (Specify R > 0 for MC test)"))
    e$estimate <- estimate
    return(e)
  }
  
  ## implement test using asymptotic distribution for p-value
  if (!is.numeric(x) || (!is.vector(x) && NCOL(x) > 1)) {
    warning("x must be a numeric vector")
    return (NA)
  } else {
    x <- as.vector(x, mode="numeric")
  }

  n <- length(x)
  t0 <- normal.e(x)
  names(t0) <- "statistic"
  
  ## load pre-computed eigenvalues
  EVnormal <- NULL
  load("./data/EVnormal.rda")
  ev <- EVnormal[, "Case4"]

  if (requireNamespace("CompQuadForm", quietly=TRUE)) {
    p <- CompQuadForm::imhof(t0, ev)$Qq
      } else {
        warning("limit distribution method requires CompQuadForm package for p-value")
        p <- NA
      }
  estimate <- c(mean(x), sd(x))
  names(estimate) <- c("mean", "sd")
  e <- list(statistic = t0, p.value = p,
            method = paste("Energy test of normality: limit distribution"),
            estimate = estimate,
            data.name = "Case 4: composite hypothesis, estimated parameters")
  class(e) <- "htest"
  e
}

