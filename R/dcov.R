dcov.test <-
function(x, y, index = 1.0, R = NULL, test = c("perm","limit")) {
    dataname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    test <- match.arg(test)
    if (test == "perm") {
      ## check for valid number of replicates R
      method <- "Specify the number of replicates R (R > 0) for a permutation test"
      if (! is.null(R)) {
        R <- floor(R)
        if (R < 1) R <- 0
        if (R > 0) 
          method <- "dCov permutation test of independence"
      } else {
        R <- 0
      }
    }
    # distance covariance test for multivariate independence
    if (!(class(x) == "dist")) x <- dist(x)
    if (!(class(y) == "dist")) y <- dist(y)
    x <- as.matrix(x)
    y <- as.matrix(y)
    dst <- TRUE
    n <- nrow(x)
    m <- nrow(y)
    if (n != m) stop("Sample sizes must agree")
    if (! (all(is.finite(c(x, y)))))
      stop("Data contains missing or infinite values")
    
    if (test == "limit") {
      return (.dcov.limit.test(x, y, index = index))
    }

    if (test == "perm") {    
      stat <- dcorr <- reps <- 0
      dcov <- rep(0, 4)
      if (R > 0) reps <- rep(0, R)
      pval <- 1
      dims <- c(n, ncol(x), ncol(y), dst, R)
  
      # dcov = [dCov,dCor,dVar(x),dVar(y)]
      a <- .C("dCOVtest",
              x = as.double(t(x)),
              y = as.double(t(y)),
              byrow = as.integer(TRUE),
              dims = as.integer(dims),
              index = as.double(index),
              reps = as.double(reps),
              DCOV = as.double(dcov),
              pval = as.double(pval),
              PACKAGE = "energy")
      # test statistic is n times the square of dCov statistic
      stat <- n * a$DCOV[1]^2
      dcorr <- a$DCOV
      V <- dcorr[[1]]
      pval <- ifelse (R < 1, NA, a$pval)
      estimates <- dcorr
      replicates <- n* a$reps^2
    }
    
    names(stat) <- "nV^2"
    names(V) <- "dCov"
    if (!isTRUE(all.equal(index, 1.0)))
      dataname <- paste(dataname, "index=", index)
    e <- list(call = match.call(),
        statistic = stat,
        method = method,
        estimate = V,
        estimates = estimates,
        p.value = pval,
        replicates = replicates,
        n = n,
        data.name = dataname)
    class(e) <- "htest"
    return(e)
}

.dcov.limit.test <- 
function(Dx, Dy, index = 1.0, 
         method = "dCov test of independence using limit distribution",
         max.eigen = 500, tol = .01, 
         show.warning = TRUE) {
  ## Dx and Dy are the distance matrices of the complete sample
  ## dcov test of independence by estimate the limit distribution of 
  ## the test statistic n V_n^2(x, y) and return htest object
  ## requires package CompQuadForm
  n <- nrow(Dx)
  if (!isTRUE(all.equal(index, 1.0))) {
    Dx <- Dx^index
    Dy <- Dy^index
  } 
  
  A <- D_center(Dx)
  B <- D_center(Dy)
  stat <- sum(A * B) / n  #n V_n^2, the test statistic
  V <- sqrt(stat / n)     #dcov
  dvx <- sqrt(mean(A * A))
  dvy <- sqrt(mean(B * B))
  dv <- dvx * dvy
  dCor <- ifelse (dv > 0, V / dv, 0) 
  estimates <- c(V, dCor, dvx, dvy)
  
  if (n <= max.eigen) {
    H <- A * B / n
  } else {
    I <- sample(1:n, size = max.eigen, replace = FALSE)
    A <- A[I, I]
    B <- B[I, I]
    H <- A * B / max.eigen
  }
  
  lambda <- eigen(H, symmetric = TRUE, only.values = TRUE)$values

  err <- tail(lambda, 1) / sum(lambda)
  if (err > tol) {
    if (show.warning == TRUE)
      warning(paste("Relative error in limit:", err, "Use permutation test."))
  }
  pval <- NA
  abserr <- NA
  if (requireNamespace("CompQuadForm", quietly = TRUE)) {
    pcalc <- CompQuadForm::imhof(stat, lambda)
    pval <- pcalc$Qq
    abserr <- pcalc$abserr
  } else {
    warning("package CompQuadForm required to compute probability")
  }
  names(stat) <- "nV^2"
  names(V) <- "dCov"
  dataname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  e <- list(call = match.call(),
            statistic = stat,
            method = method,
            estimate = V,
            estimates = estimates,
            p.value = pval,
            replicates = NA,
            n = n,
            err = err,
            abserr = abserr,
            eigenvalues = lambda,
            data.name = dataname)
  class(e) <- "test"
  return(e)
}

dcor.test <-
  function(x, y, index=1.0, R, test = c("perm", "limit")) {
    # distance correlation test for multivariate independence
    # like dcov.test but using dcor as the test statistic
    dataname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    test = match.arg(test)
    if (is.null(R)) R <- 0
    R <- ifelse(R > 0, floor(R), 0)
    RESULT <- dcov.test(x, y, index=1.0, R, test = test)
    
    # this test statistic is n times the square of dCov statistic
    DCOVteststat <- RESULT$statistic
    DCOVreplicates <- RESULT$replicates
    
    # RESULT$estimates = [dCov,dCor,dVar(x),dVar(y)]
    # dVar are invariant under permutation of sample indices
    
    DCORteststat <- RESULT$estimates[2]
    dvarX <- RESULT$estimates[3]
    dvarY <- RESULT$estimates[4]
    n <- RESULT$n
    DCORreps <- sqrt(DCOVreplicates / n) / sqrt(dvarX * dvarY)
    
    p.value <- NA
    if (test == "perm" && R > 0)
      p.value <- (1 + sum(DCORreps >= DCORteststat)) / (1 + R) 
    if (test == "limit")
      p.value <- RESULT$p.value
    
    names(DCORteststat) <- "dCor"
    method <- "dCor test of independence"
    if (test == "limit") 
      paste(method, "using limit distribution")
    
    e <- list(
      method = method,
      statistic = DCORteststat,
      estimates = RESULT$estimates,
      p.value = p.value,
      replicates = DCORreps,
      n = n,
      data.name = dataname)
    class(e) <- "htest"
    return(e)
  }


.dcov <-
function(x, y, index=1.0) {
    # distance covariance statistic for independence
    # dcov = [dCov,dCor,dVar(x),dVar(y)]   (vector)
    # this function calls C function for computing dCov
    # it is called by the dcov and dcor functions
    if (!(class(x) == "dist")) x <- dist(x)
    if (!(class(y) == "dist")) y <- dist(y)
    x <- as.matrix(x)
    y <- as.matrix(y)
    dst <- TRUE
    n <- nrow(x)
    m <- nrow(y)
    if (n != m) stop("Sample sizes must agree")
    if (! (all(is.finite(c(x, y)))))
        stop("Data contains missing or infinite values")
    dims <- c(n, NCOL(x), NCOL(y), dst)
    idx <- 1:dims[1]
    DCOV <- numeric(4)
    a <- .C("dCOV",
            x = as.double(t(x)),
            y = as.double(t(y)),
            byrow = as.integer(TRUE),
            dims = as.integer(dims),
            index = as.double(index),
            idx = as.double(idx),
            DCOV = as.double(DCOV),
            PACKAGE = "energy")
    return(a$DCOV)
}

dcov <-
function(x, y, index=1.0) {
    # distance covariance statistic for independence
    return(.dcov(x, y, index)[1])
}

dcor <-
function(x, y, index=1.0) {
    # distance correlation statistic for independence
    return(.dcov(x, y, index)[2])
}


DCOR <-
function(x, y, index=1.0) {
    # distance covariance and correlation statistics
    # alternate method, implemented in R without .C call
    # this method is usually slower than the C version
    if (!(class(x) == "dist")) x <- dist(x)
    if (!(class(y) == "dist")) y <- dist(y)
    x <- as.matrix(x)
    y <- as.matrix(y)
    n <- nrow(x)
    m <- nrow(y)
    if (n != m) stop("Sample sizes must agree")
    if (! (all(is.finite(c(x, y)))))
        stop("Data contains missing or infinite values")
    if (index < 0 || index > 2) {
        warning("index must be in [0,2), using default index=1")
        index=1.0}

    stat <- 0
    dims <- c(n, ncol(x), ncol(y))

    Akl <- function(x) {
        d <- as.matrix(x)^index
        m <- rowMeans(d)
        M <- mean(d)
        a <- sweep(d, 1, m)
        b <- sweep(a, 2, m)
        return(b + M)
    }

    A <- Akl(x)
    B <- Akl(y)
    dCov <- sqrt(mean(A * B))
    dVarX <- sqrt(mean(A * A))
    dVarY <- sqrt(mean(B * B))
    V <- sqrt(dVarX * dVarY)
    if (V > 0)
      dCor <- dCov / V else dCor <- 0
    return(list(dCov=dCov, dCor=dCor, dVarX=dVarX, dVarY=dVarY))
}
