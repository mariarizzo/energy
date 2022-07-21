pdcov.test <- function(x, y, z, R) {
  if (missing(R)) R <- 0
  
  Dx <- .arg2dist.matrix(x)
  Dy <- .arg2dist.matrix(y)
  Dz <- .arg2dist.matrix(z)
  n <- nrow(Dx)
  Pxz <- projection(Dx, Dz)  #U-center and compute projections
  Pyz <- projection(Dy, Dz)

  #PxzU <- U_center(Pxz)  #not necessary, because of invariance
  #PyzU <- U_center(Pyz)

  teststat <- n * U_product(Pxz, Pyz)
  ## calc. pdcor
  den <- sqrt(U_product(Pxz, Pxz) * U_product(Pyz, Pyz))
  if (den > 0.0) {
    estimate <- teststat / (n * den)
  } else estimate <- 0.0
  bootfn <- function(Pxz, i, Pyz) {
    # generate the permutation replicates of dcovU(Pxz, Pyz)
    # PxzU and PyzU are the U-centered matrices
    U_product(Pxz[i, i], Pyz)  #RcppExports
  }

  if (R > 0 && den > 0.0) {
    reps <- replicate(R, expr= {
      i <- sample(1:n)
      bootfn(Pxz, i, Pyz=Pyz)
    })

    replicates <- n * reps
    pval <- (1 + sum(replicates > teststat)) / (1 + R)
    #df <- n * (n-3) / 2 - 2
  } else {
    pval <- NA
    replicates <- NA
  }
  
  dataname <- paste("replicates ", R, sep="")
  if (! R>0)
    dataname <- "Specify R>0 replicates for a test"
  condition <- (den > 0.0)
  names(estimate) <- "pdcor"
  names(teststat) <- "n V^*"
  e <- list(
    call = match.call(),
    method = paste("pdcov test", sep = ""),
    statistic = teststat,
    estimate = estimate,
    p.value = pval,
    n = n,
    replicates = replicates,
    condition = condition,
    data.name = dataname)
  class(e) <- "htest"
  return(e)
}


pdcor.test <- function(x, y, z, R) {
  ## x, y, z must be dist. objects or data matrices (no dist matrix)
  ## all required calc. done in pdcov.test
  if (missing(R)) R <- 0
  result <- pdcov.test(x, y, z, R=R)

  if (result$condition) {
    ## if (A*A)(B*B) > 0
    nRootV <- result$statistic / result$estimate
    pdcor_reps <- result$replicates / nRootV
  } else pdcor_reps <- NA
  
  e <- list(
    call = match.call(),
    method = paste("pdcor test", sep = ""),
    statistic = result$estimate,
    estimate = result$estimate,
    p.value = result$p.value,
    n = result$n,
    replicates = pdcor_reps,
    condition = result$condition,
    data.name = result$data.name)
  class(e) <- "htest"
  return(e)
}

