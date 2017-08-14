pdcov.test <- function(x, y, z, R) {
  ## x, y, z must be dist. objects or data matrices (no dist matrix)
  if (!(class(x) == "dist")) x <- dist(x)
  if (!(class(y) == "dist")) y <- dist(y)
  if (!(class(z) == "dist")) z <- dist(z)
  Dx <- as.matrix(x)
  Dy <- as.matrix(y)
  Dz <- as.matrix(z)
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

  reps <- replicate(R, expr= {
    i <- sample(1:n)
    bootfn(Pxz, i, Pyz=Pyz)
  })

  replicates <- n * reps
  pval <- (1 + sum(replicates > teststat)) / (1 + R)
  #df <- n * (n-3) / 2 - 2
  dataname <- paste("replicates ", R, sep="")

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
    data.name = dataname)
  class(e) <- "htest"
  return(e)
}


pdcor.test <- function(x, y, z, R) {
  ## x, y, z must be dist. objects or data matrices (no dist matrix)
  ## all required calc. done in pdcov.test

  result <- pdcov.test(x, y, z, R)

  reps <- result$replicates
  teststat <- result$estimate
  estimate <- result$estimate
  n <- result$n
  if (estimate > 0.0) {
    nRootV <- result$statistic / result$estimate
    pdcor_reps <- reps / nRootV
  } else {
    pdcor_reps <- reps
  }

  names(estimate) <- names(teststat) <- "pdcor"
  pval <- (1 + sum(pdcor_reps > teststat)) / (1 + R)

  e <- list(
    call = match.call(),
    method = paste("pdcor test", sep = ""),
    statistic = teststat,
    estimate = estimate,
    p.value = pval,
    n = n,
    replicates = pdcor_reps,
    data.name = result$dataname)
  class(e) <- "htest"
  return(e)
}

