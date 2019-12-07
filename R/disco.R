 ### disco tests - implementation of DIStance COmponents methods in:
 ###
 ### Rizzo, M.L. and Szekely, G.J. (2010) "DISCO Analysis: A Nonparametric
 ### Extension of Analysis of Variance, Annals of Applied Statistics
 ###  Vol. 4, No. 2, 1034-1055.
 ###
 ### Sept 2010 parts of disco package merged into energy package
 ### this release supports one way models
 ### this version does not use the C library
 ###
 ### disco: computes the decomposition and test using F ratio
 ### disco.between: statistic and test using between component
 ### .disco1: internal computations for one factor
 ### .disco1stat, .disco1Bstat: internal for boot function
 ###
 ###



disco <- function(x, factors, distance = FALSE, index = 1, R,
                  method = c("disco", "discoB", "discoF")) {
  ## x is response or Euclidean distance matrix or dist() object factors
  ## is a matrix or data frame of group labels distance=TRUE if x is
  ## distance, otherwise FALSE index is the exponent on distance, in (0,2]
  ## R is number of replicates for test method: use F ratio (default) or
  ## between component (discoB) disco method is currently alias for discoF


  method <- match.arg(method)
  factors <- data.frame(factors)
  if (inherits(x, "dist")) distance <- TRUE
  if (method == "discoB")
    return(disco.between(x, factors = factors, distance = distance,
                         index = index, R = R))
  nfactors <- NCOL(factors)
  if (distance || inherits(x, "dist"))
    dst <- as.matrix(x) else dst <- as.matrix(dist(x))
  N <- NROW(dst)
  if (NCOL(dst) != N)
    stop("distance==TRUE but first argument is not distance")
  if (!isTRUE(all.equal(index, 1)))
    dst <- dst^index

  stats <- matrix(0, nfactors, 6)
  colnames(stats) <- c("Trt", "Within", "df1", "df2", "Stat", "p-value")

  for (j in 1:nfactors) {
    trt <- factors[, j]
    stats[j, 1:4] <- .disco1(trt = trt, dst = dst)
    if (R > 0) {
      b <- boot::boot(data = dst, statistic = .disco1stat, sim = "permutation",
                      R = R, trt = trt)
      stats[j, 5] <- b$t0
      stats[j, 6] <- (sum(b$t > b$t0) + 1)/(R + 1)
    } else {
      stats[j, 5] <- .disco1stat(dst, i = 1:nrow(dst), trt = trt)
      stats[j, 6] <- NA
    }
  }

  methodname <- "DISCO (F ratio)"
  dataname <- deparse(substitute(x))
  total <- sum(stats[1, 1:2])
  within <- total - sum(stats[, 1])
  Df.trt <- stats[, 3]
  factor.names <- names(factors)
  factor.levels <- sapply(factors, nlevels)
  sizes <- sapply(factors, tabulate)
  e <- list(call = match.call(), method = methodname,
            statistic = stats[, 5],
            p.value = stats[, 6],
            k = nfactors,
            N = N,
            between = stats[, 1],
            withins = stats[, 2],
            within = within,
            total = total,
            Df.trt = Df.trt,
            Df.e = nrow(dst) - sum(Df.trt) - 1,
            index = index, factor.names = factor.names,
            factor.levels = factor.levels,
            sample.sizes = sizes, stats = stats)
  class(e) <- "disco"
  e
}

disco.between <- function(x, factors, distance = FALSE, index = 1, R) {
  ## disco test based on the between-sample component similar to disco
  ## except that 'disco' test is based on the F ratio disco.between test
  ## for one factor (balanced) is asymptotically equivalent to k-sample E
  ## test (test statistics are proportional in that case but not in
  ## general).  x is response or Euclidean distance matrix or dist()
  ## object factors is a matrix or data frame of group labels
  ## distance=TRUE if x is distance, otherwise FALSE index is the exponent
  ## on distance, in (0,2]

  factors <- data.frame(factors)
  nfactors <- NCOL(factors)
  if (nfactors > 1)
    stop("More than one factor is not implemented in disco.between")
  if (distance || inherits(x, "dist"))
    dst <- as.matrix(x) else dst <- as.matrix(dist(x))
  N <- NROW(dst)
  if (NCOL(dst) != N)
    stop("distance==TRUE but first argument is not distance")
  if (!isTRUE(all.equal(index, 1)))
    dst <- dst^index

  trt <- factors[, 1]
  if (R > 0) {
    b <- boot::boot(data = dst, statistic = .disco1Bstat, sim = "permutation",
                    R = R, trt = trt)
    between <- b$t0
    reps <- b$t
    pval <- mean(reps >= between)
  } else {
    between <- .disco1Bstat(dst, i = 1:nrow(dst), trt = trt)
    pval <- NA
  }
  if (R == 0)
    return(between)

  methodname <- "DISCO (Between-sample)"
  dataname <- deparse(substitute(x))

  names(between) <- "DISCO between statistic"
  e <- list(call = match.call(), method = methodname, statistic = between,
            p.value = pval, data.name = dataname)

  class(e) <- "htest"
  e
}

.disco1 <- function(trt, dst) {
  ## dst is Euclidean distance matrix or power of it trt is the treatment,
  ## a factor

  trt <- factor(trt)
  k <- nlevels(trt)
  n <- tabulate(trt)
  N <- sum(n)
  total <- sum(dst)/(2 * N)
  y <- as.vector(dst[, 1])
  M <- model.matrix(y ~ 0 + trt)
  G <- t(M) %*% dst %*% M
  withins <- diag(G)/(2 * n)
  W <- sum(withins)
  B <- total - W
  c(B, W, k - 1, N - k)
}

.disco1stat <- function(dst, i, trt) {
  ## i is permuation vector supplied by bootstrap dst is Euclidean
  ## distance matrix or power of it trt is the treatment, a factor returns
  ## the disco 'F' ratio
  idx <- 1:nrow(dst)
  d <- .disco1(trt = trt[idx[i]], dst = dst)
  statistic <- (d[1]/d[3])/(d[2]/d[4])
}

.disco1Bstat <- function(dst, i, trt) {
  ## i is permuation vector supplied by bootstrap dst is Euclidean
  ## distance matrix or power of it trt is the treatment, a factor returns
  ## the between-sample component (for one factor)
  idx <- 1:nrow(dst)
  .disco1(trt = trt[idx[i]], dst = dst)[1]
}

print.disco <- function(x, ...) {
  k <- x$k
  md1 <- x$between/x$Df.trt
  md2 <- x$within/x$Df.e
  f0 <- x$statistic
  print(x$call)
  cat(sprintf("\nDistance Components: index %5.2f\n", x$index))
  cat(sprintf("%-20s %4s %10s %10s %10s %10s\n", "Source", "Df", "Sum Dist",
              "Mean Dist", "F-ratio", "p-value"))
  for (i in 1:k) {
    fname <- x$factor.names[i]
    cat(sprintf("%-20s %4d %10.5f %10.5f %10.3f %10s\n", fname, x$Df.trt[i],
                x$between[i], md1[i], f0[i], format.pval(x$p.value[i])))
  }
  cat(sprintf("%-20s %4d %10.5f %10.5f\n", "Within", x$Df.e, x$within,
              md2))
  cat(sprintf("%-20s %4d %10.5f\n", "Total", x$N - 1, x$total))
}


