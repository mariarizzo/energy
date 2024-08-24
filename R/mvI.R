indep.test<-
function(x, y, method = c("dcov","mvI"), index = 1, R) {
    # two energy tests for multivariate independence
    type <- match.arg(method)
    if (type == "dcov")
        return(dcov.test(x, y, index, R)) else
    if (type == "mvI")
        return(mvI.test(x, y, R))
}

mvI <-
function(x, y) {
    # energy statistic for multivariate independence
    # returns dependence coefficient I_n
    n <- NROW(x)
    m <- NROW(y)
    Dx <- .arg2dist.matrix(x)
    Dy <- .arg2dist.matrix(y)
    return(Istat(Dx, Dy)) #Rcpp
}

mvI.test<-
  function(x, y, R) {
    # an energy test for multivariate independence
    # not based on dCov or dCor
    n <- NROW(x)
    m <- NROW(y)
    if (n != m || n < 2) stop("Sample sizes must agree")
    Dx <- .arg2dist.matrix(x)
    Dy <- .arg2dist.matrix(y)
    stats <- Istats(Dx, Dy, R)
    stat <- n * stats[1]^2
    est <- stats[1]
    names(est) <- "I"
    names(stat) <- "n I^2"
    dataname <- paste("x (",n," by ",ncol(x), "), y(",n," by ", ncol(y), "), replicates ", R, sep="")
    if (R > 0) {
      p.value = (1 + sum(stats[-1] > est)) / (R+1)
    } else {
      p.value = NA
    }
    e <- list(
      method = "mvI energy test of independence",
      statistic = stat,
      estimate = est,
      replicates = stats[-1],
      p.value = p.value,
      data.name = dataname)
    class(e) <- "htest"
    e
  }
