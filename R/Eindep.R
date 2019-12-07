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
    x <- as.matrix(x)
    y <- as.matrix(y)
    n <- nrow(x)
    m <- nrow(y)
    if (n != m || n < 2) stop("Sample sizes must agree")
    if (! (all(is.finite(c(x, y)))))
        stop("Data contains missing or infinite values")

    stat <- 0
    dims <- c(n, ncol(x), ncol(y))

    e <- .C("indepE",
            x = as.double(t(x)),
            y = as.double(t(y)),
            byrow = as.integer(TRUE),
            dims = as.integer(dims),
            stat = as.double(stat),
            PACKAGE = "energy")
    sqrt(e$stat)
}

mvI.test<-
function(x, y, R) {
    # energy test for multivariate independence
    x <- as.matrix(x)
    y <- as.matrix(y)
    n <- nrow(x)
    m <- nrow(y)
    if (n != m || n < 2) stop("Sample sizes must agree")
    if (! (all(is.finite(c(x, y)))))
        stop("Data contains missing or infinite values")

    stat <- reps <- 0
    if (R > 0) reps <- rep(0, R)
    pval <- 1
    dims <- c(n, ncol(x), ncol(y), R)

    a <- .C("indepEtest",
            x = as.double(t(x)),
            y = as.double(t(y)),
            byrow = as.integer(TRUE),
            dims = as.integer(dims),
            stat = as.double(stat),
            reps = as.double(reps),
            pval = as.double(pval),
            PACKAGE = "energy")

    stat <- n*a$stat
    est <- sqrt(a$stat)
    names(est) <- "I"
    names(stat) <- "nI^2"
    dataname <- paste("x (",n," by ",ncol(x), "), y(",n," by ", ncol(y), "), replicates ", R, sep="")
    if (R > 0)
      p.value = a$pval else p.value = NA
    e <- list(
        method = "mvI energy test of independence",
        statistic = stat,
        estimate = est,
        replicates = n*reps,
        p.value = p.value,
        data.name = dataname)
    class(e) <- "htest"
    e
}


