dcov.test <-
function(x, y, index=1.0, R) {
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
    names(stat) <- "nV^2"
    names(V) <- "dCov"
    dataname <- paste("index ", index, ", replicates ", R, sep="")
    pval <- ifelse (R < 1, NA, a$pval)
    e <- list(
        method = paste("dCov test of independence", sep = ""),
        statistic = stat,
        estimate = V,
        estimates = dcorr,
        p.value = pval,
        replicates = n* a$reps^2,
        data.name = dataname)
    class(e) <- "htest"
    return(e)
}

.dcov <-
function(x, y, index=1.0) {
    # distance covariance statistic for independence
    # dcov = [dCov,dCor,dVar(x),dVar(y)]   (vector)
    # this function provides the fast method for computing dCov
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
    # distance correlation statistic for independence
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

