eqdist.e <-
function(x, sizes, distance = FALSE, method = c("original","discoB","discoF"))
{
    ## multivariate E-statistic for testing equal distributions
    ##   x:          matrix of pooled sample or distance matrix
    ##   sizes:      vector of sample sizes
    ##   distance:   logical, TRUE if x is a distance matrix, otherwise false
    ##   method:     original (default) or disco between components, or disco F ratio

    method <-match.arg(method)
    if (method=="discoB") {
        g <- as.factor(rep(1:length(sizes), sizes))
        RVAL <- disco(x, factors=g, distance=distance, R=0, method=method)
        } else {
			RVAL <- eqdist.etest(x, sizes, distance = distance, R=0, method=method)$statistic
		}
	RVAL
}

eqdist.etest <-
function(x, sizes, distance = FALSE, method = c("original","discoB","discoF"), R)
{
    ## multivariate E-test of the multisample hypothesis of equal distributions
    ##   x:          matrix of pooled sample or distance matrix
    ##   sizes:      vector of sample sizes
    ##   distance:   logical, TRUE if x is a distance matrix, otherwise false
    ##   method:     original (default) or disco components
    ##   R:          number of replicates
    ##

    method <-match.arg(method)

    if (method=="discoB" || method=="discoF") {
      g <- as.factor(rep(1:length(sizes), sizes))
	  # for other index use disco() function directly
      return(disco(x, factors=g, distance=distance, index=1.0, R=R, method=method))
      }

    nsamples <- length(sizes)
    if (nsamples < 2) return (NA)
    if (min(sizes) < 1) return (NA)
    if (!is.null(attr(x, "Size"))) distance <- TRUE

    x <- as.matrix(x)
    if (NROW(x) != sum(sizes)) stop("nrow(x) should equal sum(sizes)")
    if (distance == FALSE && nrow(x) == ncol(x))
        warning("square data matrix with distance==FALSE")
    d <- NCOL(x)
    if (distance == TRUE) d <- 0
    str <- "Multivariate "
    if (d == 1) str <- "Univariate "
    if (d == 0) str <- ""

    e0 <- 0.0
    repl <- rep(0, R)
    pval <- 1.0
    b <- .C("ksampleEtest",
        x = as.double(t(x)),
        byrow = as.integer(1),
        nsamples = as.integer(nsamples),
        sizes = as.integer(sizes),
        dim = as.integer(d),
        R = as.integer(R),
        e0 = as.double(e0),
        e = as.double(repl),
        pval = as.double(pval),
        PACKAGE = "energy")

    names(b$e0) <- "E-statistic"
    sz <- paste(sizes, collapse = " ", sep = "")
    methodname <- paste(str, length(sizes),
                  "-sample E-test of equal distributions", sep = "")
    dataname <- paste("sample sizes ", sz, ", replicates ", R, sep="")
    e <- list(
	    call = match.call(),
        method = methodname,
        statistic = b$e0,
        p.value = b$pval,
        data.name = dataname)

    class(e) <- "htest"
    e
}

ksample.e <-
function(x, sizes, distance = FALSE, method = c("original","discoB","discoF"),
         ix = 1:sum(sizes))
{
    ## computes k-sample E-statistics for equal distributions
    ## retained for backward compatibility or use with boot
    ## (this function simply passes arguments to eqdist.e)
    ##
    ##   x:          pooled sample or distance matrix
    ##   sizes:      vector of sample sizes
    ##   distance:   TRUE if x is a distance matrix, otherwise FALSE
    ##   method:     default (original) or disco between components or disco F ratio
    ##   ix:         a permutation of row indices of x
    ##
    x <- as.matrix(x)
    method <- match.arg(method)
    eqdist.e(x[ix,], sizes=sizes, distance=distance, method=method)
}

