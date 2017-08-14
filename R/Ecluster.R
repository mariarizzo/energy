
energy.hclust <-
function(dst, alpha = 1) {
    if (!(class(dst) == "dist")) 
      stop("The first argument must be a dist object.")
    d <- dst
    n <- attr(d, "Size")
    if (!isTRUE(all.equal(alpha, 1))) {
    	if (alpha > 2)
    	    warning("Exponent alpha should be in (0,2]")
      if (alpha < 0)
        stop("Cannot use negative exponent on distance.")
    	d <- d^alpha
    }
    labels <- attr(d, "Labels")
    if (is.null(labels))
        labels <- paste(1:n)
    merge <- integer(2 * (n - 1))
    height <- double(n - 1)
    order <- integer(n)
    ecl <- .C("Emin_hclust",
              diss = as.double(d),
              en = as.integer(n),
              merge = as.integer(merge),
              height = as.double(height),
              order = as.integer(order),
              PACKAGE = "energy")
    merge <- matrix(ecl$merge, nrow = n - 1, ncol = 2)
    e <- list(merge = merge,
              height = ecl$height,
              order = ecl$order,
              labels = labels,
              method = "e-distance",
              call = match.call(),
              dist.method = attr(dst, "method"))
    class(e) <- "hclust"
    e
}
