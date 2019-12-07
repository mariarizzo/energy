
energy.hclust <-
function(dst, alpha = 1) {
    if (!inherits(dst, "dist"))
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
    ## heights of hclust are half of energy; otherwise equivalent
    return(hclust(d, method = "ward.D"))
}

